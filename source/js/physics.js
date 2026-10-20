// Event generation, simplified particle propagation, light/PMT response and approximate reconstruction.
// Pure JavaScript, no DOM. Everything is seeded so the same seed reproduces the same event.
import { GEO, PMT, OPT, region, ttWalls } from './detector.js';

// ------------------------------------------------------------------ random numbers
export function makeRng(seed) {
  let a = (seed >>> 0) || 1;
  const u = () => {
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
  const r = () => { let x; do x = u(); while (x === 0); return x; };
  r.gauss = () => Math.sqrt(-2 * Math.log(r())) * Math.cos(2 * Math.PI * r());
  r.exp = (mean) => -Math.log(r()) * mean;
  r.poisson = (mu) => {
    if (mu <= 0) return 0;
    if (mu > 40) return Math.max(0, Math.round(mu + Math.sqrt(mu) * r.gauss()));
    const L = Math.exp(-mu); let k = 0, p = 1;
    do { k++; p *= r(); } while (p > L);
    return k - 1;
  };
  r.iso = () => { const z = 2 * r() - 1, s = Math.sqrt(1 - z * z), p = 2 * Math.PI * r(); return [s * Math.cos(p), s * Math.sin(p), z]; };
  return r;
}

// ------------------------------------------------------------------ particle table
const T = (name, m, q, cls) => ({ name, m, q, cls });
export const PDG = {
  11: T('e⁻', 0.000511, -1, 'e'), [-11]: T('e⁺', 0.000511, 1, 'e'),
  13: T('μ⁻', 0.105658, -1, 'mu'), [-13]: T('μ⁺', 0.105658, 1, 'mu'),
  15: T('τ⁻', 1.77686, -1, 'mu'), [-15]: T('τ⁺', 1.77686, 1, 'mu'),
  12: T('νe', 0, 0, 'nu'), [-12]: T('ν̄e', 0, 0, 'nu'), 14: T('νμ', 0, 0, 'nu'), [-14]: T('ν̄μ', 0, 0, 'nu'),
  16: T('ντ', 0, 0, 'nu'), [-16]: T('ν̄τ', 0, 0, 'nu'),
  22: T('γ', 0, 0, 'gamma'), 111: T('π⁰', 0.134977, 0, 'pi0'), 221: T('η', 0.547862, 0, 'pi0'),
  211: T('π⁺', 0.139570, 1, 'chad'), [-211]: T('π⁻', 0.139570, -1, 'chad'),
  321: T('K⁺', 0.493677, 1, 'chad'), [-321]: T('K⁻', 0.493677, -1, 'chad'),
  130: T('K⁰L', 0.497611, 0, 'nhad'), 310: T('K⁰S', 0.497611, 0, 'nhad'), 311: T('K⁰', 0.497611, 0, 'nhad'), [-311]: T('K̄⁰', 0.497611, 0, 'nhad'),
  2212: T('p', 0.938272, 1, 'p'), [-2212]: T('p̄', 0.938272, -1, 'p'),
  2112: T('n', 0.939565, 0, 'n'), [-2112]: T('n̄', 0.939565, 0, 'n'),
  3122: T('Λ', 1.115683, 0, 'nhad'), 3222: T('Σ⁺', 1.18937, 1, 'chad'), 3112: T('Σ⁻', 1.197449, -1, 'chad'), 3212: T('Σ⁰', 1.192642, 0, 'nhad'),
};
export function pinfo(pdg) {
  if (PDG[pdg]) return PDG[pdg];
  if (Math.abs(pdg) > 1e9) return T('nucleus ' + pdg, 0, 0, 'nucleus');
  return T('pdg ' + pdg, 0, 0, 'other');
}
export const CLS_COLOR_KEY = { mu: 'mu', e: 'e', gamma: 'gamma', pi0: 'gamma', chad: 'chad', p: 'p', n: 'n', nhad: 'n', nu: 'nu', nucleus: 'n', other: 'n' };

// ------------------------------------------------------------------ vectors & kinematics
const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const mul = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const len = (a) => Math.sqrt(dot(a, a));
const norm = (a) => { const l = len(a) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
export const V = { add, sub, mul, dot, len, norm, cross };

function perp(u) { const a = Math.abs(u[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]; return norm(cross(u, a)); }
export function rotateAround(u, theta, phi) { // direction at polar angle theta, azimuth phi around u
  const e1 = perp(u), e2 = cross(u, e1), c = Math.cos(theta), s = Math.sin(theta);
  return norm(add(add(mul(u, c), mul(e1, s * Math.cos(phi))), mul(e2, s * Math.sin(phi))));
}
function boost(p4, b) { // p4=[E,px,py,pz], b = velocity vector
  const b2 = dot(b, b); if (b2 < 1e-14) return p4.slice();
  const g = 1 / Math.sqrt(1 - b2), bp = b[0] * p4[1] + b[1] * p4[2] + b[2] * p4[3];
  const g2 = (g - 1) / b2;
  return [g * (p4[0] + bp), p4[1] + g2 * bp * b[0] + g * b[0] * p4[0], p4[2] + g2 * bp * b[1] + g * b[1] * p4[0], p4[3] + g2 * bp * b[2] + g * b[2] * p4[0]];
}
const mass4 = (p) => Math.sqrt(Math.max(0, p[0] * p[0] - p[1] * p[1] - p[2] * p[2] - p[3] * p[3]));
function twoBody(P, m1, m2, rng) {
  const M = mass4(P);
  const pst = Math.sqrt(Math.max(0, (M * M - (m1 + m2) ** 2) * (M * M - (m1 - m2) ** 2))) / (2 * M);
  const d = rng.iso();
  const a = [Math.hypot(pst, m1), ...mul(d, pst)], b = [Math.hypot(pst, m2), ...mul(d, -pst)];
  const bv = [P[1] / P[0], P[2] / P[0], P[3] / P[0]];
  return [boost(a, bv), boost(b, bv)];
}
function nBody(P, masses, rng) { // sequential two-body decays; exact 4-momentum conservation
  if (masses.length === 1) return [P];
  const M = mass4(P), rest = masses.slice(1), restMin = rest.reduce((s, m) => s + m, 0);
  const Mrest = rest.length === 1 ? rest[0] : restMin + (M - masses[0] - restMin) * Math.pow(rng(), 0.6) * 0.999;
  const [a, b] = twoBody(P, masses[0], Mrest, rng);
  return [a, ...nBody(b, rest, rng)];
}

// ------------------------------------------------------------------ vertex & direction helpers
function sampleVertex(rng, o) {
  const mode = o.vertex || 'fv';
  if (mode === 'centre') return [0, 0, 0];
  if (mode === 'fixed') return [+o.vx || 0, +o.vy || 0, +o.vz || 0];
  const R = mode === 'ls' ? GEO.R_LS - 0.05 : GEO.FV_R;
  for (;;) {
    const v = [(2 * rng() - 1) * R, (2 * rng() - 1) * R, (2 * rng() - 1) * R];
    if (len(v) < R && (mode !== 'fv' || Math.abs(v[2]) < GEO.FV_Z)) return v;
  }
}
// zenith convention: cosZ = cosine of the zenith angle of the direction the particle COMES FROM.
// cosZ > 0: down-going (from the sky); cosZ < 0: up-going (through the Earth). Azimuth phi from +x.
export function dirFromZenith(cosZ, phiDeg) {
  const s = Math.sqrt(Math.max(0, 1 - cosZ * cosZ)), p = (phiDeg * Math.PI) / 180;
  return [-s * Math.cos(p), -s * Math.sin(p), -cosZ];
}
export function zenithOf(dir) {
  const u = norm(dir);
  return { cosZ: -u[2], phi: ((Math.atan2(-u[1], -u[0]) * 180) / Math.PI + 360) % 360 };
}

// ------------------------------------------------------------------ Honda (HAKKM) flux table
export function parseHonda(text) {
  const blocks = []; let cur = null;
  for (const line of text.split(/\r?\n/)) {
    const m = line.match(/cosZ\s*=\s*([-\d.]+)\s*--\s*([-\d.]+)(?:.*?phi_Az\s*=\s*([-\d.]+)\s*--\s*([-\d.]+))?/);
    if (m) { cur = { c0: +m[1], c1: +m[2], p0: m[3] ? +m[3] : 0, p1: m[4] ? +m[4] : 360, rows: [] }; blocks.push(cur); continue; }
    if (!cur) continue;
    const f = line.trim().split(/\s+/).map(Number);
    if (f.length >= 5 && f.every(Number.isFinite)) cur.rows.push(f.slice(0, 5));
  }
  if (!blocks.length || !blocks[0].rows.length) throw new Error('Could not find Honda flux blocks ("average flux in [cosZ = ...]").');
  return blocks;
}
const HONDA_COL = { 14: 1, [-14]: 2, 12: 3, [-12]: 4 };
export function buildFluxSampler(blocks, flavours, Emin, Emax, xsecWeight = true) {
  const cells = []; let tot = 0;
  for (const b of blocks) {
    const E = b.rows.map(r => r[0]);
    for (let i = 0; i < E.length; i++) {
      const lo = i ? Math.sqrt(E[i] * E[i - 1]) : E[i], hi = i < E.length - 1 ? Math.sqrt(E[i] * E[i + 1]) : E[i];
      const a = Math.max(lo, Emin), c = Math.min(hi, Emax); if (c <= a) continue;
      for (const f of flavours) {
        const w = b.rows[i][HONDA_COL[f]] * (c - a) * (xsecWeight ? E[i] : 1) * (b.c1 - b.c0) * (b.p1 - b.p0);
        if (w > 0) { tot += w; cells.push({ cum: tot, f, a, c, b }); }
      }
    }
  }
  if (!cells.length) throw new Error('No flux in the chosen energy range.');
  return (rng) => {
    const x = rng() * tot; let lo = 0, hi = cells.length - 1;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (cells[mid].cum < x) lo = mid + 1; else hi = mid; }
    const c = cells[lo];
    return { pdg: c.f, E: c.a * Math.pow(c.c / c.a, rng()), cosZ: c.b.c0 + (c.b.c1 - c.b.c0) * rng(), phi: c.b.p0 + (c.b.p1 - c.b.p0) * rng() };
  };
}

// ------------------------------------------------------------------ generators
const MN = 0.9389;

function sampleE(rng, o, def) {
  if (o.emode === 'fixed') return +o.E || def;
  const g = +o.gamma || 2.7, a = +o.Emin || 0.1, b = +o.Emax || 100;
  if (Math.abs(g - 1) < 1e-6) return a * Math.pow(b / a, rng());
  const k = 1 - g; return Math.pow(Math.pow(a, k) + rng() * (Math.pow(b, k) - Math.pow(a, k)), 1 / k);
}

// Approximate shares of the CC cross section per process vs energy (shapes follow the usual
// QE / 2p2h / resonance / DIS / coherent picture; illustrative, not a GENIE tune) [D]
export const MODES_ALL = ['QE', 'MEC', 'RES', 'DIS', 'COH'];
export function processShares(E) {
  const qe = 0.95 * (1 - Math.exp(-E / 0.25)) / (1 + 0.04 * E);
  const mec = 0.28 * qe * (1 - Math.exp(-E / 0.5));
  const res = E > 0.3 ? 0.55 * (1 - Math.exp(-(E - 0.3) / 0.6)) : 0;
  const dis = E > 1.0 ? 0.68 * E * (1 - Math.exp(-(E - 1.0) / 1.5)) : 0;
  const coh = E > 0.25 ? 0.03 * (1 - Math.exp(-(E - 0.25) / 1.5)) : 0;
  const t = qe + mec + res + dis + coh;
  return { QE: qe / t, MEC: mec / t, RES: res / t, DIS: dis / t, COH: coh / t };
}
function chooseMode(rng, E, req) {
  if (req && req !== 'auto') return req;
  const f = processShares(E); let x = rng();
  for (const m of MODES_ALL) { x -= f[m]; if (x <= 0) return m; }
  return 'QE';
}

function hadronCharges(Qhad, nPi, rng) {
  // nucleon + nPi pions with total charge Qhad
  for (let tries = 0; tries < 200; tries++) {
    const nuc = rng() < 0.5 ? 2212 : 2112; let q = nuc === 2212 ? 1 : 0;
    const pis = [];
    for (let i = 0; i < nPi; i++) { const r = rng(); const c = r < 0.33 ? 211 : r < 0.66 ? 111 : -211; pis.push(c); q += PDG[c].q; }
    let i = 0;
    while (q !== Qhad && i < 3 * nPi) {
      const k = i % nPi; const cur = PDG[pis[k]].q;
      if (q < Qhad && cur < 1) { pis[k] = cur === 0 ? 211 : 111; q++; } else if (q > Qhad && cur > -1) { pis[k] = cur === 0 ? -211 : 111; q--; }
      i++;
    }
    if (q === Qhad) return [nuc, ...pis];
  }
  return [2212, ...Array(nPi).fill(111)];
}

function nuNucleon(rng, nuPdg, E, cc, modeReq, dirNu) {
  const isNu = nuPdg > 0;
  const lep = cc ? (Math.abs(nuPdg) - 1) * Math.sign(nuPdg) : nuPdg; // 14->13, -14->-13, 12->11
  const ml = PDG[lep].m;
  let mode = chooseMode(rng, E, modeReq);
  if (mode === 'DIS' && E < 1.2) mode = 'RES';
  if (mode === 'RES' && E < 0.35) mode = 'QE';
  if (mode === 'COH') { const k = cohPion(rng, nuPdg, E, cc, dirNu, lep, ml); if (k) return k; mode = 'QE'; }
  if (mode === 'MEC') { const k = mecPair(rng, nuPdg, E, cc, dirNu, lep, ml); if (k) return k; mode = 'QE'; }
  let target = rng() < 0.5 ? 2212 : 2112;
  if (cc && mode === 'QE') target = isNu ? 2112 : 2212;
  const Qhad = PDG[target].q + (cc ? (isNu ? 1 : -1) : 0);
  for (let tries = 0; tries < 400; tries++) {
    let W, Q2, x;
    if (mode === 'QE') { W = isNu && cc ? PDG[2212].m : cc ? PDG[2112].m : PDG[target].m; Q2 = rng.exp(0.35); }
    else if (mode === 'RES') { W = 1.232 + 0.06 * Math.tan(Math.PI * (rng() - 0.5)); if (W < 1.09 || W > 1.75) continue; Q2 = rng.exp(0.45); }
    else { const Wmax = Math.sqrt(MN * MN + 2 * MN * E) - 0.15; if (Wmax < 1.5) { mode = 'RES'; continue; } W = 1.5 + (Wmax - 1.5) * Math.pow(rng(), 1.6); x = 0.02 + 0.75 * Math.pow(rng(), 1.4); Q2 = (x * (W * W - MN * MN)) / (1 - x); }
    const nu = (W * W - MN * MN + Q2) / (2 * MN), El = E - nu;
    if (El <= ml * 1.001) continue;
    const pl = Math.sqrt(El * El - ml * ml), cosT = (2 * E * El - Q2 - ml * ml) / (2 * E * pl);
    if (cosT < -1 || cosT > 1) continue;
    const ul = rotateAround(dirNu, Math.acos(cosT), 2 * Math.PI * rng());
    const plv = mul(ul, pl), q = sub(mul(dirNu, E), plv);
    const Phad = [nu + MN, q[0], q[1], q[2]], Wreal = mass4(Phad);
    let pdgs;
    if (mode === 'QE') pdgs = [Qhad === 1 ? 2212 : Qhad === 0 ? 2112 : null];
    else if (mode === 'RES') {
      const opts = { 2: [[2212, 211]], 1: [[2212, 111], [2112, 211]], 0: [[2112, 111], [2212, -211]], [-1]: [[2112, -211]] }[Qhad];
      pdgs = opts[Math.floor(rng() * opts.length)];
    } else {
      let nPi = Math.max(1, rng.poisson(0.5 + 1.1 * Math.log(Wreal * Wreal)));
      while (nPi > 1 && MN + nPi * 0.14 > Wreal - 0.02) nPi--;
      pdgs = hadronCharges(Qhad, nPi, rng);
    }
    if (pdgs.includes(null)) { mode = 'RES'; continue; }
    const masses = pdgs.map(p => PDG[p].m);
    if (masses.length > 1 && masses.reduce((s, m) => s + m, 0) > Wreal - 0.002) continue;
    const had = nBody(Phad, masses, rng);
    const y = nu / E; if (x === undefined) x = Q2 / (2 * MN * nu);
    return { lep, El, plv, mode, target, Q2, W: Wreal, x, y, hadrons: pdgs.map((p, i) => ({ pdg: p, p4: had[i] })) };
  }
  return null;
}

function leptonFor(rng, E, nu, Q2, ml, dirNu) {
  const El = E - nu; if (El <= ml * 1.001) return null;
  const pl = Math.sqrt(El * El - ml * ml), cosT = (2 * E * El - Q2 - ml * ml) / (2 * E * pl);
  if (cosT < -1 || cosT > 1) return null;
  const ul = rotateAround(dirNu, Math.acos(cosT), 2 * Math.PI * rng());
  const plv = mul(ul, pl); return { El, plv, q: sub(mul(dirNu, E), plv) };
}
// 2p2h / MEC: lepton + nucleon pair, knocked out of a correlated pair at rest [D]
function mecPair(rng, nuPdg, E, cc, dirNu, lep, ml) {
  const M2 = 2 * MN, isNu = nuPdg > 0;
  const Q = 1 + (cc ? (isNu ? 1 : -1) : 0); // target pair n p
  const pair = Q === 2 ? [2212, 2212] : Q === 0 ? [2112, 2112] : [2212, 2112];
  for (let t = 0; t < 300; t++) {
    const W = M2 + 0.02 + 0.28 * rng(), Q2 = rng.exp(0.4), nu = (W * W - M2 * M2 + Q2) / (2 * M2);
    const L = leptonFor(rng, E, nu, Q2, ml, dirNu); if (!L) continue;
    const Phad = [nu + M2, ...L.q], had = twoBody(Phad, PDG[pair[0]].m, PDG[pair[1]].m, rng);
    if (!had[0].every(Number.isFinite)) continue;
    return { lep, El: L.El, plv: L.plv, mode: 'MEC', target: 1000060120, Q2, W: mass4(Phad), x: Q2 / (2 * MN * nu), y: nu / E, hadrons: pair.map((p, i) => ({ pdg: p, p4: had[i] })) };
  }
  return null;
}
// coherent pion production on 12C: forward lepton + pion, nucleus left intact (recoil neglected) [D]
function cohPion(rng, nuPdg, E, cc, dirNu, lep, ml) {
  const pi = cc ? (nuPdg > 0 ? 211 : -211) : 111, mpi = PDG[pi].m;
  for (let t = 0; t < 300; t++) {
    const Epi = mpi + 0.02 + (E - mpi - ml - 0.04) * (0.05 + 0.9 * rng()), Q2 = rng.exp(0.05);
    if (Epi <= mpi) continue;
    const L = leptonFor(rng, E, Epi, Q2, ml, dirNu); if (!L) continue;
    const ppi = Math.sqrt(Epi * Epi - mpi * mpi), dir = norm(L.q);
    return { lep, El: L.El, plv: L.plv, mode: 'COH', target: 1000060120, Q2, W: mpi, x: null, y: Epi / E, hadrons: [{ pdg: pi, p4: [Epi, ...mul(dir, ppi)] }] };
  }
  return null;
}

export const CATEGORIES = {
  custom: { label: 'Custom ν · aim & energy', kind: 'atm', flav: 14, cur: 'CC', custom: true },
  atm_numu: { label: 'Atmospheric νμ CC', kind: 'atm', flav: 14, cur: 'CC' },
  atm_nue: { label: 'Atmospheric νe CC', kind: 'atm', flav: 12, cur: 'CC' },
  atm_nc: { label: 'Atmospheric NC', kind: 'atm', flav: 14, cur: 'NC' },
  ibd: { label: 'Reactor IBD  ν̄e p → e⁺ n', kind: 'ibd' },
  cosmic: { label: 'Cosmic muon', kind: 'cosmic' },
};

export function generate(cat, seed, o = {}) {
  const rng = makeRng(seed * 2654435761 + 12345);
  const c = CATEGORIES[cat];
  if (c.kind === 'atm') return genAtm(rng, seed, cat, o);
  if (c.kind === 'ibd') return genIbd(rng, seed, o);
  return genCosmic(rng, seed, o);
}

function genAtm(rng, seed, cat, o) {
  const c = CATEGORIES[cat];
  let flav = o.flav && o.flav !== 'cat' ? +o.flav : c.flav;
  let E, cosZ, phi;
  if (c.custom) { o = { ...o, emode: 'fixed', E: o.cE, dmode: 'fixed', cosZ: Math.cos((o.cZen * Math.PI) / 180), phi: o.cAz, mode: o.cMode, cur: o.cCur }; flav = +o.cFlav; }
  if (o.dmode === 'flux' && o.fluxSampler) {
    const s = o.fluxSampler(rng); flav = s.pdg; E = s.E; cosZ = s.cosZ; phi = s.phi;
  } else {
    E = sampleE(rng, o, 4);
    if (o.dmode === 'fixed') { cosZ = +o.cosZ; phi = +o.phi || 0; } else { cosZ = 2 * rng() - 1; phi = 360 * rng(); }
  }
  let cc = o.cur && o.cur !== 'cat' ? o.cur === 'CC' : c.cur === 'CC';
  if (o.cur === 'auto') cc = rng() < 0.72;
  const dir = dirFromZenith(cosZ, phi), vtx = sampleVertex(rng, o);
  let k = null;
  for (let t = 0; t < 30 && !k; t++) k = nuNucleon(rng, flav, E, cc, o.mode, dir);
  if (!k) { E = Math.max(E, 1); k = nuNucleon(rng, flav, E, cc, 'QE', dir); }
  const parts = [{ pdg: k.lep, p4: [k.El, ...k.plv] }, ...k.hadrons];
  // illustrative nuclear de-excitation / FSI neutrons (not energy-conserving, flagged as 'nuclear')
  const nExtra = k.mode === 'COH' ? 0 : rng.poisson(0.4 + 0.25 * Math.log(1 + E)); // coherent: nucleus stays intact
  for (let i = 0; i < nExtra; i++) { const T = 0.003 + rng.exp(0.012); parts.push({ pdg: 2112, p4: [MN + T, ...mul(rng.iso(), Math.sqrt(T * T + 2 * T * MN))], origin: 'nuclear' }); }
  return {
    source: 'demo', category: cat, seed,
    nu: { pdg: flav, E, dir, cosZ, phi },
    interaction: { cc, mode: k.mode, target: 1000060120, hitnuc: k.target, Q2: k.Q2, W: k.W, x: k.x, y: k.y },
    vertex: vtx, t0: 0, primaries: parts,
  };
}

function reactorEnu(rng) {
  for (;;) {
    const E = 1.81 + rng() * 8.2, Ee = E - 1.293, pe = Math.sqrt(Math.max(0, Ee * Ee - 0.511 * 0.511));
    const f = Ee * pe * Math.exp(0.87 - 0.16 * E - 0.091 * E * E);
    if (rng() * 1.5 < f) return E;
  }
}

function genIbd(rng, seed, o) {
  const E = o.ibdMode === 'fixed' ? Math.max(1.81, +o.ibdE || 4) : reactorEnu(rng); // MeV
  let cosZ, phi;
  if (o.dmode === 'fixed') { cosZ = +o.cosZ; phi = +o.phi || 0; } else { cosZ = 0.05 * rng(); phi = 360 * rng(); } // reactors ~horizontal
  const dir = dirFromZenith(cosZ, phi), vtx = sampleVertex(rng, o);
  const Ee = (E - 1.293) / 1000, me = 0.000511, pe = Math.sqrt(Math.max(0, Ee * Ee - me * me));
  const ue = rotateAround(dir, Math.acos(Math.max(-1, Math.min(1, -0.1 + 1.1 * (2 * rng() - 1)))), 2 * Math.PI * rng());
  const pn = sub(mul(dir, E / 1000), mul(ue, pe)), mn = PDG[2112].m;
  return {
    source: 'demo', category: 'ibd', seed,
    nu: { pdg: -12, E: E / 1000, dir, cosZ, phi },
    interaction: { cc: true, mode: 'IBD', target: 2212, hitnuc: 2212, Q2: null, W: null, x: null, y: null },
    vertex: vtx, t0: 0,
    primaries: [{ pdg: -11, p4: [Ee, ...mul(ue, pe)] }, { pdg: 2112, p4: [Math.hypot(len(pn), mn), ...pn] }],
  };
}

function genCosmic(rng, seed, o) {
  const E = o.muMode === 'fixed' ? +o.muE || 200 : 15 + rng.exp(190) * (0.6 + 0.8 * rng()); // GeV, mean ~200 GeV underground [D]
  let cosZ, phi;
  if (o.dmode === 'fixed') { cosZ = Math.max(0.05, +o.cosZ); phi = +o.phi || 0; } else { cosZ = Math.max(0.12, Math.pow(rng(), 1 / 3.5)); phi = 360 * rng(); }
  const dir = dirFromZenith(cosZ, phi);
  const bmax = o.through === false ? 24 : 16.8;
  const e1 = perp(dir), e2 = cross(dir, e1), r = bmax * Math.sqrt(rng()), a = 2 * Math.PI * rng();
  const P = add(mul(e1, r * Math.cos(a)), mul(e2, r * Math.sin(a)));
  let s = (29 - P[2]) / -dir[2]; // back-propagate to z = 29 m, but stay inside the 76 m world box
  let start = add(P, mul(dir, -s));
  while ((Math.abs(start[0]) > 37 || Math.abs(start[1]) > 37) && s > 0) { s -= 0.5; start = add(P, mul(dir, -s)); }
  const m = PDG[13].m, p = Math.sqrt(E * E - m * m), q = rng() < 0.56 ? -13 : 13; // mu+/mu- ~1.27 [D]
  return {
    source: 'demo', category: 'cosmic', seed,
    nu: null, muon: { E, dir, cosZ, phi, impact: r },
    interaction: { cc: null, mode: 'cosmic μ', target: null },
    vertex: start, t0: 0, cosmic: true,
    primaries: [{ pdg: q, p4: [E, ...mul(dir, p)] }],
    spallation: rng.poisson(o.spall === false ? 0 : 1.2),
  };
}

// GENIE gst row -> event
export function fromGst(row, idx, fileName, unitScale = 1) {
  const nu = { pdg: row.neu, E: row.Ev, dir: norm([row.pxv, row.pyv, row.pzv]) };
  Object.assign(nu, zenithOf(nu.dir));
  const parts = [];
  const lep = row.fspl;
  if (lep && row.El > 0) parts.push({ pdg: lep, p4: [row.El, row.pxl, row.pyl, row.pzl] });
  const n = row.nf || 0;
  for (let i = 0; i < n; i++) parts.push({ pdg: row.pdgf[i], p4: [row.Ef[i], row.pxf[i], row.pyf[i], row.pzf[i]] });
  const mode = row.qel ? (row.mec ? 'MEC' : 'QE') : row.mec ? 'MEC' : row.res ? 'RES' : row.dis ? 'DIS' : row.coh ? 'COH' : 'other';
  return {
    source: 'genie', category: row.cc ? (Math.abs(row.neu) === 12 ? 'atm_nue' : 'atm_numu') : 'atm_nc', seed: idx, file: fileName, iev: row.iev,
    nu, interaction: { cc: !!row.cc, mode, target: row.tgt, Q2: row.Q2, W: row.W, x: row.x, y: row.y },
    vertex: [row.vtxx * unitScale, row.vtxy * unitScale, row.vtxz * unitScale], t0: 0, primaries: parts,
  };
}
export function guessUnitScale(rows) {
  let m = 0; for (const r of rows) m = Math.max(m, Math.hypot(r.vtxx, r.vtxy, r.vtxz));
  return m > 200 ? 0.01 : m > 1e4 ? 0.001 : 1; // gst vertices are SI (m); fall back if cm/mm detected
}

// ------------------------------------------------------------------ propagation
const X0 = 0.52, EC = 0.09, LAMBDA_I = 0.95; // [D] LS radiation length (m), critical energy (GeV), interaction length (m)
const THETA_C = Math.acos(1 / OPT.N_LS);

function makeCtx(ev, rng) {
  return { ev, rng, particles: [], tracks: [], deps: [], wcd: [], tt: [], windows: [{ label: 'Prompt', t0: 0 }], nid: 0 };
}
function windowFor(ctx, tAbs, label) {
  for (let i = 0; i < ctx.windows.length; i++) if (Math.abs(tAbs - ctx.windows[i].t0) < 900) return i;
  ctx.windows.push({ label, t0: tAbs }); return ctx.windows.length - 1;
}
function inWorld(x) { return Math.abs(x[0]) < 40 && Math.abs(x[1]) < 40 && x[2] > -26 && x[2] < 32; }

function deposit(ctx, P, x, t, edepMeV, q = 1, cher = 0, u = null, sig = 8) {
  if (edepMeV <= 0) return;
  const reg = region(x[0], x[1], x[2]);
  P.edepAll += edepMeV;
  if (reg !== 'LS') return;
  P.edepLS += edepMeV; P.evis += edepMeV * q;
  ctx.deps.push({ x: x[0], y: x[1], z: x[2], t, edep: edepMeV, evis: edepMeV * q, cher, u, sig, pid: P.id, cls: P.cls });
}

export function addParticle(ctx, d) {
  const info = pinfo(d.pdg);
  const p3 = [d.p4[1], d.p4[2], d.p4[3]];
  const P = {
    id: ctx.nid++, pdg: d.pdg, name: info.name, cls: info.cls, q: info.q, m: info.m,
    E: d.p4[0], p: p3, pmag: len(p3), x0: d.x0.slice(), t0: d.t0, parent: d.parent ?? -1,
    origin: d.origin || (d.parent === undefined || d.parent < 0 ? 'primary' : 'secondary'),
    edepLS: 0, edepAll: 0, evis: 0, lenLS: 0, end: null, endProcess: '', track: null,
  };
  ctx.particles.push(P);
  if (ctx.particles.length > 1500) return P;
  const fn = PROP[info.cls]; if (fn) fn(ctx, P);
  return P;
}
function track(ctx, P, pts, style) { const tr = { pid: P.id, cls: P.cls, pts, style }; ctx.tracks.push(tr); if (!P.track) P.track = tr; return tr; }

function stepCharged(ctx, P, dedx, opts = {}) {
  // straight charged track with continuous losses; returns {x, t, T, stopped}
  let x = P.x0.slice(), t = P.t0, T = P.E - P.m; const u = norm(P.p);
  const pts = [[...x, t]]; const ds = opts.ds || 0.3; let seenWater = false;
  const q = opts.quench ?? 0.97;
  for (let i = 0; i < 2000; i++) {
    const reg = region(x[0], x[1], x[2]);
    if (reg === 'ROCK' && (seenWater || !(ctx.ev.cosmic && P.parent < 0))) { P.endProcess = 'exits'; break; }
    if (!inWorld(x)) { P.endProcess = 'exits'; break; }
    if (reg !== 'ROCK' && reg !== 'AIR') seenWater = true;
    const k = reg === 'LS' ? 1 : reg === 'BUF' || reg === 'WP' ? 1.14 : reg === 'AIR' ? 0.001 : 1.5;
    const rad = P.cls === 'mu' ? 1 + 0.12 * Math.log10(Math.max(1, T)) : 1;
    let step = ds, dE = dedx * k * rad * step;
    let stop = false;
    if (dE >= T) { step = T / (dedx * k * rad); dE = T; stop = true; }
    const gam = (T + P.m) / P.m, beta = Math.sqrt(1 - 1 / (gam * gam));
    const mid = add(x, mul(u, step / 2));
    const cher = beta > 1 / OPT.N_LS ? OPT.CHER_FRAC * dE * 1000 : 0;
    deposit(ctx, P, mid, t + step / 2 / (beta * OPT.C), dE * 1000, q, cher, u, 6);
    if (reg === 'LS') P.lenLS += step;
    if ((reg === 'WP') && beta > 0.75) ctx.wcd.push({ x: mid, t, L: step, u, pid: P.id });
    if (opts.onStep) { const r = opts.onStep(mid, t, T, step); if (r) { x = mid; T -= dE / 2; pts.push([...x, t]); P.endProcess = r; return { x, t, T, stopped: false, pts }; } }
    x = add(x, mul(u, step)); t += step / (beta * OPT.C); T -= dE;
    if (opts.tt) crossTT(ctx, P, pts[pts.length - 1], [...x, t]);
    pts.push([...x, t]);
    if (stop || T <= 1e-5) { P.endProcess = 'stops'; return { x, t, T: 0, stopped: true, pts }; }
  }
  return { x, t, T, stopped: false, pts };
}

const TT_WALLS = ttWalls();
function crossTT(ctx, P, a, b) {
  for (const z of GEO.TT_Z) {
    if ((a[2] - z) * (b[2] - z) > 0) continue;
    const f = (a[2] - z) / (a[2] - b[2]), x = a[0] + f * (b[0] - a[0]), y = a[1] + f * (b[1] - a[1]);
    const w = TT_WALLS.find(w => w.z === z && Math.abs(x - w.x) < GEO.TT_WALL / 2 && Math.abs(y - w.y) < GEO.TT_WALL / 2);
    if (w) ctx.tt.push({ x, y, z, layer: w.layer, wall: `${w.ix},${w.iy}`, pid: P.id });
  }
}

function muCapture(ctx, P, x, t) {
  const dt = ctx.rng.exp(2026);
  P.endProcess = `stops, captured by ¹²C after ${(dt / 1000).toFixed(2)} µs`;
  windowFor(ctx, t + dt, 'μ⁻ nuclear capture');
  const nN = ctx.rng() < 0.6 ? 1 : 2;
  for (let i = 0; i < nN; i++) {
    const Tn = (2 + ctx.rng.exp(12)) / 1000;
    addParticle(ctx, { pdg: 2112, p4: [PDG[2112].m + Tn, ...mul(ctx.rng.iso(), Math.sqrt(Tn * Tn + 2 * Tn * PDG[2112].m))], x0: x, t0: t + dt, parent: P.id, origin: 'μ⁻ capture' });
  }
}
function michel(ctx, P, x, t, positive, tau = 2197) {
  const dt = ctx.rng.exp(tau);
  P.endProcess = `stops, decays after ${(dt / 1000).toFixed(2)} µs`;
  let e; do { e = ctx.rng(); } while (ctx.rng() > e * e * (3 - 2 * e)); // Michel spectrum
  const E = 0.0528 * e, w = windowFor(ctx, t + dt, 'Michel e (decay)');
  const pmag = Math.sqrt(Math.max(0, E * E - 0.000511 ** 2));
  addParticle(ctx, { pdg: positive ? -11 : 11, p4: [E, ...mul(ctx.rng.iso(), pmag)], x0: x, t0: t + dt, parent: P.id, origin: 'decay', win: w });
}

const PROP = {
  mu(ctx, P) {
    const r = stepCharged(ctx, P, 0.175, { tt: true, ds: P.E > 5 ? 0.4 : 0.25, quench: 0.98 });
    track(ctx, P, r.pts, 'solid'); P.end = r.x;
    if (r.stopped && region(r.x[0], r.x[1], r.x[2]) !== 'ROCK' && Math.abs(P.pdg) === 13) {
      // a stopped muon sits at rest (μ⁻ bound in a carbon atom) until it decays or, for μ⁻, is captured by the nucleus
      if (P.pdg === 13 && ctx.rng() < 0.078) muCapture(ctx, P, r.x, r.t);
      else michel(ctx, P, r.x, r.t, P.pdg < 0, P.pdg === 13 ? 2026 : 2197);
    }
    if (ctx.ev.cosmic && P.parent < 0) cosmicExtras(ctx, P);
  },
  e(ctx, P) {
    const T = P.E - P.m;
    if (T > 0.05) { shower(ctx, P, P.x0, norm(P.p), P.E, P.t0, false); }
    else {
      const R = (T * 1000) / 171, n = Math.max(2, Math.min(8, Math.ceil(R / 0.03)));
      let x = P.x0.slice(), u = norm(P.p); const pts = [[...x, P.t0]];
      for (let i = 0; i < n; i++) {
        const step = R / n, nx = add(x, mul(u, step));
        deposit(ctx, P, add(x, mul(u, step / 2)), P.t0, (T * 1000) / n, 1, T * 1000 > 0.4 ? 0.04 * (T * 1000) / n : 0, u, 25);
        x = nx; pts.push([...x, P.t0]); u = rotateAround(u, 0.35 * Math.abs(ctx.rng.gauss()), 2 * Math.PI * ctx.rng());
      }
      track(ctx, P, pts, 'solid'); P.end = x; P.lenLS = region(...P.x0) === 'LS' ? R : 0;
      if (P.pdg === -11) { // annihilation
        const d = ctx.rng.iso();
        for (const s of [1, -1]) addParticle(ctx, { pdg: 22, p4: [0.000511, ...mul(d, 0.000511 * s)], x0: x, t0: P.t0 + 0.2, parent: P.id, origin: 'annihilation' });
        P.endProcess = 'annihilates';
      } else P.endProcess = 'stops';
    }
  },
  gamma(ctx, P) {
    const u = norm(P.p);
    if (P.E > 0.01) {
      const d = ctx.rng.exp(X0 * 9 / 7), x = add(P.x0, mul(u, d));
      track(ctx, P, [[...P.x0, P.t0], [...x, P.t0 + d / OPT.C]], 'dashed');
      if (region(...x) !== 'ROCK') shower(ctx, P, x, u, P.E, P.t0 + d / OPT.C, true, true);
      P.end = x; P.endProcess = 'converts';
      return;
    }
    // MeV gamma: Compton chain
    let Eg = P.E * 1000, x = P.x0.slice(), dir = u, t = P.t0; const pts = [[...x, t]];
    for (let k = 0; k < 7 && Eg > 0.015; k++) {
      const lam = Math.max(0.08, 0.17 * Math.pow(Eg, 0.45));
      const d = ctx.rng.exp(lam); x = add(x, mul(dir, d)); t += d / OPT.C; pts.push([...x, t]);
      if (region(...x) === 'ROCK') break;
      const last = Eg < 0.12 || k === 6 || ctx.rng() < 0.08;
      const f = last ? 1 : 0.15 + 0.55 * ctx.rng();
      const dep = Eg * f; deposit(ctx, P, x, t, dep, 1, dep > 0.3 ? 0.03 * dep : 0, dir, 30);
      const E2 = Eg - dep; if (last || E2 < 0.015) { Eg = 0; break; }
      const cosT = Math.max(-1, Math.min(1, 1 - 0.511 * (1 / E2 - 1 / Eg)));
      dir = rotateAround(dir, Math.acos(cosT), 2 * Math.PI * ctx.rng()); Eg = E2;
    }
    track(ctx, P, pts, 'dashed'); P.end = x; P.endProcess = 'absorbed';
  },
  pi0(ctx, P) {
    const [a, b] = twoBody([P.E, ...P.p], 0, 0, ctx.rng);
    P.end = P.x0; P.endProcess = 'decays γγ';
    for (const g of [a, b]) addParticle(ctx, { pdg: 22, p4: g, x0: P.x0, t0: P.t0, parent: P.id, origin: 'decay' });
  },
  chad(ctx, P) {
    const DEC = { 3222: [2212, 111], 3112: [2112, -211] };
    const ctau = { 3222: 0.024, 3112: 0.0443 }[Math.abs(P.pdg)];
    let decayAt = Infinity;
    if (ctau) decayAt = ctx.rng.exp(ctau * P.pmag / P.m);
    let travelled = 0;
    const r = stepCharged(ctx, P, 0.18, {
      ds: 0.2, quench: 0.95,
      onStep: (x, t, T, step) => {
        travelled += step;
        if (travelled > decayAt) return 'decays';
        if (T > 0.08 && ctx.rng() < 1 - Math.exp(-step / LAMBDA_I)) return 'interacts';
        return null;
      },
    });
    track(ctx, P, r.pts, 'solid'); P.end = r.x;
    if (P.endProcess === 'decays' && ctau) {
      const d = DEC[Math.abs(P.pdg)].map(p => (P.pdg < 0 ? -p : p));
      const T = r.T, E = T + P.m, pm = Math.sqrt(Math.max(0, E * E - P.m * P.m));
      const prod = twoBody([E, ...mul(norm(P.p), pm)], PDG[Math.abs(d[0])].m, PDG[Math.abs(d[1])].m, ctx.rng);
      prod.forEach((p4, i) => addParticle(ctx, { pdg: d[i], p4, x0: r.x, t0: r.t, parent: P.id, origin: 'decay' }));
    } else if (P.endProcess === 'interacts') hadBlob(ctx, P, r.x, r.t, r.T);
    else if (r.stopped) {
      if (P.pdg === 211 || P.pdg === 321) michel(ctx, P, r.x, r.t, true);
      else hadBlob(ctx, P, r.x, r.t, 0.02);
    }
  },
  p(ctx, P) {
    const T = (P.E - P.m) * 1000; // MeV
    if (P.pdg < 0) { hadBlob(ctx, P, P.x0, P.t0, P.E - P.m + 1.876); P.end = P.x0; P.endProcess = 'annihilates'; track(ctx, P, [[...P.x0, P.t0]], 'solid'); return; }
    const R = (0.0022 * Math.pow(T, 1.77)) / 0.856 / 100;
    let L = R, interact = false;
    if (T > 300 && ctx.rng() < 1 - Math.exp(-R / LAMBDA_I)) { L = ctx.rng() * R; interact = true; }
    const u = norm(P.p), n = Math.max(2, Math.min(24, Math.ceil(L / 0.08)));
    const pts = [[...P.x0, P.t0]]; let Tl = T, t = P.t0;
    const w = []; for (let i = 0; i < n; i++) { const s = (i + 0.5) / n * L; w.push(Math.pow(R - s + 0.02 * R, -0.45)); }
    const Elost = interact ? T - (0.0022 > 0 ? T * Math.pow(1 - L / R, 1 / 1.77) : 0) : T;
    const ws = w.reduce((a, b) => a + b, 0);
    for (let i = 0; i < n; i++) {
      const dE = (Elost * w[i]) / ws, x = add(P.x0, mul(u, ((i + 0.5) / n) * L));
      const q = Math.max(0.25, Math.min(0.9, 0.9 - 0.65 * Math.exp(-Tl / 25)));
      deposit(ctx, P, x, t, dE, q, Tl > 480 ? 0.03 * dE : 0, u, 10);
      Tl -= dE; t += L / n / (0.3 * OPT.C); pts.push([...add(P.x0, mul(u, ((i + 1) / n) * L)), t]);
    }
    const end = add(P.x0, mul(u, L));
    if (region(...P.x0) === 'LS') P.lenLS = L;
    track(ctx, P, pts, 'solid'); P.end = end; P.endProcess = interact ? 'interacts' : 'stops';
    if (interact) hadBlob(ctx, P, end, t, Math.max(0, Tl) / 1000);
  },
  n(ctx, P) { neutronWalk(ctx, P, P.x0, P.t0, (P.E - P.m) * 1000 + (P.pdg < 0 ? 0 : 0)); if (P.pdg < 0) hadBlob(ctx, P, P.end, P.t0, 1.876); },
  nhad(ctx, P) {
    const u = norm(P.p);
    if (Math.abs(P.pdg) === 310 || P.pdg === 3122) {
      const dec = P.pdg === 310 ? [211, -211] : [2212, -211], ctau = P.pdg === 310 ? 0.0268 : 0.0789;
      const d = Math.min(3, ctx.rng.exp(ctau * P.pmag / P.m)), x = add(P.x0, mul(u, d));
      track(ctx, P, [[...P.x0, P.t0], [...x, P.t0 + d / OPT.C]], 'dashed'); P.end = x; P.endProcess = 'decays';
      twoBody([P.E, ...P.p], PDG[dec[0]].m, PDG[dec[1]].m, ctx.rng).forEach((p4, i) => addParticle(ctx, { pdg: dec[i], p4, x0: x, t0: P.t0, parent: P.id, origin: 'decay' }));
      return;
    }
    if (P.pdg === 3212) { // Σ0 -> Λ γ
      twoBody([P.E, ...P.p], PDG[3122].m, 0, ctx.rng).forEach((p4, i) => addParticle(ctx, { pdg: i ? 22 : 3122, p4, x0: P.x0, t0: P.t0, parent: P.id, origin: 'decay' }));
      P.end = P.x0; P.endProcess = 'decays'; return;
    }
    const d = ctx.rng.exp(LAMBDA_I), x = add(P.x0, mul(u, d));
    track(ctx, P, [[...P.x0, P.t0], [...x, P.t0 + d / OPT.C]], 'dashed'); P.end = x; P.endProcess = 'interacts';
    hadBlob(ctx, P, x, P.t0 + d / OPT.C, P.E - P.m);
  },
  nu(ctx, P) {
    const x = add(P.x0, mul(norm(P.p), 32));
    track(ctx, P, [[...P.x0, P.t0], [...x, P.t0 + 32 / OPT.C]], 'invisible'); P.end = x; P.endProcess = 'escapes (not detected)';
  },
};

function hadBlob(ctx, P, x, t, Tgev) {
  // illustrative hadronic interaction: quenched local deposits + evaporation neutrons [D]
  const T = Tgev * 1000; if (T <= 0) return;
  const fvis = 0.45 + 0.25 * ctx.rng(), nPts = 3 + Math.floor(ctx.rng() * 4);
  for (let i = 0; i < nPts; i++) {
    const off = mul(ctx.rng.iso(), 0.08 + 0.35 * ctx.rng());
    deposit(ctx, P, add(x, off), t, (T * fvis) / nPts / 0.7, 0.7, 0, null, 30);
  }
  const nN = Math.min(4, ctx.rng.poisson(0.3 + T / 250));
  for (let i = 0; i < nN; i++) {
    const Tn = (2 + ctx.rng.exp(10)) / 1000;
    addParticle(ctx, { pdg: 2112, p4: [PDG[2112].m + Tn, ...mul(ctx.rng.iso(), Math.sqrt(Tn * Tn + 2 * Tn * PDG[2112].m))], x0: x, t0: t, parent: P.id, origin: 'hadronic' });
  }
}

function neutronWalk(ctx, P, x0, t0, Tmev) {
  let x = x0.slice(), t = t0, T = Tmev, u = norm(P.p);
  const pts = [[...x, t]];
  for (let k = 0; k < 14 && T > 1; k++) {
    const lam = 0.04 + 0.55 * (1 - Math.exp(-T / 80));
    const d = ctx.rng.exp(lam); x = add(x, mul(u, d));
    const v = Math.sqrt(2 * T / 939.6) * OPT.C; t += d / v; pts.push([...x, t]);
    if (region(...x) === 'ROCK') { track(ctx, P, pts, 'walk'); P.end = x; P.endProcess = 'leaves detector'; return; }
    const f = T > 30 && ctx.rng() < 0.3 ? 0.6 + 0.3 * ctx.rng() : 0.2 + 0.6 * ctx.rng();
    const Tp = T * f;
    deposit(ctx, P, x, t, Tp, Math.max(0.2, Math.min(0.8, 0.85 - 0.6 * Math.exp(-Tp / 8))), 0, null, 30);
    T -= Tp; u = rotateAround(u, Math.acos(2 * ctx.rng() - 1), 2 * Math.PI * ctx.rng());
  }
  for (let k = 0; k < 4; k++) { x = add(x, mul(ctx.rng.iso(), 0.05 + 0.06 * ctx.rng())); pts.push([...x, t]); }
  const reg = region(...x);
  P.end = x;
  track(ctx, P, pts, 'walk');
  if (reg === 'ROCK' || reg === 'AIR') { P.endProcess = 'leaves detector'; return; }
  const dt = ctx.rng.exp(OPT.TAU_N || 202000);
  const onC = ctx.rng() < 0.01, Eg = onC ? 4.946 : 2.223;
  P.endProcess = `captured on ${onC ? '¹²C' : 'H'} after ${(dt / 1000).toFixed(0)} µs`;
  P.capture = { x, t: t + dt, Eg };
  const w = windowFor(ctx, t + dt, `n capture (${Eg} MeV γ)`);
  const g = addParticle(ctx, { pdg: 22, p4: [Eg / 1000, ...mul(ctx.rng.iso(), Eg / 1000)], x0: x, t0: t + dt, parent: P.id, origin: 'n capture' });
  g.win = w;
}

function shower(ctx, P, x0, u, E, t0, isGamma, isChild = false) {
  const tmax = Math.max(0.3, Math.log(E / EC) + (isGamma ? 0.5 : -0.5)), b = 0.5, a = b * tmax + 1;
  const dt = 0.25, nmax = Math.ceil((tmax + 10) / dt);
  const w = []; let ws = 0;
  for (let i = 0; i < nmax; i++) { const s = (i + 0.5) * dt; const v = Math.pow(s, a - 1) * Math.exp(-b * s); w.push(v); ws += v; }
  const pts = [[...x0, t0]];
  let cut = nmax, acc = 0;
  for (let i = 0; i < nmax; i++) { acc += w[i] / ws; if (acc > 0.97 && cut === nmax) cut = i + 1; }
  for (let i = 0; i < nmax; i++) {
    const s = (i + 0.5) * dt * X0, c = add(x0, mul(u, s)), sig = 0.03 + 0.012 * i * dt;
    const pos = add(c, [ctx.rng.gauss() * sig, ctx.rng.gauss() * sig, ctx.rng.gauss() * sig]);
    const dE = (E * 1000 * w[i]) / ws;
    deposit(ctx, P, pos, t0 + s / OPT.C, dE, 1, OPT.CHER_FRAC * dE, u, 16);
    if (region(...c) === 'LS') P.lenLS += dt * X0;
    if (i < cut) pts.push([...add(x0, mul(u, (i + 1) * dt * X0)), t0 + ((i + 1) * dt * X0) / OPT.C]);
  }
  const tr = track(ctx, P, isChild ? pts : pts, 'solid');
  tr.shower = true;
  const nb = Math.min(16, 4 + Math.floor(Math.log(E / 0.05 + 1) * 2.5));
  for (let i = 0; i < nb; i++) {
    const s = ctx.rng() * cut * dt * X0 * 0.85, st = add(x0, mul(u, s));
    const dir = rotateAround(u, Math.abs(ctx.rng.gauss()) * 0.35, 2 * Math.PI * ctx.rng());
    const l = 0.2 + ctx.rng() * 0.9, en = add(st, mul(dir, l));
    ctx.tracks.push({ pid: P.id, cls: ctx.rng() < 0.5 ? 'e' : 'gamma', pts: [[...st, t0 + s / OPT.C], [...en, t0 + (s + l) / OPT.C]], style: 'branch' });
  }
  if (!P.end) { P.end = add(x0, mul(u, cut * dt * X0)); P.endProcess = P.endProcess || 'showers'; }
}

function cosmicExtras(ctx, P) {
  // illustrative: spallation neutrons near the LS part of the track [D]
  const pts = P.track.pts.filter(p => region(p[0], p[1], p[2]) === 'LS');
  if (!pts.length) return;
  for (let i = 0; i < ctx.ev.spallation; i++) {
    const p = pts[Math.floor(ctx.rng() * pts.length)];
    const x = add([p[0], p[1], p[2]], mul(ctx.rng.iso(), Math.abs(ctx.rng.gauss()) * 0.8));
    const Tn = (1 + ctx.rng.exp(15)) / 1000;
    addParticle(ctx, { pdg: 2112, p4: [PDG[2112].m + Tn, ...mul(ctx.rng.iso(), Math.sqrt(Tn * Tn + 2 * Tn * PDG[2112].m))], x0: x, t0: p[3], parent: P.id, origin: 'spallation' });
  }
}

export function simulate(ev) {
  const rng = makeRng((ev.seed * 7919 + (ev.iev || 0) * 31 + 17) >>> 0);
  const ctx = makeCtx(ev, rng);
  for (const pr of ev.primaries) {
    const info = pinfo(pr.pdg);
    if (info.cls === 'nucleus' || info.cls === 'other') continue;
    const P = addParticle(ctx, { pdg: pr.pdg, p4: pr.p4, x0: ev.vertex, t0: ev.t0 || 0, origin: pr.origin || 'primary' });
    if (ev.cosmic) P.cosmicStart = true;
  }
  // assign deposits/tracks to readout windows by absolute time
  ctx.deps.forEach(d => { d.win = windowFor(ctx, d.t, 'Delayed'); });
  ctx.tracks.forEach(tr => { tr.win = windowFor(ctx, tr.pts[0][3], 'Delayed'); });
  ctx.wcd.forEach(s => { s.win = windowFor(ctx, s.t, 'Delayed'); });
  // sort windows by time, remap
  const order = ctx.windows.map((w, i) => i).sort((a, b) => ctx.windows[a].t0 - ctx.windows[b].t0);
  const remap = new Map(order.map((o, i) => [o, i]));
  ctx.windows = order.map(i => ctx.windows[i]);
  for (const arr of [ctx.deps, ctx.tracks, ctx.wcd]) arr.forEach(o => (o.win = remap.get(o.win)));
  ctx.windows.forEach((w, i) => {
    w.edep = 0; w.evis = 0;
    if (i) w.label = `${w.label}  Δt = ${w.t0 < 1e4 ? (w.t0 / 1000).toFixed(2) : (w.t0 / 1000).toFixed(0)} µs`;
  });
  ctx.deps.forEach(d => { ctx.windows[d.win].edep += d.edep; ctx.windows[d.win].evis += d.evis; });
  return { ...ev, particles: ctx.particles, tracks: ctx.tracks, deps: ctx.deps, wcd: ctx.wcd, tt: ctx.tt, windows: ctx.windows };
}

// ------------------------------------------------------------------ detector response
export function makeResponseModel(pm, wpm) {
  const { pos, dir, type, area, N, nL } = pm;
  let gd = 0, gux = 0, guy = 0, guz = 0; // scratch outputs of gf()
  const gf = (x, y, z, j) => {
    const dx = pos[3 * j] - x, dy = pos[3 * j + 1] - y, dz = pos[3 * j + 2] - z;
    const d2 = dx * dx + dy * dy + dz * dz, d = Math.sqrt(d2);
    gd = d; gux = dx / d; guy = dy / d; guz = dz / d;
    const ca = -(gux * dir[3 * j] + guy * dir[3 * j + 1] + guz * dir[3 * j + 2]);
    if (ca <= 0) return 0;
    const b = x * gux + y * guy + z * guz, c = x * x + y * y + z * z - GEO.R_LS * GEO.R_LS;
    const disc = b * b - c; let tl = c < 0 && disc > 0 ? -b + Math.sqrt(disc) : 0; if (tl > d) tl = d;
    return (area[j] * (0.35 + 0.65 * ca) * Math.exp(-tl / OPT.ATT_LS - (d - tl) / OPT.ATT_W)) / (4 * Math.PI * d2);
  };
  let s0 = 0; for (let j = 0; j < nL; j++) s0 += gf(0, 0, 0, j);
  const scale = OPT.LY_CENTRE / s0;
  const coarse = []; for (let j = 0; j < nL; j += 25) coarse.push(j);
  const coarseSum = (x, y, z) => { let s = 0; for (const j of coarse) s += gf(x, y, z, j); return s * 25; };
  // charge-centroid bias factor for the simple vertex estimate
  const zc = 8; let cw = 0, cz = 0;
  for (const j of coarse) { const g = gf(0, 0, zc, j); cw += g; cz += g * pos[3 * j + 2]; }
  const kCentroid = zc / (cz / cw);
  const SPMT_EFF = 0.75; // [D]
  const sigC = (s) => (s * Math.PI) / 180;

  function mergeDeps(list) {
    // combine nearby deposits (same window) into voxels to keep the PMT loop fast
    let cell = 0.35, out;
    for (;;) {
      const m = new Map();
      for (const d of list) {
        const k = `${Math.round(d.x / cell)},${Math.round(d.y / cell)},${Math.round(d.z / cell)},${d.cher > 0 ? d.sig : 0}`;
        let o = m.get(k);
        if (!o) { o = { x: 0, y: 0, z: 0, t: d.t, evis: 0, cher: 0, u: [0, 0, 0], sig: d.sig, w: 0 }; m.set(k, o); }
        const wt = d.evis + d.cher + 1e-9;
        o.x += d.x * wt; o.y += d.y * wt; o.z += d.z * wt; o.w += wt; o.evis += d.evis; o.cher += d.cher;
        if (d.t < o.t) o.t = d.t;
        if (d.cher > 0 && d.u) { o.u[0] += d.u[0] * d.cher; o.u[1] += d.u[1] * d.cher; o.u[2] += d.u[2] * d.cher; }
      }
      out = [...m.values()].map(o => ({ x: o.x / o.w, y: o.y / o.w, z: o.z / o.w, t: o.t, evis: o.evis, cher: o.cher, u: o.cher > 0 ? norm(o.u) : null, sig: o.sig }));
      if (out.length <= 110 || cell > 3) return out;
      cell *= 1.6;
    }
  }

  function windowResponse(sim, w, rng, opts = {}) {
    const deps = mergeDeps(sim.deps.filter(d => d.win === w));
    const t0 = sim.windows[w].t0;
    const mu = new Float64Array(N), tarr = new Float64Array(N).fill(Infinity);
    for (const dp of deps) {
      const lx = dp.x, ly = dp.y, lz = dp.z, tt = dp.t - t0;
      const hasC = dp.cher > 0 && dp.u;
      const sg = sigC(dp.sig || 8), norm0 = 1 / (2 * Math.PI * Math.sin(THETA_C) * sg * Math.sqrt(2 * Math.PI));
      const cLo = Math.cos(Math.min(Math.PI, THETA_C + 4 * sg)), cHi = Math.cos(Math.max(0, THETA_C - 4 * sg));
      for (let j = 0; j < N; j++) {
        const g = gf(lx, ly, lz, j);
        if (g <= 0) continue;
        const d = gd, ux = gux, uy = guy, uz = guz;
        let m = dp.evis * g;
        if (hasC) {
          const cp = ux * dp.u[0] + uy * dp.u[1] + uz * dp.u[2];
          if (cp > cLo && cp < cHi) { const psi = Math.acos(cp); m += dp.cher * g * 4 * Math.PI * norm0 * Math.exp(-((psi - THETA_C) ** 2) / (2 * sg * sg)); }
        }
        m *= scale * (type[j] === 2 ? SPMT_EFF : 1);
        mu[j] += m;
        if (m > 0.002) { const ta = tt + (d * OPT.N_EFF) / OPT.C; if (ta < tarr[j]) tarr[j] = ta; }
      }
    }
    const idx = [], npe = [], tim = [], dark = [];
    let npeL = 0, npeS = 0, firedL = 0, muL = 0;
    const tFirst = deps.length ? Math.min(...deps.map(d => d.t - t0)) : 0;
    for (let j = 0; j < N; j++) {
      if (j < nL) muL += mu[j];
      let k = rng.poisson(mu[j]); let isDark = false, t;
      if (k > 0) {
        const tt = Number.isFinite(tarr[j]) ? tarr[j] : tFirst + 60;
        t = tt + rng.exp(OPT.TAU_FAST) / Math.max(1, Math.sqrt(k)) + rng.gauss() * PMT.TTS[type[j]];
      } else if (opts.dark && rng() < PMT.DCR[type[j]] * PMT.GATE_NS * 1e-9) { k = 1; isDark = true; t = -200 + rng() * PMT.GATE_NS; }
      if (!k) continue;
      idx.push(j); npe.push(k); tim.push(t); dark.push(isDark ? 1 : 0);
      if (j < nL) { npeL += k; firedL++; } else npeS += k;
    }
    // multiplicity trigger
    const tl = []; for (let i = 0; i < idx.length; i++) if (idx[i] < nL) tl.push(tim[i]);
    tl.sort((a, b) => a - b); let maxIn = 0;
    for (let i = 0, k0 = 0; i < tl.length; i++) { while (tl[i] - tl[k0] > PMT.TRIG_WIN) k0++; maxIn = Math.max(maxIn, i - k0 + 1); }
    const res = { win: w, idx: Int32Array.from(idx), npe: Float32Array.from(npe), t: Float32Array.from(tim), dark: Uint8Array.from(dark), npeL, npeS, firedL, muL, maxIn, triggered: maxIn >= PMT.TRIG_N };
    res.wcd = wcdResponse(sim, w, rng);
    res.reco = reconstruct(res);
    res.evis = sim.windows[w].evis; res.edep = sim.windows[w].edep;
    return res;
  }

  function wcdResponse(sim, w, rng) {
    const segs = sim.wcd.filter(s => s.win === w);
    const n = wpm.N, mu = new Float64Array(n), out = { idx: [], npe: [], fired: 0 };
    if (!segs.length) return out;
    const thw = Math.acos(1 / 1.33), sg = sigC(14), nrm = 1 / (2 * Math.PI * Math.sin(thw) * sg * Math.sqrt(2 * Math.PI));
    for (const s of segs) {
      for (let j = 0; j < n; j++) {
        const dx = wpm.pos[3 * j] - s.x[0], dy = wpm.pos[3 * j + 1] - s.x[1], dz = wpm.pos[3 * j + 2] - s.x[2];
        const d = Math.hypot(dx, dy, dz), ca = -(dx * wpm.dir[3 * j] + dy * wpm.dir[3 * j + 1] + dz * wpm.dir[3 * j + 2]) / d;
        if (ca <= 0) continue;
        const psi = Math.acos(Math.max(-1, Math.min(1, (dx * s.u[0] + dy * s.u[1] + dz * s.u[2]) / d)));
        const frac = (0.2 * (0.4 + 0.6 * ca)) / (d * d) * (4 * Math.PI * nrm * Math.exp(-((psi - thw) ** 2) / (2 * sg * sg)) + 0.6); // + Tyvek diffuse [D]
        mu[j] += 900 * s.L * frac * Math.exp(-d / OPT.ATT_W);
      }
    }
    for (let j = 0; j < n; j++) { const k = rng.poisson(mu[j]); if (k) { out.idx.push(j); out.npe.push(k); } }
    out.fired = out.idx.length; out.tagged = out.fired >= 40;
    return out;
  }

  function reconstruct(r) {
    // approximate: charge centroid vertex with bias factor, energy from NPE with position correction
    let w = 0, cx = 0, cy = 0, cz = 0;
    for (let i = 0; i < r.idx.length; i++) { const j = r.idx[i]; if (j >= nL || r.dark[i]) continue; const q = r.npe[i]; w += q; cx += q * pos[3 * j]; cy += q * pos[3 * j + 1]; cz += q * pos[3 * j + 2]; }
    if (!w) return null;
    let v = [(cx / w) * kCentroid, (cy / w) * kCentroid, (cz / w) * kCentroid];
    if (len(v) > GEO.R_LS) v = mul(norm(v), GEO.R_LS - 0.1);
    const G = coarseSum(...v) / coarseSum(0, 0, 0);
    const E = w / (OPT.LY_CENTRE * Math.max(0.3, G)); // dark-noise hits excluded
    const reco = { vertex: v, E };
    // track-like: two separated high-charge clusters (entry/exit)
    if (r.npeL > 4e5) {
      const L = []; for (let i = 0; i < r.idx.length; i++) if (r.idx[i] < nL) L.push(i);
      const P3 = (i) => [pos[3 * r.idx[i]], pos[3 * r.idx[i] + 1], pos[3 * r.idx[i] + 2]];
      const top = L.reduce((a, b) => (r.npe[b] > r.npe[a] ? b : a), L[0]);
      const clus = (c) => { let s = 0, p = [0, 0, 0], tt = 0; for (const i of L) { const x = P3(i); if (len(sub(x, c)) < 3.5) { s += r.npe[i]; p = add(p, mul(x, r.npe[i])); tt += r.t[i] * r.npe[i]; } } return { q: s, p: mul(p, 1 / s), t: tt / s }; };
      const A = clus(P3(top));
      let best = -1; for (const i of L) if (len(sub(P3(i), A.p)) > 12 && (best < 0 || r.npe[i] > r.npe[best])) best = i;
      if (best >= 0) {
        const B = clus(P3(best));
        if (B.q > 0.03 * A.q) { const [en, ex] = A.t <= B.t ? [A, B] : [B, A]; reco.track = { entry: en.p, exit: ex.p, dir: norm(sub(ex.p, en.p)) }; }
      }
    }
    return reco;
  }

  function quickNpe(sim) { // expected LPMT p.e. in the prompt window, coarse estimate for batch export
    let s = 0; for (const d of sim.deps) if (d.win === 0) s += d.evis * coarseSum(d.x, d.y, d.z);
    return s * scale;
  }

  return { windowResponse, quickNpe, scale, kCentroid };
}

export function summarize(sim) {
  const prim = sim.particles.filter(p => p.parent < 0);
  const count = {}; for (const p of prim) count[p.cls] = (count[p.cls] || 0) + 1;
  const edepLS = sim.deps.reduce((s, d) => s + d.edep, 0), evisLS = sim.deps.reduce((s, d) => s + d.evis, 0);
  const R = len(sim.vertex);
  return {
    count, edepLS, evisLS, R, inFV: R < GEO.FV_R && Math.abs(sim.vertex[2]) < GEO.FV_Z,
    missing: sim.nu ? sim.nu.E * 1000 - edepLS : null,
    ttHits: sim.tt.length,
  };
}
