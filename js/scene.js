// Three.js detector model and event visuals. No DOM access here.
import * as THREE from 'three';
import { Line2 } from 'three/addons/lines/Line2.js';
import { LineMaterial } from 'three/addons/lines/LineMaterial.js';
import { LineGeometry } from 'three/addons/lines/LineGeometry.js';
import { GEO, PMT, ttWalls } from './detector.js';

export const COLORS = {
  mu: 0xff4d5e, e: 0x3ddc84, gamma: 0xffd34d, chad: 0x4da3ff, p: 0xb07cff, n: 0x9aa7b4, nu: 0xe6f0fa, photon: 0xbff3ff,
};
export const CLS_KEY = { mu: 'mu', e: 'e', gamma: 'gamma', pi0: 'gamma', chad: 'chad', p: 'p', n: 'n', nhad: 'n', nu: 'nu' };

// ------------------------------------------------------------------ detector
export function buildDetector(pm, wpm) {
  const root = new THREE.Group();
  const groups = {};
  const mats = []; // {m, base, kind}
  const reg = (m, base, kind = 'shell') => { m.transparent = true; m.opacity = base; mats.push({ m, base, kind }); return m; };
  const G = (name, parent = root) => { const g = new THREE.Group(); g.name = name; parent.add(g); groups[name] = g; return g; };
  const lineMat = (color, op, kind = 'line') => reg(new THREE.LineBasicMaterial({ color, depthWrite: false }), op, kind);
  const segs = (arr, mat) => { const g = new THREE.BufferGeometry(); g.setAttribute('position', new THREE.Float32BufferAttribute(arr, 3)); return new THREE.LineSegments(g, mat); };
  const ring = (out, r, z, n = 96, cx = 0, cy = 0) => { for (let i = 0; i < n; i++) { const a = (2 * Math.PI * i) / n, b = (2 * Math.PI * (i + 1)) / n; out.push(cx + r * Math.cos(a), cy + r * Math.sin(a), z, cx + r * Math.cos(b), cy + r * Math.sin(b), z); } };
  const sph = (R, th, ph) => [R * Math.sin(th) * Math.cos(ph), R * Math.sin(th) * Math.sin(ph), R * Math.cos(th)];

  // groups used by the exploded view
  const gOuter = G('outer'), gSss = G('sssAll'), gInner = G('inner'), gTop = G('top');

  // ---- liquid scintillator + acrylic vessel
  const ls = G('ls', gInner);
  ls.add(new THREE.Mesh(new THREE.SphereGeometry(GEO.R_LS, 72, 48), reg(new THREE.MeshBasicMaterial({ color: 0x2b62e8, depthWrite: false, side: THREE.BackSide, blending: THREE.AdditiveBlending }), 0.11, 'glow')));
  ls.add(new THREE.Mesh(new THREE.SphereGeometry(GEO.R_LS * 0.995, 72, 48), reg(new THREE.MeshBasicMaterial({ color: 0x163f9a, depthWrite: false, side: THREE.FrontSide, blending: THREE.AdditiveBlending }), 0.05, 'glow')));
  const av = G('av', gInner);
  const Rav = GEO.R_LS + GEO.AV_THICK;
  av.add(new THREE.Mesh(new THREE.SphereGeometry(Rav, 72, 48), reg(new THREE.MeshPhongMaterial({ color: 0x9fe7ff, specular: 0x88ccee, shininess: 80, depthWrite: false, side: THREE.DoubleSide }), 0.09, 'shell')));
  { // panel seams: ~263 panels
    const a = []; const rows = 19;
    for (let i = 1; i < rows; i++) { const th = (Math.PI * i) / rows; ring(a, Rav * Math.sin(th), Rav * Math.cos(th), 96); }
    for (let i = 0; i < rows; i++) {
      const th0 = (Math.PI * i) / rows, th1 = (Math.PI * (i + 1)) / rows;
      const n = Math.max(4, Math.round(26 * Math.sin((th0 + th1) / 2)));
      for (let k = 0; k < n; k++) { const ph = (2 * Math.PI * (k + 0.5 * (i % 2))) / n; a.push(...sph(Rav, th0, ph), ...sph(Rav, th1, ph)); }
    }
    av.add(segs(a, lineMat(0x8fdcff, 0.22)));
  }

  // ---- connecting bars (590) and nodes
  const bars = G('bars', gSss);
  {
    const a = [], n = GEO.N_AV_BARS, ga = Math.PI * (3 - Math.sqrt(5));
    const node = new THREE.InstancedMesh(new THREE.BoxGeometry(0.35, 0.35, 0.08), reg(new THREE.MeshPhongMaterial({ color: 0x9aa9b8, depthWrite: true }), 0.7, 'solidish'), n);
    const m4 = new THREE.Matrix4(), q = new THREE.Quaternion(), zAx = new THREE.Vector3(0, 0, 1);
    for (let i = 0; i < n; i++) {
      const z = 1 - (2 * (i + 0.5)) / n, r = Math.sqrt(1 - z * z), p = i * ga + 0.4;
      const v = new THREE.Vector3(r * Math.cos(p), r * Math.sin(p), z);
      a.push(...v.clone().multiplyScalar(Rav).toArray(), ...v.clone().multiplyScalar(GEO.R_SSS).toArray());
      q.setFromUnitVectors(zAx, v); m4.compose(v.clone().multiplyScalar(Rav + 0.05), q, new THREE.Vector3(1, 1, 1)); node.setMatrixAt(i, m4);
    }
    bars.add(segs(a, lineMat(0xb8c4d0, 0.35)));
    bars.add(node);
  }

  // ---- stainless steel latticed shell + legs
  const sss = G('sss', gSss);
  {
    const R = GEO.R_SSS, a = [], nTh = 26, nPh = 56;
    for (let i = 1; i < nTh; i++) { const th = (Math.PI * i) / nTh; ring(a, R * Math.sin(th), R * Math.cos(th), nPh * 2); }
    for (let i = 0; i < nTh; i++) for (let j = 0; j < nPh; j++) {
      const th0 = (Math.PI * i) / nTh, th1 = (Math.PI * (i + 1)) / nTh, p0 = (2 * Math.PI * j) / nPh, p1 = (2 * Math.PI * (j + 1)) / nPh;
      a.push(...sph(R, th0, p0), ...sph(R, th1, p0));
      a.push(...sph(R, th0, i % 2 ? p0 : p1), ...sph(R, th1, i % 2 ? p1 : p0));
    }
    sss.add(segs(a, lineMat(0x8d99a6, 0.3)));
    const L = [], nLeg = 30, th = (104 * Math.PI) / 180, rL = R * Math.sin(th), zTop = R * Math.cos(th);
    for (let k = 0; k < nLeg; k++) {
      const p = (2 * Math.PI * k) / nLeg, p2 = (2 * Math.PI * (k + 1)) / nLeg;
      const x = rL * Math.cos(p), y = rL * Math.sin(p), x2 = rL * Math.cos(p2), y2 = rL * Math.sin(p2);
      L.push(x, y, zTop, x, y, GEO.WP_ZMIN);
      for (const [za, zb] of [[zTop, -10.5], [-10.5, -16], [-16, GEO.WP_ZMIN]]) { L.push(x, y, za, x2, y2, zb, x2, y2, za, x, y, zb); }
    }
    for (const z of [-10.5, -16, GEO.WP_ZMIN + 0.05]) ring(L, rL, z, nLeg);
    sss.add(segs(L, lineMat(0x7d8894, 0.32)));
  }

  // ---- Earth magnetic field compensation coils
  const coils = G('coils', gSss);
  { const a = []; for (let i = 0; i < 24; i++) { const th = (Math.PI * (i + 1.5)) / 27; ring(a, GEO.R_COIL * Math.sin(th), GEO.R_COIL * Math.cos(th), 128); } coils.add(segs(a, lineMat(0xc4573f, 0.42))); }

  // ---- 20-inch and 3-inch PMTs (instanced; faces toward the centre)
  const q = new THREE.Quaternion(), m4 = new THREE.Matrix4(), one = new THREE.Vector3(1, 1, 1), yAx = new THREE.Vector3(0, 1, 0), zAx = new THREE.Vector3(0, 0, 1);
  const lpmtGeo = new THREE.SphereGeometry(PMT.R_LPMT, 10, 5, 0, Math.PI * 2, 0, Math.PI / 2);
  const lpmtMat = reg(new THREE.MeshPhongMaterial({ color: 0xa47a32, specular: 0xffe0a0, shininess: 60, emissive: 0x1a1206 }), 1, 'pmt');
  const lpmt = new THREE.InstancedMesh(lpmtGeo, lpmtMat, pm.nL);
  const spmt = new THREE.InstancedMesh(new THREE.SphereGeometry(PMT.R_SPMT * 1.2, 6, 3, 0, Math.PI * 2, 0, Math.PI / 2), reg(new THREE.MeshPhongMaterial({ color: 0xc8a050, shininess: 30 }), 1, 'pmt'), pm.nS);
  const dyn = new THREE.Color(0xb48a3e), mcp = new THREE.Color(0x8f6a2a);
  const v = new THREE.Vector3();
  for (let j = 0; j < pm.N; j++) {
    v.set(pm.pos[3 * j], pm.pos[3 * j + 1], pm.pos[3 * j + 2]);
    q.setFromUnitVectors(yAx, new THREE.Vector3(pm.dir[3 * j], pm.dir[3 * j + 1], pm.dir[3 * j + 2]));
    m4.compose(v, q, one);
    if (j < pm.nL) { lpmt.setMatrixAt(j, m4); lpmt.setColorAt(j, pm.type[j] === 0 ? dyn : mcp); } else spmt.setMatrixAt(j - pm.nL, m4);
  }
  const gL = G('lpmt', gSss); gL.add(lpmt);
  const gS = G('spmt', gSss); gS.add(spmt);
  // hit glows: one instanced disc per PMT that fired, in time order
  const glowMat = reg(new THREE.MeshBasicMaterial({ color: 0xffffff, depthWrite: false, toneMapped: false, transparent: true, opacity: 0.92, side: THREE.DoubleSide }), 0.92, 'event');
  const hitGlow = new THREE.InstancedMesh(new THREE.CircleGeometry(1, 14), glowMat, pm.N);
  hitGlow.count = 0; hitGlow.frustumCulled = false;
  hitGlow.instanceColor = new THREE.InstancedBufferAttribute(new Float32Array(pm.N * 3), 3);
  gSss.add(hitGlow);
  // precomputed glow transforms
  const glowMatrix = (j) => {
    const R = j < pm.nL ? PMT.R_LPMT * 1.25 : PMT.R_SPMT * 3;
    v.set(pm.pos[3 * j], pm.pos[3 * j + 1], pm.pos[3 * j + 2]).multiplyScalar(0.985);
    q.setFromUnitVectors(zAx, new THREE.Vector3(pm.dir[3 * j], pm.dir[3 * j + 1], pm.dir[3 * j + 2]));
    return m4.clone().compose(v.clone(), q.clone(), new THREE.Vector3(R, R, R));
  };

  // ---- water pool (Water Cherenkov Detector) and its PMTs
  const pool = G('pool', gOuter);
  const H = GEO.WP_ZMAX - GEO.WP_ZMIN, zc = (GEO.WP_ZMAX + GEO.WP_ZMIN) / 2;
  const wall = new THREE.Mesh(new THREE.CylinderGeometry(GEO.WP_R, GEO.WP_R, H, 96, 1, true), reg(new THREE.MeshPhongMaterial({ color: 0x0c3d6e, side: THREE.DoubleSide, depthWrite: false, shininess: 10 }), 0.13, 'shell'));
  wall.rotation.x = Math.PI / 2; wall.position.z = zc; pool.add(wall);
  const floor = new THREE.Mesh(new THREE.CircleGeometry(GEO.WP_R, 96), reg(new THREE.MeshPhongMaterial({ color: 0x0a2c4c, depthWrite: false, side: THREE.DoubleSide }), 0.22, 'shell'));
  floor.position.z = GEO.WP_ZMIN; pool.add(floor);
  const water = new THREE.Mesh(new THREE.CircleGeometry(GEO.WP_R, 96), reg(new THREE.MeshPhongMaterial({ color: 0x1e5c96, depthWrite: false, side: THREE.DoubleSide }), 0.08, 'shell'));
  water.position.z = GEO.WATER_TOP; pool.add(water);
  { const a = []; ring(a, GEO.WP_R, GEO.WP_ZMIN, 128); ring(a, GEO.WP_R, GEO.WP_ZMAX, 128); ring(a, GEO.WP_R, GEO.WATER_TOP, 128);
    for (let i = 0; i < 48; i++) { const p = (2 * Math.PI * i) / 48; a.push(GEO.WP_R * Math.cos(p), GEO.WP_R * Math.sin(p), GEO.WP_ZMIN, GEO.WP_R * Math.cos(p), GEO.WP_R * Math.sin(p), GEO.WP_ZMAX); }
    pool.add(segs(a, lineMat(0xd8e6f2, 0.16))); } // Tyvek lining seams
  const wGeo = new THREE.SphereGeometry(0.22, 8, 4, 0, Math.PI * 2, 0, Math.PI / 2);
  const wcd = new THREE.InstancedMesh(wGeo, reg(new THREE.MeshPhongMaterial({ color: 0x7d6a3e, shininess: 40 }), 1, 'pmt'), wpm.N);
  for (let j = 0; j < wpm.N; j++) {
    v.set(wpm.pos[3 * j], wpm.pos[3 * j + 1], wpm.pos[3 * j + 2]);
    q.setFromUnitVectors(yAx, new THREE.Vector3(wpm.dir[3 * j], wpm.dir[3 * j + 1], wpm.dir[3 * j + 2]));
    const s = wpm.kind[j] === 2 ? 0.45 : 1; m4.compose(v, q, new THREE.Vector3(s, s, s)); wcd.setMatrixAt(j, m4);
  }
  const gW = G('wcdpmt', gOuter); gW.add(wcd);
  const wcdGlow = new THREE.InstancedMesh(new THREE.CircleGeometry(0.4, 12), glowMat, wpm.N);
  wcdGlow.count = 0; wcdGlow.frustumCulled = false; wcdGlow.instanceColor = new THREE.InstancedBufferAttribute(new Float32Array(wpm.N * 3), 3);
  gW.add(wcdGlow);
  const wcdGlowMatrix = (j) => {
    v.set(wpm.pos[3 * j], wpm.pos[3 * j + 1], wpm.pos[3 * j + 2]).addScaledVector(new THREE.Vector3(wpm.dir[3 * j], wpm.dir[3 * j + 1], wpm.dir[3 * j + 2]), 0.25);
    q.setFromUnitVectors(zAx, new THREE.Vector3(wpm.dir[3 * j], wpm.dir[3 * j + 1], wpm.dir[3 * j + 2]));
    return m4.clone().compose(v.clone(), q.clone(), one);
  };

  // ---- chimney, cover, calibration house
  const top = G('chimney', gTop);
  const ch = new THREE.Mesh(new THREE.CylinderGeometry(GEO.CHIMNEY_R, GEO.CHIMNEY_R, 27.2 - Rav, 20, 1, true), reg(new THREE.MeshPhongMaterial({ color: 0x9fb4c6, side: THREE.DoubleSide, depthWrite: false }), 0.4, 'solidish'));
  ch.rotation.x = Math.PI / 2; ch.position.z = (27.2 + Rav) / 2; top.add(ch);
  const house = G('house', gTop);
  const hb = new THREE.Mesh(new THREE.BoxGeometry(6, 6, 2.6), reg(new THREE.MeshPhongMaterial({ color: 0x8e9eb0, depthWrite: false }), 0.25, 'solidish'));
  hb.position.z = 28.5; house.add(hb);
  { const e = new THREE.LineSegments(new THREE.EdgesGeometry(hb.geometry), lineMat(0xc8d6e4, 0.5)); e.position.copy(hb.position); house.add(e); }
  const cover = G('cover', gOuter);
  const cv = new THREE.Mesh(new THREE.RingGeometry(1.2, GEO.WP_R, 96), reg(new THREE.MeshPhongMaterial({ color: 0x141c24, side: THREE.DoubleSide, depthWrite: false }), 0.35, 'shell'));
  cv.position.z = GEO.WP_ZMAX; cover.add(cv);

  // ---- Top Tracker: 3 layers x 7x3 walls of plastic scintillator strips
  const tt = G('tt', gTop);
  {
    const walls = ttWalls(), W = GEO.TT_WALL;
    const box = new THREE.InstancedMesh(new THREE.BoxGeometry(W * 0.985, W * 0.985, 0.06), reg(new THREE.MeshPhongMaterial({ color: 0x6f9a7e, depthWrite: false, side: THREE.DoubleSide }), 0.16, 'shell'), walls.length);
    const a = [];
    walls.forEach((w, i) => {
      m4.compose(new THREE.Vector3(w.x, w.y, w.z), new THREE.Quaternion(), one); box.setMatrixAt(i, m4);
      const h = (W * 0.985) / 2;
      a.push(w.x - h, w.y - h, w.z, w.x + h, w.y - h, w.z, w.x + h, w.y - h, w.z, w.x + h, w.y + h, w.z, w.x + h, w.y + h, w.z, w.x - h, w.y + h, w.z, w.x - h, w.y + h, w.z, w.x - h, w.y - h, w.z);
      for (let k = 1; k < 8; k++) { const o = -h + (2 * h * k) / 8; if (w.layer % 2) a.push(w.x + o, w.y - h, w.z, w.x + o, w.y + h, w.z); else a.push(w.x - h, w.y + o, w.z, w.x + h, w.y + o, w.z); }
    });
    tt.add(box); tt.add(segs(a, lineMat(0x9fd0b0, 0.35)));
    const b = []; const X = 25.5, Y = 11.2, z0 = GEO.WP_ZMAX + 0.4;
    for (const y of [-Y, Y]) for (const z of [z0, GEO.TT_Z[2] + 0.6]) b.push(-X, y, z, X, y, z);
    for (let x = -X; x <= X + 0.1; x += 6.375) for (const y of [-Y, Y]) b.push(x, y, z0, x, y, GEO.TT_Z[2] + 0.6);
    tt.add(segs(b, lineMat(0xe0b040, 0.45))); // bridge frame
  }

  // ---- cavern outline
  const cavern = G('cavern', gOuter);
  { const a = []; const Rc = 27; ring(a, Rc, GEO.WP_ZMIN - 0.5, 128); ring(a, Rc, 31, 128);
    for (let i = 0; i < 32; i++) { const p = (2 * Math.PI * i) / 32; a.push(Rc * Math.cos(p), Rc * Math.sin(p), GEO.WP_ZMIN - 0.5, Rc * Math.cos(p), Rc * Math.sin(p), 31); }
    for (let i = 1; i < 10; i++) { const th = (Math.PI / 2) * (i / 10); ring(a, Rc * Math.sin(th), 31 + 14 * Math.cos(th), 96); }
    for (let i = 0; i < 24; i++) { const p = (2 * Math.PI * i) / 24; for (let k = 0; k < 12; k++) { const t0 = (Math.PI / 2) * (k / 12), t1 = (Math.PI / 2) * ((k + 1) / 12); a.push(Rc * Math.sin(t0) * Math.cos(p), Rc * Math.sin(t0) * Math.sin(p), 31 + 14 * Math.cos(t0), Rc * Math.sin(t1) * Math.cos(p), Rc * Math.sin(t1) * Math.sin(p), 31 + 14 * Math.cos(t1)); } }
    cavern.add(segs(a, lineMat(0x5a6b7c, 0.2))); }
  cavern.visible = false;

  // ---- calibration systems (simplified)
  const calib = G('calib', gInner);
  { const a = [0, 0, GEO.R_LS, 0, 0, -GEO.R_LS]; // ACU along the central axis
    for (const phd of [123.4, 267.4]) { const p = (phd * Math.PI) / 180; for (let k = 0; k < 48; k++) { const t0 = (Math.PI * k) / 48, t1 = (Math.PI * (k + 1)) / 48; a.push(...sph(Rav + 0.15, t0, p), ...sph(Rav + 0.15, t1, p)); } } // Guide Tube
    const anchor = [0, -12.4, 12.4], src = [0, -5.5, 1.5]; a.push(0, 0, GEO.R_LS, ...src, ...anchor, ...src); // CLS cables
    calib.add(segs(a, lineMat(0xffb547, 0.55)));
    const s = new THREE.Mesh(new THREE.SphereGeometry(0.18, 12, 8), reg(new THREE.MeshBasicMaterial({ color: 0xffb547 }), 0.9, 'solidish')); s.position.set(...src); calib.add(s);
    const rov = new THREE.Mesh(new THREE.CylinderGeometry(0.35, 0.35, 0.9, 16), reg(new THREE.MeshPhongMaterial({ color: 0xffb547 }), 0.8, 'solidish')); rov.position.set(6, 7, -5); calib.add(rov); }
  calib.visible = false;

  const L = (-131.5 * Math.PI) / 180, Rt = (48.5 * Math.PI) / 180; // silhouette azimuths for the default view
  const labels = [
    { text: 'Acrylic vessel · R 17.7 m', g: 'av', p: sph(Rav, 0.62, L) },
    { text: 'Liquid scintillator · 20 kt', g: 'ls', p: sph(12, 2.35, L) },
    { text: 'SS structure · Ø 40.1 m', g: 'sss', p: sph(GEO.R_SSS, 1.95, Rt) },
    { text: '20" + 3" PMTs', g: 'lpmt', p: sph(GEO.R_PMT, 1.15, Rt) },
    { text: 'Compensation coils', g: 'coils', p: sph(GEO.R_COIL, 2.45, L) },
    { text: '590 support bars', g: 'bars', p: sph(GEO.R_SSS - 1, 2.25, Rt) },
    { text: 'Water pool · Ø 43.5 m', g: 'pool', p: [GEO.WP_R * Math.cos(L), GEO.WP_R * Math.sin(L), -15] },
    { text: 'Top Tracker · 3 layers', g: 'tt', p: [-22, 9, GEO.TT_Z[2] + 0.5] },
    { text: 'Chimney', g: 'chimney', p: [0, 0, 21] },
    { text: 'Calibration house', g: 'house', p: [0, 0, 30.2] },
  ];

  return { root, groups, mats, lpmt, spmt, wcd, hitGlow, wcdGlow, glowMatrix, wcdGlowMatrix, labels };
}

// geometry-mode presentation: opacity per material kind and clipping wedge
export const MODES = ['cutaway', 'transparent', 'solid', 'exploded'];
const MODE_F = {
  cutaway: { shell: 1, line: 1, pmt: 1, solidish: 1, glow: 1 },
  transparent: { shell: 0.7, line: 0.6, pmt: 0.22, solidish: 0.5, glow: 1 },
  solid: { shell: 3.6, line: 1.5, pmt: 1, solidish: 2, glow: 0.7 },
  exploded: { shell: 1, line: 1, pmt: 1, solidish: 1, glow: 1 },
};
export function applyMode(det, mode, layerOpacity, planes) {
  const f = MODE_F[mode];
  const cut = mode !== 'transparent';
  for (const { m, base, kind } of det.mats) {
    if (kind === 'event') { m.clippingPlanes = cut ? planes : null; m.clipIntersection = true; m.needsUpdate = true; continue; }
    const op = Math.min(1, base * (f[kind] ?? 1) * (kind === 'pmt' ? 1 : layerOpacity * 2));
    m.opacity = op; m.transparent = op < 0.999 || kind === 'glow' || kind === 'line';
    m.depthWrite = kind === 'pmt' ? op > 0.5 : kind === 'solidish' && op > 0.6;
    m.clippingPlanes = cut ? planes : null; m.clipIntersection = true; m.needsUpdate = true;
  }
  const ex = mode === 'exploded';
  det.groups.outer.scale.setScalar(ex ? 1.3 : 1);
  det.groups.sssAll.scale.setScalar(ex ? 1.16 : 1);
  det.groups.top.position.z = ex ? 12 : 0;
}
export function wedgePlanes(phi0Deg, angleDeg) {
  const a1 = ((phi0Deg - angleDeg / 2) * Math.PI) / 180, a2 = ((phi0Deg + angleDeg / 2) * Math.PI) / 180;
  // clipped where angle in (a1, a2): n1.p > 0 and n2.p > 0; three clips where distance < 0, so negate
  const p1 = new THREE.Plane(new THREE.Vector3(Math.sin(a1), -Math.cos(a1), 0), 0);
  const p2 = new THREE.Plane(new THREE.Vector3(-Math.sin(a2), Math.cos(a2), 0), 0);
  return [p1, p2];
}

// ------------------------------------------------------------------ event visuals
const WIDTH = { mu: 3.4, e: 2.6, gamma: 1.8, chad: 2.6, p: 2.6, n: 1.4, nu: 1.4 };
export function buildEventVisuals(sim, resolution) {
  const group = new THREE.Group();
  const items = []; // track visuals {line, mat, tr, times, key}
  for (const tr of sim.tracks) {
    if (tr.pts.length < 2) continue;
    const key = CLS_KEY[tr.cls] || 'n';
    const pos = []; tr.pts.forEach(p => pos.push(p[0], p[1], p[2]));
    const geo = new LineGeometry(); geo.setPositions(pos);
    const base = new THREE.Color(COLORS[key]);
    const branch = tr.style === 'branch';
    const boost = branch ? 1.3 : key === 'nu' ? 0.9 : 2.2;
    const mat = new LineMaterial({
      color: base.clone().multiplyScalar(boost), linewidth: branch ? 1.1 : WIDTH[key] || 2,
      dashed: tr.style === 'dashed' || tr.style === 'invisible' || tr.style === 'walk',
      dashSize: tr.style === 'walk' ? 0.12 : 0.45, gapSize: tr.style === 'walk' ? 0.12 : 0.3,
      transparent: true, opacity: branch ? 0.55 : tr.style === 'invisible' ? 0.45 : 1, depthWrite: false, toneMapped: false,
    });
    mat.resolution.copy(resolution);
    const line = new Line2(geo, mat);
    if (mat.dashed) line.computeLineDistances();
    line.frustumCulled = false; line.userData = { pid: tr.pid };
    geo.instanceCount = 0;
    group.add(line);
    items.push({ line, mat, tr, key, branch, baseOpacity: mat.opacity, times: tr.pts.map(p => p[3]), n: tr.pts.length - 1 });
  }
  // energy deposits as glowing points
  const dg = new THREE.BufferGeometry();
  const dp = new Float32Array(sim.deps.length * 3), dc = new Float32Array(sim.deps.length * 3);
  sim.deps.forEach((d, i) => {
    dp.set([d.x, d.y, d.z], 3 * i);
    const c = new THREE.Color(COLORS[CLS_KEY[d.cls] || 'n']).multiplyScalar(1.6 * Math.min(1.6, 0.4 + Math.log10(1 + d.evis) / 2));
    dc.set([c.r, c.g, c.b], 3 * i);
  });
  dg.setAttribute('position', new THREE.BufferAttribute(dp, 3));
  dg.setAttribute('color', new THREE.BufferAttribute(dc, 3));
  dg.setDrawRange(0, 0);
  const deps = new THREE.Points(dg, new THREE.PointsMaterial({ size: 0.28, vertexColors: true, transparent: true, opacity: 0.9, blending: THREE.AdditiveBlending, depthWrite: false, toneMapped: false, sizeAttenuation: true }));
  deps.frustumCulled = false;
  group.add(deps);
  // order deposits by time for reveal
  const depOrder = sim.deps.map((d, i) => i).sort((a, b) => sim.deps[a].t - sim.deps[b].t);
  if (depOrder.length) {
    const p2 = new Float32Array(dp.length), c2 = new Float32Array(dc.length);
    depOrder.forEach((k, i) => { p2.set(dp.subarray(3 * k, 3 * k + 3), 3 * i); c2.set(dc.subarray(3 * k, 3 * k + 3), 3 * i); });
    dg.attributes.position.array.set(p2); dg.attributes.color.array.set(c2);
  }
  const depTimes = depOrder.map(i => sim.deps[i].t), depWin = depOrder.map(i => sim.deps[i].win);
  // vertex flash
  const flash = new THREE.Mesh(new THREE.SphereGeometry(1, 20, 14), new THREE.MeshBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0, blending: THREE.AdditiveBlending, depthWrite: false, toneMapped: false }));
  flash.position.set(...sim.vertex); group.add(flash);
  // scintillation light front
  const front = new THREE.Mesh(new THREE.SphereGeometry(1, 40, 24), new THREE.MeshBasicMaterial({ color: 0x7fd8ff, transparent: true, opacity: 0, blending: THREE.AdditiveBlending, depthWrite: false, side: THREE.DoubleSide, toneMapped: false, wireframe: true }));
  front.visible = false; group.add(front);
  // incoming neutrino (dashed, not detected)
  let nuLine = null;
  if (sim.nu) {
    const d = sim.nu.dir, s = sim.vertex.map((v, i) => v - d[i] * 36);
    const geo = new LineGeometry(); geo.setPositions([...s, ...sim.vertex]);
    const mat = new LineMaterial({ color: 0xe6f0fa, linewidth: 1.6, dashed: true, dashSize: 0.8, gapSize: 0.5, transparent: true, opacity: 0.85, toneMapped: false });
    mat.resolution.copy(resolution);
    nuLine = new Line2(geo, mat); nuLine.computeLineDistances(); nuLine.frustumCulled = false; nuLine.visible = false;
    group.add(nuLine);
  }
  // TT crossing markers
  const ttm = new THREE.Group();
  for (const h of sim.tt) {
    const m = new THREE.Mesh(new THREE.RingGeometry(0.35, 0.55, 20), new THREE.MeshBasicMaterial({ color: 0x9dffb8, transparent: true, opacity: 0.95, side: THREE.DoubleSide, toneMapped: false }));
    m.position.set(h.x, h.y, h.z + 0.05); ttm.add(m);
  }
  ttm.visible = false; group.add(ttm);
  // photon rays
  const NR = 600;
  const rg = new THREE.BufferGeometry(); rg.setAttribute('position', new THREE.BufferAttribute(new Float32Array(NR * 6), 3));
  const rays = new THREE.LineSegments(rg, new THREE.LineBasicMaterial({ color: new THREE.Color(COLORS.photon).multiplyScalar(1.4), transparent: true, opacity: 0.75, blending: THREE.AdditiveBlending, depthWrite: false, toneMapped: false }));
  rays.frustumCulled = false; rg.setDrawRange(0, 0); group.add(rays);
  // reconstruction overlay
  const reco = new THREE.Group(); reco.visible = false; group.add(reco);

  return { group, items, deps, depTimes, depWin, flash, front, nuLine, ttm, rays, NR, reco, photonData: null };
}

export function setReco(vis, reco, resolution) {
  vis.reco.clear();
  if (!reco) return;
  const m = new THREE.Mesh(new THREE.TorusGeometry(0.7, 0.08, 8, 32), new THREE.MeshBasicMaterial({ color: 0x6e8bff, toneMapped: false }));
  m.position.set(...reco.vertex); vis.reco.add(m);
  const m2 = m.clone(); m2.rotation.x = Math.PI / 2; vis.reco.add(m2);
  if (reco.track) {
    const geo = new LineGeometry(); geo.setPositions([...reco.track.entry, ...reco.track.exit]);
    const mat = new LineMaterial({ color: 0x8ea6ff, linewidth: 1.6, dashed: true, dashSize: 0.6, gapSize: 0.4, toneMapped: false, transparent: true, opacity: 0.9 });
    mat.resolution.copy(resolution);
    const l = new Line2(geo, mat); l.computeLineDistances(); vis.reco.add(l);
  }
}

export function disposeGroup(g) {
  g.traverse(o => { if (o.geometry) o.geometry.dispose(); if (o.material) (Array.isArray(o.material) ? o.material : [o.material]).forEach(m => m.dispose()); });
  if (g.parent) g.parent.remove(g);
}
