// JUNO // Event Display: interface state, rendering loop, animation, panels.
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { EffectComposer } from 'three/addons/postprocessing/EffectComposer.js';
import { RenderPass } from 'three/addons/postprocessing/RenderPass.js';
import { UnrealBloomPass } from 'three/addons/postprocessing/UnrealBloomPass.js';
import { OutputPass } from 'three/addons/postprocessing/OutputPass.js';
import { LineGeometry } from 'three/addons/lines/LineGeometry.js';
import { GEO, PMT, OPT, buildPmts, buildWcdPmts } from './detector.js';
import { generate, simulate, makeResponseModel, makeRng, summarize, CATEGORIES, fromGst, guessUnitScale, parseHonda, buildFluxSampler, pinfo, zenithOf, V } from './physics.js';
import { buildDetector, applyMode, wedgePlanes, buildEventVisuals, setReco, disposeGroup, COLORS, CLS_KEY, MODES } from './scene.js';
import { drawPmtMap, drawLegend, drawTimeHist, drawProjection, hitScale, pmtAngles } from './plots.js';
import { readGst } from './gst.js';

const $ = (s) => document.querySelector(s);
const $$ = (s) => [...document.querySelectorAll(s)];
const reduceMotion = matchMedia('(prefers-reduced-motion: reduce)').matches;
const hex = (c) => '#' + new THREE.Color(c).getHexString();
const CSSC = Object.fromEntries(Object.entries(COLORS).map(([k, v]) => [k, hex(v)]));
const KEYNAME = { mu: 'Muons', e: 'Electrons', gamma: 'Photons γ', chad: 'Charged hadrons', p: 'Protons', n: 'Neutrons', nu: 'Neutrinos' };

// ------------------------------------------------------------------ state
const S = {
  source: 'demo', cat: 'atm_numu', seed: 2026, n: 20, k: 0,
  gen: { flav: 'cat', cur: 'cat', mode: 'auto', emode: 'fixed', E: 4, gamma: 2.7, Emin: 1, Emax: 20, dmode: 'fixed', cosZ: 0.45, phi: 35, vertex: 'fv', vx: 0, vy: 0, vz: 0,
    ibdMode: 'reactor', ibdE: 4, muMode: 'spectrum', muE: 200, muDir: 'sampled', through: true, spall: true, flux: null, fluxName: '' },
  genie: { rows: [], name: '', total: 0, unit: 1, first: 0 },
  mode: 'cutaway', cutAngle: 110, opacity: 0.5, front: false, lab: {}, truth: true, orbit: false, dark: false, spmtMap: false, autoWin: true,
  playing: false, rate: 0.5, filters: { mu: 1, e: 1, gamma: 1, chad: 1, p: 1, n: 1, nu: 1 }, pmtColor: 'charge', photons: 250,
  win: 0, sel: null, selPmt: null, plotTab: 'time',
};
const PHI0 = -41.5; // cutaway faces the default camera azimuth

// ------------------------------------------------------------------ renderer
const stage = $('#stage'), canvas = $('#gl');
let renderer;
try {
  renderer = new THREE.WebGLRenderer({ canvas, antialias: true, preserveDrawingBuffer: true, powerPreference: 'high-performance' });
} catch (e) { $('#nogl').hidden = false; throw e; }
renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 1.75));
renderer.localClippingEnabled = true;
renderer.toneMapping = THREE.ACESFilmicToneMapping; renderer.toneMappingExposure = 1.05;
const scene = new THREE.Scene();
scene.fog = new THREE.FogExp2(0x050d0c, 0.0045);
const camera = new THREE.PerspectiveCamera(38, 1, 0.1, 800);
camera.up.set(0, 0, 1);
const controls = new OrbitControls(camera, canvas);
controls.enableDamping = true; controls.dampingFactor = 0.08; controls.minDistance = 3; controls.maxDistance = 260; controls.autoRotateSpeed = 0.6;
scene.add(new THREE.AmbientLight(0x6f86a6, 0.9));
scene.add(new THREE.HemisphereLight(0xbcd6ff, 0x1a2230, 0.7));
const sun = new THREE.DirectionalLight(0xfff2dc, 1.5); sun.position.set(40, -50, 70); scene.add(sun);
const fill = new THREE.DirectionalLight(0x6f9bff, 0.6); fill.position.set(-50, 40, -20); scene.add(fill);

const composer = new EffectComposer(renderer);
composer.addPass(new RenderPass(scene, camera));
const bloom = new UnrealBloomPass(new THREE.Vector2(512, 512), 0.75, 0.45, 0.78);
composer.addPass(bloom);
composer.addPass(new OutputPass());
const resolution = new THREE.Vector2(1, 1);

// ------------------------------------------------------------------ detector
const pm = buildPmts(), wpm = buildWcdPmts();
const det = buildDetector(pm, wpm);
scene.add(det.root);
let planes = wedgePlanes(PHI0, S.cutAngle);
const model = makeResponseModel(pm, wpm); // main-thread copy for batch estimates and fallback
const selRing = new THREE.Mesh(new THREE.RingGeometry(0.42, 0.56, 28), new THREE.MeshBasicMaterial({ color: 0xffffff, toneMapped: false, side: THREE.DoubleSide, depthTest: false }));
selRing.visible = false; selRing.renderOrder = 10; det.groups.sssAll.add(selRing);

function refreshMode() { planes = wedgePlanes(PHI0, S.cutAngle); applyMode(det, S.mode, S.opacity, planes); }

// ------------------------------------------------------------------ worker
let worker = null, workerReady = false, jobId = 0;
try {
  worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
  worker.onmessage = (e) => {
    const d = e.data;
    if (d.ready) { workerReady = true; return; }
    if (!cur || d.job !== cur.job) return;
    if (d.done) { cur.respMs = performance.now() - cur.respStart; updateDiag(); return; }
    onResponse(d.w, d.r);
  };
  worker.onerror = () => { worker = null; };
} catch { worker = null; }

function requestResponse() {
  cur.job = ++jobId; cur.respStart = performance.now();
  const sim = cur.sim, payload = { deps: sim.deps, wcd: sim.wcd, windows: sim.windows.map(w => ({ t0: w.t0, evis: w.evis, edep: w.edep, label: w.label })) };
  const seed = (cur.ev.seed * 977 + (cur.ev.iev || 0)) >>> 0;
  if (worker) worker.postMessage({ job: cur.job, sim: payload, seed, dark: S.dark });
  else {
    const job = cur.job;
    const run = (w) => { if (!cur || cur.job !== job || w >= sim.windows.length) { if (cur) { cur.respMs = performance.now() - cur.respStart; updateDiag(); } return; }
      onResponse(w, model.windowResponse(payload, w, makeRng((seed * 131 + w * 977 + 7) >>> 0), { dark: S.dark })); setTimeout(() => run(w + 1), 0); };
    setTimeout(() => run(0), 30);
  }
}

// ------------------------------------------------------------------ events
let cur = null;
const logList = [];

function genOpts() {
  const g = S.gen, o = { ...g };
  if (CATEGORIES[S.cat].kind === 'cosmic') o.dmode = g.muDir === 'fixed' ? 'fixed' : 'sampled';
  if (CATEGORIES[S.cat].kind === 'ibd') o.dmode = 'sampled';
  o.fluxSampler = g.dmode === 'flux' && g.flux ? g.flux.sampler : null;
  if (g.dmode === 'flux' && !o.fluxSampler) o.dmode = 'iso';
  return o;
}

function makeEvent(k) {
  if (S.source === 'genie') { const r = S.genie.rows[k]; return r ? fromGst(r, k, S.genie.name, S.genie.unit) : null; }
  return generate(S.cat, S.seed + k, genOpts());
}

function showEvent(k, { keepCamera = true } = {}) {
  const ev = makeEvent(k);
  if (!ev) return;
  if (cur) disposeGroup(cur.vis.group);
  const t0 = performance.now();
  const sim = simulate(ev);
  const vis = buildEventVisuals(sim, resolution);
  scene.add(vis.group);
  vis.group.traverse(o => { if (o.material && o.material.isLineMaterial) o.material.clippingPlanes = null; });
  cur = { k, ev, sim, vis, res: [], summary: summarize(sim), simMs: performance.now() - t0, job: 0, respMs: null, rays: [] };
  cur.truthMu = truthMuonTrack(sim);
  S.win = 0; S.sel = null; S.selPmt = null; selRing.visible = false;
  det.hitGlow.count = 0; det.wcdGlow.count = 0;
  buildWindowsUI(); renderObjList(); renderReadout(); renderHeader(); renderSummary();
  requestResponse();
  addLog();
  startAnim(0);
  if (!keepCamera && camPreset === 4) setCamera(4);
}

function truthMuonTrack(sim) {
  const mu = sim.particles.find(p => p.cls === 'mu' && p.parent < 0);
  if (!mu || !mu.track) return null;
  const inside = mu.track.pts.filter(p => Math.hypot(p[0], p[1], p[2]) < GEO.R_LS);
  if (inside.length < 2) return null;
  const a = inside[0], b = inside[inside.length - 1];
  return [[a[0], a[1], a[2]], [b[0], b[1], b[2]]];
}

function onResponse(w, r) {
  cur.res[w] = r;
  if (w === S.win) { prepareHits(); setReco(cur.vis, S.win === w ? r.reco : null, resolution); cur.vis.reco.traverse(o => { if (o.material?.isLineMaterial) o.material.resolution.copy(resolution); }); }
  if (w === 0) { const L = logList.find(l => l.k === cur.k && l.src === S.source); if (L) { L.npe = r.npeL; renderLog(); } }
  renderSummary(); renderPlots(); renderReadout();
}

// sorted hit order for the current window, glow matrices and colours
function prepareHits() {
  const r = cur.res[S.win];
  if (!r) { cur.hitOrder = null; det.hitGlow.count = 0; det.wcdGlow.count = 0; return; }
  const order = Array.from(r.idx.keys()).sort((a, b) => r.t[a] - r.t[b]);
  cur.hitOrder = order;
  const col = hitScale(r, S.pmtColor), c = new THREE.Color();
  let qmax = 1; for (let i = 0; i < r.npe.length; i++) qmax = Math.max(qmax, r.npe[i]);
  order.forEach((i, k) => {
    det.hitGlow.setMatrixAt(k, det.glowMatrix(r.idx[i]));
    const rgb = col(i), b = 0.55 + 1.15 * Math.pow(Math.log(1 + r.npe[i]) / Math.log(1 + qmax), 3);
    c.setRGB((rgb[0] / 255) * b, (rgb[1] / 255) * b, (rgb[2] / 255) * b, THREE.SRGBColorSpace);
    det.hitGlow.setColorAt(k, c);
  });
  det.hitGlow.instanceMatrix.needsUpdate = true; det.hitGlow.instanceColor.needsUpdate = true;
  let tmin = Infinity; for (let i = 0; i < r.t.length; i++) if (!r.dark[i]) tmin = Math.min(tmin, r.t[i]);
  const ts = order.map(i => r.t[i]).filter(Number.isFinite);
  cur.hitT0 = Number.isFinite(tmin) ? tmin : 0;
  cur.hitT1 = ts.length ? Math.max(cur.hitT0 + 20, Math.min(cur.hitT0 + 160, ts[Math.floor(ts.length * 0.95)])) : cur.hitT0 + 50;
  // WCD
  const wq = r.wcd; let wmax = 1; wq.npe.forEach(v => (wmax = Math.max(wmax, v)));
  wq.idx.forEach((j, k) => { det.wcdGlow.setMatrixAt(k, det.wcdGlowMatrix(j)); const f = Math.log(1 + wq.npe[k]) / Math.log(1 + wmax); c.setRGB(0.3 + 0.9 * f, 0.6 + 0.8 * f, 0.9 + 0.7 * f); det.wcdGlow.setColorAt(k, c); });
  det.wcdGlow.instanceMatrix.needsUpdate = true; if (det.wcdGlow.instanceColor) det.wcdGlow.instanceColor.needsUpdate = true;
  cur.wcdN = wq.idx.length;
  if (anim && anim.p >= 1) applyAnim(1);
}

// ------------------------------------------------------------------ animation
let anim = null;
function windowSpan(w) {
  const sim = cur.sim, t0 = sim.windows[w].t0;
  let tEnd = 20;
  for (const it of cur.vis.items) if (it.tr.win === w && !it.branch && it.tr.style !== 'invisible') tEnd = Math.max(tEnd, Math.min(300, it.times[it.times.length - 1] - t0));
  return { t0, tEnd: Math.min(tEnd, 260) };
}
function buildRays(w) {
  const sim = cur.sim, deps = sim.deps.filter(d => d.win === w), NR = Math.min(cur.vis.NR, S.photons);
  const rays = [];
  if (!deps.length || !NR) return rays;
  const rng = makeRng(cur.ev.seed * 31 + w + 5), tot = deps.reduce((s, d) => s + d.evis + d.cher, 0);
  const { t0, tEnd } = windowSpan(w);
  for (let k = 0; k < NR; k++) {
    let x = rng() * tot, d = deps[0]; for (const dd of deps) { x -= dd.evis + dd.cher; if (x <= 0) { d = dd; break; } }
    let u;
    if (d.cher > 0 && d.u && rng() < d.cher / (d.evis + d.cher) * 3) u = rotateCone(d.u, Math.acos(1 / OPT.N_LS), rng);
    else u = rng.iso();
    const o = [d.x, d.y, d.z], b = V.dot(o, u), c = V.dot(o, o) - GEO.R_PMT * GEO.R_PMT, L = -b + Math.sqrt(Math.max(0, b * b - c));
    rays.push({ o, u, L, e: Math.min(1, (d.t - t0) / tEnd) });
  }
  return rays;
}
function rotateCone(u, th, rng) { const a = Math.abs(u[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]; const e1 = V.norm(V.cross(u, a)), e2 = V.cross(u, e1), p = 2 * Math.PI * rng(); return V.norm(V.add(V.add(V.mul(u, Math.cos(th)), V.mul(e1, Math.sin(th) * Math.cos(p))), V.mul(e2, Math.sin(th) * Math.sin(p)))); }

function startAnim(w, banner) {
  S.win = w;
  cur.rays = buildRays(w);
  const span = windowSpan(w);
  const front = cur.sim.deps.filter(d => d.win === w).sort((a, b) => a.t - b.t)[0];
  cur.frontC = front ? [front.x, front.y, front.z] : cur.sim.vertex;
  const dur = Math.max(450, (1000 / S.rate) * (w === 0 ? 0.85 : 0.6));
  anim = { start: performance.now(), dur, w, span, p: 0, done: false };
  setReco(cur.vis, cur.res[w]?.reco || null, resolution);
  prepareHits(); buildWindowsUI(); renderSummary(); renderPlots();
  const b = $('#banner');
  if (banner) { b.textContent = banner; b.hidden = false; } else b.hidden = true;
  if (reduceMotion) { anim.p = 1; applyAnim(1); finishAnim(); }
}

function applyAnim(p) {
  const v = cur.vis, sim = cur.sim, w = anim.w, { t0, tEnd } = anim.span;
  const truth = S.truth;
  // neutrino line
  if (v.nuLine) {
    v.nuLine.visible = truth && w === 0 && S.filters.nu;
    if (v.nuLine.visible) {
      const f = Math.min(1, p / 0.18), d = sim.nu.dir, s = sim.vertex.map((x, i) => x - d[i] * 36), e = s.map((x, i) => x + (sim.vertex[i] - x) * f);
      v.nuLine.geometry.setPositions([...s, ...e]); v.nuLine.computeLineDistances();
    }
  }
  // vertex flash
  const fp = (p - 0.17) / 0.16;
  v.flash.material.opacity = fp > 0 && fp < 1 && truth ? 0.9 * (1 - fp) : 0;
  v.flash.scale.setScalar(0.3 + 2.6 * Math.max(0, fp));
  if (w > 0) v.flash.position.set(...cur.frontC); else v.flash.position.set(...sim.vertex);
  // tracks
  const tau = Math.max(0, Math.min(1, (p - (w === 0 && sim.nu ? 0.18 : 0.04)) / 0.42)) * tEnd + t0;
  const selOn = S.sel != null;
  for (const it of v.items) {
    const key = it.key, show = truth && S.filters[key] && it.tr.win <= w;
    it.line.visible = show;
    if (!show) continue;
    let n;
    if (it.tr.win < w) n = it.n; else { n = 0; while (n < it.n && it.times[n + 1] <= tau) n++; if (tau >= t0 + tEnd - 1e-6) n = it.n; }
    it.line.geometry.instanceCount = n;
    const dim = it.tr.win < w ? 0.25 : 1;
    it.mat.opacity = it.baseOpacity * dim * (selOn && it.tr.pid !== S.sel ? 0.2 : 1);
    it.mat.linewidth = (it.branch ? 1.1 : ({ mu: 3.4, e: 2.6, gamma: 1.8, chad: 2.6, p: 2.6, n: 1.4, nu: 1.4 }[key] || 2)) * (S.sel === it.tr.pid ? 1.7 : 1);
  }
  // deposits
  let nd = 0; while (nd < v.depTimes.length && v.depTimes[nd] <= tau) nd++;
  v.deps.geometry.setDrawRange(0, truth ? nd : 0);
  // light front
  const lf = (p - 0.24) / 0.5;
  v.front.visible = S.front && lf > 0 && lf < 1 && cur.sim.deps.some(d => d.win === w);
  if (v.front.visible) { v.front.position.set(...cur.frontC); v.front.scale.setScalar(0.5 + lf * 24); v.front.material.opacity = 0.16 * (1 - lf); }
  // photon rays
  const rg = v.rays.geometry, arr = rg.attributes.position.array; let nr = 0;
  if (p > 0.22 && p < 0.95) for (const r of cur.rays) {
    const emit = 0.22 + 0.25 * r.e, head = ((p - emit) / 0.32) * r.L;
    if (head <= 0 || head - 2 > r.L) continue;
    const a = Math.max(0, head - 2.2), b = Math.min(r.L, head);
    arr.set([r.o[0] + r.u[0] * a, r.o[1] + r.u[1] * a, r.o[2] + r.u[2] * a, r.o[0] + r.u[0] * b, r.o[1] + r.u[1] * b, r.o[2] + r.u[2] * b], 6 * nr); nr++;
  }
  rg.attributes.position.needsUpdate = true; rg.setDrawRange(0, nr * 2);
  // PMT hits in time order
  if (cur.hitOrder) {
    const r = cur.res[w], f = Math.max(0, Math.min(1, (p - 0.4) / 0.45)), tcut = cur.hitT0 + f * (cur.hitT1 - cur.hitT0) + (f >= 1 ? 1e9 : 0);
    let n = 0; const ord = cur.hitOrder; while (n < ord.length && r.t[ord[n]] <= tcut) n++;
    det.hitGlow.count = p >= 0.4 ? n : 0; cur.revealN = det.hitGlow.count;
    det.wcdGlow.count = p > 0.6 ? cur.wcdN : 0;
  } else { det.hitGlow.count = 0; det.wcdGlow.count = 0; cur.revealN = 0; }
  v.ttm.visible = p > 0.3 && w === 0;
  v.reco.visible = p > 0.86;
  // sequence steps
  const step = p < 0.17 ? 0 : p < 0.33 ? 1 : p < 0.6 ? 2 : p < 0.86 ? 3 : 4;
  $$('#steps .step').forEach((el, i) => { el.classList.toggle('on', i === step && p < 1); el.classList.toggle('done', i < step || p >= 1); });
}

let lastMapDraw = 0;
function tickAnim(now) {
  if (!anim || !cur) return;
  if (!anim.done) {
    anim.p = Math.min(1, (now - anim.start) / anim.dur);
    applyAnim(anim.p);
    if (now - lastMapDraw > 90) { renderMap(); lastMapDraw = now; }
    if (anim.p >= 1) finishAnim();
  } else if (S.playing && now > anim.nextAt) nextEvent();
}
function finishAnim() {
  anim.done = true; renderMap();
  const nW = cur.sim.windows.length, w = anim.w;
  if (S.autoWin && w + 1 < Math.min(nW, 4)) {
    const nw = cur.sim.windows[w + 1];
    const dt = nw.t0 - cur.sim.windows[w].t0;
    anim.nextAt = Infinity;
    setTimeout(() => { if (cur && anim && anim.w === w) startAnim(w + 1, `${nw.label.split('Δt')[0].trim()} · Δt = ${fmtTime(dt)} later (compressed)`); }, reduceMotion ? 0 : 380);
    return;
  }
  setTimeout(() => { $('#banner').hidden = true; }, 1200);
  anim.nextAt = performance.now() + (1000 / S.rate) * 0.15 + 150;
}
function fmtTime(ns) { return ns < 1e3 ? `${ns.toFixed(0)} ns` : ns < 1e6 ? `${(ns / 1e3).toFixed(ns < 1e4 ? 2 : 0)} µs` : `${(ns / 1e6).toFixed(2)} ms`; }

function nextEvent() { S.k = (S.k + 1) % runLength(); showEvent(S.k); }
function runLength() { return S.source === 'genie' ? Math.max(1, S.genie.rows.length) : S.n; }

// ------------------------------------------------------------------ camera
let camPreset = 1, camTween = null;
const PRESETS = {
  1: () => ({ p: [62, -55, 27], t: [0, 0, 1] }),
  2: () => ({ p: [0, -0.5, 82], t: [0, 0, 0] }),
  3: () => ({ p: [0, -84, 3], t: [0, 0, 1] }),
  4: () => { const v = cur ? (cur.sim.cosmic ? [0, 0, 0] : cur.sim.vertex) : [0, 0, 0]; return { p: [v[0] + 9, v[1] - 9, v[2] + 6], t: v }; },
  5: () => ({ p: [84, -76, 52], t: [0, 0, 4] }),
  6: () => ({ p: [2.5, -18.6, -3.5], t: [9, -15.8, 2.5] }),
};
function setCamera(i, instant = false) {
  camPreset = i; let { p, t } = PRESETS[i]();
  if (i !== 4 && i !== 6) { // pull back on narrow screens so the whole detector fits
    const a = stage.clientWidth / Math.max(1, stage.clientHeight), k = (a < 1.25 ? Math.pow(1.25 / a, 0.85) : 1) * (stage.clientWidth >= 1100 ? 1.12 : 1);
    p = p.map((x, j) => t[j] + (x - t[j]) * k);
  }
  $('#s-cam').value = String(i);
  if (instant || reduceMotion) { camera.position.set(...p); controls.target.set(...t); controls.update(); return; }
  camTween = { from: camera.position.clone(), to: new THREE.Vector3(...p), ft: controls.target.clone(), tt: new THREE.Vector3(...t), start: performance.now(), dur: 750 };
}
function tickCam(now) {
  if (!camTween) return;
  const f = Math.min(1, (now - camTween.start) / camTween.dur), e = f < 0.5 ? 2 * f * f : 1 - Math.pow(-2 * f + 2, 2) / 2;
  camera.position.lerpVectors(camTween.from, camTween.to, e); controls.target.lerpVectors(camTween.ft, camTween.tt, e);
  if (f >= 1) camTween = null;
}

// ------------------------------------------------------------------ labels
const labelEls = det.labels.map(l => { const d = document.createElement('div'); d.className = 'lbl'; d.textContent = l.text; $('#labels').appendChild(d); return d; });
const tmpV = new THREE.Vector3();
function tickLabels() {
  const W = canvas.clientWidth, H = canvas.clientHeight;
  det.labels.forEach((l, i) => {
    const g = det.groups[l.g], el = labelEls[i];
    let vis = !!S.lab[l.g] && g && (g.parent ? isVisible(g.parent) : true);
    if (vis) { tmpV.set(...l.p).applyMatrix4(g.matrixWorld).project(camera); vis = tmpV.z < 1 && Math.abs(tmpV.x) < 1.05 && Math.abs(tmpV.y) < 1.05; }
    el.style.display = vis ? '' : 'none';
    if (vis) el.style.transform = `translate(${((tmpV.x + 1) / 2) * W + 10}px, ${((1 - tmpV.y) / 2) * H}px) translateY(-50%)`;
  });
}
function isVisible(o) { for (let x = o; x; x = x.parent) if (!x.visible) return false; return true; }

// ------------------------------------------------------------------ main loop
let frames = 0, fpsT = performance.now(), fps = 0;
function loop(now) {
  requestAnimationFrame(loop);
  tickCam(now); controls.autoRotate = S.orbit && !reduceMotion; controls.update();
  tickAnim(now);
  composer.render();
  tickLabels();
  frames++; if (now - fpsT > 1000) { fps = (frames * 1000) / (now - fpsT); frames = 0; fpsT = now; updateDiag(); }
}
function resize() {
  const w = stage.clientWidth, h = stage.clientHeight;
  renderer.setSize(w, h, false); composer.setSize(w, h);
  camera.aspect = w / h;
  // keep the detector centred in the free area between the floating panels (desktop layout)
  if (w >= 1100 && !hud.classList.contains('panels-off')) { const L = $('.col-left').getBoundingClientRect().width, R = $('.col-right').getBoundingClientRect().width, B = $('#p-map').getBoundingClientRect().height + 12; camera.setViewOffset(w, h, (R - L) / 2, B / 2, w, h); }
  else camera.clearViewOffset();
  camera.updateProjectionMatrix();
  renderer.getDrawingBufferSize(resolution);
  if (cur) cur.vis.group.traverse(o => { if (o.material?.isLineMaterial) o.material.resolution.copy(resolution); });
  renderPlots();
}
new ResizeObserver(resize).observe(stage);
window.addEventListener('resize', () => renderPlots());

// ------------------------------------------------------------------ picking
let downXY = null;
canvas.addEventListener('pointerdown', (e) => { downXY = [e.clientX, e.clientY]; });
canvas.addEventListener('pointerup', (e) => {
  if (!downXY || Math.hypot(e.clientX - downXY[0], e.clientY - downXY[1]) > 5 || !cur) return;
  const rect = canvas.getBoundingClientRect(), mx = e.clientX - rect.left, my = e.clientY - rect.top;
  // tracks first: screen-space distance to segments
  let best = null, bd = 10;
  const pv = new THREE.Vector3(), W = rect.width, H = rect.height;
  const proj = (p) => { pv.set(p[0], p[1], p[2]).project(camera); return pv.z < 1 ? [((pv.x + 1) / 2) * W, ((1 - pv.y) / 2) * H] : null; };
  for (const it of cur.vis.items) {
    if (!it.line.visible || it.branch || it.tr.style === 'invisible') continue;
    const n = it.line.geometry.instanceCount; let prev = null;
    for (let i = 0; i <= Math.min(n, it.n); i++) {
      const q = proj(it.tr.pts[i]);
      if (q && prev) { const d = segDist(mx, my, prev, q); if (d < bd) { bd = d; best = it.tr.pid; } }
      prev = q;
    }
  }
  if (best != null) { selectParticle(best); return; }
  // PMTs
  const ray = new THREE.Raycaster(); ray.setFromCamera(new THREE.Vector2((mx / W) * 2 - 1, -(my / H) * 2 + 1), camera);
  const hits = ray.intersectObject(det.lpmt, false);
  const cut = S.mode !== 'transparent';
  const h = hits.find(x => !cut || !(planes[0].distanceToPoint(x.point) < 0 && planes[1].distanceToPoint(x.point) < 0));
  if (h && h.instanceId != null) { selectPmt(h.instanceId); return; }
  clearSelection();
});
function segDist(px, py, a, b) { const dx = b[0] - a[0], dy = b[1] - a[1], l2 = dx * dx + dy * dy; let t = l2 ? ((px - a[0]) * dx + (py - a[1]) * dy) / l2 : 0; t = Math.max(0, Math.min(1, t)); return Math.hypot(px - a[0] - t * dx, py - a[1] - t * dy); }

function selectParticle(pid) { S.sel = pid; S.selPmt = null; selRing.visible = false; if (anim) applyAnim(anim.done ? 1 : anim.p); renderObjList(); renderReadout(); renderPlots(); }
function selectPmt(j) {
  S.selPmt = j; S.sel = null;
  const p = new THREE.Vector3(pm.pos[3 * j], pm.pos[3 * j + 1], pm.pos[3 * j + 2]).multiplyScalar(0.98);
  selRing.position.copy(p); selRing.quaternion.setFromUnitVectors(new THREE.Vector3(0, 0, 1), p.clone().negate().normalize()); selRing.visible = true;
  if (anim) applyAnim(anim.done ? 1 : anim.p); renderObjList(); renderReadout(); renderPlots();
}
function clearSelection() { S.sel = null; S.selPmt = null; selRing.visible = false; if (anim && cur) applyAnim(anim.done ? 1 : anim.p); renderObjList(); renderReadout(); renderPlots(); }

// ------------------------------------------------------------------ formatting
const fE = (gev) => (gev == null || !Number.isFinite(gev) ? '—' : gev >= 1 ? `${gev.toFixed(gev >= 100 ? 1 : 3)} GeV` : `${(gev * 1000).toFixed(gev * 1000 >= 10 ? 1 : 3)} MeV`);
const fMeV = (mev) => (mev >= 1000 ? `${(mev / 1000).toFixed(3)} GeV` : `${mev.toFixed(mev >= 10 ? 1 : 3)} MeV`);
const fN = (n) => (n >= 1e6 ? `${(n / 1e6).toFixed(2)}M` : n >= 1e4 ? `${(n / 1e3).toFixed(1)}k` : `${Math.round(n)}`);
const fV = (v) => `(${v.map(x => x.toFixed(2)).join(', ')})`;
const esc = (s) => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
function catLabel(ev) {
  if (ev.source === 'genie') { const n = pinfo(ev.nu.pdg).name; return `${n} ${ev.interaction.cc ? 'CC' : 'NC'} · ${ev.interaction.mode}`; }
  return CATEGORIES[ev.category].label;
}

// ------------------------------------------------------------------ panels
function renderHeader() {
  const ev = cur.ev;
  $('#h-num').textContent = String(ev.source === 'genie' ? ev.iev ?? cur.k : cur.k + 1).padStart(6, '0');
  $('#h-cat').textContent = catLabel(ev);
  const b = $('#h-badge');
  b.textContent = ev.source === 'genie' ? 'SIMULATED · GENIE FILE' : 'SIMULATED · DEMO'; b.classList.toggle('genie', ev.source === 'genie');
  $('#h-meta').textContent = ev.source === 'genie' ? `${S.genie.name} · entry ${S.genie.first + cur.k}` : `seed ${ev.seed} · event ${cur.k + 1} of ${S.n}`;
}

function renderSummary() {
  if (!cur) return;
  const ev = cur.ev, sim = cur.sim, s = cur.summary, w = S.win, r = cur.res[w], win = sim.windows[w];
  $('#sum-win').textContent = w === 0 ? 'PROMPT WINDOW' : `WINDOW ${w + 1} / ${sim.windows.length}`;
  $('#s-big').textContent = fMeV(win.evis);
  $('#s-bigl').textContent = 'visible energy in LS (truth, this window)';
  const kv = [];
  if (ev.nu) {
    kv.push(['Neutrino', pinfo(ev.nu.pdg).name], ['E_ν', fE(ev.nu.E)], ['cos θ_z', `${ev.nu.cosZ.toFixed(2)} ${ev.nu.cosZ >= 0 ? '↓' : '↑'}`]);
    kv.push(['Current', ev.interaction.mode === 'IBD' ? 'CC' : ev.interaction.cc ? 'CC' : 'NC'], ['Mode', ev.interaction.mode], ['Q²', ev.interaction.Q2 != null ? `${ev.interaction.Q2.toFixed(3)} GeV²` : '—']);
    kv.push(['Vertex R', `${s.R.toFixed(2)} m`], ['Fiducial', s.inFV ? '<span class="pill yes">INSIDE</span>' : '<span class="pill no">OUTSIDE</span>'], ['W', ev.interaction.W ? `${ev.interaction.W.toFixed(3)} GeV` : '—']);
  } else if (ev.muon) {
    kv.push(['Muon E', fE(ev.muon.E)], ['cos θ_z', ev.muon.cosZ.toFixed(2)], ['Impact b', `${ev.muon.impact.toFixed(1)} m`]);
    const mu = sim.particles[0]; kv.push(['Track in LS', `${mu.lenLS.toFixed(1)} m`], ['TT hits', `${sim.tt.length} layers`], ['Charge', mu.q > 0 ? 'μ⁺' : 'μ⁻']);
  }
  kv.push(['E_dep LS', fMeV(win.edep)], ['NPE 20"', r ? fN(r.npeL) : '…'], ['Fired 20"', r ? `${fN(r.firedL)} / ${fN(pm.nL)}` : '…']);
  kv.push(['E_rec ≈', r && r.reco ? fMeV(r.reco.E) : '…'], ['NPE 3"', r ? fN(r.npeS) : '…'], ['Trigger', r ? (r.triggered ? `<span class="pill yes">YES</span> ${r.maxIn}` : `<span class="pill no">NO</span> ${r.maxIn}`) : '…']);
  kv.push(['WCD PMTs', r ? `${r.wcd.fired}${r.wcd.tagged ? ' · μ tag' : ''}` : '…'], ['Windows', `${sim.windows.length}`], [ev.nu ? 'Not in LS' : 'Neutrons', ev.nu ? fMeV(Math.max(0, s.missing)) : `${sim.particles.filter(p => p.cls === 'n').length}`]);
  if (r && r.reco && !sim.cosmic && w === 0) kv.push(['Vtx resid.', `${V.len(V.sub(r.reco.vertex, sim.vertex)).toFixed(2)} m`]);
  if (r && r.reco && r.reco.track && cur.truthMu && w === 0) { const td = V.norm(V.sub(cur.truthMu[1], cur.truthMu[0])); kv.push(['Track Δangle', `${((Math.acos(Math.max(-1, Math.min(1, V.dot(td, r.reco.track.dir)))) * 180) / Math.PI).toFixed(1)}°`]); }
  $('#s-kv').innerHTML = kv.map(([k, v]) => `<div><span>${esc(k)}</span><b>${v}</b></div>`).join('');
}

function renderObjList() {
  if (!cur) return;
  const prim = cur.sim.particles.filter(p => p.parent < 0 || p.origin === 'n capture' || p.origin === 'decay' && p.cls === 'e');
  $('#s-cnt').textContent = `${cur.sim.particles.length} particles`;
  $('#s-obj').innerHTML = prim.slice(0, 40).map(p => {
    const key = CLS_KEY[p.cls] || 'n';
    return `<button class="obj" data-pid="${p.id}" aria-pressed="${S.sel === p.id}"><i class="sw" style="background:${CSSC[key]}"></i><span>${esc(p.name)} <small>${esc(p.origin)}</small></span><em>${fE(p.E - p.m > 0 && p.cls !== 'nu' && p.cls !== 'gamma' ? p.E - p.m : p.E)}</em></button>`;
  }).join('');
}
$('#s-obj').addEventListener('click', (e) => { const b = e.target.closest('.obj'); if (!b) return; const pid = +b.dataset.pid; S.sel === pid ? clearSelection() : selectParticle(pid); });

function renderReadout() {
  const el = $('#s-read');
  if (!cur || (S.sel == null && S.selPmt == null)) { el.innerHTML = '<dt>Tip</dt><dd>click a track or PMT</dd>'; return; }
  const rows = [];
  if (S.sel != null) {
    const p = cur.sim.particles[S.sel], u = V.norm(p.p), z = zenithOf(u), par = p.parent >= 0 ? cur.sim.particles[p.parent] : null;
    const th = (Math.acos(u[2]) * 180) / Math.PI, ph = ((Math.atan2(u[1], u[0]) * 180) / Math.PI + 360) % 360;
    rows.push(['Particle', `${p.name} (PDG ${p.pdg})`], ['Origin', `${p.origin}${par ? ' of ' + par.name : ''}`], ['Charge', p.q ? (p.q > 0 ? '+1 e' : '−1 e') : '0']);
    rows.push(['Energy', fE(p.E)], ['Kinetic', fE(Math.max(0, p.E - p.m))], ['|p|', `${p.pmag.toFixed(4)} GeV/c`]);
    rows.push(['Direction θ, φ', `${th.toFixed(1)}°, ${ph.toFixed(1)}°`], ['Comes from cos θ_z', z.cosZ.toFixed(3)], ['Start', `${fV(p.x0)} m`]);
    rows.push(['Fate', p.endProcess || '—'], ['Path in LS', `${p.lenLS.toFixed(2)} m`], ['E deposited in LS', fMeV(p.edepLS)], ['Visible (quenched)', fMeV(p.evis)]);
    rows.push(['Daughters', `${cur.sim.particles.filter(q => q.parent === p.id).length}`]);
  } else {
    const j = S.selPmt, [th, ph] = pmtAngles(pm, j), r = cur.res[S.win];
    let hit = -1; if (r) for (let i = 0; i < r.idx.length; i++) if (r.idx[i] === j) { hit = i; break; }
    rows.push(['PMT', `#${j} · ${['20" dynode (HPK)', '20" MCP (NNVT)', '3"'][pm.type[j]]}`], ['θ, φ', `${((th * 180) / Math.PI).toFixed(1)}°, ${((ph * 180) / Math.PI).toFixed(1)}°`], ['Position', `${fV([pm.pos[3 * j], pm.pos[3 * j + 1], pm.pos[3 * j + 2]])} m`]);
    rows.push(['Charge', hit >= 0 ? `${r.npe[hit]} p.e.${r.dark[hit] ? ' (dark noise)' : ''}` : 'not hit'], ['First hit', hit >= 0 ? `${r.t[hit].toFixed(1)} ns` : '—'], ['TTS', `${PMT.TTS[pm.type[j]]} ns`]);
  }
  el.innerHTML = rows.map(([k, v]) => `<dt>${esc(k)}</dt><dd>${esc(v)}</dd>`).join('');
}

function buildWindowsUI() {
  if (!cur) return;
  $('#wins').innerHTML = cur.sim.windows.map((w, i) => `<button data-w="${i}" aria-pressed="${i === S.win}">${i === 0 ? 'Prompt' : esc(w.label.replace(/\s+Δt.*/, '').replace(' (', ' · ').replace(')', ''))}${i ? ' +' + fmtTime(w.t0) : ''}</button>`).join('');
}
$('#wins').addEventListener('click', (e) => { const b = e.target.closest('button'); if (!b || !cur) return; const w = +b.dataset.w; startAnim(w, w ? `${cur.sim.windows[w].label.replace(/\s+Δt.*/, '')} · Δt = ${fmtTime(cur.sim.windows[w].t0)} (compressed)` : null); });

function addLog() {
  const ev = cur.ev;
  const E = ev.nu ? fE(ev.nu.E) : fE(ev.muon.E);
  logList.unshift({ k: cur.k, src: S.source, text: `#${String(ev.source === 'genie' ? ev.iev ?? cur.k : cur.k + 1).padStart(5, '0')}  ${catLabel(ev)}  ${E}`, npe: null });
  if (logList.length > 30) logList.pop();
  renderLog();
}
function renderLog() { $('#log').innerHTML = logList.map((l, i) => `<div data-i="${i}" title="Show this event again">${esc(l.text)}${l.npe != null ? '  · ' + fN(l.npe) + ' p.e.' : ''}</div>`).join(''); }
$('#log').addEventListener('click', (e) => { const d = e.target.closest('div[data-i]'); if (!d) return; const L = logList[+d.dataset.i]; if (L.src !== S.source) return; S.k = L.k; showEvent(L.k); });

function updateDiag() {
  $('#diag').textContent = `APP DIAGNOSTICS · ${fps.toFixed(0)} FPS · SIM ${cur ? cur.simMs.toFixed(0) : '—'} MS · PMT RESPONSE ${cur && cur.respMs != null ? cur.respMs.toFixed(0) + ' MS' : '…'}`;
}

// plots
let mapGeom = null;
function renderMap() {
  if (!cur) return;
  const r = cur.res[S.win];
  const mu = S.win === 0 ? cur.truthMu : null;
  mapGeom = drawPmtMap($('#cv-map'), pm, r, {
    mode: S.pmtColor, waiting: !r, spmt: S.spmtMap, order: cur.hitOrder, reveal: anim && !anim.done ? cur.revealN : undefined,
    truthTrack: S.truth ? mu : null, recoTrack: r && r.reco && anim?.done ? r.reco.track : null, selPmt: S.selPmt,
  });
}
function renderPlots() {
  if (!cur) return;
  renderMap();
  drawLegend($('#cv-cbar'), S.pmtColor);
  $('#cb-mode').textContent = S.pmtColor === 'time' ? 'first-hit time (0–60 ns after first)' : 'charge (log p.e.)';
  $('#cb-lo').textContent = S.pmtColor === 'time' ? 'early' : '1 p.e.';
  const r = cur.res[S.win];
  $('#cb-hi').textContent = S.pmtColor === 'time' ? 'late' : r ? `${fN(r.npe.reduce((a, b) => Math.max(a, b), 0))} p.e.` : 'high';
  const cv = $('#cv-plot');
  if (S.plotTab === 'time') {
    let selT = null; if (S.selPmt != null && r) for (let i = 0; i < r.idx.length; i++) if (r.idx[i] === S.selPmt) selT = r.t[i];
    drawTimeHist(cv, pm, r, { waiting: !r, selT });
  } else drawProjection(cv, cur.sim, { view: S.plotTab, win: S.win, colors: CSSC, key: (c) => CLS_KEY[c] || 'n', filter: (k) => S.truth && S.filters[k], sel: S.sel, reco: r && r.reco && !cur.sim.cosmic ? r.reco : null });
}
$('#cv-map').addEventListener('click', (e) => {
  if (!mapGeom || !mapGeom.hits.length) return;
  const rect = e.target.getBoundingClientRect(), x = e.clientX - rect.left, y = e.clientY - rect.top;
  let best = null, bd = 64; for (const [hx, hy, i] of mapGeom.hits) { const d = (hx - x) ** 2 + (hy - y) ** 2; if (d < bd) { bd = d; best = i; } }
  if (best != null) selectPmt(cur.res[S.win].idx[best]);
});

// ------------------------------------------------------------------ generator panel
function field(label, html) { return `<label class="row"><span>${label}</span>${html}</label>`; }
const sel = (id, opts, val) => `<select id="${id}" style="width:150px">${opts.map(([v, t]) => `<option value="${v}"${String(v) === String(val) ? ' selected' : ''}>${t}</option>`).join('')}</select>`;
const num = (id, val, step = 'any', min = '', max = '') => `<input type="number" id="${id}" value="${val}" step="${step}"${min !== '' ? ` min="${min}"` : ''}${max !== '' ? ` max="${max}"` : ''}>`;
function renderGenFields() {
  const g = S.gen, kind = CATEGORIES[S.cat].kind;
  let h = '';
  const vtx = field('Vertex', sel('g-vertex', [['fv', 'Fiducial (R<16.5 m)'], ['ls', 'Whole LS'], ['centre', 'Detector centre'], ['fixed', 'Fixed x, y, z']], g.vertex)) +
    (g.vertex === 'fixed' ? `<div class="grid2" style="grid-template-columns:1fr 1fr 1fr">${['vx', 'vy', 'vz'].map(k => `<input type="number" id="g-${k}" value="${g[k]}" step="0.5" aria-label="${k} (m)" style="width:100%">`).join('')}</div>` : '');
  if (kind === 'atm') {
    h += field('Flavour', sel('g-flav', [['cat', 'From category'], ['14', 'νμ'], ['-14', 'ν̄μ'], ['12', 'νe'], ['-12', 'ν̄e']], g.flav));
    h += field('Current', sel('g-cur', [['cat', 'From category'], ['CC', 'CC'], ['NC', 'NC'], ['auto', 'Mixed (72% CC)']], g.cur));
    h += field('Process', sel('g-mode', [['auto', 'Auto by energy'], ['QE', 'QE'], ['RES', 'RES (Δ)'], ['DIS', 'DIS']], g.mode));
    h += field('Energy', sel('g-emode', [['fixed', 'Fixed'], ['power', 'Power law E^−γ'], ['flux', 'Honda flux file']], g.dmode === 'flux' ? 'flux' : g.emode));
    if (g.dmode !== 'flux' && g.emode === 'fixed') h += field('E_ν (GeV)', num('g-E', g.E, 0.1, 0.1, 1000));
    if (g.dmode !== 'flux' && g.emode === 'power') h += `<div class="grid2" style="grid-template-columns:1fr 1fr 1fr">${field('γ', num('g-gamma', g.gamma, 0.1))}</div><div class="grid2">${field('E min', num('g-Emin', g.Emin, 0.1))}${field('E max', num('g-Emax', g.Emax, 1))}</div>`;
    if (g.dmode === 'flux') {
      h += `<div class="grid2">${field('E min', num('g-Emin', g.Emin, 0.1))}${field('E max', num('g-Emax', g.Emax, 1))}</div>`;
      h += `<label class="btn" for="f-flux">${g.flux ? 'Replace' : 'Open'} Honda .d file</label><input type="file" id="f-flux" accept=".d,.txt,.dat" hidden><div class="status ${g.flux ? 'ok' : ''}" id="flux-status">${g.flux ? esc(`${g.fluxName}: ${g.flux.blocks.length} angular bins`) : 'Use the same file as gevgen_atmo (e.g. juno_solmax.d). Energy, flavour and direction are sampled from flux × E.'}</div>`;
    } else {
      h += field('Direction', sel('g-dmode', [['fixed', 'Fixed zenith'], ['iso', 'Isotropic']], g.dmode));
      if (g.dmode === 'fixed') h += `<div class="grid2">${field('cos θ_z', num('g-cosZ', g.cosZ, 0.05, -1, 1))}${field('φ (°)', num('g-phi', g.phi, 5, 0, 360))}</div><p class="note">cos θ_z &gt; 0: down-going from the sky. &lt; 0: up-going through the Earth.</p>`;
    }
    h += vtx;
  } else if (kind === 'ibd') {
    h += field('Energy', sel('g-ibdMode', [['reactor', 'Reactor spectrum'], ['fixed', 'Fixed']], g.ibdMode));
    if (g.ibdMode === 'fixed') h += field('E_ν (MeV)', num('g-ibdE', g.ibdE, 0.1, 1.81, 12));
    h += vtx + '<p class="note">Prompt: e⁺ (E_prompt ≈ E_ν − 0.784 MeV). Delayed: 2.223 MeV γ from n capture on H, τ ≈ 200 µs.</p>';
  } else {
    h += field('Energy', sel('g-muMode', [['spectrum', 'Underground spectrum'], ['fixed', 'Fixed']], g.muMode));
    if (g.muMode === 'fixed') h += field('E_μ (GeV)', num('g-muE', g.muE, 10, 1, 10000));
    h += field('Direction', sel('g-muDir', [['sampled', 'Sampled (steep)'], ['fixed', 'Fixed zenith']], g.muDir));
    if (g.muDir === 'fixed') h += `<div class="grid2">${field('cos θ_z', num('g-cosZ', g.cosZ, 0.05, 0.05, 1))}${field('φ (°)', num('g-phi', g.phi, 5, 0, 360))}</div>`;
    h += `<div class="grid2"><label class="chk">Through CD<input type="checkbox" id="g-through"${g.through ? ' checked' : ''}></label><label class="chk">Spallation n<input type="checkbox" id="g-spall"${g.spall ? ' checked' : ''}></label></div>`;
  }
  $('#gen-fields').innerHTML = `<div class="pb" style="padding:0;gap:7px">${h}</div>`;
}
$('#gen-fields').addEventListener('change', async (e) => {
  const t = e.target, id = t.id.replace(/^g-/, ''), g = S.gen;
  if (t.id === 'f-flux') { await loadFlux(t.files[0]); return; }
  if (t.type === 'checkbox') g[id] = t.checked;
  else if (t.type === 'number') g[id] = +t.value;
  else if (id === 'emode') { if (t.value === 'flux') { g.dmode = 'flux'; } else { g.emode = t.value; if (g.dmode === 'flux') g.dmode = 'iso'; } }
  else g[id] = t.value;
  if (['emode', 'dmode', 'vertex', 'ibdMode', 'muMode', 'muDir'].includes(id)) renderGenFields();
  if (g.flux && ['Emin', 'Emax', 'flav'].includes(id)) rebuildFlux();
});
async function loadFlux(file) {
  if (!file) return;
  try { const blocks = parseHonda(await file.text()); S.gen.flux = { blocks }; S.gen.fluxName = file.name; rebuildFlux(); }
  catch (err) { S.gen.flux = null; renderGenFields(); $('#flux-status').textContent = err.message; $('#flux-status').className = 'status err'; return; }
  renderGenFields();
}
function rebuildFlux() {
  const g = S.gen, fl = g.flav !== 'cat' ? [+g.flav] : [12, -12, 14, -14];
  try { g.flux.sampler = buildFluxSampler(g.flux.blocks, fl, g.Emin, g.Emax); } catch (err) { g.flux.sampler = null; }
}

// GENIE files
async function loadGst(file) {
  const st = $('#gst-status');
  if (!file) return;
  st.className = 'status'; st.textContent = `Reading ${file.name}…`;
  try {
    const first = Math.max(0, +$('#i-first').value || 0), max = Math.min(1000, Math.max(1, +$('#i-max').value || 200));
    const out = await readGst(file, { first, max });
    if (!out.rows.length) throw new Error('No entries in the chosen range.');
    S.genie = { rows: out.rows, name: file.name, total: out.total, unit: guessUnitScale(out.rows), first };
    st.className = 'status ok';
    st.textContent = `Read ${out.rows.length} of ${out.total.toLocaleString()} entries from ${file.name}${S.genie.unit !== 1 ? ` (vertex units converted ×${S.genie.unit})` : ''}.`;
    setSource('genie'); S.k = 0; showEvent(0);
  } catch (err) {
    st.className = 'status err';
    st.textContent = `Could not read the file: ${err.message || err}. It must be a GENIE gst ROOT file (gntpc -f gst). Check your internet connection: the ROOT reader loads from cdn.jsdelivr.net.`;
  }
}
$('#f-gst').addEventListener('change', (e) => loadGst(e.target.files[0]));
let dragN = 0;
stage.addEventListener('dragenter', (e) => { e.preventDefault(); dragN++; $('#drop').hidden = false; });
stage.addEventListener('dragleave', () => { if (--dragN <= 0) { dragN = 0; $('#drop').hidden = true; } });
stage.addEventListener('dragover', (e) => e.preventDefault());
stage.addEventListener('drop', (e) => { e.preventDefault(); dragN = 0; $('#drop').hidden = true; const f = e.dataTransfer.files[0]; if (f) loadGst(f); });

function setSource(src) {
  S.source = src;
  $$('#seg-src button').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.src === src)));
  $('#gen-demo').hidden = src !== 'demo'; $('#gen-genie').hidden = src !== 'genie';
  $('#s-cat').disabled = src === 'genie';
  $('#b-gen').textContent = src === 'genie' ? 'Show first event' : 'Generate run';
}
$('#seg-src').addEventListener('click', (e) => { const b = e.target.closest('button'); if (!b) return; setSource(b.dataset.src); if (b.dataset.src === 'demo') { S.k = 0; showEvent(0); } else if (S.genie.rows.length) { S.k = 0; showEvent(0); } });
$('#b-gen').addEventListener('click', () => {
  if (S.source === 'genie') { if (S.genie.rows.length) { S.k = 0; showEvent(0); } else $('#f-gst').click(); return; }
  S.n = Math.max(1, Math.min(1000, +$('#i-n').value || 20)); S.seed = Math.max(1, +$('#i-seed').value || 1);
  S.k = 0; showEvent(0);
});

// ------------------------------------------------------------------ export
function download(name, text, type) {
  try {
    const url = URL.createObjectURL(new Blob([text], { type }));
    const a = Object.assign(document.createElement('a'), { href: url, download: name }); document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 4000);
    $('#x-status').textContent = `Saved ${name}. If nothing downloaded, this preview blocks downloads; use your GitHub Pages site.`;
  } catch (err) { $('#x-status').textContent = `Download failed: ${err.message}`; }
}
function particleRows(sim) {
  return sim.particles.map(p => ({ id: p.id, pdg: p.pdg, name: p.name, parent: p.parent, origin: p.origin, E_GeV: +p.E.toFixed(6), p_GeV: p.p.map(x => +x.toFixed(6)), start_m: p.x0.map(x => +x.toFixed(4)), t0_ns: +p.t0.toFixed(2), end_m: p.end ? p.end.map(x => +x.toFixed(4)) : null, fate: p.endProcess, lenLS_m: +p.lenLS.toFixed(3), EdepLS_MeV: +p.edepLS.toFixed(4), EvisLS_MeV: +p.evis.toFixed(4) }));
}
function eventHeader(ev) {
  return { source: ev.source, file: ev.file || null, gst_entry: ev.source === 'genie' ? S.genie.first + ev.seed : null, iev: ev.iev ?? null, seed: ev.source === 'demo' ? ev.seed : null, category: ev.category,
    neutrino: ev.nu ? { pdg: ev.nu.pdg, E_GeV: ev.nu.E, dir: ev.nu.dir, cosZ: ev.nu.cosZ, phi_deg: ev.nu.phi } : null, muon: ev.muon || null, interaction: ev.interaction, vertex_m: ev.vertex };
}
$('#b-xjson').addEventListener('click', () => {
  if (!cur) return;
  const sim = cur.sim, out = { schema: 'juno-event-display/1', units: { length: 'm', time: 'ns', energy: 'GeV unless _MeV' }, ...eventHeader(cur.ev), particles: particleRows(sim),
    windows: sim.windows.map((w, i) => { const r = cur.res[i]; return { label: w.label, t0_ns: w.t0, Edep_LS_MeV: w.edep, Evis_LS_MeV: w.evis, ...(r ? { npe_lpmt: r.npeL, npe_spmt: r.npeS, fired_lpmt: r.firedL, triggered: r.triggered, max_in_304ns: r.maxIn, wcd_fired: r.wcd.fired, reco: r.reco, hits: Array.from(r.idx, (j, k) => ({ pmt: j, type: ['dyn', 'mcp', 'spmt'][pm.type[j]], npe: r.npe[k], t_ns: +r.t[k].toFixed(2), dark: !!r.dark[k] })) } : { pending: true }) }; }) };
  download(`juno_event_${cur.ev.source}_${String(cur.k + 1).padStart(4, '0')}.json`, JSON.stringify(out), 'application/json');
});
$('#b-xcsv').addEventListener('click', () => {
  if (!cur) return;
  const lines = ['window,pmt_id,type,theta_deg,phi_deg,x_m,y_m,z_m,npe,t_ns,dark'];
  cur.res.forEach((r, w) => { if (!r) return; for (let k = 0; k < r.idx.length; k++) { const j = r.idx[k], [th, ph] = pmtAngles(pm, j); lines.push([w, j, ['dyn', 'mcp', 'spmt'][pm.type[j]], ((th * 180) / Math.PI).toFixed(2), ((ph * 180) / Math.PI).toFixed(2), pm.pos[3 * j].toFixed(3), pm.pos[3 * j + 1].toFixed(3), pm.pos[3 * j + 2].toFixed(3), r.npe[k], r.t[k].toFixed(2), r.dark[k]].join(',')); } });
  download(`juno_hits_${String(cur.k + 1).padStart(4, '0')}.csv`, lines.join('\n'), 'text/csv');
});
function runSummaries(withParticles) {
  const N = runLength(), out = [];
  for (let k = 0; k < N; k++) {
    const ev = makeEvent(k); if (!ev) continue;
    const sim = simulate(ev), s = summarize(sim);
    const row = { event: k + 1, ...eventHeader(ev), Edep_LS_MeV: s.edepLS, Evis_LS_MeV: s.evisLS, npe_prompt_expected: model.quickNpe(sim), windows: sim.windows.length, n_captures: sim.particles.filter(p => p.capture).length, inFV: s.inFV, primaries: ev.primaries.map(p => p.pdg).join(' ') };
    if (withParticles) row.particles = particleRows(sim);
    out.push(row);
  }
  return out;
}
$('#b-bcsv').addEventListener('click', () => {
  $('#x-status').textContent = 'Building run summary…';
  setTimeout(() => {
    const rows = runSummaries(false);
    const cols = ['event', 'source', 'seed', 'iev', 'category', 'nu_pdg', 'E_GeV', 'cosZ', 'phi_deg', 'cc', 'mode', 'Q2', 'W', 'vx_m', 'vy_m', 'vz_m', 'R_m', 'inFV', 'Edep_LS_MeV', 'Evis_LS_MeV', 'npe_prompt_expected', 'windows', 'n_captures', 'primaries'];
    const lines = [cols.join(',')];
    for (const r of rows) {
      const nu = r.neutrino || {}, it = r.interaction || {}, v = r.vertex_m;
      lines.push([r.event, r.source, r.seed ?? '', r.iev ?? '', r.category, nu.pdg ?? 13, (nu.E_GeV ?? r.muon?.E ?? '').toString(), (nu.cosZ ?? r.muon?.cosZ ?? '').toString(), (nu.phi_deg ?? r.muon?.phi ?? '').toString(), it.cc ?? '', it.mode ?? '', it.Q2 ?? '', it.W ?? '', ...v.map(x => x.toFixed(3)), Math.hypot(...v).toFixed(3), r.inFV, r.Edep_LS_MeV.toFixed(3), r.Evis_LS_MeV.toFixed(3), r.npe_prompt_expected.toFixed(0), r.windows, r.n_captures, `"${r.primaries}"`].join(','));
    }
    download(`juno_run_${S.source}_${rows.length}ev.csv`, lines.join('\n'), 'text/csv');
  }, 30);
});
$('#b-bjson').addEventListener('click', () => {
  $('#x-status').textContent = 'Building run file…';
  setTimeout(() => download(`juno_run_${S.source}_${runLength()}ev.json`, JSON.stringify({ schema: 'juno-event-display-run/1', generator: S.source === 'demo' ? { category: S.cat, seed: S.seed, options: { ...S.gen, flux: S.gen.flux ? S.gen.fluxName : null } } : { file: S.genie.name }, events: runSummaries(true) }), 'application/json'), 30);
});
$('#b-shot').addEventListener('click', () => {
  composer.render();
  canvas.toBlob((b) => { if (!b) return; const url = URL.createObjectURL(b); const a = Object.assign(document.createElement('a'), { href: url, download: `juno_display_${String((cur?.k ?? 0) + 1).padStart(4, '0')}.png` }); document.body.appendChild(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(url), 4000); }, 'image/png');
});

// ------------------------------------------------------------------ controls wiring
const VIS = [
  ['ls', 'Liquid scintillator', '#2b62e8'], ['av', 'Acrylic vessel', '#9fe7ff'], ['bars', 'Support bars + nodes', '#b8c4d0'], ['sss', 'SS structure + legs', '#8d99a6'],
  ['coils', 'Compensation coils', '#c4573f'], ['lpmt', '20" PMTs (17,596)', '#b48a3e'], ['spmt', '3" PMTs (25,587)', '#c8a050'], ['pool', 'Water pool + Tyvek', '#2a6eb0'],
  ['wcdpmt', 'Water-pool PMTs (3,352)', '#7d6a3e'], ['tt', 'Top Tracker', '#6f9a7e'], ['chimney', 'Chimney', '#9fb4c6'], ['house', 'Calibration house', '#8e9eb0'],
  ['cover', 'Pool cover', '#3a4654'], ['calib', 'Calibration systems', '#ffb547'], ['cavern', 'Cavern outline', '#5a6b7c'],
];
const LABELLED = new Set(det.labels.map(l => l.g));
function renderVisList() {
  $('#vis-list').innerHTML = `<div class="vis-head"><span></span><span>Detector part</span><button data-all="show" title="Show or hide all parts">Show</button><button data-all="label" title="Turn all labels on or off">Label</button></div>` +
    VIS.map(([k, t, c]) => `<div class="vis-row"><i class="sw" style="background:${c}"></i><span>${t}</span><input type="checkbox" data-vis="${k}" aria-label="Show ${t}"${det.groups[k].visible ? ' checked' : ''}><input type="checkbox" data-lab="${k}" aria-label="Label for ${t}"${LABELLED.has(k) ? '' : ' disabled'}${S.lab[k] ? ' checked' : ''}></div>`).join('');
}
renderVisList();
$('#vis-list').addEventListener('change', (e) => {
  const k = e.target.dataset.vis, L = e.target.dataset.lab;
  if (k) det.groups[k].visible = e.target.checked;
  if (L) S.lab[L] = e.target.checked;
});
$('#vis-list').addEventListener('click', (e) => {
  const b = e.target.closest('[data-all]'); if (!b) return;
  if (b.dataset.all === 'show') { const on = !VIS.every(([k]) => det.groups[k].visible); VIS.forEach(([k]) => { if (k !== 'cavern' && k !== 'calib' || !on) det.groups[k].visible = on; }); }
  else setAllLabels(!VIS.some(([k]) => S.lab[k]));
  renderVisList();
});
function setAllLabels(on) { LABELLED.forEach(k => { S.lab[k] = on; }); renderVisList(); }
$('#legend').innerHTML = Object.keys(KEYNAME).map(k => `<label><input type="checkbox" data-f="${k}" checked><i class="sw" style="background:${CSSC[k]}"></i>${KEYNAME[k]}</label>`).join('');
$('#legend').addEventListener('change', (e) => { const k = e.target.dataset.f; S.filters[k] = e.target.checked ? 1 : 0; if (anim && cur) applyAnim(anim.done ? 1 : anim.p); renderPlots(); });
$('#r-op').addEventListener('input', (e) => { S.opacity = e.target.value / 100; $('#v-op').textContent = `${e.target.value}%`; refreshMode(); });
$('#r-cut').addEventListener('input', (e) => { S.cutAngle = +e.target.value; $('#v-cut').textContent = `${e.target.value}°`; refreshMode(); });
$('#r-ph').addEventListener('input', (e) => { S.photons = +e.target.value; $('#v-ph').textContent = e.target.value; if (cur) cur.rays = buildRays(S.win); });
const bindChk = (id, key, fn) => $(id).addEventListener('change', (e) => { S[key] = e.target.checked; fn && fn(); });
bindChk('#c-front', 'front', () => { if (anim && cur) applyAnim(anim.done ? 1 : anim.p); }); bindChk('#c-orbit', 'orbit'); bindChk('#c-autowin', 'autoWin');
bindChk('#c-truth', 'truth', () => { if (anim && cur) applyAnim(anim.done ? 1 : anim.p); renderPlots(); });
bindChk('#c-dark', 'dark', () => { if (cur) { cur.res = []; requestResponse(); startAnim(S.win); } });
bindChk('#c-spmt', 'spmtMap', () => renderPlots());
$$('[data-pc]').forEach(b => b.addEventListener('click', () => { S.pmtColor = b.dataset.pc; $$('[data-pc]').forEach(x => x.setAttribute('aria-pressed', String(x === b))); if (cur) { prepareHits(); if (anim) applyAnim(anim.done ? 1 : anim.p); renderPlots(); } }));
$('#seg-rate').addEventListener('click', (e) => { const b = e.target.closest('button'); if (!b) return; S.rate = +b.dataset.rate; $$('#seg-rate button').forEach(x => x.setAttribute('aria-pressed', String(x === b))); });
$('#seg-mode').addEventListener('click', (e) => { const b = e.target.closest('button'); if (b) setMode(b.dataset.mode); });
function setMode(m) { S.mode = m; $$('#seg-mode button').forEach(x => x.setAttribute('aria-pressed', String(x.dataset.mode === m))); refreshMode(); }
$('#s-cam').addEventListener('change', (e) => setCamera(+e.target.value));
$('#plot-tabs').addEventListener('click', (e) => { const b = e.target.closest('button'); if (!b) return; S.plotTab = b.dataset.pt; $$('#plot-tabs button').forEach(x => x.setAttribute('aria-pressed', String(x === b))); renderPlots(); });
$('#s-cat').innerHTML = Object.entries(CATEGORIES).map(([k, c]) => `<option value="${k}">${c.label}</option>`).join('');
$('#s-cat').addEventListener('change', (e) => { S.cat = e.target.value; if (CATEGORIES[S.cat].kind === 'cosmic') S.gen.cosZ = Math.max(0.05, S.gen.cosZ); renderGenFields(); S.k = 0; showEvent(0); });
function setPlaying(p) { S.playing = p; $('#b-play').innerHTML = p ? '❚❚ Pause <kbd>SPACE</kbd>' : '▶ Play <kbd>SPACE</kbd>'; $('#run-state').textContent = p ? 'PLAYING' : 'PAUSED'; if (p && anim && anim.done && !Number.isFinite(anim.nextAt)) anim.nextAt = performance.now(); }
$('#b-play').addEventListener('click', () => setPlaying(!S.playing));
$('#b-next').addEventListener('click', () => nextEvent());
$('#b-replay').addEventListener('click', () => { if (cur) startAnim(0); });
$('#b-reset').addEventListener('click', () => setCamera(1));
$$('[data-col]').forEach(b => b.addEventListener('click', () => { const p = b.closest('.panel'); p.classList.toggle('collapsed'); b.textContent = p.classList.contains('collapsed') ? '+' : '–'; b.setAttribute('aria-expanded', String(!p.classList.contains('collapsed'))); }));
window.addEventListener('keydown', (e) => {
  if (e.target.closest('input, select, textarea') || e.metaKey || e.ctrlKey || e.altKey) return;
  const k = e.key.toLowerCase();
  if (k === ' ') { e.preventDefault(); setPlaying(!S.playing); }
  else if (k === 'n') nextEvent();
  else if (k === 'r') setCamera(1);
  else if (k >= '1' && k <= '6') setCamera(+k);
  else if (k === 'c') setMode(MODES[(MODES.indexOf(S.mode) + 1) % MODES.length]);
  else if (k === 'l') setAllLabels(!VIS.some(([g]) => S.lab[g]));
  else if (k === 'h') togglePanels();
  else if (k === 'escape') clearSelection();
});

// ------------------------------------------------------------------ panel layout: hide all, drag to resize
const hud = $('#hud');
const LAYOUT_DEFAULT = { left: 296, right: 318, bottom: 200 };
let layout = { ...LAYOUT_DEFAULT };
try { Object.assign(layout, JSON.parse(localStorage.getItem('juno-layout') || '{}')); } catch {}
function applyLayout(save) {
  const W = window.innerWidth, H = window.innerHeight;
  layout.left = Math.max(220, Math.min(layout.left, W * 0.4));
  layout.right = Math.max(240, Math.min(layout.right, W * 0.42));
  layout.bottom = Math.max(130, Math.min(layout.bottom, H * 0.55));
  hud.style.setProperty('--left-w', `${layout.left}px`); hud.style.setProperty('--right-w', `${layout.right}px`); hud.style.setProperty('--bottom-h', `${layout.bottom}px`);
  stage.style.setProperty('--bottom-h', `${layout.bottom}px`);
  if (save) try { localStorage.setItem('juno-layout', JSON.stringify(layout)); } catch {}
  resize();
}
$$('.gutter').forEach(g => g.addEventListener('pointerdown', (e) => {
  e.preventDefault(); g.setPointerCapture(e.pointerId); g.classList.add('drag');
  const which = g.dataset.g, x0 = e.clientX, y0 = e.clientY, start = { ...layout };
  const move = (ev) => {
    if (which === 'left') layout.left = start.left + (ev.clientX - x0);
    else if (which === 'right') layout.right = start.right - (ev.clientX - x0);
    else layout.bottom = start.bottom - (ev.clientY - y0);
    applyLayout(false);
  };
  const up = () => { g.classList.remove('drag'); g.removeEventListener('pointermove', move); g.removeEventListener('pointerup', up); applyLayout(true); };
  g.addEventListener('pointermove', move); g.addEventListener('pointerup', up);
}));
function togglePanels() {
  const off = hud.classList.toggle('panels-off');
  $('#b-panels').textContent = off ? '◨ Show panels' : '◧ Hide panels';
  stage.style.setProperty('--bottom-h', off ? '-14px' : `${layout.bottom}px`);
  resize();
}
$('#b-panels').addEventListener('click', togglePanels);
$('#b-layout').addEventListener('click', () => { layout = { ...LAYOUT_DEFAULT }; applyLayout(true); });

// ------------------------------------------------------------------ boot
refreshMode();
applyLayout(false);
renderGenFields();
setCamera(1, true);
resize();
showEvent(0);
requestAnimationFrame(loop);
if (window.innerWidth < 1100) { $('#p-det').classList.add('collapsed'); $('#p-det [data-col]').textContent = '+'; }
