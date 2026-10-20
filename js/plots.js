// 2D canvas plots: PMT Mollweide map, first-hit time histogram, side/top projections.
import { GEO } from './detector.js';

const css = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
export function setupCanvas(cv) {
  const r = cv.getBoundingClientRect(), dpr = Math.min(2, window.devicePixelRatio || 1);
  const w = Math.max(10, Math.round(r.width)), h = Math.max(10, Math.round(r.height));
  if (cv.width !== w * dpr || cv.height !== h * dpr) { cv.width = w * dpr; cv.height = h * dpr; }
  const g = cv.getContext('2d'); g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, w, h);
  return { g, w, h };
}

// colour maps
const STOPS_Q = [[0, [27, 42, 107]], [0.35, [47, 140, 255]], [0.6, [47, 211, 255]], [0.82, [255, 224, 102]], [1, [255, 90, 60]]];
const STOPS_T = [[0, [90, 40, 160]], [0.35, [60, 110, 220]], [0.6, [32, 190, 190]], [0.85, [150, 225, 90]], [1, [250, 240, 80]]];
function ramp(stops, x) {
  x = Math.max(0, Math.min(1, x));
  for (let i = 1; i < stops.length; i++) if (x <= stops[i][0]) {
    const [a, ca] = stops[i - 1], [b, cb] = stops[i], f = (x - a) / (b - a);
    return [ca[0] + f * (cb[0] - ca[0]), ca[1] + f * (cb[1] - ca[1]), ca[2] + f * (cb[2] - ca[2])];
  }
  return stops[stops.length - 1][1];
}
export const cmapQ = (x) => ramp(STOPS_Q, x), cmapT = (x) => ramp(STOPS_T, x);
export function hitScale(res, mode) {
  let qmax = 1, tmin = Infinity;
  for (let i = 0; i < res.idx.length; i++) { if (res.npe[i] > qmax) qmax = res.npe[i]; if (!res.dark[i] && res.t[i] < tmin) tmin = res.t[i]; }
  if (!Number.isFinite(tmin)) tmin = 0;
  return mode === 'time' ? (i) => cmapT((res.t[i] - tmin) / 60) : (i) => cmapQ(Math.log(1 + res.npe[i]) / Math.log(1 + qmax));
}

// Mollweide projection of (theta from +z, phi)
function mollweide(theta, phi) {
  const lat = Math.PI / 2 - theta; let t = lat;
  for (let k = 0; k < 8; k++) { const f = 2 * t + Math.sin(2 * t) - Math.PI * Math.sin(lat), d = 2 + 2 * Math.cos(2 * t); if (Math.abs(d) < 1e-9) break; t -= f / d; }
  let lon = phi; if (lon > Math.PI) lon -= 2 * Math.PI;
  return [(2 * Math.SQRT2 / Math.PI) * lon * Math.cos(t), Math.SQRT2 * Math.sin(t)];
}
export function pmtAngles(pm, j) {
  const x = pm.pos[3 * j], y = pm.pos[3 * j + 1], z = pm.pos[3 * j + 2], r = Math.hypot(x, y, z);
  return [Math.acos(z / r), (Math.atan2(y, x) + 2 * Math.PI) % (2 * Math.PI)];
}
export function sphereAngles(p) { const r = Math.hypot(...p) || 1; return [Math.acos(p[2] / r), (Math.atan2(p[1], p[0]) + 2 * Math.PI) % (2 * Math.PI)]; }

export function drawPmtMap(cv, pm, res, opts) {
  const { g, w, h } = setupCanvas(cv);
  if (w < 60 || h < 40) return { P: () => [0, 0], hits: [] }; // panel hidden or collapsed
  const sx = (w - 16) / (4 * Math.SQRT2), sy = (h - 14) / (2 * Math.SQRT2), s = Math.min(sx, sy), cx = w / 2, cy = h / 2;
  const P = (th, ph) => { const [X, Y] = mollweide(th, ph); return [cx + X * s, cy - Y * s]; };
  g.fillStyle = css('--plot-bg'); g.strokeStyle = css('--line');
  g.beginPath(); g.ellipse(cx, cy, 2 * Math.SQRT2 * s, Math.SQRT2 * s, 0, 0, 2 * Math.PI); g.fill(); g.stroke();
  g.strokeStyle = css('--grid'); g.lineWidth = 0.6;
  for (const lat of [-60, -30, 0, 30, 60]) { g.beginPath(); for (let ph = 0; ph <= 360; ph += 6) { const [x, y] = P((90 - lat) * Math.PI / 180, ((ph - 180) * Math.PI) / 180); ph ? g.lineTo(x, y) : g.moveTo(x, y); } g.stroke(); }
  for (const lon of [-120, -60, 0, 60, 120]) { g.beginPath(); for (let th = 0; th <= 180; th += 6) { const [x, y] = P((th * Math.PI) / 180, (lon * Math.PI) / 180 + (lon < 0 ? 2 * Math.PI : 0)); th ? g.lineTo(x, y) : g.moveTo(x, y); } g.stroke(); }
  const out = { P, hits: [] };
  if (!res) { g.fillStyle = css('--fg2'); g.font = `11px ${css('--mono')}`; g.textAlign = 'center'; g.fillText(opts.waiting ? 'computing PMT response…' : 'no hits', cx, cy); return out; }
  const col = hitScale(res, opts.mode);
  const nshow = opts.reveal ?? res.idx.length;
  const order = opts.order || null;
  const pxs = Math.max(1.4, s / 45);
  for (let k = 0; k < nshow; k++) {
    const i = order ? order[k] : k, j = res.idx[i];
    if (j >= pm.nL && !opts.spmt) continue;
    const [th, ph] = pmtAngles(pm, j), [x, y] = P(th, ph), c = col(i);
    g.fillStyle = `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0})`;
    const r = j < pm.nL ? pxs : pxs * 0.5;
    g.fillRect(x - r / 2, y - r / 2, r, r);
    out.hits.push([x, y, i]);
  }
  const mark = (p, shape, color) => {
    const [th, ph] = sphereAngles(p), [x, y] = P(th, ph);
    g.strokeStyle = color; g.lineWidth = 1.6; g.beginPath();
    if (shape === 'star') { for (let k = 0; k < 10; k++) { const a = -Math.PI / 2 + (k * Math.PI) / 5, rr = k % 2 ? 3 : 7; g.lineTo(x + rr * Math.cos(a), y + rr * Math.sin(a)); } g.closePath(); }
    else g.arc(x, y, 6, 0, 2 * Math.PI);
    g.stroke();
  };
  if (opts.truthTrack) { mark(opts.truthTrack[0], 'star', '#ff4d5e'); mark(opts.truthTrack[1], 'star', '#ffd34d'); }
  if (opts.recoTrack) { mark(opts.recoTrack.entry, 'circ', '#ff4d5e'); mark(opts.recoTrack.exit, 'circ', '#ffd34d'); }
  if (opts.selPmt != null) { const [th, ph] = pmtAngles(pm, opts.selPmt), [x, y] = P(th, ph); g.strokeStyle = '#ffffff'; g.lineWidth = 1.4; g.strokeRect(x - 5, y - 5, 10, 10); }
  return out;
}

export function drawLegend(cv, mode, res) {
  const { g, w, h } = setupCanvas(cv);
  for (let x = 0; x < w; x++) { const c = mode === 'time' ? cmapT(x / w) : cmapQ(x / w); g.fillStyle = `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0})`; g.fillRect(x, 0, 1, h); }
}

function axes(g, w, h, m, xr, yr, xl, yl, log = false) {
  g.strokeStyle = css('--grid'); g.fillStyle = css('--fg2'); g.lineWidth = 1; g.font = `10px ${css('--mono')}`;
  const X = (v) => m.l + ((v - xr[0]) / (xr[1] - xr[0])) * (w - m.l - m.r);
  const Y = (v) => { const f = log ? (Math.log10(Math.max(v, yr[0])) - Math.log10(yr[0])) / (Math.log10(yr[1]) - Math.log10(yr[0])) : (v - yr[0]) / (yr[1] - yr[0]); return h - m.b - f * (h - m.t - m.b); };
  const nice = (a, b, n) => { const st = Math.pow(10, Math.floor(Math.log10((b - a) / n))); const k = [1, 2, 5, 10].find(k => (b - a) / (k * st) <= n) * st; const out = []; for (let v = Math.ceil(a / k) * k; v <= b + 1e-9; v += k) out.push(+v.toFixed(6)); return out; };
  g.textAlign = 'center'; g.textBaseline = 'top';
  for (const v of nice(xr[0], xr[1], 5)) { const x = X(v); g.beginPath(); g.moveTo(x, m.t); g.lineTo(x, h - m.b); g.stroke(); g.fillText(String(v), x, h - m.b + 3); }
  g.textAlign = 'right'; g.textBaseline = 'middle';
  const yt = log ? [1, 10, 100, 1e3, 1e4].filter(v => v >= yr[0] && v <= yr[1]) : nice(yr[0], yr[1], 4);
  for (const v of yt) { const y = Y(v); g.beginPath(); g.moveTo(m.l, y); g.lineTo(w - m.r, y); g.stroke(); g.fillText(v >= 1000 ? `${v / 1000}k` : String(v), m.l - 4, y); }
  g.textAlign = 'center'; g.textBaseline = 'bottom'; g.fillText(xl, (m.l + w - m.r) / 2, h - 1);
  g.save(); g.translate(10, (m.t + h - m.b) / 2); g.rotate(-Math.PI / 2); g.textBaseline = 'middle'; g.fillText(yl, 0, 0); g.restore();
  return { X, Y };
}

export function drawTimeHist(cv, pm, res, opts = {}) {
  const { g, w, h } = setupCanvas(cv);
  if (w < 60 || h < 40) return;
  if (!res || !res.idx.length) { g.fillStyle = css('--fg2'); g.font = `11px ${css('--mono')}`; g.textAlign = 'center'; g.fillText(opts.waiting ? 'computing…' : 'no hits', w / 2, h / 2); return; }
  let tmin = Infinity; for (let i = 0; i < res.idx.length; i++) if (!res.dark[i] && res.idx[i] < pm.nL) tmin = Math.min(tmin, res.t[i]);
  if (!Number.isFinite(tmin)) tmin = 0;
  const t0 = Math.floor(tmin / 10) * 10 - 20, span = 220, bw = 2, nb = span / bw;
  const H = [new Float64Array(nb), new Float64Array(nb)];
  for (let i = 0; i < res.idx.length; i++) { const j = res.idx[i]; if (j >= pm.nL) continue; const b = Math.floor((res.t[i] - t0) / bw); if (b >= 0 && b < nb) H[pm.type[j]][b]++; }
  const ymax = Math.max(10, ...H[0], ...H[1]) * 1.4;
  const m = { l: 34, r: 6, t: 6, b: 26 };
  const { X, Y } = axes(g, w, h, m, [t0, t0 + span], [0.8, ymax], 'first-hit time (ns)', 'PMTs / 2 ns', true);
  const draw = (arr, color) => { g.strokeStyle = color; g.lineWidth = 1.4; g.beginPath(); for (let b = 0; b < nb; b++) { const y = Y(Math.max(0.8, arr[b])); g.lineTo(X(t0 + b * bw), y); g.lineTo(X(t0 + (b + 1) * bw), y); } g.stroke(); };
  draw(H[1], '#ff8a5c'); draw(H[0], '#6e8bff');
  g.font = `10px ${css('--mono')}`; g.textAlign = 'left'; g.textBaseline = 'top';
  g.fillStyle = '#6e8bff'; g.fillText('dynode (TTS 1.9 ns)', m.l + 6, m.t + 2);
  g.fillStyle = '#ff8a5c'; g.fillText('MCP (TTS 7.5 ns)', m.l + 6, m.t + 14);
  if (opts.selT != null) { g.strokeStyle = '#fff'; g.beginPath(); g.moveTo(X(opts.selT), m.t); g.lineTo(X(opts.selT), h - m.b); g.stroke(); }
}

export function drawProjection(cv, sim, opts) {
  const { g, w, h } = setupCanvas(cv);
  if (w < 60 || h < 40) return;
  const side = opts.view === 'side'; // side: x-z ; top: x-y
  const m = { l: 30, r: 6, t: 6, b: 24 };
  const span = side ? [-24, 24] : [-24, 24];
  const sc = Math.min((w - m.l - m.r) / 48, (h - m.t - m.b) / 48), ox = m.l + (w - m.l - m.r) / 2, oy = m.t + (h - m.t - m.b) / 2;
  const P = (p) => [ox + p[0] * sc, oy - (side ? p[2] : p[1]) * sc];
  g.strokeStyle = css('--grid'); g.lineWidth = 1; g.font = `10px ${css('--mono')}`; g.fillStyle = css('--fg2');
  g.textAlign = 'center'; g.textBaseline = 'top';
  for (const v of [-20, -10, 0, 10, 20]) { const [x] = P([v, 0, 0]); g.beginPath(); g.moveTo(x, m.t); g.lineTo(x, h - m.b); g.stroke(); g.fillText(v, x, h - m.b + 2); }
  g.textAlign = 'right'; g.textBaseline = 'middle';
  for (const v of [-20, -10, 0, 10, 20]) { const [, y] = P([0, v, v]); g.beginPath(); g.moveTo(m.l, y); g.lineTo(w - m.r, y); g.stroke(); g.fillText(v, m.l - 3, y); }
  g.textAlign = 'center'; g.textBaseline = 'bottom'; g.fillText(side ? 'x (m)        side view, z vertical' : 'x (m)        top view, y vertical', ox, h - 0);
  g.strokeStyle = css('--line'); g.lineWidth = 1.2;
  g.beginPath(); g.arc(ox, oy, GEO.R_LS * sc, 0, 2 * Math.PI); g.stroke();
  g.setLineDash([2, 3]); g.beginPath(); g.arc(ox, oy, GEO.R_PMT * sc, 0, 2 * Math.PI); g.stroke(); g.setLineDash([]);
  if (side) { const [x0, y0] = P([-GEO.WP_R, 0, GEO.WP_ZMAX]), [x1, y1] = P([GEO.WP_R, 0, GEO.WP_ZMIN]); g.strokeRect(x0, y0, x1 - x0, y1 - y0); }
  else { g.beginPath(); g.arc(ox, oy, GEO.WP_R * sc, 0, 2 * Math.PI); g.stroke(); }
  if (!sim) return;
  const colors = opts.colors;
  for (const tr of sim.tracks) {
    if (tr.win !== opts.win && tr.style !== 'invisible') continue;
    const key = opts.key(tr.cls); if (!opts.filter(key)) continue;
    g.strokeStyle = colors[key]; g.globalAlpha = opts.sel != null && tr.pid !== opts.sel ? 0.25 : tr.style === 'branch' ? 0.5 : 1;
    g.lineWidth = tr.style === 'branch' ? 0.8 : key === 'mu' ? 2 : 1.4;
    g.setLineDash(tr.style === 'solid' || tr.style === 'branch' ? [] : [3, 3]);
    g.beginPath(); tr.pts.forEach((p, i) => { const [x, y] = P(p); i ? g.lineTo(x, y) : g.moveTo(x, y); }); g.stroke();
  }
  g.setLineDash([]); g.globalAlpha = 1;
  for (const d of sim.deps) { if (d.win !== opts.win) continue; const [x, y] = P([d.x, d.y, d.z]); const r = Math.min(4, 0.8 + Math.log10(1 + d.evis)); g.fillStyle = colors[opts.key(d.cls)]; g.globalAlpha = 0.55; g.beginPath(); g.arc(x, y, r, 0, 2 * Math.PI); g.fill(); }
  g.globalAlpha = 1;
  if (!sim.cosmic) { const [x, y] = P(sim.vertex); g.strokeStyle = '#ffffff'; g.lineWidth = 1.5; g.beginPath(); g.moveTo(x - 4, y - 4); g.lineTo(x + 4, y + 4); g.moveTo(x + 4, y - 4); g.lineTo(x - 4, y + 4); g.stroke(); }
  if (opts.reco) { const [x, y] = P(opts.reco.vertex); g.strokeStyle = '#6e8bff'; g.lineWidth = 1.5; g.beginPath(); g.arc(x, y, 5, 0, 2 * Math.PI); g.stroke(); }
}
