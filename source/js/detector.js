// JUNO geometry and detector constants.
// Coordinates: origin at the Central Detector (CD) centre, z vertical (up), metres.
// [P] = value from the papers (arXiv:2511.14593, arXiv:2511.14590).
// [D] = approximate display value chosen here; not engineering geometry.

export const GEO = {
  R_LS: 17.7,            // [P] AV inner radius (35.4 m inner diameter)
  AV_THICK: 0.12,        // [P] acrylic panel thickness
  N_AV_PANELS: 263,      // [P]
  N_AV_BARS: 590,        // [P] connecting bars AV -> SS structure
  R_PMT: 19.5,           // [D] LPMT photocathode sphere (gives ~75% coverage with 17,596 x 20")
  R_SSS: 20.05,          // [P] SS structure inner diameter 40.1 m
  R_COIL: 20.6,          // [D] compensation coils outside the SSS
  WP_R: 21.75,           // [P] water pool diameter 43.5 m
  WP_ZMIN: -21.6,        // [D] pool floor (pool height 44 m [P]; vertical placement approximate)
  WP_ZMAX: 22.4,         // [D] pool top
  WATER_TOP: 21.9,       // [D] water level (43.5 m filled [P])
  CHIMNEY_R: 0.45,       // [D]
  TT_Z: [23.6, 25.1, 26.6], // [D] three TT layers, ~1.5 m spacing [P]
  TT_WALL: 6.8,          // [P] wall size 6.8 x 6.8 m
  TT_NX: 7, TT_NY: 3,    // [P] 7 x 3 walls per layer, 63 walls total
  FV_R: 16.5,            // [P] fiducial radius used in the papers
  FV_Z: 15.5,            // [P] |z| cut
};

export const PMT = {
  N_DYN: 4939, N_MCP: 12657,   // [P] 20-inch CD PMTs
  N_LPMT: 17596,               // [P]
  N_SPMT: 25587,               // [P] 3-inch PMTs
  R_LPMT: 0.254, R_SPMT: 0.038,
  N_WCD_SSS: 2404,             // [P] 20-inch MCP on SSS, facing outward
  WCD_RING20: [20.66, 18.98],  // [P] two rings, 348 x 20-inch MCP in total
  N_WCD_RING20: 348,
  WCD_RING8: [17.05, 16.05],   // [P] two rings, 600 x 8-inch PMTs in total
  N_WCD_RING8: 600,
  TTS: [1.93, 7.53, 1.6],      // [P] ns; dynode, MCP; [D] SPMT
  DCR: [20.6e3, 22.7e3, 0.4e3],// [P] Hz dark rate; SPMT [D]
  GATE_NS: 1008,               // [P] readout window
  TRIG_N: 350, TRIG_WIN: 304,  // [P] multiplicity trigger
};

export const OPT = {
  LY_CENTRE: 1785,     // [P] p.e./MeV at CD centre (2.223 MeV n-H gamma)
  ATT_LS: 20.6,        // [P] m at 430 nm
  ATT_W: 60,           // [P] water attenuation (60-70 m measured)
  N_EFF: 1.60,         // [P] effective refractive index used in JVertex
  N_LS: 1.50,          // [D] for Cherenkov angle
  TAU_FAST: 4.6,       // [D] ns scintillation fast component
  CHER_FRAC: 0.06,     // [D] fraction of light emitted as directional Cherenkov
  C: 0.299792458,      // m/ns
};

// ---------- PMT layout (deterministic, quasi-uniform; NOT the engineering map) ----------
function fib(n, i, phase = 0) {
  const ga = Math.PI * (3 - Math.sqrt(5));
  const z = 1 - (2 * (i + 0.5)) / n;
  const r = Math.sqrt(1 - z * z);
  const p = i * ga + phase;
  return [r * Math.cos(p), r * Math.sin(p), z];
}

function sphereLayout(target, capCos, phase) {
  // fibonacci points with a hole around the top chimney (z > capCos), exact count
  let m = target;
  for (let it = 0; it < 60; it++) {
    let k = 0;
    for (let i = 0; i < m; i++) if (fib(m, i, phase)[2] < capCos) k++;
    if (k === target) break;
    m += target - k;
  }
  const out = [];
  for (let i = 0; i < m && out.length < target; i++) {
    const v = fib(m, i, phase);
    if (v[2] < capCos) out.push(v);
  }
  return out;
}

function hash01(i) { // deterministic pseudo-random in [0,1)
  let x = (i + 0x9e3779b9) | 0;
  x = Math.imul(x ^ (x >>> 16), 0x85ebca6b);
  x = Math.imul(x ^ (x >>> 13), 0xc2b2ae35);
  return ((x ^ (x >>> 16)) >>> 0) / 4294967296;
}

export function buildPmts() {
  const capCos = Math.cos(2.2 / GEO.R_PMT); // ~2.2 m arc hole at the chimney
  const L = sphereLayout(PMT.N_LPMT, capCos, 0);
  const S = sphereLayout(PMT.N_SPMT, Math.cos(1.8 / GEO.R_PMT), 1.3);
  const N = L.length + S.length;
  const pos = new Float32Array(N * 3), dir = new Float32Array(N * 3);
  const type = new Uint8Array(N), area = new Float32Array(N);
  // choose exactly N_DYN dynode PMTs, spread over the sphere
  const order = L.map((_, i) => i).sort((a, b) => hash01(a) - hash01(b));
  const isDyn = new Uint8Array(L.length);
  for (let k = 0; k < PMT.N_DYN; k++) isDyn[order[k]] = 1;
  const put = (j, v, R, t, rad) => {
    pos[3 * j] = v[0] * R; pos[3 * j + 1] = v[1] * R; pos[3 * j + 2] = v[2] * R;
    dir[3 * j] = -v[0]; dir[3 * j + 1] = -v[1]; dir[3 * j + 2] = -v[2];
    type[j] = t; area[j] = Math.PI * rad * rad;
  };
  L.forEach((v, i) => put(i, v, GEO.R_PMT, isDyn[i] ? 0 : 1, PMT.R_LPMT));
  S.forEach((v, i) => put(L.length + i, v, GEO.R_PMT + 0.05, 2, PMT.R_SPMT));
  return { N, nL: L.length, nS: S.length, pos, dir, type, area };
}

export function buildWcdPmts() {
  const list = [];
  const sss = sphereLayout(PMT.N_WCD_SSS, Math.cos(2.5 / GEO.R_SSS), 0.7);
  sss.forEach(v => list.push({ p: v.map(c => c * (GEO.R_SSS + 0.5)), d: v, kind: 0 }));
  const ring = (z, n, kind) => {
    for (let i = 0; i < n; i++) {
      const a = (2 * Math.PI * (i + 0.5 * (kind === 2))) / n;
      const c = Math.cos(a), s = Math.sin(a), r = GEO.WP_R - 0.35;
      list.push({ p: [r * c, r * s, z], d: [-c, -s, 0], kind });
    }
  };
  PMT.WCD_RING20.forEach(z => ring(z, PMT.N_WCD_RING20 / 2, 1));
  PMT.WCD_RING8.forEach(z => ring(z, PMT.N_WCD_RING8 / 2, 2));
  const N = list.length;
  const pos = new Float32Array(N * 3), dir = new Float32Array(N * 3), kind = new Uint8Array(N);
  list.forEach((o, i) => { pos.set(o.p, 3 * i); dir.set(o.d, 3 * i); kind[i] = o.kind; });
  return { N, pos, dir, kind };
}

export function ttWalls() {
  const walls = [];
  const W = GEO.TT_WALL;
  GEO.TT_Z.forEach((z, layer) => {
    for (let ix = 0; ix < GEO.TT_NX; ix++)
      for (let iy = 0; iy < GEO.TT_NY; iy++)
        walls.push({ layer, ix, iy, x: (ix - (GEO.TT_NX - 1) / 2) * W, y: (iy - (GEO.TT_NY - 1) / 2) * W, z });
  });
  return walls;
}

export function region(x, y, z) {
  const r2 = x * x + y * y + z * z;
  if (r2 < GEO.R_LS * GEO.R_LS) return 'LS';
  if (r2 < GEO.R_SSS * GEO.R_SSS) return 'BUF';
  if (x * x + y * y < GEO.WP_R * GEO.WP_R && z > GEO.WP_ZMIN && z < GEO.WATER_TOP) return 'WP';
  if (z > GEO.WATER_TOP && z < 30 && Math.abs(x) < 28 && Math.abs(y) < 28) return 'AIR';
  return 'ROCK';
}
