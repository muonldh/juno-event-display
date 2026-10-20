// Web Worker: computes the PMT response off the main thread so the 3D view stays smooth.
import { buildPmts, buildWcdPmts } from './detector.js';
import { makeResponseModel, makeRng } from './physics.js';

const pm = buildPmts(), wpm = buildWcdPmts();
const model = makeResponseModel(pm, wpm);

self.onmessage = (e) => {
  const { job, sim, seed, dark } = e.data;
  for (let w = 0; w < sim.windows.length; w++) {
    const r = model.windowResponse(sim, w, makeRng((seed * 131 + w * 977 + 7) >>> 0), { dark });
    self.postMessage({ job, w, r }, [r.idx.buffer, r.npe.buffer, r.t.buffer, r.dark.buffer]);
  }
  self.postMessage({ job, done: true });
};
self.postMessage({ ready: true });
