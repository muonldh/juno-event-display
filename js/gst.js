// Read GENIE gst trees directly in the browser with JSROOT (no server, no ROOT install).
const JSROOT_BASE = 'https://cdn.jsdelivr.net/npm/jsroot@7.10.3/modules/';

const SCALAR = ['iev','neu','fspl','tgt','cc','nc','qel','res','dis','coh','mec',
  'Ev','pxv','pyv','pzv','El','pxl','pyl','pzl','Q2','W','x','y','vtxx','vtxy','vtxz','vtxt'];
const ARRAYS = ['pdgf','Ef','pxf','pyf','pzf'];

export async function readGst(source, { first = 0, max = 500, base = JSROOT_BASE } = {}) {
  const io = await import(base + 'io.mjs');
  const tr = await import(base + 'tree.mjs');
  const file = await io.openFile(source);
  const tree = await file.readObject('gst');
  if (!tree) throw new Error('No "gst" tree found in this file. Convert GHEP with: gntpc -f gst');
  const have = new Set(tree.fBranches.arr.map(b => b.fName));
  const missing = ['neu', 'Ev', 'vtxx', 'nf', 'pdgf'].filter(b => !have.has(b));
  if (missing.length) throw new Error('gst tree is missing branches: ' + missing.join(', '));
  const sel = new tr.TSelector();
  const use = [...SCALAR, 'nf', ...ARRAYS].filter(b => have.has(b));
  use.forEach(b => sel.addBranch(b));
  const rows = [];
  sel.Process = function () {
    const o = this.tgtobj, r = {};
    for (const b of use) {
      const v = o[b];
      r[b] = (v && typeof v === 'object' && 'length' in v) ? Array.from(v) : v;
    }
    rows.push(r);
  };
  const total = tree.fEntries;
  const n = Math.max(0, Math.min(max, total - first));
  await tr.treeProcess(tree, sel, { firstentry: first, numentries: n });
  return { rows, total, name: typeof source === 'string' ? source : source.name };
}
