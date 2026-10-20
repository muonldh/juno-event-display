# Single-file build: everything in one index.html (three.js and JSROOT still load from the CDN).
import subprocess, json, pathlib
E = './node_modules/.bin/esbuild'
w = subprocess.run([E, 'site/js/worker.js', '--bundle', '--format=iife', '--minify'], capture_output=True, text=True, check=True).stdout
a = subprocess.run([E, 'site/js/app.js', '--bundle', '--format=esm', '--minify', '--external:three', '--external:three/addons/*'], capture_output=True, text=True, check=True).stdout
html = pathlib.Path('dist/index.html').read_text()
tag = '<script type="module" src="./js/app.js"></script>'
assert tag in html
inline = '<script>window.__JUNO_WORKER_SRC = ' + json.dumps(w).replace('</', '<\\/') + ';</script>\n<script type="module">\n' + a.replace('</script', '<\\/script') + '\n</script>'
out = pathlib.Path('single'); out.mkdir(exist_ok=True)
(out / 'index.html').write_text(html.replace(tag, inline))
print('single index.html', len((out / 'index.html').read_text()) // 1024, 'KB')
