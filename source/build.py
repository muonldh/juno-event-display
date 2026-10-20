# Builds dist/ for GitHub Pages: wraps index.html in a full HTML document, copies js/.
import shutil, pathlib, re
src = pathlib.Path('site'); dst = pathlib.Path('dist')
if dst.exists(): shutil.rmtree(dst)
shutil.copytree(src / 'js', dst / 'js')
body = (src / 'index.html').read_text()
title = re.search(r'<title>.*?</title>', body).group(0)
head_end = body.index('</style>') + len('</style>')
head, rest = body[:head_end], body[head_end:]
head = head.replace(title, '').replace('<meta name="claude:viewport-fit-sides" content="cover">\n', '')
html = f'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
{title}
{head}
</head>
<body>
{rest.strip()}
</body>
</html>
'''
(dst / 'index.html').write_text(html)
(dst / '.nojekyll').write_text('')
print('built', sum(1 for _ in dst.rglob('*')), 'files')
