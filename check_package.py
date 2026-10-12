from pathlib import Path
import json,hashlib
root=Path(__file__).resolve().parent
m=json.loads((root/'MANIFEST.json').read_text())
for f in m['files']:
 p=root/f['path']
 if not p.is_file() or p.stat().st_size!=f['size_bytes'] or hashlib.sha256(p.read_bytes()).hexdigest()!=f['sha256']:
  raise ValueError('manifest mismatch: '+f['path'])
print('MANIFEST VERIFIED',len(m['files']))
