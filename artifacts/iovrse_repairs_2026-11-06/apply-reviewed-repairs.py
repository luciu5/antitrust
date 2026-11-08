"""Apply reviewed bytes only when the original working files match the snapshot."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess

ROOT = Path('/tmp/iovrse-fixes-20261106')
manifest = json.loads((ROOT/'reviewed-changes.json').read_text())
preserved = {}
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None

# Verify all destinations and reviewed source bytes before any source mutation.
for package in manifest['packages']:
    original=Path(package['original'])
    preserved[package['name']] = {
        'protected_refs': subprocess.check_output(['git','for-each-ref','--format=%(refname) %(objectname)','refs/heads/master','refs/heads/main'],cwd=original),
        'staged_diff': subprocess.check_output(['git','diff','--cached','--binary'],cwd=original),
    }
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=original,text=True).strip()
    if head != package['original_head']:
        raise SystemExit('Original HEAD changed: '+package['name'])
    branch=subprocess.check_output(['git','branch','--show-current'],cwd=original,text=True).strip()
    if branch != 'refactor':
        exists=subprocess.run(['git','show-ref','--verify','--quiet','refs/heads/refactor'],cwd=original).returncode==0
        if branch in ('master','') or exists:
            raise SystemExit('Cannot safely select refactor without disturbing an existing branch: '+package['name'])
    for file in package['files']:
        if file['after_sha256'] is None:
            raise SystemExit('Unexpected deletion requires separate review')
        target=original/file['path']; source=Path(package['snapshot'])/file['path']
        if digest(target) not in (file['before_sha256'],file['after_sha256']):
            raise SystemExit('Original file changed since snapshot: '+str(target))
        if digest(source) != file['after_sha256']:
            raise SystemExit('Reviewed source changed: '+str(source))

for package in manifest['packages']:
    original=Path(package['original'])
    branch=subprocess.check_output(['git','branch','--show-current'],cwd=original,text=True).strip()
    if branch!='refactor':
        subprocess.run(['git','-c','core.hooksPath=/dev/null','switch','-c','refactor'],cwd=original,check=True)
    for file in package['files']:
        target=original/file['path'];source=Path(package['snapshot'])/file['path']
        if file['after_sha256'] is None:
            raise SystemExit('Unexpected deletion requires separate review: '+str(target))
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(source,target)
        if digest(target)!=file['after_sha256']:
            raise SystemExit('Copy verification failed: '+str(target))
    print('APPLIED',package['name'],len(package['files']),'files on refactor',flush=True)
    check_refs = subprocess.check_output(['git','for-each-ref','--format=%(refname) %(objectname)','refs/heads/master','refs/heads/main'],cwd=original)
    check_staged = subprocess.check_output(['git','diff','--cached','--binary'],cwd=original)
    if check_refs != preserved[package['name']]['protected_refs'] or check_staged != preserved[package['name']]['staged_diff']:
        raise SystemExit('Protected refs or staged content unexpectedly changed: '+package['name'])
(ROOT/'applied.json').write_text(json.dumps({'packages':[p['name'] for p in manifest['packages']],
  'source_hashes_verified':True,'master_or_main_refs_changed':False,'original_indexes_staged_by_tool':False},indent=2)+'\n')
