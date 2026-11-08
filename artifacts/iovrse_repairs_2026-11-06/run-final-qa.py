import concurrent.futures
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path('/tmp/iovrse-fixes-20261106')
NAMES = ['antitrust', 'coordination', 'trade', 'vertical', 'financial', 'competitiontoolbox', 'antitrustBayes']
LIB = Path(os.environ.get('IOVRSE_FINAL_LIB', str(ROOT / 'final-library')))
OUT = ROOT / 'final-evidence'
LIB.mkdir(parents=True, exist_ok=True)
OUT.mkdir(exist_ok=True)
env = os.environ.copy()
env.update(R_LIBS=str(LIB), ANTITRUST_TEST_TIER='extended', TRADE_TEST_TIER='extended',
           FINANCIAL_EXTENDED_TESTS='1', ANTITRUSTBAYES_RUN_STAN_TESTS='false',
           COMPETITIONTOOLBOX_TEST_TIER='EXTENDED',
           _R_CHECK_FORCE_SUGGESTS_='false')

def execute(args, cwd, log):
    with log.open('w') as f:
        result = subprocess.run(args, cwd=cwd, env=env, stdout=f, stderr=subprocess.STDOUT)
    return result.returncode

if sys.argv[1:] == ['install']:
    cubature = ROOT / 'financial-lib' / 'cubature'
    if not (LIB / 'cubature').exists():
        shutil.copytree(cubature, LIB / 'cubature', symlinks=True)
    for name in NAMES:
        out = OUT / name
        out.mkdir(exist_ok=True)
        rc = execute(['R','CMD','INSTALL','--library='+str(LIB),'--no-docs','--no-html','--no-demo',
                      str(ROOT/'sources'/name)], ROOT, out/'install.log')
        print('INSTALL', name, rc, flush=True)
        if rc:
            raise SystemExit(rc)
elif sys.argv[1:] == ['tests']:
    def test(name):
        out = OUT / name
        rc = execute(['Rscript','--vanilla',str(ROOT/'final-tests.R'),name], ROOT, out/'tests.log')
        return {'package':name,'exit':rc}
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        results=[]
        for f in concurrent.futures.as_completed([pool.submit(test,n) for n in NAMES]):
            result=f.result();results.append(result);print('TEST',result,flush=True)
    (OUT/'test-exits.json').write_text(json.dumps(results,indent=2)+'\n')
    raise SystemExit(any(r['exit'] for r in results))
elif sys.argv[1:] == ['checks']:
    def check(name):
        out = OUT / name
        build = execute(['R','CMD','build',str(ROOT/'sources'/name)],out,out/'build.log')
        if build:
            return {'package':name,'build_exit':build,'check_exit':None}
        archives=list(out.glob(name+'_*.tar.gz'))
        if len(archives)!=1:
            raise RuntimeError('Expected exactly one archive: '+str(out))
        rc=execute(['R','CMD','check','--no-manual',str(archives[0])],out,out/'check.log')
        return {'package':name,'build_exit':0,'check_exit':rc}
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        results=[]
        for f in concurrent.futures.as_completed([pool.submit(check,n) for n in NAMES]):
            result=f.result();results.append(result);print('CHECK',result,flush=True)
    (OUT/'check-exits.json').write_text(json.dumps(results,indent=2)+'\n')
    raise SystemExit(any(r['build_exit'] or r['check_exit'] for r in results))
else:
    raise SystemExit('Choose install, tests or checks')
