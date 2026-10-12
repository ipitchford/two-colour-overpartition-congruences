"""Run all declared proof-critical computations; discovery scripts are excluded."""
import importlib.util
from pathlib import Path
import subprocess
import sys
R=Path(__file__).resolve().parent
def module(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
if __name__=='__main__':
    audit=module('audit',R/'checks/independent_bsc_audit.py')
    v=module('certificate',R/'papers/certificates/code/verify.py')
    for c in audit.CASES:
        F={}
        for d,e in [(1,-2),(2,1),(c['ell'],-2),(2*c['ell'],1),*c['B'].items()]:F[d]=F.get(d,0)+e
        v.main(c['ell'],c['M'],c['r'],c['j'],c['N'],F)
    for script in ['checks/independent_bsc_audit.py','checks/families.py']:
        subprocess.run([sys.executable,str(R/script)],check=True)
    print('PASS: all declared certificates and family corroboration')
