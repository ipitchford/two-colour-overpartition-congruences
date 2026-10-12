"""Require plausible certificate mutations to fail, also under optimized Python."""
import contextlib
import io
from verify import module,R
v=module('mutant_verifier',R/'papers/certificates/code/verify.py')
F={1:-2,2:1,8:-2,16:1,32:69,64:-29,128:-2}
tests=[('wrong residue',8,32,27,5,128,F),('doubled modulus',8,32,28,6,128,F),('invalid U level',8,31,28,5,128,F),('missing eta index',8,32,28,5,128,{k:e for k,e in F.items() if k!=1})]
for label,*args in tests:
    try:
        with contextlib.redirect_stdout(io.StringIO()):v.main(*args)
    except ValueError:print('REJECTED',label)
    else:raise RuntimeError('Mutation escaped: '+label)
print('PASS: all semantic mutations rejected')
