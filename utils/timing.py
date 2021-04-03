from functools import wraps
from time import time

def timeit(f):
    @wraps(f)
    def wrap(*args, **kw):
        ts = time()
        result = f(*args, **kw)
        te = time()
        #print('Function: %s with args:[%s, %s] took: %2.4f sec' % \
        #      (f.__name__, args, kw, te-ts))
        print('Function: %s took: %2.4f sec' % (f.__name__, te-ts))
        return result
    return wrap
