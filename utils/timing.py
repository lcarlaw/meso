from functools import wraps
from time import time
import logging as log

def timeit(f):
    @wraps(f)
    def wrap(*args, **kwargs):
        ts = time()
        result = f(*args, **kwargs)
        te = time()
        log.info('Function: %s took: %2.4f sec' % (f.__name__, te-ts))
        return result
    return wrap
