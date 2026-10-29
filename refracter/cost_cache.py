"""Scoped, memory-bounded memoization of immutable optical cost blocks."""
from contextlib import contextmanager
from contextvars import ContextVar
import numpy as np

ACTIVE_COST_CACHE = ContextVar('optical_cost_cache',default=None)


class CostBlockCache:
    def __init__(self,x,y,max_bytes,block_rows):
        if not isinstance(max_bytes,int) or max_bytes<0 or block_rows<1:
            raise ValueError('Cost cache budget must be nonnegative and block_rows positive')
        # Private read-only snapshots prevent cached costs outliving mutable inputs.
        self.x=np.array(x,dtype=np.float64,order='C',copy=True)
        self.y=np.array(y,dtype=np.float64,order='C',copy=True)
        self.x.flags.writeable=self.y.flags.writeable=False
        self.max_bytes=max_bytes;self.block_rows=block_rows
        source=int(max_bytes*.8);target=int(max_bytes*.1)
        self.limits=dict(source_self=source,target_self=target,cross=max_bytes-source-target)
        self.used={kind:0 for kind in self.limits}
        self.blocks={};self.hits=0;self.misses=0;self.rejected=0;self.peak_bytes=0

    @staticmethod
    def _same_array(a,b):
        return a.shape==b.shape and a.strides==b.strides and a.ctypes.data==b.ctypes.data

    @staticmethod
    def _row_start(a,full):
        if a.ndim!=2 or a.shape[1:]!=full.shape[1:] or a.strides!=full.strides:
            return None
        offset=a.ctypes.data-full.ctypes.data
        if offset<0 or offset%full.strides[0]: return None
        start=offset//full.strides[0]
        return start if start+len(a)<=len(full) else None

    def key(self,x,y,kappa):
        if not self.max_bytes or not 0<len(x)<=self.block_rows: return None
        if self._same_array(y,self.x):
            start=self._row_start(x,self.x);kind='source_self'
        elif self._same_array(y,self.y):
            start=self._row_start(x,self.y);kind='target_self'
            if start is None:
                start=self._row_start(x,self.x);kind='cross'
        else:
            return None
        return (kind,float(kappa),start,len(x)) if start is not None else None

    def get(self,key):
        result=self.blocks.get(key)
        if result is not None: self.hits+=1
        else: self.misses+=1
        return result

    def put(self,key,value):
        kind=key[0]
        if value.nbytes>self.limits[kind]-self.used[kind]:
            self.rejected+=1
            return
        # No LRU eviction: repeated scans retain useful prefixes instead of thrashing.
        value.flags.writeable=False
        self.blocks[key]=value;self.used[kind]+=value.nbytes
        self.peak_bytes=max(self.peak_bytes,sum(self.used.values()))

    def stats(self):
        return dict(budget_bytes=self.max_bytes,peak_retained_bytes=self.peak_bytes,
            hits=self.hits,misses=self.misses,admission_rejections=self.rejected,
            retained_bytes_by_kind=dict(self.used),block_rows=self.block_rows)


@contextmanager
def bounded_cost_cache(x,y,max_bytes,block_rows=256):
    cache=CostBlockCache(x,y,max_bytes,block_rows)
    token=ACTIVE_COST_CACHE.set(cache)
    try:
        yield cache
    finally:
        ACTIVE_COST_CACHE.reset(token)
        cache.blocks.clear()
        cache.used={kind:0 for kind in cache.used}
