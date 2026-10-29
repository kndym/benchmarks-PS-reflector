"""Reuse identical optical solves across map consumers, without changing math."""
import ctypes
import hashlib
import json
import os
from pathlib import Path
import time
import numpy as np
from atomic_io import atomic_replace


def process_alive(pid):
    if os.name == 'nt':
        kernel = ctypes.WinDLL('kernel32', use_last_error=True)
        kernel.OpenProcess.restype = ctypes.c_void_p
        kernel.OpenProcess.argtypes = [ctypes.c_uint32,ctypes.c_bool,ctypes.c_uint32]
        kernel.CloseHandle.argtypes = [ctypes.c_void_p]
        kernel.GetExitCodeProcess.argtypes = [ctypes.c_void_p,ctypes.POINTER(ctypes.c_uint32)]
        handle = kernel.OpenProcess(0x1000,False,pid)
        if not handle:
            return ctypes.get_last_error() == 5  # Unknown access: don't steal a live lock.
        code = ctypes.c_uint32()
        try:
            return bool(kernel.GetExitCodeProcess(handle,ctypes.byref(code))) and code.value == 259
        finally:
            kernel.CloseHandle(handle)
    try:
        os.kill(pid,0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True


def problem_key(arrays, settings):
    h=hashlib.sha256(json.dumps(settings,sort_keys=True,allow_nan=False).encode())
    for array in arrays:
        value=np.ascontiguousarray(array,dtype='<f8')
        h.update(str(value.shape).encode());h.update(value.tobytes())
    return h.hexdigest()


def get_or_solve(directory, key, solve, poll_seconds=1):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    cache=directory/(key+'.npz');lock=directory/(key+'.lock')
    while True:
        if cache.exists():
            with np.load(cache,allow_pickle=False) as data:
                if str(data['_identity']) != key:
                    raise ValueError('Shared optical cache identity mismatch')
                solution={name:(data[name].item() if data[name].ndim==0 else data[name].copy())
                          for name in data.files if name!='_identity'}
            if not all(np.isfinite(value).all() for value in solution.values()):
                raise ValueError('Nonfinite shared optical solution')
            print('Reused shared optical solve '+key[:12],flush=True)
            return solution
        try:
            with lock.open('x') as stream:
                json.dump(dict(pid=os.getpid(),created_at=time.time()),stream)
        except FileExistsError:
            try:
                owner=json.loads(lock.read_text())
                if not process_alive(int(owner['pid'])):
                    # Recheck unchanged ownership before removing a stale lock.
                    if json.loads(lock.read_text()) == owner:
                        lock.unlink(missing_ok=True)
                        continue
            except FileNotFoundError:
                pass
            except json.JSONDecodeError:
                # A killed writer can leave an empty newly-created lock.
                try:
                    if time.time()-lock.stat().st_mtime>60:
                        lock.unlink(missing_ok=True)
                except FileNotFoundError:
                    pass
            time.sleep(poll_seconds)
            continue
        try:
            # A writer may have published between our existence check and lock.
            if cache.exists():
                continue
            solution=solve()
            if not all(np.isfinite(value).all() for value in solution.values()):
                raise ValueError('Refusing nonfinite shared optical solution')
            temporary=directory/(key+'.tmp.npz')
            np.savez(temporary,_identity=np.array(key),**solution)
            atomic_replace(temporary,cache)
            print('Saved shared optical solve '+key[:12],flush=True)
            return solution
        finally:
            lock.unlink(missing_ok=True)
