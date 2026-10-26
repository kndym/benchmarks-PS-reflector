"""Atomic publication that tolerates briefly open Windows readers."""
import time


def atomic_replace(temporary, destination):
    for attempt in range(100):
        try:
            temporary.replace(destination)
            return
        except PermissionError:
            if attempt == 99:
                raise
            time.sleep(.1)
