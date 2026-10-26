"""Atomic publication that tolerates briefly open Windows readers."""
import time
import json


def atomic_replace(temporary, destination):
    for attempt in range(100):
        try:
            temporary.replace(destination)
            return
        except PermissionError:
            if attempt == 99:
                raise
            time.sleep(.1)


def atomic_write_json(destination, value):
    temporary = destination.with_suffix('.tmp.json')
    temporary.write_text(json.dumps(value, default=str, indent=2)+'\n', encoding='utf-8')
    atomic_replace(temporary, destination)
