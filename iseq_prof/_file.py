import os
from contextlib import contextmanager
from pathlib import Path
from tempfile import NamedTemporaryFile

__all__ = ["inplace", "assert_file_exist"]


@contextmanager
def inplace(filepath: Path):
    folder = filepath.parents[0]
    name = filepath.name
    success = True
    try:
        with NamedTemporaryFile(dir=folder, prefix="." + name, delete=False) as tmp:
            yield Path(tmp.name)
        success = False
    finally:
        if success:
            os.replace(tmp.name, filepath)


def assert_file_exist(filepath: Path):
    if not filepath.exists():
        raise RuntimeError(f"{filepath} does not exist.")
