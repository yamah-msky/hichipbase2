import gzip
from pathlib import Path

def is_gzip(path):
    """Check whether the file given is gzipped"""
    path = Path(path)
    GZIP_MAGIC = b'\x1f\x8b'
    with open(path, 'rb') as f:
        return f.read(2) == GZIP_MAGIC

def open_maybe_gzip(path, mode: str):
    """Open an ambiguous file - whether is gzipped is not sure.
    Do not forget to close the file opened"""
    path = Path(path)
    if is_gzip(path):
        return gzip.open(path, mode)
    else:
        return open(path, mode)