import typing as t
from pyfaidx import Fasta


def read_fa(path: str) -> t.Iterable[t.Tuple[str, str]]:
    fa = Fasta(path)
    for key in fa.keys():
        yield key, fa[key][:].seq


