import typing as t
import re
from probelib.filters.acc import KeepDist
from probelib.gen import gen_random_seqs
from bramin import P, END
from bramin.curry import curry as c
from Bio.SeqUtils import MeltingTemp as mt


def primer_gen(seq_len: int, d: int, gc_r=(0.3, 0.7), mt_r=(0, 100)):
    return (
        gen_random_seqs(seq_len) | P
        | KeepDist(d).filter
        | c(filter, lambda s: mt_r[0] <= mt.Tm_NN(s) <= mt_r[1])
        | c(filter, lambda s: gc_r[0] <= len(re.findall("C|G", s))/len(s) <= gc_r[1])
        | END
    )


def pair(g: t.Iterable[t.Any]) -> t.Iterable[t.Tuple]:
    while True:
        e1 = next(g)
        e2 = next(g)
        yield e1, e2

