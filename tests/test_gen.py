from os.path import abspath, dirname, join
import sys
sys.path.insert(0, ".")
from probelib.gen import gen_random_seqs, slide_through, slide_through_fasta
from probelib.io.fasta import read_fa

HERE = dirname(abspath(__file__))
fapath = join(HERE, "data/lib.fa")


def test_gen_random_seqs():
    g = gen_random_seqs(10)
    for i in range(10):
        s = next(g)
        assert len(s) == 10


def test_slide_through():
    def check(seq, slen, ov):
        seq_pos = list(slide_through(seq, slen, ov))
        assert len(seq_pos) <= ((len(seq) - slen) / (slen - ov)) + 1
        assert len(seq_pos) >= ((len(seq) - slen) / (slen - ov))
        assert all([len(s) == slen for s, _, _ in seq_pos])
    for _, seq in read_fa(fapath):
        for i in range(10):
            check(seq, 10, i)


def test_slide_through_fasta():
    slen = 10; ov = 2
    old = None
    for sub, name, s, e in slide_through_fasta(fapath, slen, ov):
        assert type(sub) is type(name) is str
        assert type(s) is type(e) is int
        assert e - s == slen == len(sub)
        if (old is not None) and (old[0] == name):
            assert s - old[1] == slen - ov
        old = (name, s, e)
