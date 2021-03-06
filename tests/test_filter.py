from os.path import abspath, dirname, join
import sys
sys.path.insert(0, ".")
from probelib.filters.avoid_otp import AvoidOTP
from probelib.filters.acc import Unique, KeepDist

HERE = dirname(abspath(__file__))
fa_path = join(HERE, "data/lib.fa")


def mock_align_blocks():
    blocks = [
        ('q1', "", [
            ('chr1', 100, 110),
            ('chr1', 110, 120),
            ('chr1', 120, 130),
        ]),
        ('q2', "", [
            ('chr1', 130, 140),
            ('chr2', 120, 130),
        ]),
        ('q3', "", [
            ('chr1', 140, 150),
            ('chr2', 130, 140),
            ('chr2',  90, 100),
        ]),
        ('q4', "", [
            ('chr1', 100, 110),
        ]),
    ]
    return blocks


def test_AvoidOTP():
    align_blocks = mock_align_blocks()
    aotp = AvoidOTP(target_regions=[("chr1", 100, 200)],
                    density_thresh=0.01,
                    search_range=(-100, 100))
    filtered = list(aotp.filter(align_blocks))
    assert len(filtered) == 2
    for name, seq, aligns in filtered:
        assert name != "q3"
        assert name != "q4"


def test_Unique():
    seqs = ["AAA", "ATA", "AAA", "TTT"]
    acc = Unique()
    filtered = list(acc.filter(seqs))
    assert filtered == ["AAA", "ATA", "TTT"]


def test_keepDist():
    seqs = ["AAA", "AAT", "GGG", "GGC"]
    acc = KeepDist(2)
    filtered = list(acc.filter(seqs))
    assert filtered == ["AAA", "GGG"]

