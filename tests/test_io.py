from os.path import abspath, dirname, join
import sys
sys.path.insert(0, ".")
from probelib.io.fasta import read_fa
from probelib.io.fastq import write_fq, read_fq, fq2fa
from probelib.io.align import read_align_blocks

HERE = dirname(abspath(__file__))
fa_path = join(HERE, "data/lib.fa")


def test_read_fa():
    for name, seq in read_fa(fa_path):
        assert type(name) is str
        assert type(seq) is str


def test_write_fq():
    tmp_fq_path = "/tmp/test_write.fq"

    def gen():
        for i in range(5):
            yield "AAATTTCCC", "chr1", i, i + 10
    g = gen()
    write_fq(tmp_fq_path, g)

    with open(tmp_fq_path) as f:
        for idx, line in enumerate(f):
            if idx % 4 == 0:
                assert line.startswith("@")
            elif idx % 4 == 1:
                assert line == "AAATTTCCC\n"
            elif idx % 4 == 2:
                assert line == "+\n"
            else:
                assert line == "~" * 9 + "\n"


def test_read_fq():
    tmp_fq_path = "/tmp/test_read.fq"

    def gen():
        for i in range(5):
            yield "AAATTTCCC", "chr1", i, i + 10

    g = gen()
    write_fq(tmp_fq_path, g)
    recs = list(read_fq(tmp_fq_path))
    assert len(recs) == 5
    for seqname, seq, qual in recs:
        assert seq == "AAATTTCCC"
        assert type(seqname) is type(qual) is str


def test_fq2fa():
    tmp_fq_path = "/tmp/test_fq2fa.fq"

    def gen():
        for i in range(5):
            yield "AAATTTCCC", "chr1", i, i + 10

    g = gen()
    write_fq(tmp_fq_path, g)
    tmp_fa_path = "/tmp/test_fq2fa.fa"
    fq2fa(tmp_fq_path, tmp_fa_path)
    for name, seq in read_fa(tmp_fa_path):
        assert name.startswith("chr")
        assert seq == "AAATTTCCC"


def test_read_align_blocks():
    tmp_sam = "/tmp/test_read_aln_block.sam"
    from probelib.gen import slide_through_fasta
    from probelib.io.fastq import write_fq
    from probelib.align.bowtie2 import align_se_sen as align_bt
    from itertools import islice
    gen = slide_through_fasta(fa_path, 40, 30)
    gen = islice(gen, 0, 3)
    tmp_fq_path = "/tmp/test_read_align_blocks.fq"
    write_fq(tmp_fq_path, gen)
    bowtie_index = join(HERE, "data/bowtie2-index/lib")
    align_bt(tmp_fq_path, bowtie_index, tmp_sam, log=None)
    blocks = []
    for q_name, block in read_align_blocks(tmp_sam):
        assert type(q_name) is str
        assert type(block) is list
        assert len(block) > 0
        blocks.append(block)
    assert len(blocks) == 3
