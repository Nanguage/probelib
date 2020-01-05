from os.path import abspath, dirname, join
import sys
sys.path.insert(0, ".")
from probelib.io.fasta import read_fa
from probelib.io.fastq import write_fq, read_fq, fq2fa
from probelib.align.block import read_align_blocks

HERE = dirname(abspath(__file__))
fa_path = join(HERE, "data/lib.fa")


def test_read_fa():
    for name, seq in read_fa(fa_path):
        assert type(name) is str
        assert type(seq) is str


def gen():
    for i in range(5):
        yield f"chr1_{i}_{i+10}", "AAATTTCCC", "~~~~~~~~~"


def test_write_fq():
    tmp_fq_path = "/tmp/test_write.fq"

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

    g = gen()
    write_fq(tmp_fq_path, g)
    recs = list(read_fq(tmp_fq_path))
    assert len(recs) == 5
    for seqname, seq, qual in recs:
        assert seq == "AAATTTCCC"
        assert type(seqname) is type(qual) is str


def test_fq2fa():
    tmp_fq_path = "/tmp/test_fq2fa.fq"

    g = gen()
    write_fq(tmp_fq_path, g)
    tmp_fa_path = "/tmp/test_fq2fa.fa"
    fq2fa(tmp_fq_path, tmp_fa_path)
    for name, seq in read_fa(tmp_fa_path):
        assert name.startswith("chr")
        assert seq == "AAATTTCCC"


def mimic_sam_file(path):
    contents = """@HD     VN:1.0  SO:unsorted
@SQ     SN:chr1 LN:100
@SQ     SN:chr3 LN:300
@PG     ID:bowtie2      PN:bowtie2      VN:2.3.5        CL:"/opt/anaconda/envs/probe_design/bin/bowtie2-align-s --wrapper basic-0 -x /home/wzxu/Projects/probelib/tests/data/bowtie2-index/lib -t -k 100 --very-sensitive-local -p 10 -S /tmp/test_bowtie.sam -U /tmp/test_bowtie.fq"
chr1:0_40       0       chr1    1       1       40M     *       0       0       GGAGGACTGTGCTAAAAGAAGTAACCCATCATCCTTATCA     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        AS:i:80 XS:i:80      XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:40 YT:Z:UU
chr1:0_40       256     chr3    101     1       40M     *       0       0       GGAGGACTGTGCTAAAAGAAGTAACCCATCATCCTTATCA     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        AS:i:80 XS:i:80      XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:40 YT:Z:UU
chr1:0_40       256     chr3    151     1       40M     *       0       0       GGAGGACTGTGCTAAAAGAAGTAACCCATCATCCTTATCA     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        AS:i:56 XS:i:80      XN:i:0  XM:i:3  XO:i:0  XG:i:0  NM:i:3  MD:Z:9G1G0G27   YT:Z:UU
chr1:10_50      0       chr1    11      1       40M     *       0       0       GCTAAAAGAAGTAACCCATCATCCTTATCACTGAGTGTAA     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        AS:i:80 XS:i:80      XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:40 YT:Z:UU
chr1:10_50      256     chr3    111     1       40M     *       0       0       GCTAAAAGAAGTAACCCATCATCCTTATCACTGAGTGTAA     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        AS:i:80 XS:i:80      XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:40 YT:Z:UU
"""
    import re
    contents = re.sub(r" +", "\t", contents)
    with open(path, 'w') as f:
        f.write(contents)


def test_read_align_blocks():
    tmp_sam = "/tmp/test_read_aln_block.sam"
    mimic_sam_file(tmp_sam)
    blocks = []
    for q_name, seq, alns in read_align_blocks(tmp_sam):
        assert type(q_name) is type(seq) is str
        assert type(alns) is list
        assert len(alns) > 0
        blocks.append(alns)
        if q_name == "chr1:0_40":
            assert len(alns) == 3
        if q_name == "chr1:10_50":
            assert len(alns) == 2
    assert len(blocks) == 2
