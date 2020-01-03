from os.path import abspath, dirname, join
import sys
sys.path.insert(0, ".")
from probelib.align.bowtie2 import align_se_sen as align_bt
from probelib.align.blastn import align as align_bl

HERE = dirname(abspath(__file__))
fapath = join(HERE, "data/lib.fa")


def test_bowtie2_align():
    from probelib.gen import slide_through_fasta
    from probelib.io.fastq import write_fq
    from itertools import islice
    fa_gen = slide_through_fasta(fapath, 40, 30)
    gen = islice(fa_gen, 0, 3)
    tmp_fq_path = "/tmp/test_bowtie.fq"
    write_fq(tmp_fq_path, gen)
    tmp_sam_path = "/tmp/test_bowtie.sam"
    bowtie_index = join(HERE, "data/bowtie2-index/lib")
    assert align_bt(tmp_fq_path, bowtie_index, tmp_sam_path) == tmp_sam_path
    with open(tmp_sam_path) as f:
        for line in f:
            items = line.strip().split("\t")
            assert items[1].isdigit()


def test_blastn_align():
    from probelib.gen import slide_through_fasta
    from probelib.io.fastq import write_fq
    from itertools import islice
    fa_gen = slide_through_fasta(fapath, 40, 30)
    gen = islice(fa_gen, 0, 3)
    tmp_fq_path = "/tmp/test_blastn.fq"
    write_fq(tmp_fq_path, gen)
    tmp_sam_path = "/tmp/test_blastn.sam"
    blastn_db = join(HERE, "data/blastdb/lib")
    assert align_bl(tmp_fq_path, blastn_db, tmp_sam_path) == tmp_sam_path
    with open(tmp_sam_path) as f:
        for line in f:
            items = line.strip().split("\t")
            assert items[1].isdigit()
