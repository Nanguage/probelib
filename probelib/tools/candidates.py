from probelib.gen import slide_through_fasta, to_fq_rec
from probelib.io.fastq import write_fq


def candidates(in_file: str,
               out_file: str,
               length: int,
               overlap: int):
    gen = slide_through_fasta(in_file, length, overlap)
    gen = map(to_fq_rec, gen)
    write_fq(out_file, gen)

