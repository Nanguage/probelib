import os.path as osp
import subprocess as subp


def tab2sam(tab: str,
            sam: str) -> str:
    """Not work Yet"""
    cmd = f"blast2sam.pl {tab} > {sam}"
    subp.check_call(cmd, shell=True)
    return sam


def align(fq_path: str,
          index: str,
          sam_path: str) -> str:
    """ Alignment with blastn

    :param fq_path:
    :param index:
    :param sam_path:
    :return:
    """
    from probelib.io.fastq import fq2fa
    fa_path = osp.splitext(sam_path)[0] + '.tmp.fa'
    fq2fa(fq_path, fa_path)
    tab_path = osp.splitext(sam_path)[0] + '.tmp.tab'
    outfmt = "6 qseqid sseqid pident nident length mismatch positive gapopen gaps ppos qframe sframe sstrand qcovs qstart qend qseq sstart send sseq evalue bitscore score"
    cmd = ["blastn", "-query", fa_path, "-db", index,
           "-outfmt", outfmt, "-out", tab_path]
    subp.check_call(cmd)
    tab2sam(tab_path, sam_path)
    subp.check_call(["rm", tab_path, fa_path])
    return sam_path

