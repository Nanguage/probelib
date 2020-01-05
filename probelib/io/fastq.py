import typing as t


def write_fq(
        path: str,
        g: t.Iterable[t.Tuple[str, str, str]]) -> str:
    """Write contents in sub-sequence generator to fastq file.

    :param path: Target fastq file.
    :param g: sub-sequences generator.
    :return: output fastq path.
    """
    with open(path, 'w') as f:
        for seqname, seq, qualstr in g:
            f.write("@"+seqname+"\n")
            f.write(seq+"\n")
            f.write("+\n")
            f.write(qualstr+"\n")
    return path


def read_fq(path: str) -> t.Iterable[t.Tuple[str, str, str]]:
    """Read fastq file.

    :param path: Input fastq path
    :return: Generator of (seqname, seq, qualstr)
    """
    with open(path) as f:
        seqname, seq, qualstr = None, None, None
        for idx, line in enumerate(f):
            i = idx % 4
            if i == 0:
                if seqname is not None:
                    yield seqname, seq, qualstr
                seqname = line.strip()[1:]
            elif i == 1:
                seq = line.strip()
            elif i == 2:
                continue
            else:
                qualstr = line.strip()
        if seqname is not None:
            yield seqname, seq, qualstr


def fq2fa(fq: str,
          fa: str) -> str:
    """Convert fastq file to fasta file.

    :param fq: Path to input fastq.
    :param fa: Path to output fasta.
    :return: Path to output fasta file.
    """
    with open(fa, 'w') as f:
        for seqname, seq, _ in read_fq(fq):
            f.write(f">{seqname}\n")
            f.write(seq+"\n")
    return fa

