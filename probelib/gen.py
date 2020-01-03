import random
import typing as t


def random_seq(length: int,
               charset: str = "ATCG") -> str:
    """Get a random sequence.

    :param length: Length of sequence.
    :param charset: Char population.
    :return: Random sequence
    """
    return "".join(random.choices(charset, k=length))


def gen_random_seqs(seq_len: int,
                    charset: str = "ATCG") -> t.Iterable[str]:
    """Yield a series of random sequences.

    :param seq_len: Length of sequence.
    :param charset: Char population.
    :return: Generator of random sequences.
    """
    while True:
        yield random_seq(seq_len, charset)


def seq_permutation(seq_len: int,
                    charset: str = "ATCG") -> t.Iterable[str]:
    """Yield permutation of given length sequences.

    :param seq_len: Length of sequence.
    :param charset: Possible bases.
    :return: Generator of permutations.
    """
    if seq_len <= 0:
        yield ""
    else:
        for seq in seq_permutation(seq_len-1, charset):
            for c in charset:
                yield seq + c


def slide_through(seq: str,
                  slen: int,
                  overlap: int) -> t.Iterable[t.Tuple[str, int, int]]:
    """Slide through a sequence, generate it's subsequences.

    :param seq: Input sequence.
    :param slen: Length of sub-sequence.
    :param overlap: Overlap size between two sub-sequences.
    :return: Generator of a tuple of sub-sequences,
    and it's start and end position in original sequence.
    """
    assert slen > 0, "sub-sequence length must large than zero."
    assert overlap >= 0, 'overlap length must large or equal to zero.'
    assert overlap < slen, 'overlap length must less than sub-seq length.'
    step = slen - overlap
    tlen = len(seq)
    s = 0
    while s + slen <= tlen:
        e = s + slen
        sub = seq[s:e]
        yield sub, s, e
        s += step


def slide_through_fasta(path: str,
                        slen: int,
                        overlap: int) -> t.Iterable[t.Tuple[str, str, int, int]]:
    """Slide through all sequences in fasta file.

    :param path: Input fasta file.
    :param slen: Length of sub-sequence.
    :param overlap: Overlap size between two sub-sequences.
    :return: Generator of a tuple of sub-sequences,
    and it's seq-name and start, end position in original sequence.
    """
    from probelib.io.fasta import read_fa
    for name, seq in read_fa(path):
        for sub, s, e in slide_through(seq, slen, overlap):
            yield sub, name, s, e

