import typing as t
import pysam

Aln = t.Tuple[str, int, int]
Block = t.Tuple[str, str, t.List[Aln]]


def read_align_blocks(
        sam_path: str
        ) -> t.Iterable[Block]:
    def yield_cond(old, rec, block, end=False):
        res = (old is not None) and (len(block) > 0)
        if res and not end:
            res &= rec.query_name != old.query_name
        return res
    with pysam.AlignmentFile(sam_path, mode='r') as sam:
        alns = []
        old = None
        for rec in sam.fetch():
            if yield_cond(old, rec, alns):
                yield old.query_name, old.query_sequence, alns
                alns = []
            else:
                if rec.reference_name is not None:
                    aln = rec.reference_name, rec.reference_start, rec.reference_end
                    # is correct use this to represent the align position?
                    alns.append(aln)
            old = rec
        if yield_cond(old, rec, alns, end=True):
            yield old.query_name, old.query_sequence, alns


GenomicRegion = t.Tuple[str, int, int]


def parse_region(region: str) -> GenomicRegion:
    chr_, r_ = region.split(":")
    s, e = r_.split("-")
    s, e = int(s), int(e)
    return chr_, s, e


def is_in_region(target_region: GenomicRegion,
                 r2: GenomicRegion) -> bool:
    r = target_region
    return (r[0] == r2[0]) & (r2[1] >= r[1]) & (r2[2] <= r[2])


def is_overlap(r1: GenomicRegion, r2: GenomicRegion) -> bool:
    if r1[0] != r2[0]:
        return False
    return (
        (r1[1] <= r2[1] <= r1[2]) |
        (r1[1] <= r2[2] <= r1[2]) |
        (r2[1] <= r1[1] <= r2[2]) |
        (r2[1] <= r1[2] <= r2[2])
    )


def count_ovlap_with_region(
        block: Block,
        target_region: GenomicRegion) -> t.Tuple[int, int]:
    _, _, alns = block
    n_in = 0
    for n, s, e in alns:
        if is_overlap(target_region, (n, s, e)):
            n_in += 1
    return n_in, len(alns) - n_in
