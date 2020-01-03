import typing as t
import pysam

Aln = t.Tuple[str, int, int]
Block = t.Tuple[str, t.List[Aln]]


def read_align_blocks(
        sam_path: str
        ) -> t.Iterable[Block]:
    def yield_cond(old, rec, block, end=False):
        res = (old is not None) and (len(block) > 0)
        if res and not end:
            res &= rec.query_name != old.query_name
        return res
    with pysam.AlignmentFile(sam_path, mode='r') as sam:
        block = []
        old = None
        for rec in sam.fetch():
            if yield_cond(old, rec, block):
                yield old.query_name, block
                block = []
            else:
                if rec.reference_name is not None:
                    aln = rec.reference_name, rec.reference_start, rec.reference_end
                    block.append(aln)
            old = rec
        if yield_cond(old, rec, block, end=True):
            yield old.query_name, block


GenomicRegion = t.Tuple[str, int, int]


def in_region(target_region: GenomicRegion,
              name: str, start: int, end: int) -> bool:
    r = target_region
    return (r[0] == name) & (start >= r[1]) & (end >= r[2])


def specificity(
        block: Block,
        target_region: GenomicRegion) -> float:
    qname, alns = block
    in_range = [in_region(target_region, n, s, e) for (n, s, e) in alns]
    return len(in_range) / len(alns)
