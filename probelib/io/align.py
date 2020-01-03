import typing as t
import pysam


aln = t.Tuple[str, int, int]


def read_align_blocks(
        sam_path: str
        ) -> t.Iterable[t.Tuple[str, t.List[aln]]]:
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

