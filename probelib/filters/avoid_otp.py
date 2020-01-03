import typing as t
from collections import defaultdict
from intervaltree import IntervalTree

from probelib.align.block import Block, Aln, in_region


class AvoidOTP(object):
    """Avoid Out of Target Peak."""
    def __init__(self,
                 target_region: t.Tuple[str, int, int],
                 density_thresh: float = 1e-3,
                 search_range: t.Tuple[int, int] = (-10**6, 10**6)):
        self.trees = defaultdict(lambda: IntervalTree)
        self.target_region = target_region
        self.density_thresh = density_thresh
        self.search_range = search_range

    def add(self, alns: t.Iterable[Aln]):
        for aln in alns:
            rname, start, end = aln
            tree = self.trees[rname]
            tree[start:end] = rname

    def filter(self, g: t.Iterable[Block]) -> t.Iterable[Block]:
        for qname, alns in g:
            for aln in alns:
                rname, start, end = aln
                tree = self.trees[rname]
                if in_region(self.target_region, rname, start, end):
                    if tree[start:end]:
                        # avoid overlap in target region
                        break
                    else:
                        continue
                search = (start-self.search_range[0], start+self.search_range[1])
                in_range = tree[search[0]:search[1]]
                density = len(in_range) / (self.search_range[1] - self.search_range[1])
                if density > self.density_thresh:
                    # avoid peak
                    break
            else:
                yield qname, alns
                self.add(alns)
