import typing as t
from collections import defaultdict
from intervaltree import IntervalTree, Interval

from probelib.align.block import Block, Aln, in_region


class AvoidOTP(object):
    """Avoid Out of Target Peak."""
    def __init__(self,
                 target_region: t.Tuple[str, int, int],
                 density_thresh: float = 1e-3,
                 search_range: t.Tuple[int, int] = (-10**6, 10**6)):
        self.trees = defaultdict(IntervalTree)
        self.target_region = target_region
        self.density_thresh = density_thresh
        self.search_range = search_range

    def add(self, alns: t.Iterable[Aln]):
        for aln in alns:
            rname, start, end = aln
            tree = self.trees[rname]
            tree[start:end] = rname

    def remove_from_tree(self, inserted):
        for rname, itv in inserted:
            self.trees[rname].remove(itv)

    def filter(self, g: t.Iterable[Block]) -> t.Iterable[Block]:
        for qname, alns in g:
            inserted = []
            for aln in alns:
                search = True
                rname, start, end = aln
                tree = self.trees[rname]
                if in_region(self.target_region, rname, start, end):
                    if tree[start:end]:
                        # avoid overlap in target region
                        self.remove_from_tree(inserted)
                        break
                    else:
                        search = False
                if search:
                    range_ = (start+self.search_range[0], end+self.search_range[1])
                    in_range = tree[range_[0]:range_[1]]
                    density = len(in_range) / (self.search_range[1] - self.search_range[0])
                    if density >= self.density_thresh:
                        # avoid peak
                        self.remove_from_tree(inserted)
                        break
                itv = Interval(start, end, qname)
                tree.add(itv)
                inserted.append((rname, itv))
            else:
                yield qname, alns
