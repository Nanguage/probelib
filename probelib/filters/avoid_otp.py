import typing as t
from collections import defaultdict
from intervaltree import IntervalTree, Interval

from probelib.align.block import Block, Aln, is_overlap, GenomicRegion


class AvoidOTP(object):
    """Avoid Out of Target Peak."""
    def __init__(self,
                 target_regions: t.List[t.Tuple[str, int, int]],
                 density_thresh: float = 1e-3,
                 search_range: t.Tuple[int, int] = (-10**6, 10**6),
                 avoid_target_overlap: bool = True):
        self.trees = defaultdict(IntervalTree)
        self.target_region = target_regions
        self.density_thresh = density_thresh
        self.search_range = search_range
        self.avoid_target_overlap = avoid_target_overlap

    def add(self, alns: t.Iterable[Aln]):
        for aln in alns:
            rname, start, end = aln
            tree = self.trees[rname]
            tree[start:end] = rname

    def remove_from_tree(self, inserted: t.Iterable[t.Tuple[str, Interval]]):
        for rname, itv in list(set(inserted)):
            self.trees[rname].remove(itv)

    def overlap_with_targets(self, r: GenomicRegion):
        return [tr for tr in self.target_region if is_overlap(tr, r)]

    def filter(self, g: t.Iterable[Block]) -> t.Iterable[Block]:
        for qname, seq, alns in g:
            inserted = []
            for aln in alns:
                search = True
                rname, start, end = aln
                tree = self.trees[rname]

                if self.overlap_with_targets((rname, start, end)):
                    if self.avoid_target_overlap and tree[start:end]:
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
                yield qname, seq, alns
