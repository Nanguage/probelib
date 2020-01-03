import typing as t

import editdistance.bycython
edit_dist = editdistance.bycython.eval


class Accumulator(object):

    def __init__(self):
        self.set = dict()

    def add(self, e):
        if e in self.set:
            self.set[e] += 1
        else:
            self.set[e] = 0

    def accumulate(self, up_stream: t.Iterable) -> t.Iterable:
        """Accumulate a stream."""
        for e in up_stream:
            self.add(e)
            yield e

    def unique(self, up_stream: t.Iterable) -> t.Iterable:
        """Let the elements in stream occurrence only once."""
        for e in up_stream:
            if e not in self.set:
                self.add(e)
                yield e
            else:
                continue

    def keep_distance(self,
                      up_stream: t.Iterable,
                      d: float,
                      dist_func: t.Callable = edit_dist) -> t.Iterable:
        """Let elements in stream keep distance with existing elements."""
        for e in up_stream:
            for e2 in self.set:
                if e is e2:
                    continue
                if dist_func(e, e2) <= d:
                    continue
                else:
                    yield e
                    break
