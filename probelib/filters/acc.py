import typing as t

import editdistance.bycython
edit_dist = editdistance.bycython.eval


class Unique(object):

    def __init__(self):
        self.set = dict()

    def add(self, e):
        if e in self.set:
            self.set[e] += 1
        else:
            self.set[e] = 0

    def filter(self, up_stream: t.Iterable) -> t.Iterable:
        """Let the elements in stream occurrence only once."""
        for e in up_stream:
            if e not in self.set:
                self.add(e)
                yield e
            else:
                continue


class KeepDist(object):

    def __init__(self,
                 d: float,
                 dist_func: t.Callable = edit_dist
                 ):
        self.set = []
        self.d = d
        self.dist_func = dist_func

    def add(self, e):
        self.set.append(e)

    def filter(self, up_stream: t.Iterable) -> t.Iterable:
        """Let elements in stream keep distance with existing elements."""
        for e in up_stream:
            for e2 in self.set:
                if self.dist_func(e, e2) <= self.d:
                    break
            else:
                self.add(e)
                yield e
