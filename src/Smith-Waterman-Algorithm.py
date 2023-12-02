"""Smith Waterman Algorithm (DNA Alignment)."""

# Written by CoolCat467 07/02/2020

from __future__ import annotations

__title__ = "Smith Waterman Algorithm"
__version__ = "0.0.1"


from typing import TYPE_CHECKING, TypeVar

if TYPE_CHECKING:
    from collections.abc import Sequence

T = TypeVar("T")


class SMAlgorythm:
    def __init__(self, sequence1: str, sequence2: str) -> None:
        self.seqA = str(sequence1).upper()
        self.seqB = str(sequence2).upper()
        self.n = len(self.seqA)
        self.m = len(self.seqB)
        self.ntides = ["A", "T", "G", "C"]
        self.sub_matrix = self.get_allignment_dict()
        # Linear Width
        self.w1 = 2
        # Affine Width
        self.v = 5  # Opening gap penalty
        self.u = 1  # Gap Extension penalty
        ##        self.W = lambda k: (self.u*k) + self.v if self.v > 0 and self.u > 0 else 999999
        self.score_matrix = [
            [0 for i in range(self.m + 1)] for i in range(self.n + 1)
        ]

    def __repr__(self) -> str:
        return "SMAlgorythm()"

    def get_allignment_dict(self) -> dict[tuple[str, str], int]:
        def shift(x: list[str]) -> list[str]:
            return x[-1:] + x[:-1]

        x = list(self.ntides)
        tides = list(x)
        for _i in range(len(self.ntides) - 1):
            x = shift(x)
            tides += x

        def is_same(x: T, y: T) -> int:
            return 3 if x == y else -3

        return {
            t: is_same(*t) for t in zip(tides, self.ntides * len(self.ntides))
        }

    def score(self) -> None:
        a = self.seqA
        b = self.seqB
        H = list(self.score_matrix)
        for k in range(self.n):
            for l in range(self.m):
                H[k][0] = 0
                H[0][l] = 0
        # main
        for i in range(1, self.n):
            for j in range(1, self.m):
                H[i][j] = max(
                    H[i - 1][j - 1]
                    + self.sub_matrix[(a[i], b[j])],  # Aligning Score
                    H[i - k][j]
                    - k
                    * self.w1,  # Score if A is at the end of a gap of length k
                    H[i][j - l]
                    - self.w1,  # Score if B is at the end of a gap of length l
                    0,
                )
        self.score_matrix = H
        self.k = k
        self.l = l

    def calculate_start_index(self) -> None:
        cmax = 0
        for i in range(self.n + 1):
            for ii in range(self.m + 1):
                if self.score_matrix[i][ii] > cmax:
                    cmax = self.score_matrix[i][ii]
                    idx = [i, ii]
        self.start_index = tuple(idx)

    def calculate_path(self) -> None:
        i, j = tuple(self.start_index)
        path = [tuple(self.start_index)]
        while self.score_matrix[i][j] != 0:
            frm = (
                self.score_matrix[i - 1][j - 1]
                + self.sub_matrix[
                    (self.seqA[i], self.seqB[j])
                ],  # Aligning Score
                self.score_matrix[i - self.k][j]
                - self.k
                * self.w1,  # Score if A is at the end of a gap of length k
                self.score_matrix[i][j - self.l]
                - self.l
                * self.w1,  # Score if B is at the end of a gap of length l
                0,
            )
            idxs = {
                frm[0]: (i - 1, j - 1),
                frm[1]: (i - self.k, j),
                frm[2]: (i, j - self.l),
            }
            m = max(frm)
            if m != 0:
                i, j = idxs[m]
                path.append(idxs[m])
                continue
            break
        self.path = [[self.seqA[i], self.seqB[j]] for i, j in reversed(path)]

    @classmethod
    def path_to_string(cls, path: list[list[str]]) -> str:
        relative = []
        for i, ii in path:
            relative.append([i, "|" if i == ii else " ", ii])
        data = ["", "", ""]
        for r in range(len(relative)):
            i, rel, ii = relative[r]
            data[0] += i
            data[1] += rel
            data[2] += ii
        return "".join(["".join(data[i]) + "\n" for i in range(3)])[:-1]

    def align(self, to_string: bool = False) -> Sequence[Sequence[str]]:
        self.score()
        self.calculate_start_index()
        self.calculate_path()
        if to_string:
            return self.__class__.path_to_string(self.path)
        return self.path


if __name__ == "__main__":
    seq1 = "catatgcgcattatgat"
    ##    seq1 = seq1 + ''.join([i for i in reversed(seq1)])
    seq2 = "gatatgatcattatgct"
    print(seq1)
    print(seq2)
    print(SMAlgorythm(seq1, seq2).align(True))
