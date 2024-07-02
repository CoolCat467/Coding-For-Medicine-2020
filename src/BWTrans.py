"""Burrows-Wheeler transforming."""

# Programmed by Samuel Davenport
from __future__ import annotations

from random import randint
from typing import TypeVar

from Bio.Seq import Seq

T = TypeVar("T")


__author__ = "BWTrans"
__version__ = "0.0.0"


def get_frames(dna_seq: Seq, n: int = 3) -> list[str]:
    """Return proteins with a length at least <l> long from frame shifts <n>."""
    # For each frame shift,
    data = []
    for shift in range(abs(n) % 4):
        leng = len(dna_seq[shift:])
        mod = leng % 3
        # print(leng, mod)
        # Get the frame shift of the DNA sequence and
        # translate data to proteins
        proteins_data = dna_seq[shift:-mod].translate()
        # Split the proteins by their end codons
        proteins = str(proteins_data).split("*")
        # For each protein we found,
        ##        for protein in proteins:
        ##            # If the protein has at least 100 amino acids,
        ##            if len(protein) >= l:
        ##                # Write the protein to the file
        ##                file.write(protein+'\n')
        data += proteins
    return [protein for protein in data if len(protein) >= 5]


def gen_random_proper_seq(
    length: int,
    a: float,
    t: float,
    g: float,
    c: float,
) -> Seq:
    """Generate random proper sequence."""
    if sum([a, t, g, c]) != 1:
        raise ArithmeticError(
            "Sum of percentages of A, T, G, and C is not equal to 100 percent!",
        )
    a_copies = "A" * round(length * a)
    t_copies = "T" * round(length * t)
    g_copies = "G" * round(length * g)
    c_copies = "C" * round(length * c)
    unrand = list(f"{a}{t}{g}{c}")
    # Help free up memory
    del a_copies, t_copies, g_copies, c_copies
    seq = []
    for _i in range(length):
        idx = randint(0, len(unrand) - 1)  # noqa: S311
        seq.append(unrand[idx])
        del unrand[idx]
    return Seq("".join(seq))


def recouple(lst: list[str], scan_length: int = 8, start: str = "") -> str:
    """Program that will re-couple text from a list that has been broken into many overlapping parts."""
    # Get all the separate chunks of data in a list, no repeats.
    data = sorted(set(lst))
    if scan_length < 4:
        scan_length = min(round(min([len(i) for i in data]) / 2), 4)
        print(f"{scan_length = }")
    # Get all the words into a dictionary
    words_list = {}
    beginnings = {}
    # For each thing to test in the dictionary
    for test in data:
        # If we haven't seen this word before, set it as a possible
        # beginning value
        if test[0:scan_length] not in beginnings:
            beginnings[test[0:scan_length]] = True
        # For each index point our test string,
        for i in range(len(test) - scan_length):
            # Get the word it could be
            word = test[i : i + scan_length]
            # Add the word to the words dict, pointing to the last letter of
            # our test string.
            words_list[word] = test[i + scan_length]
            # If it's not the beginning
            if i > 0:
                beginnings[word] = False
    # If the starting word was not defined,
    if not start:
        # Set it to a random valid word
        for word in beginnings:
            if beginnings[word]:
                test = str(word)
                break
    else:
        # Otherwise, try to use it
        test = str(start)
        yes = False
        # If it's invalid, break
        if test not in words_list:
            for word in words_list:
                if test in word:
                    test = word
                    yes = True
                    break
        if not yes:
            raise ValueError(
                "Start word argument is invalid, does not exist in word list with scan_length. Try having start word len be equal to scan_length argument number.",
            )
    # Free up some memory
    del beginnings
    # Get data that points to each other
    string = test
    count = len(words_list.keys()) * scan_length
    while count:
        # If the test word is valid
        if test in words_list:
            # Add it to the string
            string += words_list[test]
            test = test[1:scan_length] + words_list[test]
            count -= 1
        else:  # Otherwise, decrement scan by one
            count = 0
    ##    # Once we have data that points to other bits, put them together properly
    ##    data = str(string[0])
    ##    # For each bit of data to add,
    ##    for i in string[1:]:
    ##        # For all possible index positions in the new bit of data,
    ##        for ii in range(len(i)):
    ##            # Get the current data's tail with ii
    ##            piece = data[-ii:]
    ##            # If the tail of this index is the start of this bit of data,
    ##            if piece == i[:len(piece)]:
    ##                # Combine the two chucks properly and stop checking indexes for this new chunk.
    ##                data = data[:len(data)-ii] + i
    ##                break
    return string


def get_nucleotide(seq: list[T], position: int) -> T | None:
    """Get nucleotide at specific position or None if out of range (just use get lol)."""
    if position < len(seq):
        return seq[position - 1]
    return None


def bwtrans(word: str, end: str = "$") -> str:
    """Return Burrows-Wheeler transform of string."""
    word += end
    words = sorted([word[i:] + word[:i] for i in range(len(word))])
    return "".join([words[i][-1] for i in range(len(words))])


def revbwtrans(string: str, end: str = "$") -> str:
    """Return inverse Burrows-Wheeler transform of string."""
    s = len(string)
    table = [""] * s
    for i in range(s):
        table = sorted(string[i] + table[i] for i in range(s))
    return next(row for row in table if row.endswith(end))
