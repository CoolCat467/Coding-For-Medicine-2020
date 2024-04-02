"""Gets specific nucleotides in DNA from a file."""

import os

from Bio import SeqIO
from Bio.Seq import Seq

# from threading import Thread

__title__ = "Get Nucleotide"
__version__ = "0.0.1"


def get_nucleotide(seq, position):
    """Get nucleotide at position or None if out of range."""
    if position < len(seq):
        return seq[position - 1]
    return None


def run():
    """Run test."""
    file = "/share/Human/chr12.fa"
    file_type = "fasta"

    file = "SARS-2020.fasta"

    if os.sys.argv[1:]:
        if "-f" in os.sys.argv:
            idx = os.sys.argv.index("-f")
            file = os.sys.argv[idx + 1]
        if "-t" in os.sys.argv:
            idx = os.sys.argv.index("-t")
            file_type = os.sys.argv[idx + 1]

    print("Running...")
    dna_seq = SeqIO.read(file, file_type)
    dna_seq = dna_seq.seq[: -(len(dna_seq) % 3)]
    while True:
        position = input("position : ")
        if not position.isnumeric():
            break
        change = input("change : ")
        one = dna_seq.translate()
        seq2 = str(dna_seq)
        seq2[position] = change
        two = Seq(seq2).translate()
        change -= change % 3
        idx = int(change / 3)
        print(one[idx], two[idx])
        print(get_nucleotide(dna_seq, int(position)))
    print("Done")


if __name__ == "__main__":
    run()
