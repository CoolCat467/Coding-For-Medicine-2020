#!/usr/bin/env python3
# Reads a file and applies all three frame shifts to
# it for all possible proteins a DNA sequence can
# code for. Programmed by Samuel Davenport

# from Bio.Seq import Seq
# from Bio.Alphabet import IUPAC
from __future__ import annotations

import sys
from typing import TYPE_CHECKING

from Bio import SeqIO
from Bio.Blast import NCBIWWW

if TYPE_CHECKING:
    from Bio.Seq import Seq

__title__ = "Protien Shift"
__version__ = "0.0.1"


# types: note: "get_frames" defined here
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


def search_blast(protien: str, hit_count: int = 50) -> list[NCBIWWW]:
    with NCBIWWW.qblast(
        "blastp",
        "nr",
        protien,
        hitlist_size=int(hit_count),
        format_type="HTML",
    ) as result_handle:
        ##        with open("my_blast.xml", "w") as save_file:
        ##            data = result_handle.read()
        ##            # text = data.split('<Iteration>')[1].split('</Iteration_hits>')[0]
        ##            # text = ' '.join([i for i in ' '.join([i for i in text.split('\n')]).split(' ') if i != ''])
        ##
        ##            save_file.write(data)
        data = result_handle.read()
    text = [i.split("</Hit_def>\n")[0] for i in data.split("</Hit_id>\n")][1:]
    return [i.split("  <Hit_def>")[1] for i in text]


def run() -> None:
    file = "/share/SARS/SARS-2020.fasta"
    ftype = "fasta"
    wfile = "output.txt"

    if sys.argv[1:]:
        if "-f" in sys.argv:
            idx = sys.argv.index("-f")
            file = sys.argv[idx + 1]
        if "-t" in sys.argv:
            idx = sys.argv.index("-t")
            ftype = sys.argv[idx + 1]

    print(f"Opening {file = } {ftype = }")

    ##    dna_seqs = SeqIO.parse(file, ftype)
    dna_seq = SeqIO.read(file, ftype)
    ##    for dna_seq in dna_seqs:
    protiens = get_frames(dna_seq, n=3)
    names = [
        "\n".join(protiens) + "\n\n"
        for names in [search_blast(prot, 5) for prot in protiens]
    ]
    # types: assignment error: Incompatible types in assignment (expression has type "TextIOWrapper", variable has type "str")
    with open(wfile, "w", encoding="utf-8") as file_handle:
        file_handle.write("\n".join(protiens) + "\n\n")
        file_handle.write(
            "\n".join([f"{idx}: {name}" for idx, name in enumerate(names)]),
        )
    print("Done")


if __name__ != "__main__":
    run()
