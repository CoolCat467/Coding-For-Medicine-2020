"""Read proteins.

Reads a file and applies all three frame shifts to
it for all possible proteins a DNA sequence can
code for.
"""

# Programmed by CoolCat467

from __future__ import annotations

from random import randint
from threading import Thread

from Bio.Seq import Seq

NAME = "Threaded Random DNA and Shifts to Proteins"
__version__ = "0.0.2"


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
    """Get random proper sequence."""
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


class GenerateProteins(Thread):
    """Generate Proteins in thread."""

    def __init__(
        self,
        times: int,
        atgc: list[float],
        sequence_size: int,
    ) -> None:
        """Initialize self."""
        super().__init__()
        self.times = int(times)
        self.atgc = atgc
        self.sequence_size = int(sequence_size)
        self.start()

    def run(self) -> None:
        """Run thread."""
        proteins = []
        for _i in range(self.times):
            seq = gen_random_proper_seq(self.sequence_size, *self.atgc)
            frames = get_frames(seq)
            # for protein in frames:
            #    proteins[len(protein)] = protein
            proteins += [len(protein) for protein in frames]
        self.data = proteins


def run() -> None:
    """Run protein thing."""
    print("Running...")
    gcrng = (20, 65, 5)
    threads = {}
    for gcperc in range(*gcrng):
        at = (100 - gcperc) / 100
        gc = gcperc / 100
        a, t, g, c = (at / 2, at / 2, gc / 2, gc / 2)
        threads[gcperc] = GenerateProteins(500, [a, t, g, c], 10000)
        print("Thread Started.")
    ##        #proteins = {}
    ##        proteins = []
    ##        for i in range(500):
    ##            seq = gen_random_proper_seq(10000, a, t, g, c)
    ##            frames = get_frames(seq)
    ##            #for protein in frames:
    ##            #    proteins[len(protein)] = protein
    ##            proteins += [len(protein) for protein in frames]
    ##        #max(proteins.keys())
    data = {}
    print("Waiting for threads to end.")
    while threads:
        for key in list(threads.keys()):
            thread = threads[key]
            if not thread.is_alive():
                print("Thread calculating", key, "% GC is Done.")
                data[key] = sum(thread.data) / len(thread.data)
                del threads[key]
    for key in data:
        print(
            "Average Protein Size was",
            data[key],
            "with",
            key,
            "% of random DNA made up of GC.",
        )
    print("Done")


if __name__ == "__main__":
    run()
