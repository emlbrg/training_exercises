from Bio import SeqIO
from Bio.Seq import Seq

from orffinder import orffinder

sequence = list(SeqIO.parse("test.fasta", "fasta"))[0]

out = orffinder.getORFProteins(sequence, return_loci=True, minimum_length=75, remove_nested=True, trim_trailing=False)

for orf in out:

    start = min(orf["start"], orf["end"])
    end = max(orf["start"], orf["end"])
    length = end - start

    # print(orf["trailing"])

    if orf["sense"] == "+":

        print("STRAND 1:")
        print(sequence.seq[start - 1 : end - 1])

    else:
        print("STRAND -1:")
        print(sequence.seq[start - 1 : end - 1])

    print(orf["protein"])
    print("\n")

print(len(out))
