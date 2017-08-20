#!/usr/bin/python3

import sys
import pathlib
import re
from Bio import SeqIO


directory = pathlib.Path(sys.argv[1]).iterdir()

pathlib.Path(sys.argv[1]).joinpath("split").mkdir(exist_ok=True)
pathlib.Path(sys.argv[1]).joinpath("blast").mkdir(exist_ok=True)
pathlib.Path(sys.argv[1]).joinpath("blast/output").mkdir(exist_ok=True)

#Splits all genomes in a folder and writes them to a new fasta file
for file in directory:

    if file.is_file() and (file.match("*.fna") or file.match("*.fn") or file.match("*.fasta")):
        with file.parents[0].joinpath("split", file.with_suffix(".fna").name).open("w") as writeSplit,\
                file.parents[0].joinpath("blast", file.with_suffix(".fna").name).open("w") as writeCandidates:
            fastaSequences = SeqIO.parse(file.open("r"), "fasta")
            for fastaSequence in fastaSequences:
                header, sequence = fastaSequence.description, str(fastaSequence.seq)
                subsequence = re.split("N{10,}", sequence)
                count = 0
                for i in subsequence:
                    if len(i) >= 100:
                        count += 1
                        writeSplit.write(">" + str(count) + "_" + header + "\n")
                        writeSplit.write(i + "\n")
                        #Write out candidates for blast
                        if len(i) < 10000:
                            writeCandidates.write(">" + str(count) + "_" + header + "\n")
                            writeCandidates.write(i + "\n")
