#!/usr/bin/python3

import argparse
import pathlib
import re
from ete3 import NCBITaxa
from Bio import SeqIO
from collections import Counter


parser = argparse.ArgumentParser()
parser.add_argument("-l", type=str) # tax level
parser.add_argument("-p", type=str) # path fasta files (split)
parser.add_argument("-t", type=float) # threshold
parser.add_argument("-q", type=int) #min query coverage
parser.add_argument("-b", type=str) #blast files folder
parser.add_argument("-o", type=str) #outputpath
args = parser.parse_args()

level = args.l
pathFasta = pathlib.Path(args.p)
ambiThreshold = args.t
ncbi = NCBITaxa()
queryCoverage = args.q
pathBlast = pathlib.Path(args.b)
pathOutput = pathlib.Path(args.o)

pathOutput.joinpath("clean").mkdir(exist_ok=True)
pathOutput.joinpath("shrunk").mkdir(exist_ok=True)
pathOutput.joinpath("log").mkdir(exist_ok=True)



#Get genus Name
def getDesiredRanks(taxid, searchedRank):
    try:
        lineage = ncbi.get_lineage(taxid)
        name = ncbi.get_taxid_translator(lineage)
        lineage2ranks = ncbi.get_rank(name)
        ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
        return name[ranks2lineage[searchedRank]]
    except:
        return "Genus not found"


def appendToHits(encountered, hits, element):
    if element in encountered:
        hits.append(encountered[element])
    else:
        tmp = getDesiredRanks(element, level)
        hits.append(tmp)
        encountered[element] = tmp
    return encountered, hits

def addAmbiSequence(levelHits, seqToRemove, writerShrunk, bestHit, lineSplit):
    taxS = '{0!s}'.format(Counter(levelHits).most_common())
    writerShrunk.write(re.sub("\n", "", bestHit) + "\t" + taxS + "\n")
    seqToRemove.add(lineSplit[0])
    return levelHits, seqToRemove


def writeOutput(file, seqToRemove):
    with pathOutput.joinpath("clean", file.name).open("w") as writeFasta:

        fastaSequences = SeqIO.parse(file.open("r"), "fasta")

        numContigsClean = 0
        numContigsSplit = 0
        numContigsAmbi = 0

        lengthGenomeSplit = 0
        lengthGenomeClean = 0
        lengthContigsAmbi = 0

        for fastaSequence in fastaSequences:
            numContigsSplit += 1

            header, sequence = fastaSequence.description, str(fastaSequence.seq)
            lengthGenomeSplit += len(sequence)
            if header.split(" ")[0] not in seqToRemove:
                lengthGenomeClean += len(sequence)
                numContigsClean += 1
                writeFasta.write(">" + header + "\n")
                writeFasta.write(sequence + "\n")
            else:
                lengthContigsAmbi += len(sequence)
                numContigsAmbi += 1
        writeOverview.write(
            str.split(file.name, ".")[0] + "\t" + str(numContigsSplit) + "\t" + str(lengthGenomeSplit) + "\t" + str(
                numContigsClean) + "\t" + str(lengthGenomeClean) + "\t" + str(numContigsAmbi) + "\t" + str(
                lengthContigsAmbi) + "\n")


print("Running script in base folder " + str(pathFasta) + " with threshold " + str(
        ambiThreshold) + ", taxonomic level " + level + " and minimal query coverage of " + str(queryCoverage))


with pathOutput.joinpath("log", "overview.txt").open("w") as writeOverview:

    directory = pathFasta.iterdir()


    encounteredTaxIDs = {}
    for file in directory:
        print(file.name)
        seqToRemove = set()
        init = False
        queryLevel = ""
        subjectLevel = ""
        if file.is_file() and (file.match("*.fna") or file.match("*.fa") or file.match("*.fasta")):
            print(pathBlast.joinpath(file.with_suffix(".blast").name))
            #Do we have a corresponding blast file?
            if pathBlast.joinpath(file.with_suffix(".blast").name).exists():

                with pathBlast.joinpath(file.with_suffix(".blast").name).open("r") as readBlast, \
                        pathOutput.joinpath("shrunk", file.with_suffix(".blast").name).open("w") as writerShrunk:

                    lastContig = ""
                    levelHits = []
                    bestHit = ""
                    bestHitScore = 0
                    for line in readBlast.readlines():
                        lineSplit = str.split(line, "\t")
                        if not init:
                            queryLevel = getDesiredRanks((str.split(lineSplit[0], "|")[-1]), level)
                            bestHit = line
                            bestHitScore = lineSplit[3]
                            init = True

                        #We are only interested in blast hits with a coverage larger than the min querycoverage
                        if float(lineSplit[18]) >= queryCoverage:

                            #Encountering new contig -> write out old contig
                            if lineSplit[0] != lastContig:

                                #Process data for last contig
                                #Check if we have a hit at the querylevel
                                if levelHits.__contains__(queryLevel):

                                    if (levelHits.count(queryLevel) / len(levelHits)) < ambiThreshold:
                                        levelHits, seqToRemove = addAmbiSequence(levelHits, seqToRemove, writerShrunk,
                                                                                 bestHit, lineSplit)

                                #ambi contig, because list is not empty and no hit a query level
                                elif levelHits:
                                    levelHits, seqToRemove = addAmbiSequence(levelHits, seqToRemove, writerShrunk,
                                                                             bestHit, lineSplit)

                                #Set data for new Contigs
                                lastContig = lineSplit[0]
                                bestHit = line
                                bestHitScore = lineSplit[3]
                                levelHits = []
                                encounteredTaxIDs, levelHits = appendToHits(encounteredTaxIDs,
                                                                            levelHits, str.split(lineSplit[6], ";")[0])

                            else:
                                if lineSplit[3] > bestHitScore:
                                    bestHit = line
                                    bestHitScore = lineSplit[3]
                                    encounteredTaxIDs, levelHits = appendToHits(encounteredTaxIDs,
                                                                                levelHits, str.split(lineSplit[6], ";")[0])
                                else:
                                    encounteredTaxIDs, levelHits = appendToHits(encounteredTaxIDs,
                                                                                levelHits, str.split(lineSplit[6], ";")[0])

                writeOutput(file, seqToRemove)

            #No blast file -> write output like input
            else:
                print("no blast file")
                writeOutput(file, seqToRemove)
