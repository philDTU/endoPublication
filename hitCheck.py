
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

def getDesiredRank(taxid, searchedRank):
    try:
        lineage = ncbi.get_lineage(taxid)
        name = ncbi.get_taxid_translator(lineage)
        lineage2ranks = ncbi.get_rank(name)
        ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
        return name[ranks2lineage[searchedRank]], True
    except:
        return str(level) + " not found", False


def appendToHits(encountered, hits, element):
    if element in encountered:
        hits.append(encountered[element])
    else:
        queryLevel, found = getDesiredRank(element, level)
        if found:
            hits.append(queryLevel)
            encountered[element] = queryLevel
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


with pathOutput.joinpath("log", "log.txt").open("w") as writeLog, pathOutput.joinpath("log", "overview.txt").open("w") as writeOverview:

    writeLog.write("Running script in base folder " + str(pathFasta) + " with threshold " + str(ambiThreshold) +
                   ", taxonomic level " + level + " and minimal query coverage of " + str(queryCoverage) + "\n")
    writeLog.flush()

    directory = pathFasta.iterdir()
    encounteredTaxIDs = {}

    for file in directory:

        seqToRemove = set()
        init = False
        queryLevel = ""
        subjectLevel = ""
        if file.is_file() and (file.match("*.fna") or file.match("*.fa") or file.match("*.fasta")):

            writeLog.write("cleaning: " + file.name + "\n")
            writeLog.flush()
            #Do we have a corresponding blast file?
            if pathBlast.joinpath(file.with_suffix(".blast").name).exists():
                writeLog.write("found blast file: " + str(pathBlast.joinpath(file.with_suffix(".blast").name)) + "\n")
                writeLog.flush()
                with pathBlast.joinpath(file.with_suffix(".blast").name).open("r") as readBlast, \
                        pathOutput.joinpath("shrunk", file.with_suffix(".blast").name).open("w") as writerShrunk:

                    lastContig = ""
                    levelHits = []
                    bestHit = ""
                    bestHitScore = 0
                    noHitContig = ""
                    for line in readBlast.readlines():
                        lineSplit = str.split(line, "\t")
                        if not init:
                            #print(str.split(str.split(lineSplit[0], "|")[-1],"_",1)[1])
                            #print(lineSplit)
                            queryLevel, found = getDesiredRank(str.split(lineSplit[0], "|")[-1], level)
                            if found:
                                bestHit = line
                                bestHitScore = lineSplit[3]
                                init = True
                            elif lineSplit[0] != noHitContig:
                                writeLog.write("No match for contig " + lineSplit[1] + "\n")
                                writeLog.flush()
                                noHitContig = lineSplit[0]

                        #We are only interested in blast hits with a coverage larger than the min querycoverage
                        if float(lineSplit[18]) >= queryCoverage:

                            #Encountering new contig -> write out old contig
                            if lineSplit[0] != lastContig:

                                #Process data for last contig
                                #Check if we have a hit at the querylevel
                                if queryLevel in levelHits:
                                    if (levelHits.count(queryLevel) / len(levelHits)) < ambiThreshold:
                                        levelHits, seqToRemove = addAmbiSequence(levelHits, seqToRemove, writerShrunk,
                                                                                 bestHit, lineSplit)

                                #ambi contig, because list is not empty and no hit on query level
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
                writeLog.write("no blast file found" + "\n")
                writeLog.flush()
                writeOutput(file, seqToRemove)
