import argparse
import pathlib
import logging
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

taxonomic_level = args.l
path_fasta = pathlib.Path(args.p)
ambi_threshold = args.t
ncbi = NCBITaxa()
query_coverage = args.q
path_blast = pathlib.Path(args.b)
path_output = pathlib.Path(args.o)

path_output.joinpath("clean").mkdir(exist_ok=True)
path_output.joinpath("shrunk").mkdir(exist_ok=True)
path_output.joinpath("log").mkdir(exist_ok=True)

logging.basicConfig(filename=path_output.joinpath("log", "hitCheck.log"), level=logging.DEBUG)


def is_ambiguous(valid_blast_hits, query_level, ambi_threshold):
    if query_level in valid_blast_hits:

        x = valid_blast_hits.count(query_level)
        y = len(valid_blast_hits)

        if x/y < ambi_threshold:
            return True

    elif valid_blast_hits:
        return True
    else:
        return False


def get_desired_rank(new_hit, taxonomic_level):
    try:
        lineage = ncbi.get_lineage(new_hit)
        name = ncbi.get_taxid_translator(lineage)
        lineage_to_ranks = ncbi.get_rank(name)
        ranks_to_lineage = dict((rank, taxonomic_id) for (taxonomic_id, rank) in lineage_to_ranks.items())
        return name[ranks_to_lineage[taxonomic_level]], True
    except:
        return str(taxonomic_level) + " not found", False


def append_to_hits(encountered_taxonomic_ids, valid_blast_hits, new_hit):
    if new_hit in encountered_taxonomic_ids:
        valid_blast_hits.append(encountered_taxonomic_ids[new_hit])
    else:
        taxonomic_id_for_new_hit, found_taxonomic_id = get_desired_rank(new_hit, taxonomic_level)
        if found_taxonomic_id:
            valid_blast_hits.append(taxonomic_id_for_new_hit)
            encountered_taxonomic_ids[new_hit] = taxonomic_id_for_new_hit
    return encountered_taxonomic_ids, valid_blast_hits


def add_ambiguous_sequence(level_hits, sequences_to_remove, writer_shrunk, best_hit, split_line):
    taxonomic_string = '{0!s}'.format(Counter(level_hits).most_common())
    writer_shrunk.write(re.sub("\n", "", best_hit) + "\t" + taxonomic_string + "\n")
    sequences_to_remove.add(split_line[0])
    return level_hits, sequences_to_remove


def write_output(file, sequences_to_remove):

    with path_output.joinpath("clean", file.name).open("w") as write_fasta:
        fasta_sequences = SeqIO.parse(file.open("r"), "fasta")

        num_contigs_clean = 0
        num_contigs_split = 0
        num_contigs_ambiguous = 0

        length_split_genome = 0
        length_clean_genome = 0
        length_ambiguous_genome = 0

        for fasta_sequence in fasta_sequences:

            num_contigs_split += 1
            header, sequence = fasta_sequence.description, str(fasta_sequence.seq)
            length_split_genome += len(sequence)

            if header.split(" ")[0] not in sequences_to_remove:
                length_clean_genome += len(sequence)
                num_contigs_clean += 1
                write_fasta.write(">" + header + "\n")
                write_fasta.write(sequence + "\n")

            else:
                length_ambiguous_genome += len(sequence)
                num_contigs_ambiguous += 1

        write_overview.write(str.split(file.name, ".")[0] + "\t" + str(num_contigs_split) + "\t" +
                             str(length_split_genome) + "\t" + str(num_contigs_clean) + "\t" +
                             str(length_clean_genome) + "\t" + str(num_contigs_ambiguous) + "\t" +
                             str(length_ambiguous_genome) + "\n")


with path_output.joinpath("log", "overview.txt").open("w") as write_overview:

    logging.info("Running script in base folder " + str(path_fasta) + " with threshold " + str(ambi_threshold) +
                 ", taxonomic level " + taxonomic_level + " and minimal query coverage of " + str(query_coverage))

    directory = path_fasta.iterdir()
    encountered_taxonomic_ids = {}

    for file in directory:

        sequences_to_remove = set()
        init = False
        query_level = ""
        subject_level = ""

        if file.is_file() and (file.match("*.fna") or file.match("*.fa") or file.match("*.fasta")):

            logging.info("cleaning genome: " + file.name)

            #Do we have a corresponding blast file?
            if path_blast.joinpath(file.with_suffix(".blast").name).exists():

                logging.info("found blast file: " + str(path_blast.joinpath(file.with_suffix(".blast").name)))

                with path_blast.joinpath(file.with_suffix(".blast").name).open("r") as readBlast,\
                        path_output.joinpath("shrunk", file.with_suffix(".blast").name).open("w") as writer_shrunk:

                    last_contig = ""
                    valid_blast_hits = []
                    best_hit = ""
                    best_scoring_blast_hit = 0
                    no_contig_hit = ""

                    for line in readBlast.readlines():
                        split_line = str.split(line, "\t")
                        if not init:
                            query_level, found_taxonomic_id = get_desired_rank(str.split(split_line[0], "|")[1], taxonomic_level)
                            if found_taxonomic_id:
                                best_hit = line
                                best_scoring_blast_hit = split_line[3]
                                init = True
                            elif split_line[0] != no_contig_hit:
                                logging.info("No taxonomic id found for: " + split_line[1])
                                no_contig_hit = split_line[0]
                                continue

                        #We are only interested in blast hits with a coverage larger than the min query_coverage
                        if float(split_line[18]) >= query_coverage:

                            #Encountering new contig -> write out old contig
                            if split_line[0] != last_contig:

                                #Process data for last contig
                                #Check if we have a hit at the query_level
                                if is_ambiguous(valid_blast_hits, query_level, ambi_threshold):

                                    valid_blast_hits, sequences_to_remove = add_ambiguous_sequence(valid_blast_hits,
                                                                                                  sequences_to_remove,
                                                                                                  writer_shrunk,
                                                                                                  best_hit,
                                                                                                  split_line)

                                #Set data for new Contigs
                                last_contig = split_line[0]
                                best_hit = line
                                best_scoring_blast_hit = split_line[3]
                                valid_blast_hits = []

                            else:
                                if split_line[3] > best_scoring_blast_hit:
                                    best_hit = line
                                    best_scoring_blast_hit = split_line[3]

                            encountered_taxonomic_ids, valid_blast_hits = append_to_hits(encountered_taxonomic_ids,
                                                                                         valid_blast_hits, str.split(split_line[6], ";")[0])


                write_output(file, sequences_to_remove)

            #No blast file -> write output like input
            else:
                logging.warning("no blast file found")
                write_output(file, sequences_to_remove)
