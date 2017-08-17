# endoPublication

Metagenomic Approaches for Determining the Structure and Function of Complex Microbiomes

## Getting Started

These instructions will let you build a cleaned database for metagenomic classification with kraken/bracken

### Prerequisites

Both scripts are written in Python and need the BioPython and ete3 packages to run.

```
pip3 install biopython
pip3 install ete3
```

### Building the database

1. Download the desired reference genomes and manipulate the headers fitting for kraken. This can be done conveniently with the scripts from Mick Watson[https://github.com/mw55309/Kraken_db_install_scripts] or manually by following the instructions in the kraken manual [https://github.com/jenniferlu717/Bracken/blob/master/README]


2. Use the splitGenomes.py script. This will split sequences at stretches of at least 10 Ns and discards fragments shorter than 100 bases.
  
  Enter the folder with the your genomes and run the splitGenomes script.
  ```
  cd <folderWithGenomes>
  python3 splitGenomes.py <pathToGenomes>
  ```

  The script generates 3 folders in your current directory
  - split will contain the fasta files after splitting them.
  - blast will contain the fasta files with short fragments to blast against the nt database
  - blast/output will contain the output of blast from the blast against the nt database
 
 
3. Mask the genomes using ncbi-dustmasker and convert lower-case bases to N

  Enter the folder with the split genomes and dust the genomes.
  ```
  cd split
  for i in *.fna; do dustmasker -in $i -infmt fasta -outfmt fasta | sed -e '/>/!s/a\|c\|g\|t/N/g' > tmp; mv tmp $i; done
  ```


4. Use blastn with the short sequences file against the nt database. (It is important to install taxdb for blast: ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz)

  Enter the folder blast which contains the short fragments after splitting and use blast against the nt database
  ```
  cd blast
  for i in *.fna; do blastn -query $i -db <blastDB> -out output/${i//.fna/.blast} -outfmt '6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle pident length mismatch gapopen qstart qend sstart send qcovs' -evalue 1e-6;done
  ```

5. Use the hitCheck.py to filter out the sequences with a hit in a different genus by running.

  ```
  python3 hitCheck.py -l <taxaLevel> -b <pathToBlastFiles> -p <pathToSplitFastaFiles> -o <pathToWriteOutput> -t <threshold> -q <minQueryCoverage>
  ```

  hitCheck generates 3 new folders
   - shrunk will contain the best hit of the blast results together with a taxonomic string
   - log will contain a log file with information for each input file
   - clean will contain the cleaned genome files

  log file content
   - filename: name of processed file
   - numContigsSplit: number of contigs after splitting
   - lengthGenomeSplit: number of bases after splitting
   - numContigsClean: number of contigs after cleaning
   - lengthGenomeClean: number of bases after cleaning
   - numContigsAmbi: number of removed ambiguous contigs
   - lengthContigsAmbi: number of removed ambiguous bases


6. Use the fasta files in the 'clean' folder to build the kraken database as described in the manual at http://ccb.jhu.edu/software/kraken/MANUAL.html#custom-databases

7. Build the Bracken database as in the manual at https://github.com/jenniferlu717/Bracken/blob/master/README

## Author of scripts

* **Philipp Kirstahler** - [philDTU](https://github.com/philDTU)


## Acknowledgments

* Download scripts for genomes - Mike Watson - [mw55309](https://github.com/mw55309)
* Kraken - Derrick Wood - [DerrickWood](https://github.com/DerrickWood)
* Bracken - Jennifer Lu - [jenniferlu717](https://github.com/jenniferlu717)
