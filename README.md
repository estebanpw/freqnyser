# freqnyser
An frequency analysis suite for DNA sequences

Freqnyser is composed of several programs which can do a lot for you! 
Find below descriptions of each program as well as examples of the command line usage:

### consecutive-finder
This program finds all sequences of a given polymer that is repeated 3 times consecutively with a maximum distance between each of them
For instance:


../bin/consecutive-finder -query HOMSA.Chr.16.fasta -motif AAAAAAAAAAAA -pident 0.8 -out output_results


This program will find all regions which have motif "AAAAAAAAAAAA" in sequence HOMSA.Chr.16.fasta with a minimum percentage of identity of 80% and write the results to "output_results"

### find-all-consecutive
This script is a wrapper on top of consecutive-finder which iterates over a given sequence to generate all unique kmers and then run the previous algorithm for each one

./find-all-consecutive.sh <seq> <ksize>

example: ./find-all-consecutive.sh HOMSA.Chr.16.fasta 14

This will append all results from consecutive-finder into a file named <sequence>-all-motifs.txt

### match-polys-hotspots

This script will ...to be continued
