# freqnyser
An frequency analysis suite for DNA sequences

Freqnyser is composed of several programs which can do a lot for you! 
Find below descriptions of each program as well as examples of the command line usage:


### one-motif-finder
Finds a given motif (sequence pattern) in a fasta file. NOTICE: use letter 'N' for any letter match and it will not impact percentage of identity. You can also vary the minimum percentage of identity with -pident. E.g.:

>bin/one-motif-finder -query FASTA -out FILEOUTPUT -motif CCGCCGTNNCCNC -pident 0.75

### consecutive-finder
This program finds all sequences of a given polymer that is repeated 3 times consecutively with a maximum distance between each of them
For instance:


>./consecutive-finder -query HOMSA.Chr.16.fasta -motif AAAAAAAAAAAA -pident 0.8 -out output_results


This program will find all regions which have motif "AAAAAAAAAAAA" in sequence HOMSA.Chr.16.fasta with a minimum percentage of identity of 80% and write the results to "output_results"

### find-all-consecutive
This script is a wrapper on top of consecutive-finder which iterates over a given sequence to generate all unique kmers and then run the previous algorithm for each one

>./find-all-consecutive.sh seq ksize

example: 

>./find-all-consecutive.sh HOMSA.Chr.16.fasta 14

This will append all results from consecutive-finder into a file named <sequence>-all-motifs.txt

### match-polys-hotspots

This script will match hotspots coordinates to coordinates of polymer regions (any overlapped ones will be reported). Use it with the output of consecutive-finder. 

>./match-polys-hotspots.sh hotspots file consecutive-finder output

E.g.

>./match-polys-hotspots.sh chr16 polyA-3-20-HOMSA-12.fasta


### Freqgen

With freqgen you will be able to perform two tasks (1) search a particular kmer without hotstpots, in which case the algorithm will automatically detect which regions present an increase in the kmer frequency or (2) search a particular kmer within hotspots, case in which hotspots with higher frequency of the kmer will be reported.

Working mode type ONE:

>./freqgen -query sequence.fasta -kmer kmersize -word actualkmer -out outputFile

e.g.

>./freqgen -query /home/estebanpw/data/chromosomes/homo_sapiens/HOMSA.Chr.16.fasta -kmer 6 -word AAAAAA -out HOMSA-16-k6-AAAAAA.dat

Output will be:
In terminal, stdout will show the top50 results, e.g.:

[INFO] Top 50 results:

 R: 21901757, 21902325 - L 568 - S (norm) 0.096831 Ratio 7.99158 Periods 10.3273/82.5312
 
 Which is interpreted as "Coordinates between 21901757, 21902325 of length 568 have a normalised score of 0.096831 (this is the average of the kmer in the region), which means that the frequency of the kmer is 7.99158 times higher than on average, and that if on average you will find the kmer every 10.3273 base pairs, in such region you find it every 82.5312 base pairs"
 
 If there are not enough significative regions reported on stdout you can try lowering the lambda parameter with -lambda which by default is 4.
 
 Output files:
 
 1. Instead of just the top 50 results, all results are written to the file "outputFileName_density.vec"
 2. A file "outputFileName_raw.vec" which has N rows, where N is the length of the input sequence and for every row there is a zero if the kmer is NOT at the given position and a one if it is
 3. A file "outputFileName_average.vec" which has N rows as well, but instead of a one or zero it contains the average kmer frequency on the surrounding 100 base pairs
 4. A file named "outputFileName" (not format) containing (a) selected kmer, (b) total number of times found in the sequence, (c) length of the sequence, (d) probability of finding that kmer if uniform positions are assumed, (e) list of bins (histogram) of the number of kmers found throughout every bin; this is terminated by command "EOP" then (f) the same information as in "outputFileName_density.vec"
 
As a side note, you can plot the histogram contained in (4) by using "python3 simpleplot.py outputFileName"  and it will generate the plot "outputFileName.png"


 Working mode type TWO (i.e. include a hotspots file):
 
>./freqgen -query sequence.fasta -kmer kmersize -word actualkmer -out outputFile -hotspot hotspot 
 
 This will generate only one file "outputFileName_hotspots.vec" containing the frequencies of the kmer for each hotspot.
 
 ### region-average-plot.sh
 
 This script will generate a plot view of the frequency of a kmer given a region. It uses the file "outputFileName_average.vec" from freqgen working mode ONE. Just use any starting coordinates or use those reported in the Top 50 for instance. Usage:
 
 >./region-average-plot.sh start end vectorAverageFile kmer
 
 It will generate a GNU plot and show it in the X11 terminal view.
 
 ### region-cumulative-plot.sh
 
 Similar to the previous one but showing the cumulative frequencies.
 
 >./region-cumulative-plot.sh start end vectorRawFile kmer
 
 ### region-DNA-extract.sh
 
 Extracts a region from a DNA sequence fasta file (IT CAN ONLY HAVE ONE ">" sequence !!!!!!!!!!). Use it as follows:
 
 >./region-DNA-extract.sh start end DNAsequence
 
 ### power-averager.sh
 
 Re-averages a raw or average file generated by freqgen in working mode one. Run it with either file and the length of the smoothing average. This is useful if you want to plot for instance a full average file which is too long; then you can average it several times and it will reduce in size by factor TIMES given as parameter. See usage:
 
 >./power-averager.sh vectorPoint/averageFile times
 
 ### ez-runner.sh and all-ez-runner.sh
 
 You shouldnt use these scripts: they generate trees of significative kmers recursively by using the porcentual difference.
 If you were to use them, use only ez-runner.sh with parameters, e.g.:
 
 >./ez-runner.sh 30 10 sequence.fasta
 
 Redirect the stdout and use the script tonewick.sh to get a tree representation.
 
 ### massive-gen
 
 This program generates tables of kmer indices. So if you want to see the frequency of all kmers, from k=1 to k=12 for instance you just run:
 
 >./massive-gen -query sequenceFasta -kmer 12 -out fulldistr-12.kmers
 
 You can run this program for a full sequence such as a chromosome, and then also for hotspots only (use hotspot-extractor, see below) and then compare them with compare-gen (see below as well)
 
 ### hotspot-extractor
 
 Self explanatory, extracts the DNA sequence corresponding to hotspots (ONLY SEQUENCES WITH ONE HEADER ">" !!!!!!!). Usage:
 
 >./hotspot-extractor.sh sequence hotspotsFile
 
 
 ### compare-gen
 
 This program will compare two outputs from massive-gen, i.e. will compare the tables for two files. It is aimed at comparing nohotspots vs hotspots tables. It will write an output file containing which kmers are significant given a porcentual difference filter. Use:
 
 
 >./compare-gen noHotspotMers hotspotMers ksize output minDifPor

Where the two "noHotspotMers" and "hotspotMers" files are generated by massive-gen. Ksize should be the same value as used in generation. To alleviate workflow execution, use the script below which runs massive gen and everything.

### generate-all.sh

Runs massive-gen and then compare-gen easily. Just remember that the hotspots file here corresponds to the extracted DNA sequences produced by hotspot-extractor

>./generate-all.sh chromo hotspots k

### runAllSplits.sh

This script finds a given motif in (1) the whole chromosomes of human, (2) the hotspots of human and (3) the hotspots of human AGAIN but counting hotspots only once. Thus it generates a few files (.motifs) that contain lines indicating how many times the motif was found. 
Use it as follows:

>./runAllSplits.sh MOTIF
