
#To use with gnuplot

if [ $# != 4 ]; then
	echo "***ERROR*** Use: $0 <start> <end> <vector point file> <kmer>"
	exit -1
fi

start=$1
end=$2
vector=$3
kmer=$4

len=$(expr $2 - $1)

echo "Region length: $len"


tail -n +$start $vector | head -n $len | awk -v start="$start" 'BEGIN{x=0;} {print start+x, $1 ; x = x + 1; }' > $vector.plot

# Get how many kmers there are

words=$(awk 'BEGIN{x=0;} {if($2==1) x = x+1;} END{print x}' $vector.plot)

echo "The kmer is found $words times"

# For the cumulative plot

awk -v total=$words 'BEGIN{x=0;} {x=x+$2; print $1,x/total;}' $vector.plot > $vector.cuplot



gnuplot <<- EOF
		set key left top
        set ylabel "Cumulative percentage of the k-mer found in the region [0-1]"
        set xlabel "Coordinates"
        set title "Cumulative plot of kmer $kmer"   	
		
        plot "${vector}.cuplot" with lines
		pause mouse "Click any mouse button on selected data point"

EOF


echo "Closed plot"

rm $vector.plot
rm $vector.cuplot
