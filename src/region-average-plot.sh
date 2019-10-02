
#To use with gnuplot

if [ $# != 4 ]; then
	echo "***ERROR*** Use: $0 <start> <end> <vector average file> <kmer>"
	exit -1
fi

start=$1
end=$2
vector=$3
kmer=$4

len=$(expr $2 - $1)

echo "$len"

# For the density plot


tail -n +$start $vector | head -n $len | awk -v start="$start" 'BEGIN{x=0;} {print start+x, $1 ; x = x + 1; }' > $vector.plot


gnuplot <<- EOF
        set ylabel "Density [0-1]"
        set xlabel "Coordinates"
        set title "Density plot of kmer $kmer"   	
		
        plot "${vector}.plot" with lines
		pause mouse "Click any mouse button on selected data point"

EOF


echo "Closed plot"

rm $vector.plot
