
#To use with gnuplot

if [ $# != 3 ]; then
	echo "***ERROR*** Use: $0 <start> <end> <DNA sequence>"
	exit -1
fi

start=$1
end=$2
vector=$3

len=$(expr $2 - $1)

#echo ">$len"

# For the seq

grep -v ">" $vector | tr -d '\n' | tail -c +$start | head -c $len

echo "" 
