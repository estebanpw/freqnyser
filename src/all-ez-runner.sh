#!/bin/bash

if [ $# != 1 ]; then
    echo "***ERROR*** Use: $0 <max_depth>"
    exit -1
fi


add_dimension () { #$1 is nucleotide string

	echo "Entered with $1"

	for k in ${nucl[@]};
	do
	
		local _word=$1$k
		local _len=${#_word}

		if [ $_len -gt $maxDepth ];
		then
			return
		fi
		
		local _result=$(../bin/freqgen -query /home/estebanpw/data/chromosomes/homo_sapiens/HOMSA.Chr.$chromo.fasta -kmer $_len -word $_word -out $var -hotspots ../../chr$chromo | grep "Average" | sed 's/,//g')
        local _formatted=(${_result// / })
        local _distrAvg=${_formatted[3]}
        local _hotsAvg=${_formatted[7]}

		local _porcentualDiff="$(awk "BEGIN{ print ("$_hotsAvg" - "$_distrAvg")/("$_distrAvg")*100 }"  )"
        local _significative="$(awk "BEGIN{ if("$_porcentualDiff" > "${depths["$_len"]}") print 1; else print 0; }")"
        echo "$_word produces: ${_formatted[3]} ${_formatted[7]} ${_formatted[10]} $_porcentualDiff $_significative"

			
		add_dimension $_word



	done
}

maxDepth=$1

nucl=(A C G T)
depths=(0 0 0 0 0 0 0 0 0 0 0 0 0 0)

#echo ${nucl[@]}

#echo ${nucl[1]}

#make all 

# Average distribution = 0.0441287, Average hotspots = 0.0498741, Variance = 0.000664325
#a=100
#val="$(awk "BEGIN{print ("$a")/2}")"
#echo "$val"


chromo="16"
var="homsa_${chromo}_$kmer.dat"

for i in ${nucl[@]};
do	

	for j in ${nucl[@]};
	do

		str=$i$j
		len=${#str}
		result=$(../bin/freqgen -query /home/estebanpw/data/chromosomes/homo_sapiens/HOMSA.Chr.$chromo.fasta -kmer $len -word $str -out $var -hotspots ../../chr$chromo | grep "Average" | sed 's/,//g')
		formatted=(${result// / })
		distrAvg=${formatted[3]}
		hotsAvg=${formatted[7]}
		

		porcentualDiff="$(awk "BEGIN{ print ("$hotsAvg" - "$distrAvg")/("$distrAvg")*100 }"  )"
		significative="$(awk "BEGIN{if("$porcentualDiff" > 0) print 1; else print 0; }")"
		echo "$str produces: ${formatted[3]} ${formatted[7]} ${formatted[10]} $porcentualDiff $significative"


		currNucl=$str
		add_dimension $currNucl


	done

done

