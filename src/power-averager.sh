

if [ $# != 2 ]; then
	echo "***ERROR*** Use: $0 <vector point/average file> <times>"
	exit -1
fi

vector=$1
times=$2


awk -v times="$times" '

BEGIN{ n=0; }

{ 

sum = sum + $1;

if(n == times){
	print sum/n;
	n = 0;
	sum = 0;
}

n = n + 1;

}

END{

if(n > 0){

	print sum/n;

}

}


'  $vector

