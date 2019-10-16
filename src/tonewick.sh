#!/bin/bash

if [ $# != 1 ]; then
    echo "***ERROR*** Use: $0 <ez-results file>"
    exit -1
fi

FILE=$1


grep -v "Entered" $FILE | awk '{if($7==1) print $1, $6}' | awk '

BEGIN{ prevl = 0; nrow = 0; }

{
	l = length($1)
	
	if(l > prevl){

		if(nrow == 0) printf("(%s:%f", $1, $2);
		if(nrow > 0) printf(",(%s:%f", $1, $2);

	}

	if(l == prevl){

		printf(",%s:%f", $1, $2);

	}

	if(l < prevl){

		for (i = l; i < prevl; i++){
			printf(")");
		}
		printf(",%s:%f", $1, $2);

	}

	prevl = l;
	nrow = nrow + 1;
}

END{

	printf(");\n");

}

'


