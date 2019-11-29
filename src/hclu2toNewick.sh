# NOTICE

FILE=$1

awk '

BEGIN{ total = 0; }

{

	# Store values in table and go to END statement when finding a space
	if($1 == "") exit
	table[$1] = $3 " " $4
	total = total + 1

}

func recursive(array, line,    result){
	
	split(line, result, " ")

	if(result[1] > total) {
		printf("(");
		recursive(array, array[result[1]])
		printf(")");
	}else{
		printf("%d,", result[1])
	}



	if(result[2] > total) {
		printf("(");
	       	recursive(array, array[result[2]])
		printf(")");
	}else{
		printf(",%d,", result[2])
	}

}

END{
	LAST=total*2
	recursive(table, table[LAST])
	printf("\n");

}


' $FILE | sed 's/,)/)/g' | sed 's/,,/,/g' 
