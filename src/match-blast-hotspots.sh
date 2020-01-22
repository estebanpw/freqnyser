HOTSPOTS=$1
BLAST=$2

grep -v "#" $BLAST > $BLAST.temp

awk '

NR==1{ x=0; print "blastStart blastEnd hotspotStart hotspotEnd" }

func BinarySearch(array,len,target){  


 	low=0;
 	high=len;
 	while(low<high){
             	    
        mid=int((low+high)/2);
             	    
        if(array[mid]==target) return mid;             		
        
        else if(array[mid]>target)high=mid-1;
        
        else low=mid+1;
    }
    return high;
}

func calculate_overlap(x1,x2,y1,y2){

	
	lower = (x1 <= y1) ? y1 : x1;
	higher = (x2 <= y2) ? x2 : y2;

	return (higher-lower)
}

FNR==NR{
	
	h_start[x] = $2
	h_end[x] = $3
	x = x + 1;
	next

}

FNR!=NR{

	island_start = $9
	island_end = $10

	closest = BinarySearch(h_start,length(h_start),island_start);

	if(island_start <= h_end[closest] && h_start[closest] <= island_end){
		if(island_start < island_end){
			overlap = calculate_overlap(island_start,island_end,h_start[closest],h_end[closest])
			print island_start,island_end,h_start[closest],h_end[closest],overlap
		}else{
			overlap = calculate_overlap(island_end,island_start,h_start[closest],h_end[closest])
            print island_start,island_end,h_start[closest],h_end[closest],overlap
		}
	}
	closest = closest - 1;
		
	if(closest >= 0 && island_start <= h_end[closest] && h_start[closest] <= island_end){
		if(island_start < island_end){
            overlap = calculate_overlap(island_start,island_end,h_start[closest],h_end[closest])
            print island_start,island_end,h_start[closest],h_end[closest],overlap
        }else{
            overlap = calculate_overlap(island_end,island_start,h_start[closest],h_end[closest])
            print island_start,island_end,h_start[closest],h_end[closest],overlap
        }

	}
	closest = closest + 2;

	if(closest < x && island_start <= h_end[closest] && h_start[closest] <= island_end){
		if(island_start < island_end){
            overlap = calculate_overlap(island_start,island_end,h_start[closest],h_end[closest])
            print island_start,island_end,h_start[closest],h_end[closest],overlap
        }else{
            overlap = calculate_overlap(island_end,island_start,h_start[closest],h_end[closest])
            print island_start,island_end,h_start[closest],h_end[closest],overlap
        }

	}

	
	next

}

' $HOTSPOTS $BLAST.temp


rm $BLAST.temp
