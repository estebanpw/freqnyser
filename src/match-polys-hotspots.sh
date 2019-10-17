HOTSPOTS=$1
POLYS=$2

grep ">" $POLYS | sed 's/> //g' | sed 's/,/ /g' > $POLYS.temp

awk '

NR==1{ x=0; print "polyStart polyEnd hotspotStart hotspotEnd" }

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



FNR==NR{
	
	h_start[x] = $2
	h_end[x] = $3
	x = x + 1;
	next

}

FNR!=NR{

	island_start = $1
	island_end = $3

	closest = BinarySearch(h_start,length(h_start),island_start);

	if(island_start <= h_end[closest] && h_start[closest] <= island_end){
		print island_start,island_end,h_start[closest],h_end[closest]
	}
	closest = closest - 1;
		
	if(closest >= 0 && island_start <= h_end[closest] && h_start[closest] <= island_end){
		print island_start,island_end,h_start[closest],h_end[closest]
	}
	closest = closest + 2;

	if(closest < x && island_start <= h_end[closest] && h_start[closest] <= island_end){
		print island_start,island_end,h_start[closest],h_end[closest]
	}

	
	next

}

' $HOTSPOTS $POLYS.temp


rm $POLYS.temp
