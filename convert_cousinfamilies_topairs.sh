#!/bin/bash


## remove 'next' from the first line if you have no header 
awk 'NR==1{print "fid" FS "mom1" FS "kid1" FS "mom2" FS "kid2"; next}
{split($3, m, ","); split($4, k, ",");
for(i in k){
	for(j=i+1; j in k; j++){
		if(m[i]!=m[j]){
			print NR-1, m[i], k[i], m[j], k[j]
		}
	}
}
}' input_file_maternal.csv > output_file_maternal.csv

awk 'NR==1{print "fid" FS "dad1" FS "kid1" FS "dad2" FS "kid2"; next}
{split($3, m, ","); split($4, k, ",");
for(i in k){
	for(j=i+1; j in k; j++){
		if(m[i]!=m[j]){
			print NR-1, m[i], k[i], m[j], k[j]
		}
	}
}
}' input_file_paternal.csv > output_file_paternal.csv
