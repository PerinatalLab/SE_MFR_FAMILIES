#!/bin/bash


## remove 'next' from the first line if you have no header 
awk 'NR==1{print "fid" FS "mom1" FS "kid1" FS "mom2" FS "kid2"
	FS "dad1" FS "dad2" FS "birthyear1" FS "birthyear2" FS "dadage1"
	FS "dadage2" FS "momage1" FS "momage2"; next}
{split($3, m, ","); split($4, k, ","); split($5, d, ","); split($6, by, ",");
split($7, da, ","); split($8, ma, ",");
for(i in k){
	for(j=i+1; j in k; j++){
		if(m[i]!=m[j]){
			print NR-1, m[i], k[i], m[j], k[j],
			d[i], d[j], by[i], by[j],
			da[i], da[j], ma[i], ma[j]
		}
	}
}
}' input_file_maternal.csv > output_file_maternal.csv

awk 'NR==1{print "fid" FS "dad1" FS "kid1" FS "dad2" FS "kid2"
	FS "mom1" FS "mom2" FS "birthyear1" FS "birthyear2" FS "dadage1"
	FS "dadage2" FS "momage1" FS "momage2"; next}
{split($3, d, ","); split($4, k, ","); split($5, m, ","); split($6, by, ",");
split($7, da, ","); split($8, ma, ",");
for(i in k){
	for(j=i+1; j in k; j++){
		if(d[i]!=d[j]){
			print NR-1, d[i], k[i], d[j], k[j],
			m[i], m[j], by[i], by[j],
			da[i], da[j], ma[i], ma[j]
		}
	}
}
}' input_file_paternal.csv > output_file_paternal.csv
