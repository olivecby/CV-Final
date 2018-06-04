#!/bin/sh
counter=1
for f in *.jpg
 	do 
 		mv $f replace$counter.jpg;
 		counter=$counter + 1;
done