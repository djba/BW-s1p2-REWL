#!/bin/bash
for i in {1..24}
do
	../bw-s1p2-gpu $i > 36-$i.out
done
