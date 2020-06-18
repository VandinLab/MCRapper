#! /bin/sh

for ds in `ls *vars.sh`; do
	echo "### $ds"
	./experiment.sh run $ds
done
