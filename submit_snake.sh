#!/bin/bash

j=$1
if [[ -d ".snakemake" ]];
then
	echo "SNAKEMAKE <CHANNI>"
	# rm -rf .snakemake
fi

if [[ ! -n ${SSH_CONNECTION} ]];
then
  snakemake --snakefile scripts/g4_snake.smk \
                --latency-wait 5 \
                --keep-going \
                --cores 1 \
                --keep-incomplete
else
  snakemake --use-conda \
        --snakefile scripts/g4_snake.smk \
	    --keep-incomplete \
	    --rerun-triggers mtime \
	    --keep-going \
	    --latency-wait 45 \
	    --cluster-config cluster_settings.yaml \
	    --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem={cluster.mem} -c {cluster.ncpus} --nodes={cluster.nodes} --account={cluster.account} -o jobOut/{cluster.jobName}-%j.out -J {cluster.jobName} -e jobOut/{cluster.jobName}-%j.err" -j $j
fi