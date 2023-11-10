#!/bin/sh

fa=$1 # Input fasta (all sequences in the indices)
rlcsa=$2 # rlcsa index
fmd=$3 # fmd index
wd=$4 # Output directory

NR=1000 # Number of Reads
n=100 # Number of iterations

mkdir -p $wd

for i in $(seq $n)
do
    mkdir -p $wd/$i

    dwgsim -e 0 -E 0 -N $NR -1 100 -2 100 -r 0 -F 0 -R 0 -X 0 -y 0 -A 1 $fa $wd/$i/PREFIX
    zcat $wd/$i/PREFIX.bwa.read1.fastq.gz > $wd/$i/queries.fq
    rm $wd/$i/PREFIX.*

    ../PANPP exact $rlcsa $wd/$i/queries.fq > $wd/$i/rlcsa.res
    ../PANPP fmdexact $fmd $wd/$i/queries.fq > $wd/$i/ropebwt2.res

    diff $wd/$i/rlcsa.res $wd/$i/ropebwt2.res > $wd/$i.diff
done
