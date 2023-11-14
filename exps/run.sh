#!/bin/sh

set -e

FA=$1
AGC=$2
WD=$3
threads=$4

Ns=(4) # (1 2 4 8 16 32)
SEED=23
COV=(5) # 10 15 30) # coverage

WD=$WD/$SEED

mkdir -p $WD

# Simulate reads from reference
mkdir -p $WD/reads
for cov in "${COV[@]}"
do
    pbsim --model_qc models/model_qc_ccs --data-type CCS --length-min 7500 --length-max 12500 --length-mean 10000 \
	  --accuracy-min 0.98 --accuracy-max 1 --accuracy-mean 0.99 --accuracy-sd 0.01 \
	  --depth $cov --seed $SEED $FA --prefix $WD/reads/${cov}x
    mv $WD/reads/${cov}x_0001.fastq $WD/reads/${cov}x.fq 
done

# Get random list of haplotypes
mkdir -p $WD/seqs
agc listset $AGC | (RANDOM=$SEED; while read line ; do echo "$RANDOM $line" ; done ) | sort | head -n ${Ns[${#Ns[@]}-1]} | cut -f2 -d' ' > $WD/hlist.txt

# Get sequences of selected haplotypes
while read line
do
    agc getset $AGC $line > $WD/seqs/$line.fa
done < $WD/hlist.txt

# Prepare input haplotypes based on Ns
for n in "${Ns[@]}"
do
    mkdir -p $WD/$n/seqs
    head -n $n $WD/hlist.txt | while read idx 
    do
 	ln -sf $WD/seqs/$idx.fa $WD/$n/seqs/$idx.fa
    done 
    cat $WD/$n/seqs/*.fa > $WD/$n/seqs.fa
done

# Indexing
for n in "${Ns[@]}"
do
    \time -vo $WD/$n/PANPP-index.time ../PANPP index -@ $threads -i $WD/$n/index $WD/$n/seqs/*.fa
    \time -vo $WD/$n/SVDSS-index.time SVDSS index --fasta $WD/$n/seqs.fa --index $WD/$n/index.fmd
done

# Search
for n in "${Ns[@]}"
do
    for cov in "${COV[@]}"
    do
	mkdir -p $WD/$n/${cov}x/
	\time -vo $WD/$n/${cov}x/PANPP-search.time ../PANPP search $WD/$n/index $WD/reads/${cov}x.fq > $WD/$n/${cov}x/PANPP.sfs
	\time -vo $WD/$n/${cov}x/SVDSS-search.time SVDSS search --index $WD/$n/index.fmd --fastq $WD/reads/${cov}x.fq --workdir $WD/$n/${cov}x/SVDSS --threads $threads
	cat $WD/$n/${cov}x/SVDSS/*.sfs > $WD/$n/${cov}x/SVDSS.sfs

	python3 compare_sfs.py $WD/$n/${cov}x/PANPP.sfs $WD/$n/${cov}x/SVDSS.sfs > $WD/$n.diff
    done
done
