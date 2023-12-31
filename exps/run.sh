#!/bin/sh

set -e

mode=$1

SEED=23
Ns=(4) # (1 2 4 8 16 32)
COV=(5) # 10 15 30) # coverage

function log {
    echo "[$(date)] $1"
}

if [[ $mode == "simulate" ]]
then
    FA=$2
    WD=$3

    mkdir -p $WD/reads
    for cov in "${COV[@]}"
    do
	if [ ! -f $WD/reads/${cov}x.fq ]
	then
	    log "Simulating ${cov}x sample"
	    pbsim --model_qc models/model_qc_ccs --data-type CCS --length-min 7500 --length-max 12500 --length-mean 10000 \
		  --accuracy-min 0.98 --accuracy-max 1 --accuracy-mean 0.99 --accuracy-sd 0.01 \
		  --depth $cov --seed $SEED $FA --prefix $WD/reads/${cov}x
	    cat $WD/reads/${cov}x_*.fastq > $WD/reads/${cov}x.fq
	else
	    log "Skipping simulation of ${cov}x sample"
	fi
    done
elif [[ $mode == "select" ]]
then
    AGC=$2
    WD=$3
    log "Extracting haplotypes"
    mkdir -p $WD/seqs
    np1=$((${Ns[${#Ns[@]}-1]} + 1))
    agc listset $AGC | (RANDOM=$SEED; while read line ; do echo "$RANDOM $line" ; done ) | sort | head -n $np1 | cut -f2 -d' ' > $WD/hlist.txt

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

    idx=$(tail -1 $WD/hlist.txt)
    ln -sf $WD/seqs/$idx.fa $WD/$idx.fa
elif [[ $mode == "index1" ]]
then
    WD=$2
    threads=$3
    for n in "${Ns[@]}"
    do
	if [ ! -f $WD/$n/index.rlcsa.array ]
	then
	    log "Indexing $n haplotypes"
	    \time -vo $WD/$n/PANPP-index.time ../PANPP index -@ $threads -i $WD/$n/index $WD/$n/seqs/*.fa
	else
	    log "Skipping indexing $n haplotypes"
	fi
    done
elif [[ $mode == "index2" ]]
then
    WD=$2
    for n in "${Ns[@]}"
    do
	if [ "$n" -gt "4" ]
	then
	    log "NOT indexing $n haplotypes"
	    continue
	fi
	if [ ! -f $WD/$n/index.fmd ]
	then
	    log "Indexing $n haplotypes"
	    \time -vo $WD/$n/SVDSS-index.time SVDSS index --fasta $WD/$n/seqs.fa --index $WD/$n/index.fmd
	else
	    log "Skipping indexing $n haplotypes"
	fi
    done
elif [[ $mode == "search1" ]]
then
    WD=$2
    threads=$3
    for n in "${Ns[@]}"
    do
	for cov in "${COV[@]}"
	do
	    if [ ! -f $WD/$n/${cov}x/PANPP.sfs ]
	    then
		mkdir -p $WD/$n/${cov}x/
		log "Searching from ${cov}x sample against $n haplotypes"
		\time -vo $WD/$n/${cov}x/PANPP-search.time ../PANPP search $WD/$n/index $WD/reads/${cov}x.fq -@ $threads > $WD/$n/${cov}x/PANPP.sfs
	    else
	        log "Skipping searching from ${cov}x sample against $n haplotypes"
	    fi
	done
    done
elif [[ $mode == "search11" ]]
then
    WD=$2
    threads=$3
    idx=$(basename $WD/*.fa .fa)
    for n in "${Ns[@]}"
    do
	if [ ! -f $WD/$n/$idx/PANPP.sfs ]
	then
	    mkdir -p $WD/$n/$idx/
	    log "Searching from $idx against $n haplotypes"
	    \time -vo $WD/$n/$idx/PANPP-search.time ../PANPP search $WD/$n/index $WD/$idx.fa -@ $threads > $WD/$n/$idx/PANPP.sfs
	else
	    log "Skipping searching from $idx against $n haplotypes"
	fi
    done
elif [[ $mode == "search2" ]]
then
    WD=$2
    threads=$3
    for n in "${Ns[@]}"
    do
	if [ "$n" -gt "4" ]
	then
	    log "NOT searching from ${cov}x sample against $n haplotypes"
	    continue
	fi
	for cov in "${COV[@]}"
	do
	    if [ ! -f $WD/$n/${cov}x/SVDSS.sfs ]
	    then
		mkdir -p $WD/$n/${cov}x/
		log "Searching from ${cov}x sample against $n haplotypes"
		\time -vo $WD/$n/${cov}x/SVDSS-search.time SVDSS search --index $WD/$n/index.fmd --fastq $WD/reads/${cov}x.fq --workdir $WD/$n/${cov}x/SVDSS --threads $threads --assemble
		cat $WD/$n/${cov}x/SVDSS/*.sfs > $WD/$n/${cov}x/SVDSS.sfs
	    else
	        log "Skipping searching from ${cov}x sample against $n haplotypes"
	    fi
	done
    done
elif [[ $mode == "search22" ]]
then
    WD=$2
    threads=$3
    idx=$(basename $WD/*.fa .fa)
    for n in "${Ns[@]}"
    do
	if [ ! -f $WD/$n/$idx/SVDSS.sfs ]
	then
	    mkdir -p $WD/$n/$idx/
	    log "Searching from $idx against $n haplotypes"
	    \time -vo $WD/$n/$idx/SVDSS-search.time SVDSS search --index $WD/$n/index.fmd --fastq $WD/$idx.fa --workdir $WD/$n/$idx/SVDSS --threads $threads --assemble
	    cat $WD/$n/$idx/SVDSS/*.sfs > $WD/$n/$idx/SVDSS.sfs
	else
	    log "Skipping searching from $idx against $n haplotypes"
	fi
    done
elif [[ $mode == "compare" ]]
then
    WD=$2
    threads=$3
    for n in "${Ns[@]}"
    do
	for cov in "${COV[@]}"
	do
	    python3 compare_sfs.py $WD/$n/${cov}x/PANPP.sfs $WD/$n/${cov}x/SVDSS.sfs > $WD/$n.diff
	done
    done
else
    echo "Unkown command"
fi
