# Experiments

Script to:
* Simulate HiFi reads from reference genome
* Get x haplotypes from an AGC index
* Index the haplotypes (rlcsa and ropebwt2)
* Compute specific strings (rlcsa and ropebwt2)

``` sh
mamba create -c bioconda -c conda-forge -n panpp agc samtools SVDSS pbsim
mamba activate panpp
bash run.sh [.fa] [.agc] [wd] [threads]
```
