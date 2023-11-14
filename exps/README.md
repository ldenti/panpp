# Experiments

Script to:
* Simulate HiFi reads from reference genome
* Get x haplotypes from an AGC index
* Index the haplotypes (rlcsa and ropebwt2)
* Compute specific strings (rlcsa and ropebwt2)

``` sh
# T2T reference genome
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
gunzip chm13v2.0.fa.gz

# Assemblies from HPRC Year 1
wget https://zenodo.org/record/5826274/files/HPRC-yr1.agc

mamba create -c bioconda -c conda-forge -n panpp agc samtools SVDSS pbsim
mamba activate panpp

# Simulate HiFi reads
bash run.sh simulate chm13v2.0.fa OUT_DIR
# Extract random haplotypes
bash run.sh select HPRC-yr1.agc OUT_DIR
# Index with PANPP
bash run.sh index1 OUT_DIR 32
# Index with SVDSS
bash run.sh index2 OUT_DIR
# Index with PANPP
bash run.sh search1 OUT_DIR 32
# Index with SVDSS
bash run.sh search2 OUT_DIR 32
# Compare results
bash run.sh compare OUT_DIR
```
