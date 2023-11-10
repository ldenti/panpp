### check_search.sh
Take a fasta, index it with rlcsa and FMD. Then simulate 1000 reads on forward strand with no errors (using dwgsim). Run exact search on the two indices and check the outputs (diff). Iterate 10 times.

Hardcoded parameters: number of threads, reads, and iterations.

``` sh
# install dwgsim
# move to this folder (path to PANPP is hardcoded to ../PANPP)
bash check_search.sh $FA $OUTPREFIX
```