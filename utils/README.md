### check_search.sh
Take a fasta, index it with rlcsa and FMD. Then simulate 1000 reads on forward strand with no errors (using dwgsim). Run exact search on the two indices and check the outputs (diff). Iterate 100 times.

Hardcoded parameters: number of threads, reads, and iterations.

``` sh
# install dwgsim
# move to this folder (path to PANPP is hardcoded to ../PANPP)
bash check_search.sh $FA $OUTPREFIX
```

### check_search_big.sh
Same a check_search.sh but we do not index. We take indices as input.

Hardcoded parameters: number of reads and iterations.

``` sh
# install dwgsim
# move to this folder (path to PANPP is hardcoded to ../PANPP)
bash check_search_big.sh FA RLCSAPREFIX FMD OUTPREFIX
```

### compare_sfs.py
Compare this pingpong output and original pingpong output and check if are the same.

``` sh
python3 compare_sfs.py [SFS] [SFS]
```
