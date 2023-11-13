# panpp
Specific strings from multiple genomes

``` sh
git clone git@github.com:ldenti/panpp.git
cd panpp
mkdir build
cd build
cmake ..
make
cd ..
```

``` sh
./PANPP index -i INDEX -@4 data/tiny.fa
./PANPP search -v INDEX AGA
./PANPP exact INDEX data/tiny.fq

1: 4,4
2: 4,4
3: 2,4
4: 5,5
5: 14,15

./PANPP fmdindex data/tiny.fa
./PANPP fmdexact data/tiny.fa.fmd data/tiny.fq
```

### TODO
- [X] ~FMD-index construction~
- [X] queries on FMD
- [ ] start from given index
- [ ] if partial indexes are already present, don't compute them
- [ ] ping-pong
