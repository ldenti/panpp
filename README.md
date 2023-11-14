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
# Index one or more fasta
./PANPP index -i INDEX -@4 [FA...]
# Run ping-pong
./PANPP search INDEX FQ

# Search a string S against the index
./PANPP test -v INDEX S
# Search multiple strings (in FASTQ format)
./PANPP exact INDEX FQ

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
- [X] ~queries on FMD~
- [ ] start from given index
- [ ] if partial indexes are already present, don't compute them
- [X] ~ping-pong~
- [ ] parallel ping-pong
- [X] ~assembler~
- [X] ~flanking~
