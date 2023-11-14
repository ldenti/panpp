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
# Run ping-pong over input FASTA/FASTQ
./PANPP search INDEX FX
```
### Example
``` sh
./PANPP index example/tiny.fa -i example/tiny
./PANPP search example/tiny example/tiny.fq
# 6	0	5
# *	7	4
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
