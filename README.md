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

./PANPP fmdindex data/tiny.fa
./PANPP fmdexact data/tiny.fa.fmd data/tiny.fq
```

### TODO
- [ ] FMD-index construction
