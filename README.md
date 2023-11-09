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
```