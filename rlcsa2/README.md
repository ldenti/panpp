
```
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
cd ..

./RLCSA example/tiny.fa example/tiny.fq
```

```
# Against forward only
1: 4,4
2: 4,4
3: 2,4
4: 5,5
5: 14,15
6: 2,1

# Against FMD
1: 5,5
2: 5,5
3: 3,5
4: 6,6
5: 19,20
6: 3,2
```