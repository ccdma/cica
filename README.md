## compile

```
$ g++ ica.cpp -I ./include/ -fopenmp  
```
- `-O3 -mtune=native -march=native`は最適化系のオプション
- `-lblas`でblas利用（eigenのinclude前に`#define EIGEN_USE_BLAS`すること）
  - https://eigen.tuxfamily.org/dox/TopicUsingBlasLapack.html

## exec

```
$ OMP_NUM_THREADS=6 ./a.out
```

## ssh

```
$ scp -Cr calgorithm b36697@cinnamon.kudpc.kyoto-u.ac.jp:~/calgorithm
$ ssh b36697@cinnamon.kudpc.kyoto-u.ac.jp
```