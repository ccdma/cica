## compile

- `-lblas`でblas利用（eigenのinclude前に`#define EIGEN_USE_BLAS`すること）
  - [Using BLAS/LAPACK from Eigen](https://eigen.tuxfamily.org/dox/TopicUsingBlasLapack.html)


## exec
- 環境変数`OMP_NUM_THREADS=6`でスレッド数を指定可能
- `OMP_NESTED=TRUE`でOpenMPのネストが可能
  - https://stackoverflow.com/questions/15830526/c-openmp-parallel-for-inside-a-parallel-region

## スパコン関連
- [スーパーコンピュータシステムの使い方](https://web.kudpc.kyoto-u.ac.jp/manual/ja)
- [利用者ポータル](https://web.kudpc.kyoto-u.ac.jp/portal/)
