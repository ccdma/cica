## compile

- `-lblas`でblas利用（eigenのinclude前に`#define EIGEN_USE_BLAS`すること）
  - [Using BLAS/LAPACK from Eigen](https://eigen.tuxfamily.org/dox/TopicUsingBlasLapack.html)

vscodeでエラーが出る場合
```
apt-get install g++-multilib
```

## exec
- 環境変数`OMP_NUM_THREADS=6`でスレッド数を指定可能

## スパコン関連
- [スーパーコンピュータシステムの使い方](https://web.kudpc.kyoto-u.ac.jp/manual/ja)
- [利用者ポータル](https://web.kudpc.kyoto-u.ac.jp/portal/)
