#define NDEBUG
#define NPROGLESS

#ifndef COMMIT_ID
	#define COMMIT_ID "undefined"
#endif

#include <cmath>
#include <cassert>
#include <chrono>
#include <thread>
#include "cica.hpp"

struct test_report {

	double ber;
	double time;
};

/**
 * 
 * @param K ユーザー数
 * @param N 拡散符号長
 * @param stddev 
 * @return test_report 
 */
test_report test(const int K, const int N, const float stddev){
	cica::util::timer timer;
	cica::random_engine random_engine(K*N*cica::util::get_seed_by_time());
	const cica::ivector BITS = cica::random_bits(1, N, random_engine).row(0);

	const cica::vector DBITS = BITS.cast<double>(); // BPSK

	cica::matrix B(K, N);
	for (int i=0; i<K; i++){
		B.row(i) = DBITS;
	}

	cica::cmatrix S(K, N);
	for (int i=0; i<K; i++){
		std::uniform_real_distribution<double> distribution(-0.99, 0.99);
		S.row(i) = cica::const_powerd_sampling(2, 2*M_PI*distribution(random_engine), N);
	}

	const cica::matrix T = (S.array() * B.array()).matrix();
	const cica::vector A = cica::vector::Ones(N);
	const cica::vector AWGN = cica::cgauss_matrix(1, N, stddev, random_engine).row(0);
	const cica::cvector X = A * T + AWGN;

	const cica::matrix RB = (X.array() * S.conjugate().array()).matrix().rowwise().mean();
	

}

int main(){
	
	const cica::cvector S = cica::weyl_sampling(std::sqrt(0.1), 0.1, 1000);

	const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");

	for (int i=0; i<1000; i++){
		std::cout << S(i).real() << ",";
	}
	std::cout << std::endl;
	for (int i=0; i<1000; i++){
		std::cout << S(i).imag() << ",";
	}
	std::cout << std::endl;
	return 0;
}