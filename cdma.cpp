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
test_report test(const int K, const int N, const int seed, const float stddev){
	cica::util::timer timer;
	cica::random_engine random_engine(K*N*cica::util::get_seed_by_time());

	const cica::ivector BITS = cica::random_bits(1, K, random_engine).row(0);
	const cica::cvector BPSK_DATA = BITS.cast<double>(); // BPSK

	const cica::cmatrix B = BPSK_DATA.replicate(1, N); //拡散符号分の長さにする
	cica::cmatrix S(K, N);
	for (int i=0; i<K; i++){
		S.row(i) = cica::weyl_sampling((double)i/K, (double)1/(2*N), N);
	}

	const cica::cmatrix T = (S.array() * B.array()).matrix();
	const cica::vector A = cica::vector::Ones(K);
	const cica::cvector AWGN = cica::cgauss_matrix(1, N, stddev, random_engine).row(0);
	const cica::cvector X = T.transpose() * A + AWGN;

	const cica::cmatrix RB = (X.transpose().replicate(K, 1).array() * S.conjugate().array()).matrix();
	
	const cica::cvector RBPSK_DATA = RB.rowwise().mean();	// BPSKの形に戻す
	const cica::ivector RBITS = RBPSK_DATA.real().array().sign().matrix().cast<int>();

	const auto ber = cica::bit_error_rate(BITS, RBITS);
	return test_report{.ber=ber, .time=(double)timer.from_start()}; 
}

int main(){
	
	const auto trials = 1000*1000;
	const auto sep = ",";
	auto timer = new cica::util::timer();
	std::cout << "commit" << ":" << COMMIT_ID << std::endl;
	std::cout
		<< "K" << sep
		<< "N" << sep
		<< "stddev" << sep
		<< "ber" << sep
		<< "complete" << sep
		<< "time(ms)"
	<< std::endl;	// header
	// const auto N = 1000;
	// const auto K = 100;
	const auto stddev = 0.01;
	std::vector<int> v1{31}; // v1{10, 20, 30}
	std::vector<int> v2 = cica::util::range(2, 60);
	for(const auto& N : v1){
	for(const auto& K : v2){

		int complete = 0;
		double ber_sum = 0.0;
		double time = 0.0;
		#pragma omp parallel for
		for (int seed=0; seed<trials; seed++){
			try {
				const auto report = test(K, N, seed, stddev);
				#pragma omp critical
				{
					ber_sum += report.ber;
					time += report.time;
					complete += 1;
				}
			} catch (cica::exception::base e) {}
		}
		std::cout
			<< K << sep
			<< N << sep
			<< stddev << sep
			<< ber_sum/complete << sep 
			<< complete << sep
			<< time/complete
		<< std::endl;
		if (ber_sum/complete > 0.005) break;
	}} // end root for
	return 0;
}