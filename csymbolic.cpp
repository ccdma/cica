#define NDEBUG
#define NPROGLESS

#ifndef COMMIT_ID
	#define COMMIT_ID "undefined"
#endif

#include <cmath>
#include <cassert>
#include "cica.hpp"

struct test_report {

	double ber;
	double rber;
	double time;
};

test_report test(const int signals, const int samplings, const int seed, const double norm_stddev){
	cica::util::timer timer;
	cica::random_engine random_engine(seed*signals*cica::util::get_seed_by_time());
	const cica::imatrix B = cica::random_bits(signals, samplings, random_engine);
	cica::cmatrix noncenterS(signals, samplings);

	for (int i=0; i<signals; i++){
		std::uniform_real_distribution<double> distribution(-0.99, 0.99);
		noncenterS.row(i) = cica::const_powerd_sampling(2, 2*M_PI*distribution(random_engine), samplings);
	}

	const cica::cmatrix S = cica::centerize(noncenterS);

	const cica::cmatrix T = (S.array() * B.cast<double>().array()).matrix();
	const cica::matrix A = cica::random_uniform_matrix(signals, random_engine);
	const cica::cmatrix X = A * T + cica::cgauss_matrix(signals, samplings, norm_stddev, random_engine);

	const auto r_res = cica::fastica::fastica(X.real());
	const auto i_res = cica::fastica::fastica(X.imag());

	const cica::imatrix rP = cica::estimate_circulant_matrix(A, r_res.W);
	const cica::imatrix iP = cica::estimate_circulant_matrix(A, i_res.W);

	const cica::cmatrix Z = rP.cast<double>().transpose() * r_res.Y + iP.cast<double>().transpose() * i_res.Y * cica::dcomplex(0, 1);
	const cica::imatrix RB = (Z.real().array() * S.real().array() + Z.imag().array() * S.imag().array()).sign().matrix().cast<int>();
	const double ber = cica::bit_error_rate(B, RB);

	const cica::imatrix rRB = (Z.real().array() * S.real().array()).sign().matrix().cast<int>();
	const double rber = cica::bit_error_rate(B, rRB);

	const double time = (double)timer.from_start();
	return test_report{.ber=ber, .rber=rber, .time=time};
}

int main(){
	const auto trials = 1000;
	const auto sep = "\t";
	std::cout << "commit" << sep << COMMIT_ID << std::endl;
	std::cout << "trials" << sep << trials << std::endl;
	std::cout
		<< "signals" << sep
		<< "samplings" << sep
		<< "stddev" << sep
		<< "ber" << sep
		<< "rber" << sep
		<< "complete" << sep
		<< "time(ms)"
	<< std::endl;	// header
	const auto samplings = 10000;
	// const auto signals = 3;
	const auto stddev = 0.1;
	for(int i=0; i<=1; i++){
	for(double j=1; j<=100; j++){
		const auto signals = j+1;
		// const auto stddev = 0.01 * (double)(i*i);
		// const auto samplings = 1000 * i;
		int complete = 0;
		double ber_sum = 0.0;
		double rber_sum = 0.0;
		double time = 0.0;
		#pragma omp parallel for
		for (int seed=0; seed<trials; seed++){
			try {
				const auto report = test(signals, samplings, seed, stddev);
				#pragma omp critical
				{
					ber_sum += report.ber;
					rber_sum += report.rber;
					time += report.time;
					complete += 1;
				}
			} catch (cica::exception::base e) {}
		}
		std::cout
			<< signals << sep
			<< samplings << sep
			<< stddev << sep
			<< ber_sum/trials << sep 
			<< rber_sum/trials << sep 
			<< complete << sep
			<< time/trials
		<< std::endl;
		// if (ber_sum/trials > 0.001) break;
	}} // end root for
	return 0;
}