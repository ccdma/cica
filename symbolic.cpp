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
	double cte;
	double ncte;
	double mse;
	double loop_ave;
	double correlaion_mse;
	double res_correlaion_mse;
	double time;
};

test_report test(const int signals, const int samplings, const int seed, const double norm_stddev){
	cica::util::timer timer;
	cica::random_engine random_engine(cica::util::get_seed_by_time());
	const cica::imatrix B = cica::random_bits(signals, samplings, random_engine);
	cica::matrix noncenterS(signals, samplings);

	for (int i=0; i<signals; i++){
		std::uniform_real_distribution<double> distribution(-0.99, 0.99);
		noncenterS.row(i) = cica::chebyt_sampling(2, samplings, distribution(random_engine));
	}
	const cica::matrix S = cica::centerize(noncenterS);

#ifndef NDEBUG
	const cica::matrix CS = cica::correlation_matrix(S);
	const double correlaion_msea = cica::mean_squared_error(CS, cica::matrix::Identity(CS.rows(), CS.cols()));
	assert(correlaion_msea < 1.0e-2);
#endif

	const cica::matrix T = (S.array() * B.cast<double>().array()).matrix();
	const cica::matrix A = cica::random_uniform_matrix(signals, random_engine);
	const cica::matrix X = A * T + cica::gauss_matrix(signals, samplings, norm_stddev, random_engine);

	const auto res = cica::fastica::fastica(X);
	const cica::imatrix P = cica::estimate_circulant_matrix(A, res.W);

	const cica::matrix Z = P.cast<double>().transpose() * res.Y;
	const cica::imatrix RB = (Z.array() * S.array()).sign().matrix().cast<int>();
	const double ber = cica::bit_error_rate(B, RB);
	const double cte = cica::cross_talk_error(A, res.W);
	const double ncte = cica::normal_cross_talk_error(A, res.W);
	const double mse = cica::mean_squared_error(T, Z);
	const cica::matrix CT = cica::correlation_matrix(T);
	const double correlaion_mse = cica::mean_squared_error(CT, cica::matrix::Identity(CT.rows(), CT.cols()));
	const cica::matrix CY = cica::correlation_matrix(res.Y);
	const double res_correlaion_mse = cica::mean_squared_error(CY, cica::matrix::Identity(CY.rows(), CY.cols()));
	const double loop_ave = res.loop.cast<double>().mean();

	const double time = (double)timer.from_start();
	return test_report{.ber=ber, .cte=cte, .ncte=ncte, .mse=mse,
		 .loop_ave=loop_ave, .correlaion_mse=correlaion_mse, .res_correlaion_mse=res_correlaion_mse, .time=time};
}

int main(){
	const auto trials = 10;
	const auto sep = "\t";
	std::cout << "commit" << sep << COMMIT_ID << std::endl;
	std::cout << "trials" << sep << trials << std::endl;
	std::cout
		<< "signals" << sep
		<< "samplings" << sep
		<< "stddev" << sep
		<< "mse" << sep
		<< "correlaion_mse" << sep
		<< "res_correlaion_mse" << sep
		<< "loop_ave" << sep
		<< "ber" << sep
		<< "cte" << sep
		<< "ncte" << sep
		<< "complete" << sep
		<< "time(ms)"
	<< std::endl;	// header
	// const auto samplings = 10000;
	// const auto signals = 100;
	const auto stddev = 0.0;
	std::vector<int> v1{1000, 2000, 5000};
	std::vector<int> v2 = cica::util::range(2, 100);
	for(const auto& samplings : v1){
	for(const auto& j : v2){
		const int signals = j;
		if (samplings < 5000) break;
		if (signals < 70) break;
		int complete = 0;
		double ber_sum = 0.0;
		double cte_sum = 0.0;
		double ncte_sum = 0.0;
		double mse_sum = 0.0;
		double correlaion_mse_sum = 0.0;
		double res_correlaion_mse_sum = 0.0;
		double loop_ave_sum = 0.0;
		double time = 0.0;
		// #pragma omp parallel for
		for (int seed=0; seed<trials; seed++){
			try {
				const auto report = test(signals, samplings, seed, stddev);
				#pragma omp critical
				{
					ber_sum += report.ber;
					cte_sum += report.cte;
					ncte_sum += report.ncte;
					mse_sum += report.mse;
					loop_ave_sum += report.loop_ave;
					correlaion_mse_sum += report.correlaion_mse;
					res_correlaion_mse_sum += report.res_correlaion_mse;
					time += report.time;
					complete += 1;
				}
			} catch (cica::exception::base e) {}
		}
		std::cout
			<< signals << sep
			<< samplings << sep
			<< stddev << sep
			<< mse_sum/trials << sep
			<< correlaion_mse_sum/trials << sep
			<< res_correlaion_mse_sum/trials << sep
			<< loop_ave_sum/trials << sep
			<< ber_sum/trials << sep 
			<< cte_sum/trials << sep 
			<< ncte_sum/trials << sep 
			<< complete << sep 
			<< time/trials
		<< std::endl;
		if (ber_sum/trials > 0.3) break;
	}} // end root for
	return 0;
}