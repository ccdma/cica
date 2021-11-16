// #define NDEBUG
#define NPROGLESS

#ifndef COMMIT_ID
	#define COMMIT_ID "undefined"
#endif

#include <cmath>
#include <cassert>
#include "cica.cpp"

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

test_report test(const int signals, const int samplings, const int seed, const double norm_stddev, const int chebyt_n){
	cica::util::timer timer;
	cica::random_engine random_engine(seed);
	const cica::imatrix B = cica::random_bits(signals, samplings, random_engine);
	cica::matrix noncenterS(signals, samplings);
	#pragma omp parallel for
	for (int i=0; i<signals; i++){
		cica::random_engine random_engine(i*(seed+1));
		std::uniform_real_distribution<double> distribution(-0.99, 0.99);
		noncenterS.row(i) = cica::chebyt_sampling(chebyt_n, samplings, distribution(random_engine));	// 変更する場合はヘッダも変更する
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
	// const auto samplings = 10000;
	const auto signals = 100;
	const auto stddev = 0.0;
	const auto chebyt_n = 2;
	const auto trials = 10;
	std::cout << "commit" << "\t" << COMMIT_ID << std::endl;
	std::cout << "chebyt_n(fixed)" << "\t" << chebyt_n << std::endl;
	std::cout << "trials" << "\t" << trials << std::endl;
	std::cout
		<< "signals" << "\t"
		<< "samplings" << "\t"
		<< "stddev" << "\t"
		<< "mse" << "\t"
		<< "correlaion_mse" << "\t"
		<< "res_correlaion_mse" << "\t"
		<< "loop_ave" << "\t"
		<< "ber" << "\t"
		<< "cte" << "\t"
		<< "ncte" << "\t"
		<< "complete" << "\t"
		<< "time(ms)"
	<< std::endl;	// header
	// for(int i=1; i<=10; i++){
	for(double j=1; j<=20; j++){
		// const int signals = i*50;
		const int samplings = 5000 * j;
		int complete = 0;
		double ber_sum = 0.0;
		double cte_sum = 0.0;
		double ncte_sum = 0.0;
		double mse_sum = 0.0;
		double correlaion_mse_sum = 0.0;
		double res_correlaion_mse_sum = 0.0;
		double loop_ave_sum = 0.0;
		double time = 0.0;
		#pragma omp parallel for
		for (int seed=0; seed<trials; seed++){
			try {
				const auto report = test(signals, samplings, seed, stddev, chebyt_n);
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
			<< signals << "\t"
			<< samplings << "\t"
			<< stddev << "\t"
			<< mse_sum/trials << "\t"
			<< correlaion_mse_sum/trials << "\t"
			<< res_correlaion_mse_sum/trials << "\t"
			<< loop_ave_sum/trials << "\t"
			<< ber_sum/trials << "\t" 
			<< cte_sum/trials << "\t" 
			<< ncte_sum/trials << "\t" 
			<< complete << "\t" 
			<< time/trials
		<< std::endl;
	}
	return 0;
}