// #define NPARALLELIZE
#define NDEBUG
#define NPROGLESS

#ifndef COMMIT_ID
	#define COMMIT_ID "undefined"
#endif

#include <cmath>
#include "cica.cpp"

struct test_report {

	double ber;
	double cte;
	double ncte;
	double mse;
	double loop_ave;
	double correlaion_mse;
	double time;
};

test_report test(const int signals, const int samplings, const int seed, const double norm_stddev, const int chebyt_n){
	const auto start = std::chrono::system_clock::now();

	cica::random_engine random_engine(seed);
	std::uniform_real_distribution<double> distribution(-0.99, 0.99);
	const cica::imatrix B = cica::random_bits(signals, samplings, random_engine);
	cica::matrix noncenterS(signals, samplings);
	#pragma omp parallel for
	for (int i=0; i<signals; i++){
		const int additional_n = i%10;
		noncenterS.row(i) = cica::chebyt_sampling(chebyt_n+additional_n, samplings, distribution(random_engine));	// 変更する場合はヘッダも変更する
	}
	const cica::matrix S = cica::centerize(noncenterS);

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
	const double loop_ave = res.loop.cast<double>().mean();

	const double time = (double)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start).count();
	return test_report{.ber=ber, .cte=cte, .ncte=ncte, .mse=mse,
		 .loop_ave=loop_ave, .correlaion_mse=correlaion_mse, .time=time};
}

int main(){
	cica::random_engine random_engine(0);
	const auto samplings = 10000;
	const auto signals = 200;
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
		<< "loop_ave" << "\t"
		<< "ber" << "\t"
		<< "cte" << "\t"
		<< "ncte" << "\t"
		<< "complete" << "\t"
		<< "time(ms)"
	<< std::endl;	// header
	for(int signals=50; signals<500; signals+=50){
	for(double i=0.0; i<8; i+=1.0){
		const int samplings = 1000 * (int)std::pow(2, i);
		int complete = 0;
		double ber_sum = 0.0;
		double cte_sum = 0.0;
		double ncte_sum = 0.0;
		double mse_sum = 0.0;
		double correlaion_mse_sum = 0.0;
		double loop_ave_sum = 0.0;
		double time = 0.0;
		#pragma omp parallel for reduction(+:complete,ber_sum,cte_sum,ncte_sum,mse_sum,loop_ave_sum,correlaion_mse_sum,time)
		for (int seed=0; seed<trials; seed++){
			try {
				const auto report = test(signals, samplings, seed, stddev, chebyt_n);
				ber_sum += report.ber;
				cte_sum += report.cte;
				ncte_sum += report.ncte;
				mse_sum += report.mse;
				loop_ave_sum += report.loop_ave;
				correlaion_mse_sum += report.correlaion_mse;
				time += report.time;
				complete += 1;
			} catch (cica::exception::base e) {}
		}
		std::cout
			<< signals << "\t"
			<< samplings << "\t"
			<< stddev << "\t"
			<< mse_sum/trials << "\t"
			<< correlaion_mse_sum/trials << "\t"
			<< loop_ave_sum/trials << "\t"
			<< ber_sum/trials << "\t" 
			<< cte_sum/trials << "\t" 
			<< ncte_sum/trials << "\t" 
			<< complete << "\t" 
			<< time/trials
		<< std::endl;
	}}
	return 0;
}