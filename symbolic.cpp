#define NPARALLELIZE
#define NDEBUG
#define NPROGLESS

#ifndef COMMIT_ID
	#define COMMIT_ID "undefined"
#endif

#include "cica.cpp"

struct test_report {

	double ber;
	double cte;
	double ncte;
	double mse;
	double loop_ave;
	double correlaion_mse;
};

test_report test(const int signals, const int samplings, const int seed, const double norm_stddev, const int chebyt_n){
	cica::random_engine random_engine(seed);
	std::uniform_real_distribution<double> distribution(-0.99, 0.99);
	const cica::imatrix B = cica::random_bits(signals, samplings, random_engine);
	cica::matrix noncenterS(signals, samplings);
	#pragma omp parallel for
	for (int i=0; i<signals; i++){
		noncenterS.row(i) = cica::chebyt_sampling(chebyt_n, samplings, distribution(random_engine));	// 変更する場合はヘッダも変更する
	}
	const cica::matrix S = cica::centerize(noncenterS);

	const cica::matrix T = (S.array() * B.cast<double>().array()).matrix();
	const cica::matrix A = cica::random_uniform_matrix(signals, random_engine);
	const cica::matrix X = A * T + cica::gauss_matrix(signals, samplings, norm_stddev, random_engine);

	const auto res = cica::fastica(X);
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

	return test_report{.ber=ber, .cte=cte, .ncte=ncte, .mse=mse, .loop_ave=loop_ave, .correlaion_mse=correlaion_mse};
}

int main(){
	cica::random_engine random_engine(0);
	const auto samplings = 1000;
	const auto chebyt_n = 2;
	const auto trials = 100;
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
		<< "ncte"
	<< std::endl;	// header
	for(int signals=10; signals<150; signals+=20){
		for(double stddev=0; stddev<0.15; stddev+=0.005){ 
			double ber_sum = 0.0;
			double cte_sum = 0.0;
			double ncte_sum = 0.0;
			double mse_sum = 0.0;
			double correlaion_mse_sum = 0.0;
			double loop_ave_sum = 0.0;
			#pragma omp parallel for reduction(+:ber_sum,cte_sum,ncte_sum,mse_sum,loop_ave_sum,correlaion_mse_sum)
			for (int i=0; i<trials; i++){
				const auto report = test(signals, samplings, i, stddev, chebyt_n);
				ber_sum += report.ber;
				cte_sum += report.cte;
				ncte_sum += report.ncte;
				mse_sum += report.mse;
				loop_ave_sum += report.loop_ave;
				correlaion_mse_sum += report.correlaion_mse;
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
				<< ncte_sum/trials
			<< std::endl;
		}
	}
	return 0;
}