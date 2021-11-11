#define NPARALLELIZE
#define NDEBUG
#define NPROGLESS

#define COMMIT_ID "undefined"

#include <vector>
#include "cica.cpp"

struct test_report {

	double ber;
	double cte;
};

test_report test(const int signals, const int samplings, const int seed, const double norm_stddev, const int chebyt_n){
	cica::random_engine random_engine(seed);
	std::uniform_real_distribution<double> distribution(-0.9, 0.9);
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

	const cica::matrix S2 = P.cast<double>().transpose() * res.Y;
	const cica::imatrix RB = (S2.array() * S.array()).sign().matrix().cast<int>();
	const double ber = cica::bit_error_rate(B, RB);
	const double cte = cica::cross_talk_error(A, res.W);
	
	return test_report{.ber=ber, .cte=cte};
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
		<< "ber" << "\t"
		<< "cte"
	<< std::endl;	// header
	for(int signals=2; signals<50; signals+=5){
		for(double stddev=0; stddev<0.15; stddev+=0.005){ 
			double ber_sum = 0.0;
			double cte_sum = 0.0;
			#pragma omp parallel for reduction(+:ber_sum,cte_sum)
			for (int i=0; i<trials; i++){
				const auto report = test(signals, samplings, i, stddev, chebyt_n);
				ber_sum += report.ber;
				cte_sum += report.cte;
			}
			std::cout
				<< signals << "\t"
				<< samplings << "\t"
				<< stddev << "\t"
				<< ber_sum/trials << "\t" 
				<< cte_sum/trials
			<< std::endl;
		}
	}
	return 0;
}