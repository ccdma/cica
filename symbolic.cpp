#define NPARALLELIZE
#define NDEBUG
#define NPROGLESS

#include <vector>
#include "cica.cpp"

struct test_report {

	double ber;
	double cte;
};

test_report test(const int signals, const int samplings, const int seed, const double norm_stddev){
	cica::random_engine random_engine(seed);
	std::uniform_real_distribution<double> distribution(-0.9, 0.9);
	const cica::imatrix B = cica::random_bits(signals, samplings, random_engine);
	std::vector<cica::vector> s(signals);
	#pragma omp parallel for
	for (int i=0; i<signals; i++){
		s.at(i) = cica::chebyt_sampling(2, samplings, distribution(random_engine));
	}
	cica::matrix noncenterS = cica::vstack(s);
	cica::matrix S = cica::centerize(noncenterS);

	const cica::matrix T = (S.array() * B.cast<double>().array()).matrix();
	const cica::matrix A = cica::random_uniform_matrix(signals, random_engine);
	const cica::matrix X = A * T + cica::gauss_matrix(signals, samplings, norm_stddev, random_engine);

	auto res = cica::fastica(X);
	const cica::imatrix P = cica::estimate_circulant_matrix(A, res.W);

	const cica::matrix S2 = P.cast<double>().transpose() * res.Y;
	const cica::imatrix RB = (S2.array() * S.array()).sign().matrix().cast<int>();
	const double ber = cica::bit_error_rate(B, RB);
	const double cte = cica::cross_talk_error(A, res.W);
	
	return test_report{.ber=ber, .cte=cte};
}

int main(){
	cica::random_engine random_engine(0);
	const auto signals = 3;
	const auto samplings = 10000;
	const auto times = 100;
	std::cout << "samplings\t" << samplings << std::endl;
	std::cout << "times\t" << times << std::endl;
	std::cout << "stddev\tber\tcte" << std::endl;	// header
	for(int i=0; i<20; i++){
		const double stddev = (double)i/100; 
		double ber_sum = 0.0;
		double cte_sum = 0.0;
		#pragma omp parallel for reduction(+:ber_sum,cte_sum)
		for (int i=0; i<times; i++){
			const auto report = test(signals, samplings, i, stddev);
			ber_sum += report.ber;
			cte_sum += report.cte;
		}
		std::cout << stddev << "\t" << ber_sum/times << "\t" << cte_sum/times << std::endl;
	}
	return 0;
}