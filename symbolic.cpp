// #define NDEBUG

#include <vector>
#include "cica.cpp"

void test(const int signals, const int samplings, const int seed, const double norm_scale){
	cica::random_engine random_engine(seed);
	const cica::imatrix B = cica::random_bits(signals, samplings, random_engine);
	std::vector<cica::vector> s(signals);
	#pragma omp parallel for
	for (int i=0; i<signals; i++){
		s.at(i) = cica::chebyt_sampling(i+2, samplings, 0.1);
	}
	cica::matrix noncenterS = cica::vstack(s);
	cica::matrix S = cica::centerize(noncenterS);

	const cica::matrix T = S * B.cast<double>();
	const cica::matrix A = cica::random_uniform_matrix(signals, random_engine);
	const cica::matrix X = A * T + cica::gauss_matrix(signals, samplings, norm_scale, random_engine);

	auto res = cica::fastica(X);
	const cica::imatrix P = cica::estimate_circulant_matrix(A, res.W);

	const cica::matrix S2 = P.cast<double>().transpose() * res.Y;
	const cica::imatrix RB = (S2.array() * S.array()).sign().matrix().cast<int>();
	const double ber = cica::bit_error_rate(B, RB);
	const double cte = cica::cross_talk_error(A, res.W);

	std::cout << ber << std::endl;
	std::cout << cte << std::endl;
}

int main(){
	cica::random_engine random_engine(0);
	const auto signals = 3;
	const auto samplings = 1000;

	test(signals, samplings, 1, 0.1);


	return 0;
}