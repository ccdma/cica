#define NPARALLELIZE
#define NDEBUG
#define NPROGLESS

#include "cica.cpp"

struct test_report {

	double mse;
	double loop_ave;
};

test_report test(const int signals, const int samplings, const int seed, const int chebyt_start_n){
	cica::random_engine random_engine(seed);

	std::vector<cica::vector> s(signals);
	#pragma omp parallel for
	for (int i=0; i<signals; i++){
		s.at(i) = cica::chebyt_sampling(i+chebyt_start_n, samplings, 0.2);
	}
	cica::matrix noncenterS = cica::vstack(s);
	cica::matrix S = cica::centerize(noncenterS);

	cica::matrix A = cica::random_uniform_matrix(signals, random_engine);
	cica::matrix X = A * S;
	auto result = cica::fastica(X);
	cica::imatrix P = cica::estimate_circulant_matrix(A, result.W);
	cica::matrix S2 = P.cast<double>().transpose() * result.Y;

	// 平均2乗誤差
	const double mse = cica::mean_squared_error(S, S2);
	// ループ回数の平均
	const double loop_ave = result.loop.cast<double>().mean();

	return test_report{.mse=mse, .loop_ave=loop_ave};
}

int main(){
	auto samplings = 10000;
	const auto times = 100;
	const auto chebyt_start_n = 2;
	std::cout << "samplings\t" << samplings << std::endl;
	std::cout << "times\t" << times << std::endl;
	std::cout << "chebyt_start_n\t" << chebyt_start_n << std::endl;
	std::cout << "signals\tmse\tloop_ave" << std::endl;	// header
	const auto sample_max = 200;
	for(int signals=2; signals<sample_max; signals++){
		double mse_sum = 0.0;
		double loop_ave_sum = 0.0;
		#pragma omp parallel for reduction(+:mse_sum,loop_ave_sum)
		for (int i=0; i<times; i++){
			auto report = test(signals, samplings, i, chebyt_start_n);
			mse_sum += report.mse;
			loop_ave_sum += report.loop_ave;
		}
		std::cout << signals << "\t" << mse_sum/times << "\t" << loop_ave_sum/times << std::endl;
	}
	return 0;
}