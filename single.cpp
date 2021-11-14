#define NDEBUG
// #define NPROGLESS

#include "cica.cpp"

int main(){
	cica::random_engine random_engine(1);
	const auto signals = 400;
	const auto samplings = 500000;
	
	cica::matrix noncenterS(signals, samplings);
	#pragma omp parallel for
	for (int i=0; i<signals; i++){
		cica::random_engine random_engine(i);
		std::uniform_real_distribution<double> distribution(-0.99, 0.99);
		noncenterS.row(i) = cica::chebyt_sampling(2, samplings, distribution(random_engine));
	}

	const cica::matrix S = cica::centerize(noncenterS);
	const cica::matrix A = cica::random_uniform_matrix(signals, random_engine);
	const cica::matrix X = A * S;
	auto result = cica::fastica::fastica(X);
	const cica::imatrix P = cica::estimate_circulant_matrix(A, result.W);
	const cica::matrix S2 = P.cast<double>().transpose() * result.Y;

    const int plot_num = std::min(1000, samplings);
	std::stringstream ss;
	const cica::matrix SR = S.rightCols(plot_num);
	cica::util::write_matrix(ss, SR);
	ss << std::endl;
	const cica::matrix XR = X.rightCols(plot_num);
	cica::util::write_matrix(ss, XR);
	ss << std::endl;
	const cica::matrix YR = result.Y.rightCols(plot_num);
	cica::util::write_matrix(ss, YR);

	cica::util::save_stream(ss, "result.csv");
	return 0;
}