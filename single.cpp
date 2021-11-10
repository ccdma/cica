#define NDEBUG
// #define NPROGLESS

#include <vector>
#include "cica.cpp"

int main(){
	cica::random_engine random_engine(0);
	const auto signals = 3;
	const auto samplings = 1000;
	
	cica::matrix noncenterS(signals, samplings);
	#pragma omp parallel for
	for (int i=0; i<signals; i++){
		noncenterS.row(i) = cica::chebyt_sampling(i+2, samplings, 0.1);
	}
	cica::matrix S = cica::centerize(noncenterS);
	cica::matrix A = cica::random_uniform_matrix(signals, random_engine);
	cica::matrix X = A * S;
	auto result = cica::fastica(X);
	cica::imatrix P = cica::estimate_circulant_matrix(A, result.W);
	cica::matrix S2 = P.cast<double>().transpose() * result.Y;

    const int plot_num = std::min(1000, samplings);
	std::stringstream ss;
	cica::matrix SR = S.rightCols(plot_num);
	cica::write_matrix(ss, SR);
	ss << std::endl;
	cica::matrix XR = X.rightCols(plot_num);
	cica::write_matrix(ss, XR);
	ss << std::endl;
	cica::matrix YR = result.Y.rightCols(plot_num);
	cica::write_matrix(ss, YR);

	cica::save_stream(ss, "result.csv");
	return 0;
}