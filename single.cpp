#define NDEBUG
#define NPROGLESS

#include "cica.hpp"

int main(){
	const auto signals = 4;
	const auto samplings = 1000;
	
	cica::matrix noncenterS(signals, samplings);
	double prev_init = 0.7;
	for (int i=0; i<signals; i++){
		noncenterS.row(i) = cica::chebyt_sampling(2, samplings, prev_init);
		prev_init = cica::eval_chebyshev(2, prev_init);
	}
	// for (int i=2; i<signals; i++){
	// 	noncenterS.row(i) = cica::sine_sampling(std::sqrt(i-1)/10.0, samplings);
	// }
	
	cica::random_engine random_engine(2);
	// noncenterS = noncenterS + cica::gauss_matrix(signals, samplings, 0.1, random_engine);

	const cica::matrix S = noncenterS;
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
	const cica::matrix S2R = S2.rightCols(plot_num);
	cica::util::write_matrix(ss, S2R);

	cica::util::save_stream(ss, "data.csv");
	return 0;
}