#define DEBUG_PROGLESS

#include <vector>
#include "ica.cpp"

int main(){
	cica::reng reng(0);
	const auto signals = 3;
	const auto samplings = 1000;
	
	std::vector<cica::vector> s(signals);
	#pragma omp parallel for
	for (int i=0; i<signals; i++){
		s.at(i) = cica::chebyt_sampling(i+2, samplings, 0.1);
	}
	cica::matrix noncenterS = cica::vstack(s);
	cica::matrix S = cica::centerize(noncenterS);

	cica::matrix A = cica::rand_matrix(signals, reng);
	cica::matrix X = A * S;
	auto result = cica::fast_ica(X);
	cica::matrix P = cica::simple_circulant_P(A, result.W);
	cica::matrix S2 = P.transpose() * result.Y;

    const int plot_num = std::min(1000, samplings);
	std::stringstream ss;
	cica::matrix SR = S.rightCols(plot_num);
	cica::write_matrix(ss, S);
	ss << std::endl;
	cica::matrix XR = X.rightCols(plot_num);
	cica::write_matrix(ss, X);
	ss << std::endl;
	cica::matrix YR = S2.rightCols(plot_num);
	cica::write_matrix(ss, YR);

	std::ofstream outputfile("result.csv");
	outputfile << ss.rdbuf();
	outputfile.close();
	return 0;
}