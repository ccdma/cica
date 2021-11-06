#define DEBUG_PROGLESS

#include <vector>
#include "ica.cpp"

int main(){
	ICA::Reng reng(0);
	const auto signals = 20;
	const auto samplings = 1000;
	
	std::vector<ICA::Vector> s(signals);
	#pragma omp parallel for
	for (int i=0; i<signals; i++){
		s.at(i) = ICA::ChebytSeries(i+2, samplings, 0.1);
	}
	ICA::Matrix noncenterS = ICA::VStack(s);
	ICA::Matrix S = ICA::Centerize(noncenterS);

	ICA::Matrix A = ICA::RandMatrix(signals, reng);
	ICA::Matrix X = A * S;
	auto result = ICA::FastICA(X);
	ICA::Matrix P = ICA::SimpleCirculantP(A, result.W);
	ICA::Matrix S2 = P.transpose() * result.Y;

    const int plot_num = std::min(1000, samplings);
	std::stringstream ss;
	ICA::Matrix SR = S.rightCols(plot_num);
	ICA::WriteMatrix(ss, S);
	ss << std::endl;
	ICA::Matrix XR = X.rightCols(plot_num);
	ICA::WriteMatrix(ss, X);
	ss << std::endl;
	ICA::Matrix YR = S2.rightCols(plot_num);
	ICA::WriteMatrix(ss, YR);

	std::ofstream outputfile("result.csv");
	outputfile << ss.rdbuf();
	outputfile.close();
	return 0;
}