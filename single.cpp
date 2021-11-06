#define DEBUG_PROGLESS

#include <vector>
#include "ica.cpp"

int main(){
	ICA::Reng reng(1);
	const auto sample = 4;
	const auto series = 1000;
	
	std::vector<ICA::Vector> s(sample);
	#pragma omp parallel for
	for (int i=0; i<sample-2; i++){
		s.at(i) = ICA::SinSeries(1.0/((double)i+1.0), series);
	}
	ICA::Matrix noncenterS = ICA::VStack(s);
	ICA::Matrix S = ICA::Centerize(noncenterS);

	ICA::Matrix A = ICA::RandMatrix(sample, reng);
	ICA::Matrix X = A * S;
	auto result = ICA::FastICA(X);
	ICA::Matrix P = ICA::SimpleCirculantP(A, result.W);
	ICA::Matrix S2 = P.transpose() * result.Y;

    const int plot_num = std::min(1000, series);
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