// #define DEBUG_PROGLESS

#include "ica.cpp"

ICA::Matrix Pmatrix(ICA::Matrix& G){
	ICA::Matrix P = ICA::Matrix::Zero(G.rows(), G.cols());
	#pragma omp parallel for
	for(int i=0; i<G.rows(); i++){
		ICA::Vector row = G.row(i);
		ICA::Vector::Index maxId;
		row.cwiseAbs().maxCoeff(&maxId);
		const double x = row(maxId);
		P(i, maxId) = (x > 0) ? 1 : -1;
	}
	return P;
}

double test(const int sample, const int series, const int seed){
	ICA::Reng reng(seed);

	std::vector<ICA::Vector> s(sample);
	#pragma omp parallel for
	for (int i=0; i<sample; i++){
		s.at(i) = ICA::ChebytSeries(i+2, series, 0.2);
	}
	ICA::Matrix noncenterS = ICA::VStack(s);
	ICA::Matrix S = ICA::Centerize(noncenterS);

	const ICA::Matrix A = ICA::RandMatrix(sample, reng);
	ICA::Matrix X = A * S;
	auto result = ICA::FastICA(X);
	ICA::Matrix G = result.W * A;
	ICA::Matrix P = Pmatrix(G);
	ICA::Matrix S2 = P.transpose() * result.Y;

	// 平均2乗誤差
	const double mse = (S2-S).array().pow(2).mean();

	return mse;
}

int main(){
	auto series = 10000;
	const auto times = 100;
	const auto sample_max = 30;
	for(int sample=2; sample<sample_max; sample++){
		double mse_sum = 0.0;
		#pragma omp parallel for reduction(+:mse_sum)
		for (int i=0; i<times; i++){
			mse_sum += test(sample, series, i);
		}
		std::cout << sample << "\t" << mse_sum/times << std::endl;
	}
	return 0;
}