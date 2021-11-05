// #define DEBUG_PROGLESS
// #define DEBUG_MATRIX
// #define EIGEN_USE_BLAS

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>
#include <sstream>
#include <Eigen/Dense>

namespace ICA {

	using Matrix = Eigen::MatrixXd;
	using Vector = Eigen::VectorXd;
	using Reng = std::mt19937; 

	const int LOOP = 500;
	const int WRITE_LIMIT = 10000;

	Matrix RandMatrix(int size, Reng& engine){
		std::uniform_real_distribution<double> distribution(-0.5, 0.5);
		auto generator = [&] (double dummy) {return distribution(engine);};
		return Matrix::Zero(size, size).unaryExpr(generator);
	};

	/**
	 * 正方行列でなくてはいけない
	 * i番目を直行空間に射影
	 */
	void Normalize(Matrix& M, int i){
		const auto size = M.cols();
		if (i>0){
			M.col(i) = M.col(i) - M.block(0, 0, size, i) * M.block(0, 0, size, i).transpose() * M.col(i);
		}
		M.col(i) = M.col(i) / std::sqrt(M.col(i).squaredNorm());
	}

	Matrix Centerize(Matrix& M){
		return M.colwise() - M.rowwise().mean();
	}

	struct FastICAResult {
		Matrix W;	// 復元行列
		Matrix Y;	// 復元信号
	};

	/**
	 * X: 内部で中心化は行うが、すでに中心化されていることが望ましい（元信号Sの中心化ができていれば、混合されたXも自然と中心化されるはず）
	 */
	FastICAResult FastICA(Matrix& X) {

#ifdef DEBUG_PROGLESS
	std::chrono::system_clock::time_point start, prev, now;
	start = std::chrono::system_clock::now();
	prev = start;
	now = start;
	std::cout 
	<< "[PROGLESS] start fastica session"
	<< "\t" << std::chrono::duration_cast<std::chrono::milliseconds>(now-prev).count() << std::endl;
	prev = now;
#endif

	ICA::Reng reng(0);
	const auto sample = X.rows();
	const auto series = X.cols();
	
	const Matrix X_center = Centerize(X);
	const Matrix X_cov = (X_center * X_center.transpose()) / double(X_center.cols() - 1);

	Eigen::SelfAdjointEigenSolver<Matrix> es(X_cov);
	if (es.info() != Eigen::Success) abort();

	Vector lambdas = es.eigenvalues().real();
	Matrix P = es.eigenvectors().real();
	Matrix Atilda = lambdas.cwiseSqrt().asDiagonal().inverse() * P.transpose();
	Matrix X_whiten = Atilda * X_center;

#ifdef DEBUG_PROGLESS
	now = std::chrono::system_clock::now();
	std::cout 
	<< "[PROGLESS] start fixed point method"
	<< "\t" << std::chrono::duration_cast<std::chrono::milliseconds>(now-prev).count() << std::endl;
	prev = now;
#endif

#ifdef DEBUG_MATRIX
	// 単位行列であることを確認
	std::cout << (X_whiten * X_whiten.transpose()) / double(X_whiten.cols() - 1) << std::endl;
#endif

	const auto g = [](double bx) { return std::pow(bx, 3); };
	const auto g2 = [](double bx) { return 3*std::pow(bx, 2); };

	const auto I = X_whiten.rows();
	auto B = RandMatrix(I, reng);

	for(int i=0; i<I; i++){
		Normalize(B, i);
	}

#ifdef DEBUG_MATRIX
		// 単位行列であることを確認
	std::cout << B * B.transpose() << std::endl;
#endif

	for(int i=0; i<I; i++){
		for(int j=0; j<LOOP; j++){
			// 値のコピー
			const Vector prevBi = B.col(i);

			const auto collen = X_whiten.cols();
			Matrix ave(I, collen);
			#pragma omp parallel for
			for(int k=0; k<collen; k++){
				const Vector x = X_whiten.col(k);
				ave.col(k) = g(x.dot(B.col(i)))*x - g2(x.dot(B.col(i)))*B.col(i);  
			}
			B.col(i) = ave.rowwise().mean();
			Normalize(B, i);
			const auto diff = std::abs(prevBi.dot(B.col(i)));
			if (1.0 - 1.e-8 < diff && diff < 1.0 + 1.e-8) break;
			if (j==LOOP-1) printf("[WARN] loop limit exceeded\n");
		}

#ifdef DEBUG_PROGLESS
		now = std::chrono::system_clock::now();
		std::cout
		<< "[PROGLESS] end loop " << i
		<< "\t" << std::chrono::duration_cast<std::chrono::milliseconds>(now-prev).count() << std::endl;
		prev = now;
#endif
		}
		Matrix Y = B.transpose() * X_whiten;

#ifdef DEBUG_PROGLESS
		now = std::chrono::system_clock::now();
		std::cout
		<< "[PROGLESS] end fastica "
		<< "\t" << std::chrono::duration_cast<std::chrono::milliseconds>(now-prev).count()
		<< "\ttotal:" << std::chrono::duration_cast<std::chrono::milliseconds>(now-start).count() << std::endl;
		prev = now;
#endif
		return FastICAResult{.W = B.transpose()*Atilda, .Y = Y};
	};

	class EASI {
	
	public:
		ICA::Matrix B;

		EASI(const int size){
			Reng reng(0);
			this->size = size;
			B = RandMatrix(size, reng);
		}

		Vector update(Vector& x){
			Matrix y = B * x;
			Matrix V = y * y.transpose() - Matrix::Identity(size, size) + g(y) * y.transpose() - y * g(y).transpose();
			B = B - EASI_MU * V * B;
			return y.col(0);
		}

	private:

		const double EASI_MU = 0.001953125;
		int size;

		Matrix g(Matrix y){
			return -y.array().tanh().matrix();
		}
	};

	struct EasiResult {
		Matrix W;	// 復元行列
		Matrix Y;	// 復元信号
	};

	EasiResult BatchEasi(Matrix& X) {
		EASI easi(X.rows());
		Matrix Y(X.rows(), X.cols());
		for (int i=0; i<X.cols(); i++){
			Vector x = X.col(i);
			Vector y = easi.update(x);
			Y.col(i) = y;
		}
		return EasiResult{.W = easi.B, .Y = Y};
	}
    
	std::vector<double> ToStdVec(Vector& v1){
		std::vector<double> v2(v1.data(), v1.data() + v1.size());
		return v2;
	}

	void WriteMatrix(std::stringstream& ss, Matrix& mat){
		const auto sample = mat.rows();
		const auto series = mat.cols();
		for (int i=0; i<sample; i++){
			for (int j=0;j<series;j++){
				ss << mat(i,j) << ",";
				if (WRITE_LIMIT < j) break;
			}
			ss << std::endl;
		}
	}

	/**
	 * チェビシェフ多項式の値を計算
	 * オーバーフローを避けるため、漸化式を利用
	 * n: 次数(n>=0) return: T_n(x)の値
	 */
	double EvalChebyt(const double x, const int n){
		double t0 = 1.0;
		double t1 = x;
		if (n == 0) {
			return t0;
		} else if (n == 1) {
			return t1;
		}
		double t2;
		for(int i=2; i<=n; i++){
			t2 = 2*x*t1-t0;
			t0 = t1;
			t1 = t2;
		}
		return t2;
	}

	/**
	 * 縦にベクトルを積む
	 * 必ず1つ以上を渡し、複数の場合、横幅は統一すること
	 */
	Matrix VStack(std::vector<Vector>& vecs){
		const auto num = vecs.size();
		Matrix C(num, vecs.at(0).size());
		for (int i=0; i<num; i++){
			C.row(i) = vecs.at(i);
		}
		return C;
	}

	/**
	 * n:次数、a0：初期値、len：長さ
	 */
	Vector ChebytSeries(const int n, const int len, const double a0){
		Vector S(len);
		double prev = a0;
		for (int i=0; i<len; i++){
			S(i) = prev;
			prev = EvalChebyt(prev, n);
		}
		return S;
	}

	Vector SinSeries(const double w, const int len){
		Vector S(len);
		const double gap = 0.1;
		for (int i=0; i<len; i++){
			S(i) = std::sin(w*(double)i*gap);
		}
		return S;
	}
}

int sample(){
	ICA::Reng reng(0);
	const auto sample = 3;
	const auto series = 10000;
	
	std::vector<ICA::Vector> s(sample);
	#pragma omp parallel for
	for (int i=0; i<sample; i++){
		s.at(i) = ICA::ChebytSeries(i+2, series, 0.2);
	}
	ICA::Matrix S = ICA::VStack(s);

	const ICA::Matrix A = ICA::RandMatrix(sample, reng);
	ICA::Matrix X = A * S;
	auto result = ICA::BatchEasi(X);

	std::stringstream ss;
	ICA::WriteMatrix(ss, S);
	ss << std::endl;
	ICA::WriteMatrix(ss, X);
	ss << std::endl;
	ICA::Matrix YR = result.Y.rightCols(1000);
	ICA::WriteMatrix(ss, YR);

	std::ofstream outputfile("test.csv");
	outputfile << ss.rdbuf();
	outputfile.close();
	return 0;
}