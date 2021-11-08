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
#include <complex>
#include <Eigen/Dense>

namespace cica {

	using matrix = Eigen::MatrixXd;
	using vector = Eigen::VectorXd;
	using cmatrix = Eigen::MatrixXcd;
	using cvector = Eigen::VectorXcd;
	using dcomplex = std::complex<double>;
	using reng = std::mt19937; 

	const int FASTICA_LOOP_MAX = 500;
	const int WRITE_LIMIT = 10000;

	/**
	 * -0.5~0.5までの一様乱数からなる正方行列を生成
	 */ 
	matrix rand_matrix(int size, reng& engine){
		std::uniform_real_distribution<double> distribution(-0.5, 0.5);
		auto generator = [&] (double dummy) {return distribution(engine);};
		return matrix::Zero(size, size).unaryExpr(generator);
	};

	/**
	 * 正方行列でなくてはいけない
	 * i番目を0~i-1番目までの縦ベクトルの直行空間に射影＋長さの正規化
	 */
	void normalize(matrix& M, int i){
		const auto size = M.cols();
		if (i>0){
			M.col(i) = M.col(i) - M.block(0, 0, size, i) * M.block(0, 0, size, i).transpose() * M.col(i);
		}
		M.col(i) = M.col(i) / std::sqrt(M.col(i).squaredNorm());
	}

	/**
	 * 系列の中心化を行う
	 * 横に並ぶ値の平均が0になるように修正される
	 */ 
	matrix centerize(matrix& M){
		return M.colwise() - M.rowwise().mean();
	}

	struct FastICAResult {
		matrix W;	// 復元行列
		matrix Y;	// 復元信号
		Eigen::VectorXi loop; 	// 不動点法のループ回数
	};

	/**
	 * X: 内部で中心化は行うが、すでに中心化されていることが望ましい（元信号Sの中心化ができていれば、混合されたXも自然と中心化されるはず）
	 */
	FastICAResult fast_ica(matrix& X) {

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

	cica::reng reng(0);
	const auto signals = X.rows();
	const auto samplings = X.cols();
	
	const matrix X_center = centerize(X);
	const matrix X_cov = (X_center * X_center.transpose()) / double(X_center.cols() - 1);	// 分散共分散行列を作成

	Eigen::SelfAdjointEigenSolver<matrix> es(X_cov);
	if (es.info() != Eigen::Success) abort();

	vector lambdas = es.eigenvalues().real();
	matrix P = es.eigenvectors().real();
	matrix Atilda = lambdas.cwiseSqrt().asDiagonal().inverse() * P.transpose();
	matrix X_whiten = Atilda * X_center;	// 無相関化

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

	const auto g = [](double bx) { return std::pow(bx, 3); };	// 不動点法に4次キュムラントを使用する
	const auto g2 = [](double bx) { return 3*std::pow(bx, 2); };	// g()の微分

	const auto I = X_whiten.rows();
	Eigen::VectorXi loop(I);
	auto B = rand_matrix(I, reng);

	for(int i=0; i<I; i++){
		normalize(B, i);
	}

#ifdef DEBUG_MATRIX
		// 単位行列であることを確認
	std::cout << B * B.transpose() << std::endl;
#endif

	for(int i=0; i<I; i++){
		for(int j=0; j<FASTICA_LOOP_MAX; j++){
			loop(i) = j+1;	// ループ回数を記録
			const vector prevBi = B.col(i);	// 値のコピー

			const auto collen = X_whiten.cols();
			matrix ave(I, collen);
			#pragma omp parallel for
			for(int k=0; k<collen; k++){	// 不動点法による更新
				const vector x = X_whiten.col(k);
				ave.col(k) = g(x.dot(B.col(i)))*x - g2(x.dot(B.col(i)))*B.col(i);  
			}
			B.col(i) = ave.rowwise().mean();
			normalize(B, i);
			const auto diff = std::abs(prevBi.dot(B.col(i)));
			if (1.0 - 1.e-8 < diff && diff < 1.0 + 1.e-8) break;
#ifdef DEBUG_PROGLESS
			if (j==FASTICA_LOOP_MAX-1) printf("[WARN] loop limit exceeded\n");
#endif
		}

#ifdef DEBUG_PROGLESS
		now = std::chrono::system_clock::now();
		std::cout
		<< "[PROGLESS] end loop " << i
		<< "\t" << std::chrono::duration_cast<std::chrono::milliseconds>(now-prev).count() << std::endl;
		prev = now;
#endif
		}
		matrix Y = B.transpose() * X_whiten;

#ifdef DEBUG_PROGLESS
		now = std::chrono::system_clock::now();
		std::cout
		<< "[PROGLESS] end fastica "
		<< "\t" << std::chrono::duration_cast<std::chrono::milliseconds>(now-prev).count()
		<< "\ttotal:" << std::chrono::duration_cast<std::chrono::milliseconds>(now-start).count() << std::endl;
		prev = now;
#endif
		return FastICAResult{.W = B.transpose()*Atilda, .Y = Y, .loop = loop};
	};

	class EASI {
	
	public:
		cica::matrix B;

		EASI(const int size){
			reng reng(0);
			this->size = size;
			B = rand_matrix(size, reng);
		}

		vector update(vector& x){
			matrix y = B * x;
			matrix V = y * y.transpose() - matrix::Identity(size, size) + g(y) * y.transpose() - y * g(y).transpose();
			B = B - EASI_MU * V * B;
			return y.col(0);
		}

	private:
		const double EASI_MU = 0.001953125;
		int size;

		matrix g(matrix y){
			return -y.array().tanh().matrix();
		}
	};

	struct EasiResult {
		matrix W;	// 復元行列
		matrix Y;	// 復元信号
	};

	EasiResult batch_easi(matrix& X) {
		EASI easi(X.rows());
		matrix Y(X.rows(), X.cols());
		for (int i=0; i<X.cols(); i++){
			vector x = X.col(i);
			vector y = easi.update(x);
			Y.col(i) = y;
		}
		return EasiResult{.W = easi.B, .Y = Y};
	}
    
	std::vector<double> to_std_vector(vector& v1){
		std::vector<double> v2(v1.data(), v1.data() + v1.size());
		return v2;
	}

	void write_matrix(std::stringstream& ss, matrix& mat){
		const auto signals = mat.rows();
		const auto samplings = mat.cols();
		for (int i=0; i<signals; i++){
			for (int j=0;j<samplings;j++){
				ss << mat(i,j);
				if (std::min<int>(WRITE_LIMIT, samplings) == j+1) break;
				ss << ",";
			}
			ss << std::endl;
		}
	}

	enum class cheby {
		T,	// 第一種
		U	// 第二種
	};

	/**
	 * チェビシェフ多項式の値を計算
	 * オーバーフローを避けるため、漸化式を利用
	 * n: 次数(n>=0) return: T_n(x)の値
	 */
	double eval_chebyshev(const int n, const double x, const cheby type = cheby::T){
		double c0 = 1.0;
		double c1 = type == cheby::T ? x : 2*x;
		if (n == 0) {
			return c0;
		} else if (n == 1) {
			return c1;
		}
		double c2;
		for(int i=2; i<=n; i++){
			c2 = 2*x*c1-c0;
			c0 = c1;
			c1 = c2;
		}
		return c2;
	}

	/**
	 * 縦にベクトルを積む
	 * 必ず1つ以上を渡し、複数の場合、横幅は統一すること
	 */
	matrix vstack(std::vector<vector>& vecs){
		const auto num = vecs.size();
		matrix C(num, vecs.at(0).size());
		for (int i=0; i<num; i++){
			C.row(i) = vecs.at(i);
		}
		return C;
	}

	/**
	 * n:次数、a0：初期値、len：長さ
	 */
	vector chebyt_sampling(const int n, const int len, const double a0){
		vector S(len);
		double prev = a0;
		for (int i=0; i<len; i++){
			S(i) = prev;
			prev = eval_chebyshev(n, prev);
		}
		return S;
	}

	/**
	 * w: 角周波数、len: 長さ
	 */
	vector sine_sampling(const double w, const int len){
		vector S(len);
		const double gap = 0.1;
		for (int i=0; i<len; i++){
			S(i) = std::sin(w*(double)i*gap);
		}
		return S;
	}

	/**
	 * パワー一定のカオス符号
	 */ 
	cvector const_powerd_sampling(const int n, const double rad_0, const int len){
		cvector S(len);
		dcomplex a_0(std::cos(rad_0), std::sin(rad_0));
		S(0) = a_0;
		for (int i=1; i<len; i++){
			const double real = eval_chebyshev(n, a_0.real());
			const double imag = eval_chebyshev(n-1, a_0.real(), cheby::U)*a_0.imag();
			a_0 = dcomplex(real, imag);
			S(i) = a_0;
		}
		return S;
	}

	/**
	 * S ≒ P.T*Y なる循環行列Pを生成します
	 * 
	 * A: 混合行列
	 * W: 復元行列
	 */
	matrix simple_circulant_P(matrix& A, matrix& W){
		matrix G = W * A;
		matrix P = matrix::Zero(G.rows(), G.cols());
		#pragma omp parallel for
		for(int i=0; i<G.rows(); i++){
			vector row = G.row(i);
			vector::Index maxId;
			row.cwiseAbs().maxCoeff(&maxId);
			const double x = row(maxId);
			P(i, maxId) = (x > 0) ? 1 : -1;
		}
		return P;
	}
}
