// #define EIGEN_USE_BLAS
// #define NDEBUG
// #define NPROGLESS
// #define NPARALLELIZE

#ifdef NPARALLELIZE
	#define EIGEN_DONT_PARALLELIZE
#endif

#ifdef NDEBUG
	#define EIGEN_NO_DEBUG
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>
#include <sstream>
#include <complex>
#include <cassert>
#include <Eigen/Dense>

namespace cica {

	using matrix = Eigen::MatrixXd;
	using vector = Eigen::VectorXd;
	using cmatrix = Eigen::MatrixXcd;
	using cvector = Eigen::VectorXcd;
	using imatrix = Eigen::MatrixXi;
	using ivector = Eigen::VectorXi;
	using dcomplex = std::complex<double>;
	using random_engine = std::mt19937; 

	const int WRITE_LIMIT = 10000;

	/**
	 * -0.5~0.5までの一様乱数からなる正方行列を生成
	 */ 
	matrix random_uniform_matrix(const int size, random_engine& engine){
		std::uniform_real_distribution<double> distribution(-0.5, 0.5);
		auto generator = [&] (double dummy) {return distribution(engine);};
		return matrix::Zero(size, size).unaryExpr(generator);
	};

	/**
	 * ガウス分布による行列を生成
	 * stddev: 標準偏差
	 */ 
	matrix gauss_matrix(const int rows, const int cols, const double stddev, random_engine& engine){
		std::normal_distribution<double> distribution(0.0, stddev);
		auto generator = [&] (double dummy) {return distribution(engine);};
		return matrix::Zero(rows, cols).unaryExpr(generator);
	};

	/**
	 * ±1で表現されるランダムなビット行列を生成
	 * imatrixで返るのでdoubleなどを利用する場合、random_bits().cast<double>()などでcastすること
	 */ 
	imatrix random_bits(const int rows, const int cols, random_engine& engine){
		std::uniform_int_distribution<int> distribution(0, 1);
		auto generator = [&](int dummy) {
			const int bit = distribution(engine);
			return bit == 1 ? 1 : -1;
		};
		return imatrix::Zero(rows, cols).unaryExpr(generator);
	}

	/**
	 * 系列の中心化を行う
	 * 横に並ぶ値の平均が0になるように修正される
	 */ 
	matrix centerize(const matrix& M){
		return M.colwise() - M.rowwise().mean();
	}

	/**
	 * 分散共分散行列を計算する
	 * E[(M−μ)(M−μ)^⊤]である
	 */ 
	matrix cov(const matrix& M){
		const matrix M_center = centerize(M);
		return (M_center * M_center.transpose()) / double(M_center.cols() - 1);
	}

    
	std::vector<double> to_std_vector(const vector& v1){
		std::vector<double> v2(v1.data(), v1.data() + v1.size());
		return v2;
	}

	void write_matrix(std::stringstream& ss, const matrix& mat){
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

	void save_stream(std::stringstream& ss, std::string filename){
		std::ofstream outputfile(filename);
		outputfile << ss.rdbuf();
		outputfile.close();
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
	 * 第一種チェビシェフ多項式による系列（漸化式）
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
#ifndef NPARALLELIZE
		#pragma omp parallel for
#endif
		for (int i=0; i<len; i++){
			S(i) = std::sin(w*(double)i*gap);
		}
		return S;
	}

	/**
	 * パワー一定のカオス符号
	 * n: n倍角
	 * rad_0: 初期値の偏角(rad)
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
	imatrix estimate_circulant_matrix(const matrix& A, const matrix& W){
		matrix G = W * A;
		imatrix P = imatrix::Zero(G.rows(), G.cols());
#ifndef NPARALLELIZE
		#pragma omp parallel for
#endif
		for(int i=0; i<G.rows(); i++){
			vector row = G.row(i);
			vector::Index maxId;
			row.cwiseAbs().maxCoeff(&maxId);
			const double x = row(maxId);
			P(i, maxId) = (x > 0) ? 1 : -1;
		}
		return P;
	}

	/**
	 * CTE(cross talk error)
	 * A: 混合行列
	 * W: 復元行列
	 * https://www.jstage.jst.go.jp/article/elex/5/14/5_14_510/_pdf
	 */ 
	double cross_talk_error(const matrix& A, const matrix& W){
		const matrix C = W * A;
		const int size = A.cols();
		const auto row_cte = [](const vector& vec){
			const vector absvec = vec.cwiseAbs();
			const double max = absvec.maxCoeff();
			const double sum = absvec.sum();
			return sum/max-1;
		};
		double cte_sum = 0.0;
#ifndef NPARALLELIZE
		#pragma omp parallel for reduction(+:cte_sum)
#endif
		for (int i=0; i<size; i++){
			cte_sum += row_cte(C.row(i));
			cte_sum += row_cte(C.col(i));
		}
		return cte_sum;
	}

	/**
	 * CTEを異なる信号数でも比較できるように正規化したもの
	 * 
	 * size-1で割っているのは正解部分は-1で既に除かれているから
	 * 
	 * (referenceなし)
	 */
	double normal_cross_talk_error(const matrix& A, const matrix& W){
		const auto cte = cross_talk_error(A, W);
		const int size = A.cols();
		return cte/((size-1)*size);
	}

	/**
	 * BER(bit error rate)
	 * B1とB2のビット列(±1で表現)のBERを計算する(横ベクトルの順序は予め対応付けすること)
	 */
	double bit_error_rate(const imatrix& B1, const imatrix& B2){
		return (B1-B2).cwiseAbs().cast<double>().mean()/2.0;
	}

	/**
	 * MSE(mean squared error,平均2乗誤差)
	 * SとS2の2乗誤差を計算する(横ベクトルの順序は予め対応付けすること)
	 */
	double mean_squared_error(const matrix& S1, const matrix& S2){
		return (S1-S2).array().pow(2).mean();
	}

	/**
	 * 横ベクトル同士の内積を計算
	 * 
	 * normalize=true: 分散で正規化する
	 * 		→ 完全に直交する場合は単位行列になるはず
	 * 
	 * return: Xの縦の長さ分の正方行列
	 */ 
	matrix correlation_matrix(const matrix& X, const bool normalize=true){
		const int size = X.rows();
		const int n = X.cols();
		matrix M(size, size);
		for (int i=0; i<size; i++){
			for(int j=0; j<i+1; j++){
				const vector rowi = X.row(i);
				const vector rowj = X.row(j);
				double correlation = rowi.dot(rowj) / (double)n;	// 相関
				if (normalize){
					const double var = rowi.cwiseAbs().dot(rowj.cwiseAbs()) / (double)n;	// 分散
					correlation /= var;
				}
				M(i,j) = correlation;
				if(j!=i){ // 対角成分以外なら
					M(j,i) = correlation;
				}
			}
		}
		return M;
	}
}

namespace cica::fastica {

	const int LOOP_MAX = 100;

	/**
	 * 正方行列でなくてはいけない
	 * i番目を0~i-1番目までの縦ベクトルの直行空間に射影＋長さの正規化
	 */
	void _normalize(matrix& M, const int i){
		const auto size = M.cols();
		if (i>0){
			M.col(i) = M.col(i) - M.block(0, 0, size, i) * M.block(0, 0, size, i).transpose() * M.col(i);
		}
		M.col(i) = M.col(i) / std::sqrt(M.col(i).squaredNorm());
	}

	struct result {
		matrix W;	// 復元行列
		matrix Y;	// 復元信号
		ivector loop; 	// 不動点法のループ回数
	};

	/**
	 * 目的関数に使用する関数g
	 */ 
	struct objective_func {
		const std::function<double(double)> g;
		const std::function<double(double)> g2;

		objective_func(): 
			objective_func(
				[](double bx) { return std::pow(bx, 3); },
				[](double bx) { return 3*std::pow(bx, 2); }
			)
		{};

		objective_func(
			const std::function<double(double)> g,	// 不動点法におけるg()
			const std::function<double(double)> g2	// g()の微分
		): g(g), g2(g2) {};
	};

	// 4次キュムラント
	const auto kum4 = objective_func();

	const auto cpnscp1 = objective_func(
		[](double x) { return   39.7*std::pow(x, 3) - 14.5*x; },
		[](double x) { return 39.7*3*std::pow(x, 2) -   14.5; }
	);

	const auto cpnscp2 = objective_func(
		[](double x) { return   158.7*std::pow(x, 5) -   136.5*std::pow(x, 3) + 23.2*x; },
		[](double x) { return 158.7*5*std::pow(x, 4) - 136.5*3*std::pow(x, 2) +   23.2; }
	);

	/**
	 * implements of FastICA
	 * 
	 * X: 観測信号
	 * 
	 * Xについて、内部で中心化は行うが先にに中心化されていることが望ましい（元信号Sの中心化ができていれば、混合されたXも自然と中心化されるはず）
	 * 		→ 中心化されているものを扱っていれば、2乗和誤差などの計算でズレが生じない
	 * 
	 * [reference]
	 * https://ieeexplore.ieee.org/document/761722
	 * http://manabukano.brilliant-future.net/document/text-ICA.pdf
	 */
	result fastica(const matrix& X, const objective_func& func=kum4) {

#ifndef NPROGLESS
	std::chrono::system_clock::time_point start, prev, now;
	start = std::chrono::system_clock::now();
	prev = start;
	now = start;
	std::cout 
	<< "[PROGLESS] start fastica session"
	<< "\t" << std::chrono::duration_cast<std::chrono::milliseconds>(now-prev).count() << std::endl;
	prev = now;
#endif

	cica::random_engine random_engine(0);
	const auto signals = X.rows();
	const auto samplings = X.cols();
	
	const matrix X_center = centerize(X);
	const matrix X_cov = cov(X_center);	// 分散共分散行列を作成

	Eigen::SelfAdjointEigenSolver<matrix> es(X_cov);
	if (es.info() != Eigen::Success) abort();

	const vector lambdas = es.eigenvalues().real();
	const matrix P = es.eigenvectors().real();
	const matrix Atilda = lambdas.cwiseSqrt().asDiagonal().inverse() * P.transpose();
	const matrix X_whiten = Atilda * X_center;	// 無相関化

#ifndef NPROGLESS
	now = std::chrono::system_clock::now();
	std::cout 
	<< "[PROGLESS] start fixed point method"
	<< "\t" << std::chrono::duration_cast<std::chrono::milliseconds>(now-prev).count() << std::endl;
	prev = now;
#endif

	assert(cov(X_whiten).isApprox(matrix::Identity(cov(X_whiten).rows(), cov(X_whiten).cols())));

	const auto g = func.g;
	const auto g2 = func.g2;

	const auto I = X_whiten.rows();
	Eigen::VectorXi loop(I);
	auto B = random_uniform_matrix(I, random_engine);

	for(int i=0; i<I; i++){
		_normalize(B, i);
	}

	assert((B * B.transpose()).isApprox(matrix::Identity(B.rows(), B.cols())));

	for(int i=0; i<I; i++){
		for(int j=0; j<LOOP_MAX; j++){
			loop(i) = j+1;	// ループ回数を記録
			const vector prevBi = B.col(i);	// 値のコピー

			const auto collen = X_whiten.cols();
			matrix ave(I, collen);
#ifndef NPARALLELIZE
			#pragma omp parallel for
#endif
			for(int k=0; k<collen; k++){	// 不動点法による更新
				const vector x = X_whiten.col(k);
				ave.col(k) = g(x.dot(B.col(i)))*x - g2(x.dot(B.col(i)))*B.col(i);  
			}
			B.col(i) = ave.rowwise().mean();
			_normalize(B, i);
			const auto diff = std::abs(prevBi.dot(B.col(i)));
			if (1.0 - 1.e-8 < diff && diff < 1.0 + 1.e-8) break;
#ifndef NPROGLESS
			if (j==LOOP_MAX-1) printf("[WARN] loop limit exceeded\n");
#endif
		}

#ifndef NPROGLESS
		now = std::chrono::system_clock::now();
		std::cout
		<< "[PROGLESS] end loop " << i
		<< "\t" << std::chrono::duration_cast<std::chrono::milliseconds>(now-prev).count() << std::endl;
		prev = now;
#endif
		}
		matrix Y = B.transpose() * X_whiten;

#ifndef NPROGLESS
		now = std::chrono::system_clock::now();
		std::cout
		<< "[PROGLESS] end fastica "
		<< "\t" << std::chrono::duration_cast<std::chrono::milliseconds>(now-prev).count()
		<< "\ttotal:" << std::chrono::duration_cast<std::chrono::milliseconds>(now-start).count() << std::endl;
		prev = now;
#endif
		return result{.W = B.transpose()*Atilda, .Y = Y, .loop = loop};
	};
}

namespace cica::easi {

	/**
	 * implements of EASI
	 * 
	 * X: 観測信号
	 * 
	 * [reference]
	 * https://ieeexplore.ieee.org/document/553476
	 * https://www.ieice.org/ken/download/200703059ATF/
	 */ 
	class easi {
	
	public:
		cica::matrix W;

		easi(const int size): size(size){
			random_engine random_engine(0);
			W = random_uniform_matrix(size, random_engine);
		}

		vector update(const vector& x){
			matrix y = W * x;
			matrix V = y * y.transpose() - matrix::Identity(size, size) + g(y) * y.transpose() - y * g(y).transpose();
			W = W - EASI_MU * V * W;
			return y.col(0);
		}

	private:
		const double EASI_MU = 0.001953125;
		const int size;

		matrix g(const matrix& y){
			return -y.array().tanh().matrix();
		}
	};

	struct result {
		matrix W;	// 復元行列
		matrix Y;	// 復元信号
	};

	/**
	 * EASI(バッチ処理)
	 * シュミレーション時はfasticaと同じように扱える分、easiクラスよりも利用しやすい
	 * 
	 * final_recover: trueの場合、学習後の復元行列でXを再計算する。逆にfalseの場合、resultの復元信号は復元行列に対応しないので注意すること
	 */ 
	result batch_easi(const matrix& X, const bool final_recover=true) {
		easi easi(X.rows());
		matrix Y(X.rows(), X.cols());
		for (int i=0; i<X.cols(); i++){
			vector x = X.col(i);
			vector y = easi.update(x);
			Y.col(i) = y;
		}
		if (final_recover){
			Y = easi.W * X;
		}
		return result{.W = easi.W, .Y = Y};
	}
}