#define NDEBUG
// #define NPROGLESS

#include "cica.hpp"

int main(){
	
	const cica::cvector S = cica::weyl_sampling(std::sqrt(0.1), 0.1, 1000);

	const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");

	for (int i=0; i<1000; i++){
		std::cout << S(i).real() << ",";
	}
	std::cout << std::endl;
	for (int i=0; i<1000; i++){
		std::cout << S(i).imag() << ",";
	}
	std::cout << std::endl;
	return 0;
}