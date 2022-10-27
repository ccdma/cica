#define NDEBUG
#define NPROGLESS

#ifndef COMMIT_ID
	#define COMMIT_ID "undefined"
#endif

#include <cmath>
#include <cassert>
#include <chrono>
#include <thread>
#include "cica.hpp"


int main(){

    cica::cvector v1 = cica::exact_const_powerd_sampling(61, 2, 1);
    cica::cvector v2 = cica::exact_const_powerd_sampling(61, 2, 2);
	
    std::stringstream ss;
    cica::util::write_matrix(ss, v1.real());
	ss << std::endl;
    cica::util::write_matrix(ss, v1.imag());
    ss << std::endl;

    cica::util::write_matrix(ss, v2.real());
	ss << std::endl;
    cica::util::write_matrix(ss, v2.imag());
    ss << std::endl;

    cica::util::save_stream(ss, "data.csv");

	return 0;
}