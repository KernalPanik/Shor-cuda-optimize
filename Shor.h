#pragma once

#include <cstdio>
// #include "types.h" defined from unitary
#include "unitary.h"
#include <functional>

#define __int64 long long 
const double pi = acos(-1.0);

unsigned int Shor(unsigned int N, unsigned int depth_limit, __int64 &timeElapsedForQFT, __int64 &tryAmount,
	std::vector<std::vector<amp>> noisyHadamard);
