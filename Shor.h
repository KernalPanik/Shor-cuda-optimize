#pragma once

#include <cstdio>
// #include "types.h" defined from unitary
#include "unitary.h"
#include <functional>

const double pi = acos(-1.0);

unsigned int Shor(unsigned int N, unsigned int depth_limit, __int64 &timeElapsedForQFT);