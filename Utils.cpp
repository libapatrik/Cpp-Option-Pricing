#include "Utils.h"
#include<cmath>

// Standard normal CDF
double Utils::stdNormCdf(double x)
{
    return 0.5 * erfc(-x / sqrt(2.0));
}
