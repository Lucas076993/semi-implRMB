#ifndef CG_H
#define CG_H

#include <vector>

#include "algebra.h"

std::vector<double> CG(const algebra::matriz<double>&, const std::vector<double>&, const std::vector<double>&, const std::size_t, const double);

#endif
