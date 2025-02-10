#include "mfem.hpp"
#include <cmath>

#include "coeff.h"

const double PI = 3.14159265358979323846;

double xi_coeff_func(const mfem::Vector &x) {
  return std::pow(x(0), 2);
}

void vec_xi_coeff_func(const mfem::Vector &x, mfem::Vector &v) {
  v.SetSize(1);
  v[0] = std::pow(x(0), 2);
}

double theta_initial_guess(const mfem::Vector &x) {
  double xi = x[0];
  return 1-std::pow(xi/(PI), 2);
}
