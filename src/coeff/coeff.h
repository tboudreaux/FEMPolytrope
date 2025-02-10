#include "mfem.hpp"
#include <cmath>

double xi_coeff_func(const mfem::Vector &x);

void vec_xi_coeff_func(const mfem::Vector &x, mfem::Vector &v);

double theta_initial_guess(const mfem::Vector &x, double root);
