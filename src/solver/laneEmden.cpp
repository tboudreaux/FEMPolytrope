#include "mfem.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "mfemUtils.h"
#include "io.h"
#include "coeff.h"

const double epsilon = 1e-12;
const double u_safe_regularization = 0.0;
const double PI = 3.14159265358979323846;

int main(int argc, char* argv[]) {
  // Parse n from command line arguments as a float
  double n = 1.5;
  int order = 1;
  int numElements = 10;
  double root = PI;

  mfem::OptionsParser args(argc, argv);
  args.AddOption(&order, "-o", "--order", "Order Solver to use [1]");
  args.AddOption(&numElements, "-ne", "--numberOfElements", "Number of Elements to use [10]");
  args.AddOption(&n, "-n", "--polytropicIndex", "Poytropic Index to use [1.5]");
  args.AddOption(&root, "-r", "--root", "root to set second tdof at [PI]");
  args.ParseCheck();

  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(numElements, root);

  mfem::H1_FECollection fec(order, 1);
  mfem::FiniteElementSpace fes(&mesh, &fec);

  // This makes some function u over the finite element space
  mfem::GridFunction u(&fes);
  u = 0.0; // This sets the initial conditions of u (0 temp everywhere)
  mfem::FunctionCoefficient initCoeff(
    [root](const mfem::Vector &x) -> double {
      return theta_initial_guess(x, root);
    });
  u.ProjectCoefficient(initCoeff);

  mfem::Array<int> ess_tdof_list;
  mfem::Array<int> ess_bdr(mesh.bdr_attributes.Max());
  ess_bdr = 0;
  ess_bdr[0] = 1;
  ess_bdr[1] = 1;
  fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
  
  mfem::VectorFunctionCoefficient diffCoeff(1, vec_xi_coeff_func);
  mfem::FunctionCoefficient xi_coeff(xi_coeff_func);

  CompositeNonlinearIntegrator compositeIntegrator;
  compositeIntegrator.add_integrator(
      new BilinearIntegratorWrapper(new mfem::DiffusionIntegrator(diffCoeff))
  );
  compositeIntegrator.add_integrator(new NonlinearPowerIntegrator(xi_coeff, n));
   
  mfem::NonlinearForm nfl(&fes);
  nfl.AddDomainIntegrator(&compositeIntegrator);
  nfl.SetEssentialTrueDofs(ess_tdof_list);

  write_solution_to_csv(u, mesh, "laneEmdenNonlinear_presolve.csv");

  mfem::NewtonSolver newtonSolver;
  newtonSolver.SetRelTol(1e-8);
  newtonSolver.SetAbsTol(1e-10);
  newtonSolver.SetMaxIter(200);
  newtonSolver.SetPrintLevel(1);
  mfem::GMRESSolver gmresSolver;
  gmresSolver.SetRelTol(1e-8);
  gmresSolver.SetMaxIter(200);
  gmresSolver.SetPrintLevel(1);
  newtonSolver.SetSolver(gmresSolver);
  newtonSolver.SetOperator(nfl);
  newtonSolver.SetAdaptiveLinRtol();

  mfem::Vector b;
  b.SetSize(1);
  b = 0.0;
  newtonSolver.Mult(b, u);

  write_solution_to_csv(u, mesh, "laneEmdenNonlinear.csv");

}
