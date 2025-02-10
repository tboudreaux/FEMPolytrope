#include "mfem.hpp"
#include <string>
#include <iostream>
#include <cmath>

#include "mfemUtils.h"

NonlinearPowerIntegrator::NonlinearPowerIntegrator(
 mfem::FunctionCoefficient &coeff,
 double n) : coeff_(coeff), polytropicIndex(n) { 

}

void NonlinearPowerIntegrator::AssembleElementVector(
 const mfem::FiniteElement &el,
 mfem::ElementTransformation &Trans,
 const mfem::Vector &elfun,
 mfem::Vector &elvect) {
  
  const mfem::IntegrationRule *ir = &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + 3);
  int dof = el.GetDof();
  elvect.SetSize(dof);
  elvect = 0.0;

  mfem::Vector shape(dof);

  for (int iqp = 0; iqp < ir->GetNPoints(); iqp++) {
    mfem::IntegrationPoint ip = ir->IntPoint(iqp);
    Trans.SetIntPoint(&ip);
    double weight = ip.weight * Trans.Weight();

    el.CalcShape(ip, shape);

    double u_val = 0.0;
    for (int j = 0; j < dof; j++) {
      u_val += elfun(j) * shape(j);
    }
    double u_safe = std::max(u_val, 0.0);
    double u_nl = std::pow(u_safe, polytropicIndex);

    double coeff_val = coeff_.Eval(Trans, ip);
    double x2_u_nl = coeff_val * u_nl;

    for (int i = 0; i < dof; i++){
       elvect(i) += shape(i) * x2_u_nl * weight;
    }
  }
}

void NonlinearPowerIntegrator::AssembleElementGrad (
 const mfem::FiniteElement &el,
 mfem::ElementTransformation &Trans,
 const mfem::Vector &elfun,
 mfem::DenseMatrix &elmat) {

  const mfem::IntegrationRule *ir = &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + 3);
  int dof = el.GetDof();
  elmat.SetSize(dof);
  elmat = 0.0;
  mfem::Vector shape(dof);

  for (int iqp = 0; iqp < ir->GetNPoints(); iqp++) {
    mfem::IntegrationPoint ip = ir->IntPoint(iqp);
    Trans.SetIntPoint(&ip);
    double weight = ip.weight * Trans.Weight();

    el.CalcShape(ip, shape);

    double u_val = 0.0;

    for (int j = 0; j < dof; j++) {
      u_val += elfun(j) * shape(j);
    }
    double coeff_val = coeff_.Eval(Trans, ip);
    

    // Calculate the Jacobian
    double u_safe = std::max(u_val, 0.0);
    double d_u_nl = coeff_val * polytropicIndex * std::pow(u_safe, polytropicIndex - 1);
    double x2_d_u_nl = d_u_nl;

    for (int i = 0; i < dof; i++) {
      for (int j = 0; j < dof; j++) {
        elmat(i, j) += shape(i) * x2_d_u_nl * shape(j) * weight;
      }
    }

  }
}

BilinearIntegratorWrapper::BilinearIntegratorWrapper(
 mfem::BilinearFormIntegrator *integratorInput
 ) : integrator(integratorInput) { }

BilinearIntegratorWrapper::~BilinearIntegratorWrapper() {
  delete integrator;
}

void BilinearIntegratorWrapper::AssembleElementVector(
 const mfem::FiniteElement &el,
 mfem::ElementTransformation &Trans,
 const mfem::Vector &elfun,
 mfem::Vector &elvect) {
  int dof = el.GetDof();
  mfem::DenseMatrix elMat(dof);
  integrator->AssembleElementMatrix(el, Trans, elMat);
  elvect.SetSize(dof);
  elvect = 0.0;
  for (int i = 0; i < dof; i++)
  {
     double sum = 0.0;
     for (int j = 0; j < dof; j++)
     {
        sum += elMat(i, j) * elfun(j);
     }
     elvect(i) = sum;
  }
}

void BilinearIntegratorWrapper::AssembleElementGrad(const mfem::FiniteElement &el,
 mfem::ElementTransformation &Trans,
 const mfem::Vector &elfun,
 mfem::DenseMatrix &elmat) {
  int dof = el.GetDof();
  elmat.SetSize(dof, dof);
  elmat = 0.0;
  integrator->AssembleElementMatrix(el, Trans, elmat);
}

CompositeNonlinearIntegrator::CompositeNonlinearIntegrator() { }


CompositeNonlinearIntegrator::~CompositeNonlinearIntegrator() {
  for (size_t i = 0; i < integrators.size(); i++) {
    delete integrators[i];
  }
}

void CompositeNonlinearIntegrator::add_integrator(mfem::NonlinearFormIntegrator *integrator) {
  integrators.push_back(integrator);
}

void CompositeNonlinearIntegrator::AssembleElementVector(
 const mfem::FiniteElement &el,
 mfem::ElementTransformation &Trans,
 const mfem::Vector &elfun,
 mfem::Vector &elvect) {
  int dof = el.GetDof();
  elvect.SetSize(dof);
  elvect = 0.0;
  mfem::Vector temp(dof);

  for (size_t i = 0; i < integrators.size(); i++) {
    temp= 0.0;
    integrators[i]->AssembleElementVector(el, Trans, elfun, temp);
    elvect.Add(1.0, temp);
  }
}

void CompositeNonlinearIntegrator::AssembleElementGrad(
 const mfem::FiniteElement &el,
 mfem::ElementTransformation &Trans,
 const mfem::Vector &elfun,
 mfem::DenseMatrix &elmat) {
  int dof = el.GetDof();
  elmat.SetSize(dof, dof);
  elmat = 0.0;
  mfem::DenseMatrix temp(dof);
  temp.SetSize(dof, dof);
  for (size_t i = 0; i < integrators.size(); i++) {
    temp = 0.0;
    integrators[i] -> AssembleElementGrad(el, Trans, elfun, temp);
    elmat.Add(1.0, temp);
  }
}
