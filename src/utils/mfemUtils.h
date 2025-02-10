#include "mfem.hpp"
#include <string>


void write_solution_to_csv(const mfem::GridFunction &u, const mfem::Mesh &mesh, const std::string &filename);

class NonlinearPowerIntegrator: public mfem::NonlinearFormIntegrator {
private: 
  mfem::FunctionCoefficient coeff_;
  double polytropicIndex;
public:
  NonlinearPowerIntegrator(mfem::FunctionCoefficient &coeff, double n);

  virtual void AssembleElementVector(const mfem::FiniteElement &el, mfem::ElementTransformation &Trans, const mfem::Vector &elfun, mfem::Vector &elvect) override;
  virtual void AssembleElementGrad (const mfem::FiniteElement &el, mfem::ElementTransformation &Trans, const mfem::Vector &elfun, mfem::DenseMatrix &elmat) override;
};

class BilinearIntegratorWrapper : public mfem::NonlinearFormIntegrator
{
private:
   mfem::BilinearFormIntegrator *integrator;
public:
   BilinearIntegratorWrapper(mfem::BilinearFormIntegrator *integratorInput);

   virtual ~BilinearIntegratorWrapper() ;

   virtual void AssembleElementVector(const mfem::FiniteElement &el, mfem::ElementTransformation &Trans, const mfem::Vector &elfun, mfem::Vector &elvect) override;
   virtual void AssembleElementGrad(const mfem::FiniteElement &el,mfem::ElementTransformation &Trans, const mfem::Vector &elfun, mfem::DenseMatrix &elmat) override;
};

class CompositeNonlinearIntegrator: public mfem::NonlinearFormIntegrator {
  private:
    std::vector<mfem::NonlinearFormIntegrator*> integrators;
  public:
    CompositeNonlinearIntegrator();

    virtual ~CompositeNonlinearIntegrator();

    void add_integrator(mfem::NonlinearFormIntegrator *integrator);

    virtual void AssembleElementVector(const mfem::FiniteElement &el, mfem::ElementTransformation &Trans, const mfem::Vector &elfun, mfem::Vector &elvect) override;
    virtual void AssembleElementGrad(const mfem::FiniteElement &el, mfem::ElementTransformation &Trans, const mfem::Vector &elfun, mfem::DenseMatrix &elmat) override;
};
