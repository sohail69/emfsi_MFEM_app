#pragma once
#include "mfem.hpp"
#include <functional>

using namespace std;
using namespace mfem;

/*****************************************\
! Author: Sohail Rathore
! Date  : 2025-05-18
!
! Integrates a matrix contaction against
! the gradient of the shape function
! S_ij dN_im/dx_j
!
\*****************************************/
class gradContractionIntegrator : public DeltaLFIntegrator
{
private:
   Vector shape, Qvec;
   MatrixCoefficient &stress;
   DenseMatrix dshape, dshapedxt, rho_ij;

public:
   /// Constructs the domain integrator (Q, grad v)
   VectorDomainLFGradIntegrator(MatrixCoefficient &QF)
      : DeltaLFIntegrator(QF), Q(QF) {}

   bool SupportsDevice() const override { return false; }

   /** Given a particular Finite Element and a transformation (Tr)
       computes the element right hand side element vector, elvect. */
   void AssembleRHSElementVect(const FiniteElement &el,
                               ElementTransformation &Tr,
                               Vector &elvect) override;

   using LinearFormIntegrator::AssembleRHSElementVect;
};

//
// Assemble the element Vector
//
void gradContractionIntegrator::AssembleRHSElementVect(const FiniteElement &el
                                                     , ElementTransformation &Tr
                                                     , Vector &elvect)
{
  int L=0;
  int dof = el.GetDof();
  int spaceDim = Tr.GetSpaceDim();
  dshape.SetSize(dof, spaceDim);
  dshapedxt.SetSize(dof, spaceDim);

  elvect.SetSize(dof);
  elvect = 0.0;

  const IntegrationRule *ir = GetIntegrationRule(el, Tr);
  if (ir == NULL)
  {
    int intorder = 2 * el.GetOrder();
    ir = &IntRules.Get(el.GetGeomType(), intorder);
  }

  for (int i = 0; i < ir->GetNPoints(); i++)
  {
    const IntegrationPoint &ip = ir->IntPoint(i);
    Tr.SetIntPoint(&ip);
    el.CalcPhysDShape(Tr, dshape);
    Mult(dshape, Tr.AdjugateJacobian(), dshapedxt);
    stress.Eval(rho_ij, Tr, ip);

    for(int I=0; I<dim; I++){
      for(int J=0; J<dim; J++){
        for(int K=0; K<dof; K++){
          L = I*dof + K;
          elvect(L) += rho_ij(I,J) * dshape(K,J) * ip.weight * Tr.Weight();
        }
      }
    }
/*
    Q.Eval(Qvec, Tr, ip);
    Qvec *= ip.weight * Tr.Weight();
    dshape.AddMult(Qvec, elvect);
*/
  }
}
