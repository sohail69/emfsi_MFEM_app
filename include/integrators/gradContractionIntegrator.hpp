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
   VectorCoefficient &Q;
   DenseMatrix dshape;

public:
   /// Constructs the domain integrator (Q, grad v)
   VectorDomainLFGradIntegrator(VectorCoefficient &QF)
      : DeltaLFIntegrator(QF), Q(QF) { }

   bool SupportsDevice() const override { return true; }

   /// Method defining assembly on device
   void AssembleDevice(const FiniteElementSpace &fes,
                       const Array<int> &markers,
                       Vector &b) override;

   /** Given a particular Finite Element and a transformation (Tr)
       computes the element right hand side element vector, elvect. */
   void AssembleRHSElementVect(const FiniteElement &el,
                               ElementTransformation &Tr,
                               Vector &elvect) override;

   void AssembleDeltaElementVect(const FiniteElement &fe,
                                 ElementTransformation &Trans,
                                 Vector &elvect) override;

   using LinearFormIntegrator::AssembleRHSElementVect;
};

//
// Assemble the 
//
