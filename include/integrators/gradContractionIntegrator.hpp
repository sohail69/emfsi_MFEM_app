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
   MatrixCoefficient &Q;
   DenseMatrix dshape;

public:
   /// Constructs the domain integrator (Q, grad v)
   VectorDomainLFGradIntegrator(MatrixCoefficient &QF)
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
// Assemble the element Vector
//
void ElasticityIntegrator::AssembleElementMatrix(
   const FiniteElement &el, ElementTransformation &Trans, DenseMatrix &elmat)
{
   int dof = el.GetDof();
   int dim = el.GetDim();
   real_t w, L, M;

   MFEM_ASSERT(dim == Trans.GetSpaceDim(), "");

#ifdef MFEM_THREAD_SAFE
   DenseMatrix dshape(dof, dim), gshape(dof, dim), pelmat(dof);
   Vector divshape(dim*dof);
#else
   dshape.SetSize(dof, dim);
   gshape.SetSize(dof, dim);
   pelmat.SetSize(dof);
   divshape.SetSize(dim*dof);
#endif

   elmat.SetSize(dof * dim);

   const IntegrationRule *ir = GetIntegrationRule(el, Trans);
   if (ir == NULL)
   {
      int order = 2 * Trans.OrderGrad(&el); // correct order?
      ir = &IntRules.Get(el.GetGeomType(), order);
   }

   elmat = 0.0;

   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);

      el.CalcDShape(ip, dshape);

      Trans.SetIntPoint(&ip);
      w = ip.weight * Trans.Weight();
      Mult(dshape, Trans.InverseJacobian(), gshape);
      MultAAt(gshape, pelmat);
      gshape.GradToDiv (divshape);

      M = mu->Eval(Trans, ip);
      if (lambda)
      {
         L = lambda->Eval(Trans, ip);
      }
      else
      {
         L = q_lambda * M;
         M = q_mu * M;
      }

      if (L != 0.0)
      {
         AddMult_a_VVt(L * w, divshape, elmat);
      }

      if (M != 0.0)
      {
         for (int d = 0; d < dim; d++)
         {
            for (int k = 0; k < dof; k++)
               for (int l = 0; l < dof; l++)
               {
                  elmat (dof*d+k, dof*d+l) += (M * w) * pelmat(k, l);
               }
         }
         for (int ii = 0; ii < dim; ii++)
            for (int jj = 0; jj < dim; jj++)
            {
               for (int kk = 0; kk < dof; kk++)
                  for (int ll = 0; ll < dof; ll++)
                  {
                     elmat(dof*ii+kk, dof*jj+ll) +=
                        (M * w) * gshape(kk, jj) * gshape(ll, ii);
                  }
            }
      }
   }
}

