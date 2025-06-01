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

  elvect.SetSize(dof*spaceDim);
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

    for(int K=0; K<dof; K++){
      for(int J=0; J<dim; J++){
        for(int I=0; I<dim; I++){
          L = I*dof + K;
          elvect(L) += rho_ij(I,J) * dshape(K,J) * ip.weight * Tr.Weight();
        }
      }
    }
  }
}


/*****************************************\
! Author: Sohail Rathore
! Date  : 2025-05-18
!
! Integrates a matrix contaction against
! the gradient of the shape function
! v u.grad(v)
!
\*****************************************/
/** An abstract class for integrating the product of a scalar basis function and
    the inner product of a vector basis function with a vector coefficient. In
    2D the inner product can be replaced with a cross product. */
class MixedScalarVectorIntegrator: public BilinearFormIntegrator
{
public:

   void AssembleElementMatrix2(const FiniteElement &trial_fe,
                               const FiniteElement &test_fe,
                               ElementTransformation &Trans,
                               DenseMatrix &elmat) override;

   /// Support for use in BilinearForm. Can be used only when appropriate.
   /** Appropriate use cases are classes derived from
       MixedScalarVectorIntegrator where the trial and test spaces can be the
       same. Examples of such classes are: MixedVectorDivergenceIntegrator,
       MixedScalarWeakDivergenceIntegrator, etc. */
   void AssembleElementMatrix(const FiniteElement &fe,
                              ElementTransformation &Trans,
                              DenseMatrix &elmat) override
   { AssembleElementMatrix2(fe, fe, Trans, elmat); }

protected:

   MixedScalarVectorIntegrator(VectorCoefficient &vq, bool transpose_ = false,
                               bool cross_2d_ = false)
      : VQ(&vq), transpose(transpose_), cross_2d(cross_2d_) {}

   inline virtual bool VerifyFiniteElementTypes(
      4const FiniteElement & trial_fe,
      const FiniteElement & test_fe) const
   {
      return ((transpose &&
               trial_fe.GetRangeType() == mfem::FiniteElement::VECTOR &&
               test_fe.GetRangeType()  == mfem::FiniteElement::SCALAR ) ||
              (!transpose &&
               trial_fe.GetRangeType() == mfem::FiniteElement::SCALAR &&
               test_fe.GetRangeType()  == mfem::FiniteElement::VECTOR )
             );
   }

   inline virtual const char * FiniteElementTypeFailureMessage() const
   {
      if ( transpose )
      {
         return "MixedScalarVectorIntegrator:  "
                "Trial space must be a vector field "
                "and the test space must be a scalar field";
      }
      else
      {
         return "MixedScalarVectorIntegrator:  "
                "Trial space must be a scalar field "
                "and the test space must be a vector field";
      }
   }

   inline virtual int GetIntegrationOrder(const FiniteElement & trial_fe,
                                          const FiniteElement & test_fe,
                                          ElementTransformation &Trans)
   { return trial_fe.GetOrder() + test_fe.GetOrder() + Trans.OrderW(); }


   inline virtual int GetVDim(const FiniteElement & vector_fe)
   { return std::max(space_dim, vector_fe.GetRangeDim()); }

   inline virtual void CalcVShape(const FiniteElement & vector_fe,
                                  ElementTransformation &Trans,
                                  DenseMatrix & shape_)
   { vector_fe.CalcVShape(Trans, shape_); }

   inline virtual void CalcShape(const FiniteElement & scalar_fe,
                                 ElementTransformation &Trans,
                                 Vector & shape_)
   { scalar_fe.CalcPhysShape(Trans, shape_); }

   VectorCoefficient *VQ;
   int space_dim;
   bool transpose;
   bool cross_2d;  // In 2D use a cross product rather than a dot product

private:

#ifndef MFEM_THREAD_SAFE
   Vector V;
   DenseMatrix vshape;
   Vector      shape;
   Vector      vshape_tmp;
#endif
};
