#pragma once
#include "mfem.hpp"
#include <functional>

using namespace std;
using namespace mfem;

/*****************************************\
! Author: Sohail Rathore
! Date  : 2025-05-18
!
! Generate a diffusion contribution
!
\*****************************************/
class vectorDiffusionIntegrator : public DeltaLFIntegrator
{
private:
   Vector shape, Qvec;
   Coefficient       *Diff_val;
   VectorCoefficient *Diff_diag;
   MatrixCoefficient *Diff_mat;
   DenseMatrix dshape, dshapedxt;
   GridFunction *vec_gf; //vector field
public:
   /// Constructs the domain integrator (Q, grad v)
   DeltaLFIntegrator(MatrixCoefficient &QF)
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
void vectorDiffusionIntegrator::AssembleRHSElementVect(const FiniteElement &el
                                               , ElementTransformation &Tr
                                               , Vector &elvect)
{
   int nd = el.GetDof();
   dim = el.GetDim();
   int spaceDim = Tr.GetSpaceDim();
   real_t w;


#ifdef MFEM_THREAD_SAFE
   DenseMatrix dshape(nd,dim), invdfdx(dim, spaceDim), M(MQ ? spaceDim : 0);
   Vector D(VQ ? VQ->GetVDim() : 0);
#else
   dshape.SetSize(nd,dim);
   invdfdx.SetSize(dim, spaceDim);
   M.SetSize(MQ ? spaceDim : 0);
   D.SetSize(VQ ? VQ->GetVDim() : 0);
#endif
   vec.SetSize(dim);
   vecdxt.SetSize((VQ || MQ) ? spaceDim : 0);
   pointflux.SetSize(spaceDim);

  elvect.SetSize(nd);
  Vector &elfun
  vec_gf


   const IntegrationRule *ir = GetIntegrationRule(el, Tr);
   elvect = 0.0;
   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      el.CalcDShape(ip, dshape);

      Tr.SetIntPoint(&ip);
      CalcAdjugate(Tr.Jacobian(), invdfdx); // invdfdx = adj(J)
      w = ip.weight / Tr.Weight();

      if (!MQ && !VQ)
      {
         dshape.MultTranspose(elfun, vec);
         invdfdx.MultTranspose(vec, pointflux);
         if (Q)
         {
            w *= Q->Eval(Tr, ip);
         }
      }
      else
      {
         dshape.MultTranspose(elfun, vec);
         invdfdx.MultTranspose(vec, vecdxt);
         if (MQ)
         {
            MQ->Eval(M, Tr, ip);
            M.Mult(vecdxt, pointflux);
         }
         else
         {
            VQ->Eval(D, Tr, ip);
            for (int j=0; j<spaceDim; ++j)
            {
               pointflux[j] = D[j] * vecdxt[j];
            }
         }
      }
      pointflux *= w;
      invdfdx.Mult(pointflux, vec);
      dshape.AddMult(vec, elvect);

   }
};
