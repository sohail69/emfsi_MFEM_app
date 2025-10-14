#pragma once
#include "mfem.hpp"
#include <functional>
#include "../../../MFEM_STUFF/mfem-4.7/build/include/mfem/general/forall.hpp"

using namespace std;
using namespace mfem;

//
// Kroncecker delta function
//
real_t kDelta(int I, int J){return ( (I==J)? 1.00:0.00);};


// Given an ElementTransformation and IntegrationPoint in a refined mesh,
// return the ElementTransformation of the parent coarse element, and set
// coarse_ip to the location of the original ip within the coarse element.
ElementTransformation *RefinedToCoarse(
   Mesh &coarse_mesh, const ElementTransformation &T,
   const IntegrationPoint &ip, IntegrationPoint &coarse_ip)
{
   const Mesh &fine_mesh = *T.mesh;
   // Get the element transformation of the coarse element containing the
   // fine element.
   int fine_element = T.ElementNo;
   const CoarseFineTransformations &cf = fine_mesh.GetRefinementTransforms();
   int coarse_element = cf.embeddings[fine_element].parent;
   ElementTransformation *coarse_T = coarse_mesh.GetElementTransformation(
                                        coarse_element);
   // Transform the integration point from fine element coordinates to coarse
   // element coordinates.
   Geometry::Type geom = T.GetGeometryType();
   IntegrationPointTransformation fine_to_coarse;
   IsoparametricTransformation &emb_tr = fine_to_coarse.Transf;
   emb_tr.SetIdentityTransformation(geom);
   emb_tr.SetPointMat(cf.point_matrices[geom](cf.embeddings[fine_element].matrix));
   fine_to_coarse.Transform(ip, coarse_ip);
   coarse_T->SetIntPoint(&coarse_ip);
   return coarse_T;
}


//
// Copy one vector to another
//
void copyVec(const Vector & x, Vector & y){y = x;};

//
// Apply Dirchelet Boundary conditions
//
void applyDirchValues(const mfem::Vector &k, mfem::Vector &y, mfem::Array<int> dofs)
{
  if(dofs.Size() > 0){ //Only apply if there are constrained DOF's
    const bool use_dev = dofs.UseDevice() || k.UseDevice() || y.UseDevice();
    const int n = dofs.Size();
    // Use read+write access for X - we only modify some of its entries
    auto d_X = y.ReadWrite(use_dev);
    auto d_y = k.Read(use_dev);
    auto d_dofs = dofs.Read(use_dev);
    mfem::forall_switch(use_dev, n, [=] MFEM_HOST_DEVICE (int i)
    {
      const int dof_i = d_dofs[i];
      if (dof_i >= 0)   d_X[dof_i]    =  d_y[dof_i];
      if (!(dof_i >= 0))d_X[-1-dof_i] = -d_y[-1-dof_i];
    });
  }
};

//
// Apply Dirchelet Boundary conditions
//
template<typename TCoeff>
void ApplyDircheletBCs(const Array<int>     ess_BCTags
                     , const Array<TCoeff*> ess_BCs
                     , ParGridFunction      *X)
{
  if( ess_BCTags.Size() != 0){
    for(int I=0; I<ess_BCTags.Size(); I++){
      if(ess_BCTags[I] == 1){
        Array<int> tmp_BDR_tags( ess_BCTags.Size() );
        tmp_BDR_tags = 0;
        tmp_BDR_tags[I] = 1;
        X->ProjectBdrCoefficient(*(ess_BCs[I]), tmp_BDR_tags);
      }
    }
  }
};

