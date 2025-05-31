#pragma once
#include "mfem.hpp"
#include <functional>

using namespace std;
using namespace mfem;

//
// Kroncecker delta function
//
real_t kDelta(int I, int J){return ( (I==J)? 1.00:0.00);};

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
