#pragma once
#include "mfem.hpp"

using namespace std;
using namespace mfem;

/*****************************************\
!
! This adds an abstract interface for the
! single Cell models so that they can
! communicate to the rest of MFEM
!
! Effectively an ODE-Black-Box model used
! used for the Dark-Art that is single cell 
! modelling
!
/*****************************************/
class mulitCellOperator : public Operator
{
protected:
   //Has its own private FE-space
   // that never changes
   ParFiniteElementSpace *_fespace;

   //Has grid functions that auto-update
   mutable ParGridFunction *v, *gama, *phi;
public:
   //Constructs multicell Operator
   mulitCellOperator(ParFiniteElementSpace &f);

   //Solves the Poisson equation
   virtual void Mult(const Vector &x, Vector &y) const;

   //Destroys the fibre map operator
   virtual ~mulitCellOperator();
};

