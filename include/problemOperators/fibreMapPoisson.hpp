#pragma once
#include "mfem.hpp"
#include "../coefficients/fibreFieldCoeff.hpp"

using namespace std;
using namespace mfem;

int stimulus_mask(double phi){
  return ((phi <= 0.25) ? 1 : 0);
}


int CellType_mask(double phi){
  return ((phi <= 0.6) ? 5 : ((phi <= 0.9) ? 4 : 3));
}

/*****************************************\
!
! Solves the poisson equation and generates
! the necesssary fibre maps needed for other
! calculations (Active strain, orthotropic
! diffusion, etc...)
!
/*****************************************/
class fibreMapPoissonOper : public Operator
{
protected:
   Array<int>    & _ess_bcs_markers;
   Array<double> & _BC_Vals;
   ParFiniteElementSpace & _fespace;
   ParMixedBilinearForm *_K;
   ParBilinearForm *_K_uncon;
   HypreParMatrix _Kmat, _Kmat_uncon;

   CGSolver _K_solver;    // Krylov solver for inverting the mass matrix K
   HypreSmoother _K_prec; // Preconditioner for the diffusion matrix K

   //The Preconditioning objects
   HypreBoomerAMG *invM=NULL;
   mutable ParGridFunction *z; // auxiliary GFunc used for BC's

   //Variables for the actual fibre fields 
   //including the coefficients
   fibreFieldCoeff1 *cf1=NULL;
   fibreFieldCoeff2 *cf2=NULL;
   fibreFieldCoeff3 *cf3=NULL;
   fibreFieldCoeff4 *cf4=NULL;
public:
   //Constructs fibre map operator
   fibreMapPoissonOper(ParFiniteElementSpace &f, Array<int> & ess_bcs_markers, Array<double> & BC_Vals);

   //Solves the Poisson equation
   virtual void Mult(const Vector &x, Vector &y) const;

   //Destroys the fibre map operator
   virtual ~fibreMapPoissonOper();

   //Gets a grid-function array of the vector basis
   //needed to define the fibre fields
   void getFibreMapGFuncs3D(const Vector &x, Array<ParGridFunction*> & FibreFields);
   void getFibreMapGFuncs2D(const Vector &x, Array<ParGridFunction*> & FibreFields);
};

/*****************************************\
!
! Constructs the fibre map operator
! Assumes a H1 FE-space
!
/*****************************************/
fibreMapPoissonOper::fibreMapPoissonOper(ParFiniteElementSpace &f
                                       , Array<int> & ess_bcs_markers
                                       , Array<double> & BC_Vals): 
                                       Operator(f.GetTrueVSize()), _ess_bcs_markers(ess_bcs_markers)
                                       , _BC_Vals(BC_Vals), _fespace(f)
{
   MFEM_VERIFY(BC_Vals.Size() == ess_bcs_markers.Size(), "The BC arrays should be same size fibrePoisson");

  //Construct the matrix Operators
  Array<int> ess_bcs_tdofs, empty_tdofs;
  if(_ess_bcs_markers.Size() != 0) _fespace.GetEssentialTrueDofs(_ess_bcs_markers, ess_bcs_tdofs);
  ConstantCoefficient one(1.0);

  //Constrained
  _K = new ParMixedBilinearForm(&_fespace,&_fespace);
  _K->AddDomainIntegrator(new DiffusionIntegrator(one));
  _K->Assemble();
  _K->FormRectangularSystemMatrix(empty_tdofs, ess_bcs_tdofs, _Kmat);

  //Unconstrained
  _K_uncon = new ParBilinearForm(&_fespace);
  _K_uncon->AddDomainIntegrator(new DiffusionIntegrator(one));
  _K_uncon->Assemble();
  _K_uncon->FormSystemMatrix(empty_tdofs, _Kmat_uncon);

  //Build the preconditioner
  invM = new HypreBoomerAMG(_Kmat_uncon);
  invM->SetInterpolation(6);
  invM->SetCoarsening(8);
  invM->SetRelaxType(6);
  invM->SetCycleNumSweeps(2,2);
  invM->SetCycleType(2);

  //Construct the solver
  const real_t rel_tol = 1e-8;
  _K_solver.iterative_mode = false;
  _K_solver.SetRelTol(rel_tol);
  _K_solver.SetAbsTol(0.0);
  _K_solver.SetMaxIter(500);
  _K_solver.SetPrintLevel(0);
  _K_solver.SetPreconditioner(*invM);
  _K_solver.SetOperator(_Kmat);

  //Construct the temporary gridFunction
  z = new ParGridFunction(&_fespace);
};

/*****************************************\
!
! Solves the Poisson problem
!
/*****************************************/
void fibreMapPoissonOper::Mult(const Vector &x, Vector &y) const
{
  //Apply the BC's
  (*z)=0.00;
  if(_BC_Vals.Size() != 0){
    Array<int> ess_bcs_tmp(_BC_Vals.Size()), ess_tdofs;
    for(int I=0; I<_BC_Vals.Size(); I++){
      if(_ess_bcs_markers[I] != 0){
        ess_bcs_tmp = 0;
        ess_bcs_tmp[I] = 1;
        _fespace.GetEssentialTrueDofs(ess_bcs_tmp, ess_tdofs);
        z->SetSubVector(ess_tdofs,-_BC_Vals[I]);
      }
    }
  }
  y = *z;
  _Kmat.Mult( y, *z);

  //Solve the problem
  _K_solver.Mult(*z,y);

  //Repair the boundary issue (Hack)
  if(_BC_Vals.Size() != 0){
    Array<int> ess_bcs_tmp(_BC_Vals.Size()), ess_tdofs;
    for(int I=0; I<_BC_Vals.Size(); I++){
      if(_ess_bcs_markers[I] != 0){
        ess_bcs_tmp = 0;
        ess_bcs_tmp[I] = 1;
        _fespace.GetEssentialTrueDofs(ess_bcs_tmp, ess_tdofs);
        y.SetSubVector(ess_tdofs,_BC_Vals[I]);
      }
    }
  }
};

/*****************************************\
!
! Destroys the fibre map operator
!
\*****************************************/
fibreMapPoissonOper::~fibreMapPoissonOper(){
  delete _K, _K_uncon, invM;
  delete z;
};


/*****************************************\
!
! Constructs the actual fibre maps and puts
! them into a gridfunction 3D
!
\*****************************************/
void fibreMapPoissonOper::getFibreMapGFuncs3D(const Vector &x, Array<ParGridFunction*> & FibreFields){
   const int dim = 3;
   *z = x;
   real_t PI=3.14159265;
   Vector C0(dim); 
   real_t theta_epi=(60.0*PI/180.0), theta_end=(-60.0*PI/180.0);
   C0(0)=0.0;     C0(1)=1.0;   C0(2)=0.0;

   if(cf1 == NULL) cf1 = new fibreFieldCoeff1(dim, z);
   FibreFields[0]->ProjectCoefficient(*cf1);

   if(cf2 == NULL) cf2 = new fibreFieldCoeff2(dim, *z, *FibreFields[0], C0, theta_epi, theta_end);
   FibreFields[1]->ProjectCoefficient(*cf2);

   if(cf3 == NULL) cf3 = new fibreFieldCoeff3(dim, *FibreFields[1], *FibreFields[0]);
   FibreFields[2]->ProjectCoefficient(*cf3);

   delete cf1, cf2, cf3;
   cf1=NULL; cf2=NULL; cf3=NULL;
};

/*****************************************\
!
! Constructs the actual fibre maps and puts
! them into a gridfunction 2D
!
\*****************************************/
void fibreMapPoissonOper::getFibreMapGFuncs2D(const Vector &x, Array<ParGridFunction*> & FibreFields){
   const int dim = 2;
   *z = x;

   if(cf1 == NULL) cf1 = new fibreFieldCoeff1(dim, z);
   FibreFields[0]->ProjectCoefficient(*cf1);

   if(cf4 == NULL) cf4 = new fibreFieldCoeff4(dim, FibreFields[0]);
   FibreFields[1]->ProjectCoefficient(*cf4);

   delete cf1, cf4;
   cf1=NULL; cf4=NULL;
};
