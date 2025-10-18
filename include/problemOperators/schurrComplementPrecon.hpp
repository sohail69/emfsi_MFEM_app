#pragma once

#include "mfem.hpp"
#include <memory>
#include <iostream>
#include <fstream>


using namespace std;
using namespace mfem;


/*****************************************\
!
! Given a set of operator blocks solves
! a schurr complement type system, to be
! used as a preconditioner. For saddle point
! systems of the form:
!
!   [A  B][x1] = [b1]
!   [C  0][x2] = [b2]
!
! The preconditioner is thus:
!   invP = [I  invA.B][invA   0 ]
!          [0     I  ][ 0   invS]
!
! Where invS is equal to:
!   invS = C.invA.B
!
/*****************************************/
class schurrComplementPrecon : public Solver
{
  private:
    Array<int> btoffs;
    mutable Vector xtmp, rtmp;
    mutable BlockVector xBlock, rhsBlock; //auxillary vector(s)
    mutable Operator *OpA=NULL, *OpB=NULL, *OpC=NULL;


    mutable HypreParVector *Ad = NULL;
    mutable HypreParMatrix *matA=NULL, *matB=NULL, *matC=NULL, *matS=NULL;


 //   mutable IterativeSolver *invS=NULL;
    mutable Solver *invS=NULL, *invA2=NULL;
    mutable HypreBoomerAMG *invA=NULL;
    mutable BlockDiagonalPreconditioner *diagPrecon=NULL;
  public:
    //Construct the preconditioner
    schurrComplementPrecon(Operator *OpA_, Operator *OpB_, Operator *OpC_);

   //Rebuild preconditioner
   void preconRebuild() const;

   //Preconditioner Mult
   virtual void Mult(const Vector &x, Vector &y) const;

   //Set the Operator (Not used)
   virtual void SetOperator(const Operator &op){};
};


/*****************************************\
!
!     Construct the preconditioner
!
/*****************************************/
schurrComplementPrecon::schurrComplementPrecon(Operator *OpA_, Operator *OpB_, Operator *OpC_):
                     Solver(OpA_->Height()+OpB_->Height(), OpA_->Width()+OpC_->Width()),
                     btoffs(3), OpA(OpA_), OpB(OpB_), OpC(OpC_)
{
   //Set the Offsets and auxillary block Vectors
   btoffs[0] = 0;
   btoffs[1] = OpA_->Height();
   btoffs[2] = OpB_->Height();
   btoffs.PartialSum();
   rhsBlock.Update(btoffs);
   xBlock.Update(btoffs);
   xtmp.SetSize(OpB_->Height());
   rtmp.SetSize(OpA_->Height());

   if(diagPrecon != NULL ) delete diagPrecon;
   diagPrecon = NULL;
   diagPrecon = new BlockDiagonalPreconditioner(btoffs);
};

/*****************************************\
!
!         Rebuild preconditioner
!
/*****************************************/
void schurrComplementPrecon::preconRebuild() const
{
  if(Ad    != NULL) delete Ad;
  if(invA  != NULL) delete invA;
  if(invS  != NULL) delete invS;

  //Construct the a Schurr Complement
  //Gauss-Seidel block Preconditioner
  matA = static_cast<HypreParMatrix*>( OpA );
  matB = static_cast<HypreParMatrix*>( OpB );
  matC = static_cast<HypreParMatrix*>( OpC );

  invA2 = new HypreDiagScale(*matA);
  invA  = new HypreBoomerAMG(*matA);
  invA->SetInterpolation(9);
  invA->SetRelaxType(12);
  invA->SetCycleType(2);
  invA->SetCoarsening(6);

  Ad = new HypreParVector(MPI_COMM_WORLD, matA->GetGlobalNumRows(),matA->GetRowStarts());
  matA->GetDiag(*Ad);
  matC->InvScaleRows(*Ad);
  matS = ParMult(matB, matC);
  invS = new SLISolver(MPI_COMM_WORLD);
  invS->SetOperator(*matS);

  diagPrecon->SetDiagonalBlock(0,invA);
  diagPrecon->SetDiagonalBlock(1,invS);
};

/*****************************************\
!
!        Preconditioner Mult
!
/*****************************************/
void schurrComplementPrecon::Mult(const Vector &x, Vector &y) const
{
  copyVec(x, xBlock);
  diagPrecon->Mult(x,rhsBlock);
  OpB->Mult(rhsBlock.GetBlock(1),xtmp);
  invA->Mult(xtmp,rtmp);
  rhsBlock.GetBlock(0) -= rtmp;
//  rhsBlock.GetBlock(1) = 0.00;
  real_t eps = 1.1;
  real_t dotProd = InnerProduct(Ad,Ad);
  rhsBlock.GetBlock(1) = xBlock.GetBlock(1);
  rhsBlock.GetBlock(1) *= (1.0/(dotProd+eps));
  copyVec(rhsBlock, y);
};
