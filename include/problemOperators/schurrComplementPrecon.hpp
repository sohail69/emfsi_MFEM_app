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
    mutable BlockVector rhsBlock; //auxillary vector(s)
    mutable Operator *OpA=NULL, *OpB=NULL, *OpC=NULL;


    mutable HypreParVector *Ad = NULL;
    mutable HypreParMatrix *matA=NULL, *matB=NULL, *matC=NULL;

    mutable Solver *invS=NULL;
    mutable HypreBoomerAMG *invA=NULL;
    mutable BlockDiagonalPreconditioner *diagPrecon=NULL;
  public:
    //Construct the preconditioner
    schurrComplementPrecon(Operator *OpA_, Operator *OpB_, Operator *OpC_);

   //Preconditioner Mult
   virtual void Mult(const Vector &x, Vector &y) const;

   //Rebuild preconditioner
   void preconRebuild() const;
};

/*****************************************\
!
!     Construct the preconditioner
!
/*****************************************/
schurrComplementPrecon::schurrComplementPrecon(Operator *OpA_, Operator *OpB_, Operator *OpC_):
     btoffs(3), OpA(OpA_), OpB(OpB_), OpC(OpC_)
     , Solver(OpA_->Height()+OpC_->Height(), OpA_->Width()+OpB_->Width())
{
   //Set the Offsets and auxillary block Vectors
   btoffs[0] = 0;
   btoffs[1] = OpA_->Width();
   btoffs[2] = OpB_->Width();
   btoffs.PartialSum();
   rhsBlock.Update(btoffs);
   xtmp.SetSize(OpB_->Height());
   rtmp.SetSize(OpA_->Height());

   if(diagPrecon != NULL ) delete diagPrecon;
   diagPrecon = NULL;
   diagPrecon = new BlockDiagonalPreconditioner(btoffs);
};

/*****************************************\
!
!        //Rebuild preconditioner
!
/*****************************************/
void schurrComplementPrecon::preconRebuild() const
{
/*

    Array<int> btoffs;
    mutable Vector xtmp, rtmp;
    mutable BlockVector rhsBlock; //auxillary vector(s)
    mutable Operator *OpA=NULL, *OpB=NULL, *OpC=NULL;


    mutable HypreParVector *Ad = NULL;
    mutable HypreParMatrix *matA=NULL, *matB=NULL, *matC=NULL;

    mutable Solver *invS=NULL;
    mutable HypreBoomerAMG *invA=NULL;
    mutable BlockDiagonalPreconditioner *diagPrecon=NULL;

  if(Md    != NULL) delete Md;
  if(invM  != NULL) delete invM;
  if(invS  != NULL) delete invS;

  //Construct the a Schurr Complement
  //Gauss-Seidel block Preconditioner
  matM = static_cast<HypreParMatrix*>( opUU.Ptr() );
  matB = static_cast<HypreParMatrix*>( opUP.Ptr() );
  matC = static_cast<HypreParMatrix*>( opPU.Ptr() );

  Md = new HypreParVector(MPI_COMM_WORLD, matM->GetGlobalNumRows(),matM->GetRowStarts());
  matM->GetDiag(*Md);
//  matC->InvScaleRows(*Md);
//  invM = new HypreDiagScale(*matM);

  matS = ParMult(matB, matC);
  invS = new HypreBoomerAMG(*matM);
//  invS = new HypreBoomerAMG(*matS);
  invS->SetInterpolation(9);
  invS->SetRelaxType(12);
  invS->SetCycleType(2);
  invS->SetCoarsening(6);

  nssPrecon->SetDiagonalBlock(0,invS);
//  nssPrecon->SetDiagonalBlock(0,invM);
//  nssPrecon->SetDiagonalBlock(1,invS);*/
};



/*****************************************\
!
!   //Preconditioner Mult
!
/*****************************************/
void schurrComplementPrecon::Mult(const Vector &x, Vector &y) const
{
  diagPrecon->Mult(x,rhsBlock);
  OpB->Mult(rhsBlock.GetBlock(0),xtmp);
  invA->Mult(xtmp,rtmp);
  rhsBlock.GetBlock(0) -= rtmp;
  copyVec(rhsBlock, y);
};
