
/*****************************************\
!
! produces a orthotropic matrix coefficient
! based on a given set of orthotropic fibre 
! vector fields and a given vector field of
! scalar coefficients
!
! The diffusion tensor is given as
! D_ij = a_p (f_i f_j)_p
!
\*****************************************/
template<unsigned nDIM>
class orthoDiffGFuncCoeff : MatrixCoefficient
{
protected:
   int height, width;
   real_t time;
   bool symmetric;  // deprecated
   Array<GridFunction*> & fibreBasis; //Fibre GFuncs
   GridFunction         & diff;       //Diffusion GFun

public:
   /// Construct a dim x dim matrix coefficient.
   explicit orthoDiffCoeff(int dim, bool symm=false)
   { height = width = dim; time = 0.; symmetric = symm; }

   /// Construct a h x w matrix coefficient.
   orthoDiffCoeff(int h, int w, bool symm=false) :
      height(h), width(w), time(0.), symmetric(symm) { }


   /// Construct the matrix coefficients
   orthoDiffCoeff(Array<GridFunction&> & GFuncs)



   /** @brief Evaluate the matrix coefficient in the element described by @a T
       at the point @a ip, storing the result in @a K. */
   /** @note When this method is called, the caller must make sure that the
       IntegrationPoint associated with @a T is the same as @a ip. This can be
       achieved by calling T.SetIntPoint(&ip). */
   virtual void Eval(DenseMatrix &K, ElementTransformation &T,
                     const IntegrationPoint &ip) = 0;

   /// @brief Fill the QuadratureFunction @a qf by evaluating the coefficient at
   /// the quadrature points. The matrix will be transposed or not according to
   /// the boolean argument @a transpose.
   ///
   /// The @a vdim of the QuadratureFunction should be equal to the height times
   /// the width of the matrix.
   virtual void Project(QuadratureFunction &qf, bool transpose=false);
};




/*****************************************\
!
! Evaluate the diffusion tensor at the 
! integration points
!
! The diffusion tensor is given as
! D_ij = a_p (f_i f_j)_p
!
\*****************************************/
void Eval(DenseMatrix &K, ElementTransformation &T,const IntegrationPoint &ip){
/*
  for(int I=0; I<fibreBasis.Size(); I++){
    Vector f_i;
    fibreBasis[I]->;
    for(int J=0; J<; J++){
      for(int K=0; K<; K++){

      }
    }
  }*/
};

