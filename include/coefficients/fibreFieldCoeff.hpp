
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
class fibreFieldCoeff : public VectorCoefficient
{
private:
   std::function<void(const Vector &, Vector &)> Function;
   Coefficient *Q;

public:
   /// Define a time-independent vector coefficient from a std function
   /** \param dim - the size of the vector
       \param F - time-independent function
       \param q - optional scalar Coefficient to scale the vector coefficient */
   VectorFunctionCoefficient(int dim,
                             std::function<void(const Vector &, Vector &)> F,
                             Coefficient *q = nullptr)
      : VectorCoefficient(dim), Function(std::move(F)), Q(q)
   { }

   /// Define a time-dependent vector coefficient from a std function
   /** \param dim - the size of the vector
       \param TDF - time-dependent function
       \param q - optional scalar Coefficient to scale the vector coefficient */
   VectorFunctionCoefficient(int dim,
                             std::function<void(const Vector &, real_t, Vector &)> TDF,
                             Coefficient *q = nullptr)
      : VectorCoefficient(dim), TDFunction(std::move(TDF)), Q(q)
   { }

   using VectorCoefficient::Eval;
   /// Evaluate the vector coefficient at @a ip.
   void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip) override;

   virtual ~VectorFunctionCoefficient() { }
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
void fibreFieldCoeff::Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip){

  int nDIM = fibreBasis.Size();
  Vector f_i;

  for(int I=0; I<nDIM; I++){
    fibreBasis[I]->;
    for(int J=0; J<nDIM; J++){
      for(int K=0; K<nDIM; K++){

      }
    }
  }
//   Array<GridFunction*> & fibreBasis; //Fibre GFuncs
  // GridFunction         & diff;       //Diffusion GFun

   real_t x[3];
   Vector transip(x, 3);

   T.Transform(ip, transip);

   V.SetSize(vdim);
   if (Function)
   {
      Function(transip, V);
   }
   else
   {
      TDFunction(transip, GetTime(), V);
   }
   if (Q)
   {
      V *= Q->Eval(T, ip, GetTime());
   }


};

