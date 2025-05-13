#include "mfem.hpp"

#include <memory>
#include <iostream>
#include <fstream>
#include "include/dirchBCs.hpp"

using namespace std;
using namespace mfem;

void ParaViewVisualise(ParGridFunction* Field
                     , string           FieldName
					 , string           ProbName
					 , int order
					 , ParMesh *pmesh
					 , double time){
  
  ParaViewDataCollection paraview_dc(ProbName, pmesh);
  paraview_dc.SetPrefixPath("ParaView");
  paraview_dc.SetLevelsOfDetail(order);
  paraview_dc.SetDataFormat(VTKFormat::BINARY);
  paraview_dc.SetHighOrderOutput(true);
  paraview_dc.SetCycle(0);
  paraview_dc.SetTime(time);
  paraview_dc.RegisterField(FieldName,Field);
  paraview_dc.Save();
};



double f(double u) { return u; }
double df(double u) { return 1.0; }

// Define a coefficient that, given a grid function u, function func, returns func(u)
class NonlinearGridFunctionCoefficient : public Coefficient
{
   GridFunction &gf;
   std::function<double(double)> func;
public:
   NonlinearGridFunctionCoefficient(GridFunction &gf_, std::function<double(double)> func_) : gf(gf_), func(func_) { }
   double Eval(ElementTransformation &T, const IntegrationPoint &ip)
   {
      return func(gf.GetValue(T, ip));
   }
};


// Define a nonlinear integrator that computes (f(u), v) and its linearized
// operator, (u df(u), v).
//
// Note that the action (f(u), v) can be computed using DomainLFIntegrator
// and the Jacobian matrix linearized operator can be computed using
// MassIntegrator with the appropriate coefficients.
class NonlinearMassIntegrator : public NonlinearFormIntegrator
{
   FiniteElementSpace &fes;
   GridFunction gf;
   Array<int> dofs;

public:
   NonlinearMassIntegrator(FiniteElementSpace &fes_) : fes(fes_), gf(&fes) { }

   virtual void AssembleElementVector(const FiniteElement &el,
                                      ElementTransformation &Tr,
                                      const Vector &elfun, Vector &elvect)
   {
      fes.GetElementDofs(Tr.ElementNo, dofs);
      gf.SetSubVector(dofs, elfun);
      NonlinearGridFunctionCoefficient coeff(gf, f);
      DomainLFIntegrator integ(coeff);
      integ.AssembleRHSElementVect(el, Tr, elvect);
   }

   virtual void AssembleElementGrad(const FiniteElement &el,
                                    ElementTransformation &Tr,
                                    const Vector &elfun, DenseMatrix &elmat)
   {
      fes.GetElementDofs(Tr.ElementNo, dofs);
      gf.SetSubVector(dofs, elfun);
      NonlinearGridFunctionCoefficient coeff(gf, df);
      MassIntegrator integ(coeff);
      integ.AssembleElementMatrix(el, Tr, elmat);
   }
};

class MY_NLOperator : public Operator{
  private:
    //The FE-spaces and operators
    int nBlocks;
    ParFiniteElementSpace *feSpace=NULL;//FE-Spaces
    Operator              *NLForm=NULL; //Operators used for E(U)
    BilinearForm          *a=NULL;
    OperatorPtr A;

    //Forms used for Forcing vector (b)
    ParLinearForm *LForm=NULL; //Linear forms

    //Array of constrainedNodes
    Array<int> bdr_tDofs;

    //vectors
    mutable Vector b_vec, x_Dirch; //BC vectors
	mutable Vector z_vec, w_vec;   //Temporary vectors

  public:
    //Creates my block non-linear form
    MY_NLOperator(ParFiniteElementSpace *feSpace_
	            , ParNonlinearForm *NLForm_
                , ParLinearForm *LForms_
                , const MemoryType deviceMT);

    //Destroys my block linear form
    ~MY_NLOperator();

    //Sets reference configuration for a particular block
    //within the Block vector (set to zero by default)
    //used to set initial conditions as well and for
    //DirchBC's
    void setRefDirchVec(Vector xref);

    //Sets the DOF's to be eliminated for a particular block
    void setDirchBCs(Array<int> bdr_tDofs_);

    //The residual operator
    virtual void Mult(const Vector &x, Vector &y) const override;

    //The Jacobian operator (Can be used for preconditioning)
    Operator &GetGradient(const Vector &x) const override;

    //A null and void function that inherits from
    //the Operator class
    virtual void SetOperator(const Operator &op){};
};


MY_NLOperator::MY_NLOperator(ParFiniteElementSpace *feSpace_
	            , ParNonlinearForm *NLForm_
                , ParLinearForm *LForms_
                , const MemoryType deviceMT):
				Operator(feSpace_->TrueVSize())
                ,LForm(LForms_)
				,NLForm(NLForm_)
{
   feSpace = new ParFiniteElementSpace(*feSpace_);
   b_vec   = Vector(feSpace_->TrueVSize(),deviceMT); b_vec   = 0.0;
   x_Dirch = Vector(feSpace_->TrueVSize(),deviceMT); x_Dirch = 0.0;
   z_vec   = Vector(feSpace_->TrueVSize(),deviceMT); z_vec   = 0.0;
   w_vec   = Vector(feSpace_->TrueVSize(),deviceMT); w_vec   = 0.0;
   LForm->ParallelAssemble(b_vec);

   // Solve as a linear problem.
   a = new BilinearForm(feSpace);
   a->AddDomainIntegrator(new DiffusionIntegrator);
   a->AddDomainIntegrator(new MassIntegrator);
   a->Assemble();
}

MY_NLOperator::~MY_NLOperator(){}; 

void MY_NLOperator::setRefDirchVec(Vector xref){setValues(xref,x_Dirch);};

void MY_NLOperator::setDirchBCs(Array<int> bdr_tDofs_){bdr_tDofs = Array<int>(bdr_tDofs_);};

void MY_NLOperator::Mult(const Vector &x, Vector &y) const {
  //Apply the Dirchelet values to the solution vector
  setValues(x,w_vec);
  applyDirchValues(x_Dirch, w_vec, bdr_tDofs);

  //Calculate the unconstrained residual
  NLForm->Mult(w_vec,z_vec);
  y = 0.0;
  y.Add(-1.0,b_vec);
  y.Add(1.0,z_vec);

  //Apply the Dirchelet elimination on the residual
  eliminateDirchValues(y, bdr_tDofs);
};

Operator &MY_NLOperator::GetGradient(const Vector &x) const{
  return (NLForm->GetGradient(x));
};

int main(int argc, char *argv[])
{
   // 1. Initialize MPI and HYPRE.
   Mpi::Init();
   const int myid = Mpi::WorldRank();

   // 1. Parse command-line options
   const char *mesh_file = "mesh/star.mesh";
   int ref_levels = -1;
   int order = 2;
   bool visualization = true;
   double newton_rel_tol = 1e-4;
   double newton_abs_tol = 1e-6;
   int newton_iter = 500;
   double mu = 1.0;
   const char *device_config = "cpu";

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly.");
   args.AddOption(&order, "-o", "--order",
                  "Order (degree) of the finite elements.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&newton_rel_tol, "-rel", "--relative-tolerance",
                  "Relative tolerance for the Newton solve.");
   args.AddOption(&newton_abs_tol, "-abs", "--absolute-tolerance",
                  "Absolute tolerance for the Newton solve.");
   args.AddOption(&newton_iter, "-it", "--newton-iterations",
                  "Maximum iterations for the Newton solve.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");

   args.Parse();

   // 1. Parse command line options
   if (!args.Good()) args.PrintUsage(cout);
   if (!args.Good()) return 1;
   args.PrintOptions(cout);
   Device device(device_config);
   if (myid == 0) { device.Print(); }
   MemoryType mt = device.GetMemoryType();

   bool nonzero_rhs = true;

   // 2. Read the mesh from the given mesh file, and refine once uniformly.
   Mesh mesh(mesh_file);
   int dim = mesh.Dimension();
   if (ref_levels == -1) ref_levels = (int)floor(log(10000./mesh.GetNE())/log(2.)/dim);
   for (int l = 0; l < ref_levels; l++) mesh.UniformRefinement();

   ParMesh pmesh(MPI_COMM_WORLD, mesh);

   // 3. Define a finite element space on the mesh. Here we use H1 continuous
   //    high-order Lagrange finite elements of the given order.
   H1_FECollection fec(order, dim);
   ParFiniteElementSpace fespace(&pmesh, &fec);

   // 4. Extract the list of all the boundaries. These will be marked as
   //    Dirichlet in order to enforce zero boundary conditions.
   Array<int> ess_bdr(pmesh.bdr_attributes.Max());
   ess_bdr = 0;
   ess_bdr[0] = 1;
   Array<int> ess_tdof_list;
   fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

   // 5. Define the solution x as a finite element grid function in fespace. Set
   //    the initial guess to zero, which also sets the boundary conditions.
   ParGridFunction u1(&fespace), u2(&fespace);
   u1 = 0.0;
   u2 = 0.0;

   // Solve as a non-linear problem.

   // 6. Set up the nonlinear form n(u,v) = (grad u, grad v) + (f(u), v)
   ParNonlinearForm n(&fespace);
   n.AddDomainIntegrator(new NonlinearMassIntegrator(fespace));
   n.AddDomainIntegrator(new DiffusionIntegrator);

   // 7. Set up the the right-hand side. For simplicitly, we just use a zero
   //    vector. Because of the form of the nonlinear function f, it is still
   //    nontrivial to solve n(u,v) = 0.
   ParLinearForm b(&fespace);
   b = 0.0;
if (nonzero_rhs) {
      ConstantCoefficient five(5.0);
      b.AddDomainIntegrator(new DomainLFIntegrator(five));
      b.Assemble();
}


   Vector x_Dirch(fespace.TrueVSize(),mt);
   x_Dirch = 5.0;


   // 8. Get true dof vectors and set essential BCs on rhs.
   Vector X(fespace.GetTrueVSize()), B(fespace.GetTrueVSize());
   u1.GetTrueDofs(X);
   n.SetEssentialBC(ess_bdr, &B);

   MY_NLOperator myOp(&fespace, &n, &b, mt);
   myOp.setDirchBCs(ess_tdof_list);
   myOp.setRefDirchVec(x_Dirch);



   // 9. Set up the Newton solver. Each Newton iteration requires a linear
   //    solve. Here we use UMFPack as a direct solver for these systems.
   CGSolver solver(MPI_COMM_WORLD);
   NewtonSolver newton(MPI_COMM_WORLD);
   newton.SetOperator(myOp);
   newton.SetSolver(solver);
   newton.SetPrintLevel(1);
   newton.SetRelTol(1e-10);
   newton.SetMaxIter(100);

   // 10. Solve the nonlinear system.
   applyDirchValues(x_Dirch, X, ess_tdof_list); //Apply before solve
   newton.Mult(B, X);
   u1.Distribute(X);

   {
     // Solve as a linear problem.
     BilinearForm a(&fespace);
     a.AddDomainIntegrator(new DiffusionIntegrator);
     a.AddDomainIntegrator(new MassIntegrator);
     a.Assemble();


//================================
/*  setValues(x,w_vec);
  applyDirchValues(x_Dirch, w_vec, bdr_tDofs);

  //Calculate the unconstrained residual
  NLForm->Mult(w_vec,z_vec);
  y = 0.0;
  y.Add(-1.0,b_vec);
  y.Add(1.0,z_vec);

  //Apply the Dirchelet elimination on the residual
  eliminateDirchValues(y, bdr_tDofs);*/
//================================

     OperatorPtr A;
     Vector C, Y;
     a.FormLinearSystem(ess_tdof_list, u2, b, A, Y, C);
     GSSmoother M((SparseMatrix&)(*A));

     applyDirchValues(x_Dirch, Y, ess_tdof_list); //Apply before solve
     applyDirchValues(x_Dirch, b, ess_tdof_list); //Apply before solve


     PCG(*A, M, C, Y, 1, 500, 1e-12, 0.0);
     a.RecoverFEMSolution(Y, b, u2);
   }
   
   ParaViewVisualise(&u1 , "Diff", "NLPDiff", order, &pmesh, 0.0);

   u1 -= u2;

   cout << u1.Norml2() << "\n";
   return 0;
}
