//                       MFEM Example 16 - Parallel Version
//
// Compile with: make ex16p
//
// Sample runs:  mpirun -np 4 ex16p
//               mpirun -np 4 ex16p -m ../data/inline-tri.mesh
//               mpirun -np 4 ex16p -m ../data/disc-nurbs.mesh -tf 2
//               mpirun -np 4 ex16p -s 1 -a 0.0 -k 1.0
//               mpirun -np 4 ex16p -s 2 -a 1.0 -k 0.0
//               mpirun -np 8 ex16p -s 3 -a 0.5 -k 0.5 -o 4
//               mpirun -np 4 ex16p -s 14 -dt 1.0e-4 -tf 4.0e-2 -vs 40
//               mpirun -np 16 ex16p -m ../data/fichera-q2.mesh
//               mpirun -np 16 ex16p -m ../data/fichera-mixed.mesh
//               mpirun -np 16 ex16p -m ../data/escher-p2.mesh
//               mpirun -np 8 ex16p -m ../data/beam-tet.mesh -tf 10 -dt 0.1
//               mpirun -np 4 ex16p -m ../data/amr-quad.mesh -o 4 -rs 0 -rp 0
//               mpirun -np 4 ex16p -m ../data/amr-hex.mesh -o 2 -rs 0 -rp 0
//
// Description:  This example solves a time dependent nonlinear heat equation
//               problem of the form du/dt = C(u), with a non-linear diffusion
//               operator C(u) = \nabla \cdot (\kappa + \alpha u) \nabla u.
//
//               The example demonstrates the use of nonlinear operators (the
//               class ConductionOperator defining C(u)), as well as their
//               implicit time integration. Note that implementing the method
//               ConductionOperator::ImplicitSolve is the only requirement for
//               high-order implicit (SDIRK) time integration. In this example,
//               the diffusion operator is linearized by evaluating with the
//               lagged solution from the previous timestep, so there is only
//               a linear solve. Optional saving with ADIOS2
//               (adios2.readthedocs.io) is also illustrated.
//
//               We recommend viewing examples 2, 9 and 10 before viewing this
//               example.

#include "mfem.hpp"
#include <vector>
#include <fstream>
#include <iostream>

////#include "include/Visualisation.hpp"
#include "include/problemOperators/navierStokesSteady.hpp"

using namespace std;
using namespace mfem;

void topSurface(const Vector & x, Vector & u){
  unsigned nDim = x.Size();
  u.SetSize(nDim);
  u = 0.00;
  real_t  x1 = 8.00;
  real_t  x0 = 0.00;
  real_t  f_norm = 0.25*(x0 - x1)*(x1 - x0);
  real_t  f = (x(0) - x1)*(x(0) - x0);
  u(0) = 5.0*f/f_norm;
}


void InletSurface(const Vector & x, Vector & u){
  unsigned nDim = x.Size();
  u.SetSize(nDim);
  u = 0.00;
  real_t  x1 = 0.95;
  real_t  x0 = 0.35;
  real_t  f_norm = 0.25*(x0 - x1)*(x1 - x0);
  real_t  f = (x(1) - x1)*(x(1) - x0);

  u(1) = 0.00;
  if( (x(1) < x1)and(x(1) > x0) ) u(1) = 2.0*f/f_norm;
}


int main(int argc, char *argv[])
{
   // 1. Initialize MPI and HYPRE.
   Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();

   // 2. Parse command-line options.
   const char *mesh_file = "mesh/beam-quad.mesh";
   int ser_ref_levels = 3;
   int par_ref_levels = 1;
   int order = 2;
   real_t dt = 1.0e-2;
   real_t alpha = 1.0e-2;
   real_t kappa = 0.5;
   bool visualization = true;
   int vis_steps = 5;

   int precision = 8;
   cout.precision(precision);

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                  "Number of times to refine the mesh uniformly in serial.");
   args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                  "Number of times to refine the mesh uniformly in parallel.");
   args.AddOption(&order, "-o", "--order",
                  "Order (degree) of the finite elements.");
   args.AddOption(&dt, "-dt", "--time-step",
                  "Time step.");
   args.AddOption(&alpha, "-a", "--alpha",
                  "Alpha coefficient.");
   args.AddOption(&kappa, "-k", "--kappa",
                  "Kappa coefficient offset.");

   args.Parse();
   if (!args.Good()){
      args.PrintUsage(cout);
      return 1;
   }

   if (myid == 0) args.PrintOptions(cout);

   // 3. Read the serial mesh from the given mesh file on all processors. We can
   //    handle triangular, quadrilateral, tetrahedral and hexahedral meshes
   //    with the same code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();


   // 5. Refine the mesh in serial to increase the resolution. In this example
   //    we do 'ser_ref_levels' of uniform refinement, where 'ser_ref_levels' is
   //    a command-line parameter.
   for (int lev = 0; lev < ser_ref_levels; lev++) mesh->UniformRefinement();

   // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh further in parallel to increase the resolution. Once the
   //    parallel mesh is defined, the serial mesh can be deleted.
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;
   for (int lev = 0; lev < par_ref_levels; lev++) pmesh->UniformRefinement();

   // 7. Define the vector finite element space representing the velocity and the
   //    pressure.
   Array<int> btoffs(3);
   H1_FECollection fe_coll1(order+1, dim);
   H1_FECollection fe_coll2(order,   dim);
   ParFiniteElementSpace u_fes(pmesh, &fe_coll1,dim);
   ParFiniteElementSpace p_fes(pmesh, &fe_coll2);
   btoffs[0] = 0;
   btoffs[1] = u_fes.GetTrueVSize();
   btoffs[2] = p_fes.GetTrueVSize();
   btoffs.PartialSum();

   HYPRE_BigInt fe_size = u_fes.GlobalTrueVSize() +  p_fes.GlobalTrueVSize();
   if (myid == 0) cout << "Number of unknowns: " << fe_size << endl;
   ParGridFunction u_gf(&u_fes), p_gf(&p_fes);
   BlockVector xBlock(btoffs), zero_Vec(btoffs);


   /*****************************************\
   !
   ! Construct the boundary and initial
   ! conditions
   !
   \*****************************************/
   // 8. Set the initial conditions for u. All boundaries are considered
   //    natural.
   Vector ZeroVecND(dim); ZeroVecND=0.00;
   VectorFunctionCoefficient TopBC(dim, InletSurface);
   VectorConstantCoefficient OtherBCs(ZeroVecND);

   //Apply the initial conditions (zeroed)
   u_gf = 0.00;
   p_gf = 0.00;

   //The BC Tags
   int nTags = pmesh->bdr_attributes.Max();
   Array<int> totMarkers(nTags), otMarkers(nTags), tpMarkers(nTags);
   totMarkers = 1;
   totMarkers[1] = 0;
   otMarkers=1;    tpMarkers=0;
   otMarkers[0]=0; tpMarkers[0]=1;
   otMarkers[1]=0;

   //Apply the essential BC's to the
   u_gf.ProjectCoefficient(TopBC);
   u_gf.ProjectBdrCoefficient(TopBC   , tpMarkers);
   u_gf.ProjectBdrCoefficient(OtherBCs, otMarkers);
   xBlock.GetBlock(0) = u_gf;
   xBlock.GetBlock(1) = p_gf;

   /*****************************************\
   !
   ! Construct the steady-state Navier-Stokes
   !               operator
   !
   \*****************************************/
   navierStokesSteadyOper snsOper(&u_fes, &p_fes, dim, dt);
   snsOper.SetVelocityBCTags(totMarkers);
   snsOper.SetVelocityBCGfs(u_gf);

   /*****************************************\
   !
   ! Preprocess data
   !
   \*****************************************/
   ParaViewDataCollection paraview_dc("NavierStokes", pmesh);
   paraview_dc.SetPrefixPath("ParaView");
   paraview_dc.SetLevelsOfDetail(order);
   paraview_dc.SetDataFormat(VTKFormat::ASCII);
   paraview_dc.SetHighOrderOutput(true);
   paraview_dc.SetCycle(0);
   paraview_dc.SetTime(0.0);
   paraview_dc.RegisterField("Velocity",&u_gf);
   paraview_dc.RegisterField("Pressure",&p_gf);
   paraview_dc.Save();

   /*****************************************\
   !
   ! Run the Non-linear problem
   !
   \*****************************************/
   MINRESSolver solver(MPI_COMM_WORLD);
   solver.SetPreconditioner( snsOper.GetPrecon() );

   NewtonSolver newton(MPI_COMM_WORLD);
   newton.SetOperator(snsOper);
   newton.SetSolver(solver);
   newton.SetPrintLevel(1);
   newton.SetRelTol(1e-10);
   newton.SetMaxIter(60);

   // 10. Solve the nonlinear system.
   zero_Vec = 0.0;
   newton.Mult(zero_Vec, xBlock);

   /*****************************************\
   !
   ! Output the data
   !
   \*****************************************/
   u_gf.Distribute(xBlock.GetBlock(0));
   p_gf.Distribute(xBlock.GetBlock(1));

   paraview_dc.SetCycle(1);
   paraview_dc.SetTime(1.0);
   paraview_dc.RegisterField("Velocity",&u_gf);
   paraview_dc.RegisterField("Pressure",&p_gf);
   paraview_dc.Save();


   // 12. Free the used memory.
   delete pmesh;
   return 0;
}
