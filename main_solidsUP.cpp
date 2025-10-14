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

#include "include/problemOperators/fibreMapPoisson.hpp"
#include "include/problemOperators/monodomain.hpp"
#include "include/Visualisation.hpp"
#include "include/coefficients/orthoDiffCoeff.hpp"
#include "include/coefficients/PK2StressCoeff.hpp"
#include "include/problemOperators/incompActiveStrainSolidMechanics.hpp"


using namespace std;
using namespace mfem;

real_t InitialTemperature(const Vector &x);
real_t InitialTemperature(const Vector &x)
{
   if (x.Norml2() < 0.5)
   {
      return 2.0;
   }
   else
   {
      return 1.0;
   }
};

int main(int argc, char *argv[])
{
   // 1. Initialize MPI and HYPRE.
   Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();

   // 2. Parse command-line options.
   const char *mesh_file = "mesh/beam-quad.mesh";
   int ser_ref_levels = 2;
   int par_ref_levels = 1;
   int order = 2;
   int ode_solver_type = 3;
   real_t t_final = 0.5;
   real_t dt = 1.0e-2;
   real_t alpha = 1.0e-2;
   real_t kappa = 0.5;
   bool visualization = true;
   bool visit = false;
   int vis_steps = 5;
   bool adios2 = false;

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
   args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                  "ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3,\n\t"
                  "\t   11 - Forward Euler, 12 - RK2, 13 - RK3 SSP, 14 - RK4.");
   args.AddOption(&t_final, "-tf", "--t-final",
                  "Final time; start time is 0.");
   args.AddOption(&dt, "-dt", "--time-step",
                  "Time step.");
   args.AddOption(&alpha, "-a", "--alpha",
                  "Alpha coefficient.");
   args.AddOption(&kappa, "-k", "--kappa",
                  "Kappa coefficient offset.");

   args.AddOption(&visit, "-visit", "--visit-datafiles", "-no-visit",
                  "--no-visit-datafiles",
                  "Save data files for VisIt (visit.llnl.gov) visualization.");
   args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                  "Visualize every n-th timestep.");
   args.AddOption(&adios2, "-adios2", "--adios2-streams", "-no-adios2",
                  "--no-adios2-streams",
                  "Save data using adios2 streams.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }

   if (myid == 0)
   {
      args.PrintOptions(cout);
   }

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

   // 7. Define the vector finite element space representing the current and the
   //    initial temperature, u_ref.
   H1_FECollection fe_coll(order, dim);
   ParFiniteElementSpace fespace(pmesh, &fe_coll);
   ParFiniteElementSpace vecFespace(pmesh, &fe_coll,dim);

   HYPRE_BigInt fe_size = fespace.GlobalTrueVSize();
   if (myid == 0) cout << "Number of temperature unknowns: " << fe_size << endl;

   Array<ParGridFunction*> fibre(dim);
   ParGridFunction u_gf(&fespace), f_gf(&fespace);
   ParGridFunction disp_gf(&vecFespace), press_gf(&fespace), gama_gf(&fespace);

   // 8. Set the initial conditions for u. All boundaries are considered
   //    natural.
   FunctionCoefficient u_0(InitialTemperature);
   u_gf.ProjectCoefficient(u_0);
   Vector u, phi;
   u_gf.GetTrueDofs(u);
   u_gf.GetTrueDofs(phi);

   /*****************************************\
   !
   ! Construct the fibre map
   !
   \*****************************************/
   //The BC's
   Array<int>     ess_bcs_markers(pmesh->bdr_attributes.Max());
   Array<double>  BC_Vals(pmesh->bdr_attributes.Max());
   BC_Vals = 0.00;
   ess_bcs_markers = 0;
   BC_Vals[0] =  5.00;
   ess_bcs_markers[0] = 1;
   BC_Vals[1] = -5.00;
   ess_bcs_markers[1] = 1;

   //Construct the operator
   //and solve the system
   fibreMapPoissonOper fpOper(fespace, ess_bcs_markers, BC_Vals);
   f_gf = u_gf;
   fpOper.Mult(f_gf,phi);

   //Construct the fibre field
   for(int I=0; I<dim; I++) fibre[I] = new ParGridFunction(&vecFespace);
   fpOper.getFibreMapGFuncs2D(phi, fibre);

   /*****************************************\
   !
   ! Construct the PK2 stress coefficient
   ! and 
   !
   \*****************************************/
   Vector diff(dim); diff = 0.1;
   orthoDiffCoeff D_ij(&fibre, diff,  dim);
   monodomainOper mnOper(fespace, D_ij);

   /*****************************************\
   !
   ! Setup the hyperelastic solid mechanics
   ! problem using 
   !
   \*****************************************/
   NeoHookeanPK2StressCoeff solidStress(&fibre, &disp_gf, &press_gf, &gama_gf, dim);
/*

disp_gf(&vecFespace), press_gf(&fespace), fibre, gama_gf;
*/

   /*****************************************\
   !
   ! Run the Non-linear solid mechanics
   ! problem
   !
   \*****************************************/
   u_gf.SetFromTrueDofs(u);
   u_gf.SetFromTrueDofs(phi);

   ParaViewDataCollection paraview_dc("IASM_solids", pmesh);
   paraview_dc.SetPrefixPath("ParaView");
   paraview_dc.SetLevelsOfDetail(order);
   paraview_dc.SetDataFormat(VTKFormat::ASCII);
   paraview_dc.SetHighOrderOutput(true);
   paraview_dc.SetCycle(0);
   paraview_dc.SetTime(0.0);
   paraview_dc.RegisterField("Potential",&u_gf);
   paraview_dc.RegisterField("Phi",&f_gf);
   for(int I=0; I<dim; I++) paraview_dc.RegisterField("f_"+to_string(I),fibre[I]);
   paraview_dc.Save();



   // 12. Free the used memory.
   delete ode_solver;
   delete pmesh;
   for(int I=0; I<dim; I++) delete fibre[I];
   return 0;
}
