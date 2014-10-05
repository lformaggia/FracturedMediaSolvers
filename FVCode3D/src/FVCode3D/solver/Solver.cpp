/*!
 * @file solver.cpp
 * @brief These classes allow to solve a linear system (definitions).
 */

#include <FVCode3D/solver/Solver.hpp>
#ifdef FVCODE3D_HAS_UMFPACK
#include <Eigen/UmfPackSupport>
#endif // FVCODE3D_HAS_UMFPACK

#ifdef FVCODE3D_HAS_SAMG
#include <samg.h>
#endif // FVCODE3D_HAS_SAMG

namespace FVCode3D
{

void EigenCholesky::solve()
{
    Eigen::SimplicialCholesky<SpMat> chol(M_A);
    M_x = chol.solve(M_b);
} // EigenCholesky::solve

void EigenLU::solve()
{
    Eigen::SparseLU<SpMat> lu;
    lu.analyzePattern(M_A);
    lu.factorize(M_A);
    M_x = lu.solve(M_b);
} // EigenLU::solve

#ifdef FVCODE3D_HAS_UMFPACK
void EigenUmfPack::solve()
{
    Eigen::UmfPackLU<SpMat> lu( M_A );
    M_x = lu.solve( M_b );
} // EigenUmfPack::solve
#endif // FVCODE3D_HAS_UMFPACK

Real IterativeSolver::S_referenceTol = 1e-6;
UInt IterativeSolver::S_referenceMaxIter = 100;

void EigenCG::solve()
{
    Eigen::ConjugateGradient<SpMat> cg;

    cg.setMaxIterations(M_maxIter);
    cg.setTolerance(M_tol);

    cg.compute(M_A);
    M_x = cg.solve(M_b);

    M_iter = cg.iterations();
    M_res = cg.error();
} // EigenCG::solve

void EigenBiCGSTAB::solve()
{
    Eigen::BiCGSTAB<SpMat> bicgstab;

    bicgstab.setMaxIterations(M_maxIter);
    bicgstab.setTolerance(M_tol);

    bicgstab.compute(M_A);
    M_x = bicgstab.solve(M_b);

    M_iter = bicgstab.iterations();
    M_res = bicgstab.error();
} // EigenBiCGSTAB::solve

#ifdef FVCODE3D_HAS_SAMG
void Samg::solve()
{
    // set somehow the max iter
    // ... M_maxIter

    // set somehow the tolerance
    // ... M_tol

    // copy from eigen to samg

    // solve with samg

    // copy solution from samg to eigen
    // ... M_x

    // get iteration used
    // ... M_iter

    // get residual
    // ... M_res
/*
      cout << " *** Demo C++ driver which reads demo.* files and calls SAMG ***" << endl;

// ===> Set primary parameters. Others can be set by access functions as shown below.

      int    ndiu      = 1;        // dimension of (dummy) vector iu
      int    ndip      = 1;        // dimension of (dummy) vector ip

      int    nsolve    = 2;        // results in scalar approach (current system is scalar)
      int    ifirst    = 1;        // first approximation = zero
      double eps       = 1.0e-8;   // required (relative) residual reduction
      int    ncyc      = 11050;    // V-cycle as pre-conditioner for CG; at most 50 iterations

      double a_cmplx   = 2.2;      // estimated dimensioning
      double g_cmplx   = 1.7;      // estimated dimensioning
      double w_avrge   = 2.4;      // estimated dimensioning
      double p_cmplx   = 0.0;      // estimated dimensioning (irrelevant for scalar case)

      double chktol    = -1.0;     // input checking de-activated (we know it's ok!)
      int    idump     = 0;        // minimum output during setup
      int    iout      = 2;        // display residuals per iteration and work statistics

      int    n_default = 20;       // select default settings for secondary parameters
                                   // CURRENTLY AVAILABLE: 10-13, 15-18, 20-23, 25-28
                                   // NOTE: the higher the respective SECOND digit, the
                                   // more aggressive the coarsening (--> lower memory at
                                   // the expense of slower convergence)
      int    iswtch    = 5100+n_default; // complete SAMG run ....
                                   // ... memory de-allocation upon return ....
                                   // ... memory extension feature activated ....
                                   // ... residuals measured in the L2-norm .....
                                   // ... secondary parameter default setting # n_default

// ===> Secondary parameters which have to be set if n_default=0
//      (at the same time a demonstration of how to access secondary or hidden parameters)

      int intin;
      double dblin;

      if (n_default == 0) {
         intin=25;    SAMG_SET_LEVELX(&intin);
         intin=100;   SAMG_SET_NPTMN(&intin);
         intin=4;     SAMG_SET_NCG(&intin);
         intin=2;     SAMG_SET_NWT(&intin);
         intin=1;     SAMG_SET_NTR(&intin);
         intin=131;   SAMG_SET_NRD(&intin);
         intin=131;   SAMG_SET_NRU(&intin);
         intin=0;     SAMG_SET_NRC(&intin);
         intin=0;     SAMG_SET_NP_OPT(&intin);

         dblin=21.25; SAMG_SET_ECG(&dblin);
         dblin=0.20;  SAMG_SET_EWT(&dblin);
         dblin=12.20; SAMG_SET_ETR(&dblin);
      }

// ===> amg declarations

      // input:

      int npnt,nsys,matrix,nnu,nna;

      int * ia, * ja;
      int * iu     = new int[1];
      int * ip     = new int[1];
      int * iscale = new int[1];

      double * a, * u, * f;

      // output:
      int ierr,ierrl,ncyc_done;
      double res_out,res_in;

// ===> Read data from stdin

      // header data

      ifstream frmfile("demo.frm", ios::in);    // open demo.frm

      int iversion;
      char ch;

      frmfile >> ch >> iversion;
      frmfile >> nna >> nnu >> matrix >> nsys >> npnt;

      frmfile.close();

      if (ch!='f' || iversion != 4) {
          cout << "invalid file format for this test driver " << endl;
          ierr=1; return ierr;
      }

      // matrix

      int i;

      ia = new int[nnu+1];
      ja = new int[nna];
      a  = new double[nna];

      if (!(ia && ja && a)) {
          cout << " allocation failed (ia,ja,a) " << endl;
          ierr=1; return ierr;
      }

      ifstream amgfile("demo.amg", ios::in);    // open demo.amg

      for (i=0;i<nnu+1;i++) amgfile >> ia[i];
      for (i=0;i<nna;i++)   amgfile >> ja[i];
      for (i=0;i<nna;i++)   amgfile >> a[i];

      amgfile.close();

      // right hand side

      f = new double[nnu]; u = new double[nnu];
      if (!(f && u)) {
          cout << " allocation failed (f,u) " << endl;
          ierr=1; return ierr;
      }
      for (i=0;i<nnu;i++) u[i]=0.0;

      ifstream rhsfile("demo.rhs", ios::in);    // open demo.rhs

      for (i=0;i<nnu;i++) rhsfile >> f[i];

      rhsfile.close();

// ===> if, e.g., the finest-level matrix shall be dumped:
//
//    char *string ="myfilename";
//    int length=10;
//    SAMG_SET_FILNAM_DUMP(&string[0],&length,&length);
//    idump=8;

// ===> Call FORTRAN90 SAMG

//    SAMG_RESET_SECONDARY();   // necessary before a second SAMG run
                            // if secondary parameters have to be reset
                            // see manual.

      float told,tnew,tamg;
      SAMG_CTIME(&told);

      SAMG(&nnu,&nna,&nsys,
           &ia[0],&ja[0],&a[0],&f[0],&u[0],&iu[0],&ndiu,&ip[0],&ndip,&matrix,&iscale[0],
           &res_in,&res_out,&ncyc_done,&ierr,
           &nsolve,&ifirst,&eps,&ncyc,&iswtch,
           &a_cmplx,&g_cmplx,&p_cmplx,&w_avrge,
           &chktol,&idump,&iout);

      if (ierr > 0) {
          cout << endl << " SAMG terminated with error code "
               << ierr << " **** " << endl;
      }
      else if (ierr < 0) {
          cout << endl << " SAMG terminated with warning code "
               << ierr << " **** " << endl;
      }

      SAMG_CTIME(&tnew);
      tamg=tnew-told;
      cout << endl << " ***** total run time: " << tamg << " ***** " << endl;

      SAMG_LEAVE(&ierrl);

      if (ierrl != 0) {
          cout << endl << " error at samg_leave"
               << ierr << " **** " << endl;
      }

      delete[] ia,ja,a,f,u,iu,ip,iscale;

      if (ierr == 0 && ierrl == 0) return 0;
      if (ierr > 0 || ierrl > 0) return 2;
      return 1;
*/
} // Samg::solve
#endif // FVCODE3D_HAS_SAMG

} // namespace FVCode3D
