/*!
 * @file solver.cpp
 * @brief These classes allow to solve a linear system (definitions).
 */

#include <sstream>

#include <FVCode3D/solver/Solver.hpp>
#include <FVCode3D/solver/cg.hpp>
#include <FVCode3D/solver/bicgstab.hpp>
#include <FVCode3D/solver/gmres.hpp>
#include <FVCode3D/preconditioner/preconditioner.hpp>
#include <FVCode3D/preconditioner/preconHandler.hpp>
#include <Eigen/UmfPackSupport>

#ifdef FVCODE3D_HAS_SAMG
#include <samg.h>
#endif // FVCODE3D_HAS_SAMG

namespace FVCode3D
{

void EigenCholesky::solve()
{
    Eigen::SimplicialCholesky<SpMat, Eigen::Upper> chol(M_A);
    M_x = chol.solve(M_b);
} // EigenCholesky::solve


void EigenLU::solve()
{
    Eigen::SparseLU<SpMat> lu;
    lu.analyzePattern(M_A);
    lu.factorize(M_A);
    M_x = lu.solve(M_b);
} // EigenLU::solve


void EigenUmfPack::solve()
{
    Eigen::UmfPackLU<SpMat> lu( M_A );
    M_x = lu.solve( M_b );
} // EigenUmfPack::solve

Real constexpr IterativeSolver::S_referenceTol;
UInt constexpr IterativeSolver::S_referenceMaxIter;

constexpr bool imlBiCGSTAB::Default_restart;
void imlBiCGSTAB::solve()
{
	// Define the initial guess to zero
	M_x = Vector::Zero(M_A.getM().rows()+M_A.getB().rows());
	// Set the restart
//	setRestart(true);
	// Conversion needed
	int iter = (int) M_maxIter;
	// Solve the system
	int conv = BiCGSTAB(M_A, M_x, M_b, *preconPtr, restart, iter, M_res);
	// Conversion needed
	CIndex = (UInt) conv;	
	M_iter = (UInt) iter;

} // imlBiCGSTAB::solve


constexpr UInt imlGMRES::Default_m;
void imlGMRES::solve()
{
	// Define the initial guess to zero
	M_x = Vector::Zero(M_A.getM().cols()+M_A.getB().rows());
	// Set the restart level
//	set_m(60);
	// Conversion needed
	int iter = (int) M_maxIter;
	int m_int = (int) m;
	// Solve the system
	int conv = GMRES(M_A, M_x, M_b, *preconPtr, m_int, iter, M_res);
	// Conversion needed
	CIndex = (UInt) conv;	
	M_iter = (UInt) iter;
	
} // imlGMRES::solve


#ifdef FVCODE3D_HAS_SAMG
void SamgSolver::solve()
{
    // Set general SAMG parameters
    SamgSolver::SamgParameters SP;

    // Set specific SAMG parameters, max-iterations and tolerance
    setSamgParameters(SP);

    // Copy from Eigen to SAMG
    int * ia;   // i indices
    int * ja;   // j indices
    double * a; // A entries
    double * u; // solution
    double * f; // rhs

    const UInt rows = M_A.rows();
    const UInt cols = M_A.cols();
    SP.nnu     = cols;              // Number of unknowns
    SP.nna     = M_A.nonZeros();    // Number of entries
                                    // Hp: diagonal entries are non-zeros

    ia = new int[SP.nnu+1];
    ja = new int[SP.nna];
    a  = new double[SP.nna];

    UInt i=0;
    double tmp;
    for(UInt row=0, i=0; row < rows; ++row)
    {
        ia[row] = static_cast<int>(i+1);
        a[i] = static_cast<double>(M_A.coeff(row, row));
        ja[i] = static_cast<int>(row+1);
        ++i;
        for(UInt col=0; col < row; ++col)
        {
            tmp = M_A.coeff(row, col);
            if (tmp != 0.)
            {
                a[i] = tmp;
                ja[i] = static_cast<int>(col+1);
                ++i;
            }
        }
        for(UInt col=row+1; col < cols; ++col)
        {
            tmp = M_A.coeff(row, col);
            if (tmp != 0.)
            {
                a[i] = tmp;
                ja[i] = static_cast<int>(col+1);
                ++i;
            }
        }
    }
    ia[rows] = SP.nna + 1;

    f = new double[SP.nnu];
    u = new double[SP.nnu];

    for(i=0; i < cols; ++i)
    {
        f[i] = M_b[i];
    }

    for(i=0; i < cols; ++i)
    {
        u[i] = 0.;
    }

    // Solve with SAMG
    SAMG( &SP.nnu, &SP.nna, &SP.nsys,
          &ia[0], &ja[0], &a[0], &f[0], &u[0],
          &SP.iu[0], &SP.ndiu, &SP.ip[0], &SP.ndip,
          &SP.matrix, &SP.iscale[0],
          &SP.res_in, &SP.res_out, &SP.ncyc_done, &SP.ierr,
          &SP.nsolve, &SP.ifirst, &SP.eps, &SP.ncyc, &SP.iswtch,
          &SP.a_cmplx, &SP.g_cmplx, &SP.p_cmplx, &SP.w_avrge,
          &SP.chktol, &SP.idump, &SP.iout);

    if (SP.ierr > 0)
    {
        std::cout << std::endl << " Error occurred at SAMG."
                  << std::endl;
    }
    else if (SP.ierr < 0)
    {
        std::cout << std::endl << " Warning occurred at SAMG."
                  << std::endl;
    }

    // Copy solution from SAMG to Eigen
    M_x.resize( cols );
    for(i=0; i < cols; ++i)
    {
        M_x[i] = u[i];
    }

    // Save iterations done and final residual
    M_iter = SP.ncyc_done;
    M_res = SP.res_out / SP.res_in;

    // Clear SAMG parameters
    SAMG_LEAVE(&SP.ierrl);

    delete[] ia;
    delete[] ja;
    delete[] a;
    delete[] f;
    delete[] u;

    if (SP.ierrl != 0)
    {
        std::cout << std::endl << "Error at samg_leave call."
                  << std::endl;
    }
} // Samg::solve

SamgSolver::SamgParameters::SamgParameters()
{
    // Set primary parameters
    ndiu      = 1;      // Dummy for scalar system, = 1
    ndip      = 1;      // Dummy for scalar system, = 1
    nsolve    = 10;     // Variable-based approach, scalar system, = 1
                        // Smoothing strategy, plain relaxation, = 0
    ifirst    = 1;      // Initialization of the solution. In this case, = 0

    eps       = static_cast<double> // Tolerance on the residual. Res_0 = ifirst = 0. In this case normalized_residual < eps
                (IterativeSolver::S_referenceTol);

    std::stringstream ss;
    ss << 120           // Cycles: V-cycle, =1
                        // Precond: CG for sym matrix, =1
                        //          BiCGstab for not sym matrix, =2
                        // Re-start points: default, =0
       << static_cast<int>  // # of cycles: default iterations, =dafault
          (IterativeSolver::S_referenceMaxIter);
    ss >> ncyc;

    // All 0.0: initially a default dimensioning, it can be automatically extended
    // All 1.0: initially no memory, it can be automatically extended
    a_cmplx   = 0.0;    // 2.2
    g_cmplx   = 0.0;    // 1.7
    w_avrge   = 0.0;    // 2.4
    p_cmplx   = 0.0;    // 0.0

    chktol    = -1.0;   // No checking, <0
                        // Logical checking, =0
                        // Full checking, >0 (very expensive, for tests)
    idump     = -1;      // No printout (except warnings and errors), <0
                        // Minimal output =0
                        // More output, >0
    iout      = 0;      // No printout (except warnings and errors), <0
                        // Minimal output =0
                        // More output, >0
    n_default = 10;     // Coarsening strategy: no "critical" positive off-diagonal entries =10
                        //                      "critical" positive off-diagonal entries =20
    iswtch    = 5100+n_default; // No re-call of SAMG, =5
                                // Memory extension feature activted, =1
                                // Coarsening strategy, see n_default
                                // Residuals measured in the L2-norm, =0 (implicit)

    // Set secondary parameters which have to be set if n_default=0
    int intin;
    double dblin;

    if (n_default == 0)
    {
        intin = 25;    SAMG_SET_LEVELX(&intin);
        intin = 100;   SAMG_SET_NPTMN(&intin);
        intin = 4;     SAMG_SET_NCG(&intin);
        intin = 2;     SAMG_SET_NWT(&intin);
        intin = 1;     SAMG_SET_NTR(&intin);
        intin = 131;   SAMG_SET_NRD(&intin);
        intin = 131;   SAMG_SET_NRU(&intin);
        intin = 0;     SAMG_SET_NRC(&intin);
        intin = 0;     SAMG_SET_NP_OPT(&intin);
//        intin = 0;     SAMG_SET_ICOLOR_OMP(&intin);
//        intin = 0;     SAMG_SET_IRESTRICTION_OPENMP(&intin);
//        intin = 0;     SAMG_SET_NTYP_GALERKIN(&intin);
        //intin = -3;    SAMG_SET_MODE_MESS(&intin); // No output at all

        dblin = 21.25; SAMG_SET_ECG(&dblin);
        dblin = 0.20;  SAMG_SET_EWT(&dblin);
        dblin = 12.20; SAMG_SET_ETR(&dblin);
    }

    // AMG declarations
    nsys    = 1;            // Scalar system, =1
    npnt    = 0;            // Number of points
    matrix = 22;         // Sym or not: not sym, =2
                            //             sym, =1
                            // Sum zero or not: not sum zero, =2

    iu      = new int[1];   // Dummy for scalar system, = 0
    iu[0]   = 0;

    ip      = new int[1];   // Dummy for scalar system, = 0
    ip[0]   = 0;

    iscale  = new int[1];   // No scaling of solution, = 0
    iscale[0] = 0;
} // SamgSolver::SamgParameters::SamgParameters

SamgSolver::SamgParameters::~SamgParameters()
{
    delete[] iu;
    delete[] ip;
    delete[] iscale;
} // SamgSolver::SamgParameters::~SamgParameters

void SamgSym::setSamgParameters(SamgSolver::SamgParameters & SP)
{
    std::stringstream ss;
    ss << 110           // Cycles: V-cycle, =1
                        // Precond: CG for sym matrix, =1
                        //          BiCGstab for not sym matrix, =2
                        // Re-start points: default, =0
       << static_cast<int>(M_maxIter);  // # of cycles: M_maxIter iterations, =M_maxIter
    ss >> SP.ncyc;

    SP.eps = static_cast<double>(M_tol); // Tolerance on the residual. Res_0 = ifirst = 0. In this case normalized_residual < eps

    SP.matrix = 12;     // Sym or not: not sym, =2
                        //             sym, =1
                        // Sum zero or not: not sum zero, =2
} // SamgSym::setSamgParameters

void SamgNotSym::setSamgParameters(SamgSolver::SamgParameters & SP)
{
    std::stringstream ss;
    ss << 120           // Cycles: V-cycle, =1
                        // Precond: CG for sym matrix, =1
                        //          BiCGstab for not sym matrix, =2
                        // Re-start points: default, =0
       << static_cast<int>(M_maxIter);  // # of cycles: M_maxIter iterations, =M_maxIter
    ss >> SP.ncyc;

    SP.eps = static_cast<double>(M_tol); // Tolerance on the residual. Res_0 = ifirst = 0. In this case normalized_residual < eps

    SP.matrix = 22;     // Sym or not: not sym, =2
                        //             sym, =1
                        // Sum zero or not: not sum zero, =2
} // SamgNotSym::setSamgParameters
#endif // FVCODE3D_HAS_SAMG

} // namespace FVCode3D
