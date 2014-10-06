
#define SOLVE_POISSON solve_poisson_
#define SOLVE_POISSON2 solve_poisson2_

#include <iostream>

#include <HYPRE_struct_ls.h>
#include <HYPRE_struct_mv.h>

#include <_hypre_utilities.h>



extern "C" 
void SOLVE_POISSON(
    const double *adiag, // [nx+2,ny+2]
    const double *bcoefx, // [nx+1,ny+2]
    const double *bcoefy, // [nx+2,ny+1]
    const double *rhs, 
    double *sol,
    const int &nx, const int &ny,
    int *period)
{
    int num_procs = 1;
    int my_id = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    if (num_procs!=1 || my_id!= 0) {
        fprintf(stderr, "MPI error\n"); exit(1);
    }
    
    const int ndim = 2;
    const int nentry = ndim*2 + 1;
    int domlo[ndim] = { 1, 1 };
    int domhi[ndim] = { nx, ny };
    // int period[ndim] = { 0, 0 };
    
    // HYPRE stuffs
    struct hypre_StructGrid_struct *hypgrid;
    struct hypre_StructMatrix_struct *hypA;
    struct hypre_StructVector_struct *hypb, *hypx;
    
    // create HYPRE grid
    HYPRE_StructGridCreate(MPI_COMM_WORLD, ndim, &hypgrid);
    HYPRE_StructGridSetPeriodic(hypgrid, period);
    HYPRE_StructGridSetExtents(hypgrid, domlo, domhi);
    // grid ok
    HYPRE_StructGridAssemble(hypgrid);
    
    {
        // define 5-point Laplacian
        int offsets[nentry][ndim] = {
            { -1,  0 }, // 0, XLO 
            {  0, -1 }, // 1, YLO
            {  1,  0 }, // 2, XHI
            {  0,  1 }, // 3, YHI
            {  0,  0 }, // 4
        };
        
        //
        HYPRE_StructStencil stencil;
        HYPRE_StructStencilCreate(ndim, nentry, &stencil);
        for (int i=0; i<nentry; i++) {
            HYPRE_StructStencilSetElement(stencil, i, offsets[i]);
        }
        
        // create matrix
        HYPRE_StructMatrixCreate(MPI_COMM_WORLD, hypgrid, stencil, &hypA);
        HYPRE_StructMatrixInitialize(hypA);
        
        // discard stencil
        HYPRE_StructStencilDestroy(stencil);
    }
    
    {
        // define vector
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, hypgrid, &hypb);
        HYPRE_StructVectorInitialize(hypb);
        
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, hypgrid, &hypx);
        HYPRE_StructVectorInitialize(hypx);
    }
    
    // 
    double *tmp = new double[nx*ny];
    
    
    { // assign initial guess solution, sol -> x
        for (int j=1; j<=ny; j++) {
            for (int i=1; i<=nx; i++) {
                // int idx = i + j*(nx+2);
                // tmp[cnt] = sol[idx];
                tmp[(i-1)+(j-1)*nx] = 0;
                // tmp[(i-1)+(j-1)*nx] = sol[i+j*(nx+2)];
            }
        }
        
        HYPRE_StructVectorSetBoxValues(hypx, domlo, domhi, tmp);
    }
    { // assign RHS
        for (int j=1; j<=ny; j++) {
            for (int i=1; i<=nx; i++) {
                tmp[(i-1)+(j-1)*nx] = rhs[i+j*(nx+2)];
            }
        }
        HYPRE_StructVectorSetBoxValues(hypb, domlo, domhi, tmp);
    }
    
    { // assign matrix coefficient, a,bx,by -> A
        // double *mat = hypre_CTAlloc(double, nentry*nx*ny);
        double *mat = new double[nentry*nx*ny];
        
        for (int j=1; j<=ny; j++) {
            for (int i=1; i<=nx; i++) {
                int off = (i-1)*nentry + (j-1)*nx*nentry;
                { // diagonal
                    mat[off+4] = adiag[i+j*(nx+2)];
                    mat[off+4] = 0.0;
                }
                if (i != 1 || true) { //
                    mat[off+0] = -bcoefx[(i-1)+j*(nx+1)];
                    mat[off+4] += bcoefx[(i-1)+j*(nx+1)];
                } else {
                    mat[off+0] = 0;
                }
                if (i != nx || true) {
                    mat[off+2] = -bcoefx[i+j*(nx+1)];
                    mat[off+4] += bcoefx[i+j*(nx+1)];
                } else {
                    mat[off+2] = 0;
                }
                if (j != 1 || true) {
                    mat[off+1] = -bcoefy[i+(j-1)*(nx+2)];
                    mat[off+4] += bcoefy[i+(j-1)*(nx+2)];
                } else {
                    mat[off+1] = 0;
                }
                if (j != ny || true) {
                    mat[off+3] = -bcoefy[i+j*(nx+2)];
                    mat[off+4] += bcoefy[i+j*(nx+2)];
                } else {
                    mat[off+3] = 0;
                }
                
                //tmp[(i-1)+(j-1)*nx] = rhs[i+j*(nx+2)];
            }
        }
        
        //
        int indices[nentry];
        for (int i=0; i<nentry; i++) { indices[i] = i; }
        HYPRE_StructMatrixSetBoxValues(hypA, domlo, domhi, nentry, indices, mat);
        
        //
        // HYPRE_StructVectorSetBoxValues(hypb, domlo, domhi, tmp);
        
        // hypre_TFree(mat);
        delete[] mat;
    }
    
    // final pack
    HYPRE_StructMatrixAssemble(hypA);
    HYPRE_StructVectorAssemble(hypb);
    HYPRE_StructVectorAssemble(hypx);
    
    //
    HYPRE_StructSolver solver;
    HYPRE_StructSolver precond;
    
    if (0) {
    //
    HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &solver);
    //
    HYPRE_StructPFMGSetMaxIter(solver, 2000);
    HYPRE_StructPFMGSetTol(solver, 1.0e-6);
    //
    HYPRE_StructPFMGSetNumPreRelax(solver, 2);
    HYPRE_StructPFMGSetNumPostRelax(solver, 2);
    //
    HYPRE_StructPFMGSetLogging(solver, 1);
    HYPRE_StructPFMGSetPrintLevel(solver, 3);
    
    HYPRE_StructPFMGSetup(solver, hypA, hypb, hypx);
    }
    if (1) {
        //
        HYPRE_StructBiCGSTABCreate(MPI_COMM_WORLD, &solver);
        HYPRE_StructBiCGSTABSetMaxIter(solver, 2000);
        HYPRE_StructBiCGSTABSetTol(solver, 1.0e-6);
        HYPRE_StructBiCGSTABSetLogging(solver, 1);
        
        //
        HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &precond);
        HYPRE_StructPFMGSetMaxIter(precond, 1);
        HYPRE_StructPFMGSetTol(precond, 0.0);
        HYPRE_StructPFMGSetZeroGuess(precond);
        HYPRE_StructPFMGSetRelaxType(precond, 2);
        HYPRE_StructPFMGSetNumPreRelax(precond, 1);
        HYPRE_StructPFMGSetNumPostRelax(precond, 1);
        HYPRE_StructPFMGSetLogging(precond, 1);
        
        // HYPRE_StructSMGCreate(MPI_COMM_WORLD, &precond);
        // HYPRE_StructSMGSetMemoryUse(precond, 0);
        // HYPRE_StructSMGSetMaxIter(precond, 1);
        // HYPRE_StructSMGSetTol(precond, 0.0);
        // HYPRE_StructSMGSetZeroGuess(precond);
        // HYPRE_StructSMGSetNumPreRelax(precond, 1);
        // HYPRE_StructSMGSetNumPostRelax(precond, 1);
        
        // HYPRE_StructJacobiCreate(MPI_COMM_WORLD, &precond);
        // HYPRE_StructJacobiSetMaxIter(precond, 1);
        // HYPRE_StructJacobiSetTol(precond, 0.0);
        // HYPRE_StructJacobiSetZeroGuess(precond);
        
        HYPRE_StructBiCGSTABSetPrecond(solver, 
            HYPRE_StructPFMGSolve, HYPRE_StructPFMGSetup, precond);
        // HYPRE_StructBiCGSTABSetPrecond(solver, 
            // HYPRE_StructSMGSolve, HYPRE_StructSMGSetup, precond);
        // HYPRE_StructBiCGSTABSetPrecond(solver, 
            // HYPRE_StructJacobiSolve, HYPRE_StructJacobiSetup, precond);
        
        HYPRE_StructBiCGSTABSetup(solver, hypA, hypb, hypx);
    }
    
    //
    int ret;
    int niter;
    
    if (0) {
    ret = HYPRE_StructPFMGSolve(solver, hypA, hypb, hypx);
    HYPRE_StructPFMGGetNumIterations(solver, &niter);
    }
    if (1) {
        ret = HYPRE_StructBiCGSTABSolve(solver, hypA, hypb, hypx);
        HYPRE_StructBiCGSTABGetNumIterations(solver, &niter);
    }
    std::cout << "HYPRE return=" << ret << "; iter=" << niter << std::endl;
    
    // get solution
    HYPRE_StructVectorGetBoxValues(hypx, domlo, domhi, tmp);
    
    for (int j=1; j<=ny; j++) {
        for (int i=1; i<=nx; i++) {
            sol[i+j*(nx+2)] = tmp[(i-1)+(j-1)*nx];
        }
    }
    
    if (0) {
    HYPRE_StructPFMGDestroy(solver);
    }
    
    if (1) {
        HYPRE_StructBiCGSTABDestroy(solver);
        HYPRE_StructPFMGDestroy(precond);
        // HYPRE_StructSMGDestroy(precond);
        // HYPRE_StructJacobiDestroy(precond);
    }
    
    delete[] tmp;
    
    {
        // destroy HYPRE objects
        HYPRE_StructGridDestroy(hypgrid);
        HYPRE_StructMatrixDestroy(hypA);
        HYPRE_StructVectorDestroy(hypb);
        HYPRE_StructVectorDestroy(hypx);
    }
    
    // MPI_Finalize();
} // end_SOLVE_POISSON




extern "C" 
void SOLVE_POISSON2(
    double *mat,
    double *rhs, 
    double *sol,
    const int &nx, const int &ny)
{
    int num_procs = 1;
    int my_id = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    if (num_procs!=1 || my_id!= 0) {
        fprintf(stderr, "MPI error\n"); exit(1);
    }
    
    const int ndim = 2;
    const int nentry = ndim*2 + 1;
    int domlo[ndim] = { 1, 1 };
    int domhi[ndim] = { nx, ny };
    int period[ndim] = { 0, 0 };
    
    // HYPRE stuffs
    struct hypre_StructGrid_struct *hypgrid;
    struct hypre_StructMatrix_struct *hypA;
    struct hypre_StructVector_struct *hypb, *hypx;
    
    // create HYPRE grid
    HYPRE_StructGridCreate(MPI_COMM_WORLD, ndim, &hypgrid);
    HYPRE_StructGridSetPeriodic(hypgrid, period);
    HYPRE_StructGridSetExtents(hypgrid, domlo, domhi);
    // grid ok
    HYPRE_StructGridAssemble(hypgrid);
    
    {
        // define 5-point Laplacian
        int offsets[nentry][ndim] = {
            { -1,  0 }, // 0, XLO 
            {  0, -1 }, // 1, YLO
            {  1,  0 }, // 2, XHI
            {  0,  1 }, // 3, YHI
            {  0,  0 }, // 4
        };
        
        //
        HYPRE_StructStencil stencil;
        HYPRE_StructStencilCreate(ndim, nentry, &stencil);
        for (int i=0; i<nentry; i++) {
            HYPRE_StructStencilSetElement(stencil, i, offsets[i]);
        }
        
        // create matrix
        HYPRE_StructMatrixCreate(MPI_COMM_WORLD, hypgrid, stencil, &hypA);
        HYPRE_StructMatrixInitialize(hypA);
        
        // discard stencil
        HYPRE_StructStencilDestroy(stencil);
    }
    
    {
        // define vector
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, hypgrid, &hypb);
        HYPRE_StructVectorInitialize(hypb);
        
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, hypgrid, &hypx);
        HYPRE_StructVectorInitialize(hypx);
    }
    
    { // assign initial guess solution, sol -> x
        HYPRE_StructVectorSetBoxValues(hypx, domlo, domhi, sol);
    }
    { // assign RHS
        HYPRE_StructVectorSetBoxValues(hypb, domlo, domhi, rhs);
    }
    { // assign matrix coefficient, a,bx,by -> A
        int indices[nentry];
        for (int i=0; i<nentry; i++) { 
            indices[i] = i; 
        }
        HYPRE_StructMatrixSetBoxValues(hypA, domlo, domhi, nentry, indices, mat);
    }
    
    // final pack
    HYPRE_StructMatrixAssemble(hypA);
    HYPRE_StructVectorAssemble(hypb);
    HYPRE_StructVectorAssemble(hypx);
    
    //
    HYPRE_StructSolver solver;
    HYPRE_StructSolver precond;
    
    if (0) {
    //
    HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &solver);
    //
    HYPRE_StructPFMGSetMaxIter(solver, 2000);
    HYPRE_StructPFMGSetTol(solver, 1.0e-6);
    //
    HYPRE_StructPFMGSetNumPreRelax(solver, 2);
    HYPRE_StructPFMGSetNumPostRelax(solver, 2);
    //
    HYPRE_StructPFMGSetLogging(solver, 1);
    HYPRE_StructPFMGSetPrintLevel(solver, 3);
    
    HYPRE_StructPFMGSetup(solver, hypA, hypb, hypx);
    }
    if (1) {
        //
        // HYPRE_StructBiCGSTABCreate(MPI_COMM_WORLD, &solver);
        // HYPRE_StructBiCGSTABSetMaxIter(solver, 2000);
        // HYPRE_StructBiCGSTABSetTol(solver, 1.0e-5);
        // HYPRE_StructBiCGSTABSetLogging(solver, 1);
        
        HYPRE_StructGMRESCreate(MPI_COMM_WORLD, &solver);
        HYPRE_StructGMRESSetMaxIter(solver, 2000);
        HYPRE_StructGMRESSetTol(solver, 1.0e-5);
        HYPRE_StructGMRESSetLogging(solver, 1);
        
        //
        HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &precond);
        HYPRE_StructPFMGSetMaxIter(precond, 1);
        HYPRE_StructPFMGSetTol(precond, 0.0);
        // HYPRE_StructPFMGSetZeroGuess(precond);
        // HYPRE_StructPFMGSetRAPType(precond, 1);
        // HYPRE_StructPFMGSetRelaxType(precond, 3);
        HYPRE_StructPFMGSetNumPreRelax(precond, 1);
        HYPRE_StructPFMGSetNumPostRelax(precond, 1);
        HYPRE_StructPFMGSetLogging(precond, 1);
        
        // HYPRE_StructSMGCreate(MPI_COMM_WORLD, &precond);
        // HYPRE_StructSMGSetMemoryUse(precond, 0);
        // HYPRE_StructSMGSetMaxIter(precond, 1);
        // HYPRE_StructSMGSetTol(precond, 0.0);
        // HYPRE_StructSMGSetZeroGuess(precond);
        // HYPRE_StructSMGSetNumPreRelax(precond, 1);
        // HYPRE_StructSMGSetNumPostRelax(precond, 1);
        
        // HYPRE_StructJacobiCreate(MPI_COMM_WORLD, &precond);
        // HYPRE_StructJacobiSetMaxIter(precond, 1);
        // HYPRE_StructJacobiSetTol(precond, 0.0);
        // HYPRE_StructJacobiSetZeroGuess(precond);
        
        // HYPRE_StructBiCGSTABSetPrecond(solver, 
            // HYPRE_StructPFMGSolve, HYPRE_StructPFMGSetup, precond);
        // HYPRE_StructBiCGSTABSetPrecond(solver, 
            // HYPRE_StructSMGSolve, HYPRE_StructSMGSetup, precond);
        // HYPRE_StructBiCGSTABSetPrecond(solver, 
            // HYPRE_StructJacobiSolve, HYPRE_StructJacobiSetup, precond);
        
        // HYPRE_StructBiCGSTABSetup(solver, hypA, hypb, hypx);
        
        HYPRE_StructGMRESSetPrecond(solver,
            HYPRE_StructPFMGSolve, HYPRE_StructPFMGSetup, precond);
        
        HYPRE_StructGMRESSetup(solver, hypA, hypb, hypx);
    }
    
    //
    int ret;
    int niter;
    
    if (0) {
    ret = HYPRE_StructPFMGSolve(solver, hypA, hypb, hypx);
    HYPRE_StructPFMGGetNumIterations(solver, &niter);
    }
    if (1) {
        // ret = HYPRE_StructBiCGSTABSolve(solver, hypA, hypb, hypx);
        // HYPRE_StructBiCGSTABGetNumIterations(solver, &niter);
        ret = HYPRE_StructGMRESSolve(solver, hypA, hypb, hypx);
        HYPRE_StructGMRESGetNumIterations(solver, &niter);
    }
    std::cout << "HYPRE return=" << ret << "; iter=" << niter << std::endl;
    
    // get solution
    HYPRE_StructVectorGetBoxValues(hypx, domlo, domhi, sol);
        
    if (0) {
    HYPRE_StructPFMGDestroy(solver);
    }
    
    if (1) {
        // HYPRE_StructBiCGSTABDestroy(solver);
        HYPRE_StructGMRESDestroy(solver);
        HYPRE_StructPFMGDestroy(precond);
        // HYPRE_StructSMGDestroy(precond);
        // HYPRE_StructJacobiDestroy(precond);
    }
    
    {
        // destroy HYPRE objects
        HYPRE_StructGridDestroy(hypgrid);
        HYPRE_StructMatrixDestroy(hypA);
        HYPRE_StructVectorDestroy(hypb);
        HYPRE_StructVectorDestroy(hypx);
    }
    
    // MPI_Finalize();
} // end_SOLVE_POISSON2


























