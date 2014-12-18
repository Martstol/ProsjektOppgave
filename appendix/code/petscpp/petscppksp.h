#ifndef PETSC_CPP_KSP_H
#define PETSC_CPP_KSP_H

#include <petscksp.h>

#include "petscppmatrix.h"
#include "petscppvector.h"

namespace petscpp {

    class KspSolver {
    private:
        KSP ksp;
    public:
        KspSolver();
        ~KspSolver();
		
        KspSolver & operator=(const KspSolver&) = delete;
        KspSolver(const KspSolver&) = delete;
		
        void setTolerances(PetscReal rtol, PetscReal abstol, PetscReal dtol, PetscInt maxits);
        void setOperators(Matrix &A, Matrix &P);
        KSP getPetscKSP();
        void setup();
        void solve(Vector &b, Vector &x);
        void createPreconditioner();
		
		void printDebug();
		void printConvergenceReason();
		void printIterationNumber();
    };

}

#endif
