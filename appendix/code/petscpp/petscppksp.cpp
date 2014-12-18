// This file is created by students of the HPC-Lab at the Norwegian University of Science and Technology (NTNU)
// and is part of the HPC-Lab Snow Simulator application distributed under the GPL license. 
// Copyright (c) 2006-2013 High Performance Lab at the Norwegian University of Science and Technology (NTNU)
// Department of Computer and Information Science (IDI). All rights reserved.
// See the file README.md for more information.

#include <stdexcept>

#include <petscsys.h>

#include "petscpp.h"
#include "petscppksp.h"
#include "petscpperror.h"

petscpp::KspSolver::KspSolver() {
    PetscErrorCode err = KSPCreate(petscpp::worldComm(), &ksp);
    petscpp::handleError(err);
}

petscpp::KspSolver::~KspSolver() {
    PetscErrorCode err = KSPDestroy(&ksp);
    petscpp::handleError(err);
}

void petscpp::KspSolver::setTolerances(PetscReal rtol, PetscReal abstol, PetscReal dtol, PetscInt maxits) {
    PetscErrorCode err = KSPSetTolerances(ksp, rtol, abstol, dtol, maxits);
    petscpp::handleError(err);
}

void petscpp::KspSolver::setOperators(Matrix &A, Matrix &P) {
    PetscErrorCode err = KSPSetOperators(ksp, A.getPetscMatrix(), P.getPetscMatrix());
    petscpp::handleError(err);
}

KSP petscpp::KspSolver::getPetscKSP() {
    return ksp;
}

void petscpp::KspSolver::setup() {
    PetscErrorCode err = KSPSetFromOptions(ksp);
    petscpp::handleError(err);

    err = KSPSetUp(ksp);
    petscpp::handleError(err);
}

void petscpp::KspSolver::solve(Vector &b, Vector &x) {
    PetscErrorCode err = KSPSolve(ksp, b.getPetscVector(), x.getPetscVector());
    petscpp::handleError(err);
}

void petscpp::KspSolver::createPreconditioner() {
	PC pc;
    PetscErrorCode err = KSPGetPC(ksp, &pc);
    petscpp::handleError(err);

    err = PCSetFromOptions(pc);
    petscpp::handleError(err);

    err = PCSetUp(pc);
    petscpp::handleError(err);
}

void petscpp::KspSolver::printDebug() {
	printConvergenceReason();
	printIterationNumber();
}

void petscpp::KspSolver::printConvergenceReason() {
	KSPConvergedReason reason;
	
	PetscErrorCode err = KSPGetConvergedReason(ksp, &reason);
	petscpp::handleError(err);
	
	PetscPrintf(petscpp::worldComm(), "Convergence reason: %d\n", reason);
}

void petscpp::KspSolver::printIterationNumber() {
	PetscInt its;
	KSPGetIterationNumber(ksp, &its);
	PetscPrintf(petscpp::worldComm(), "Iterations: %d\n", its);
}

