// This file is created by students of the HPC-Lab at the Norwegian University of Science and Technology (NTNU)
// and is part of the HPC-Lab Snow Simulator application distributed under the GPL license. 
// Copyright (c) 2006-2013 High Performance Lab at the Norwegian University of Science and Technology (NTNU)
// Department of Computer and Information Science (IDI). All rights reserved.
// See the file README.md for more information.

#include <stdexcept>
#include <iostream>

#include <petscsys.h>

#include "petscpp.h"
#include "petscppmatrix.h"
#include "petscpperror.h"

petscpp::Matrix::Matrix(PetscInt rows, PetscInt columns) {
    PetscErrorCode err = MatCreate(petscpp::worldComm(), &matrix);
    petscpp::handleError(err);

    err = MatSetSizes(matrix, PETSC_DECIDE, PETSC_DECIDE, rows, columns);
    petscpp::handleError(err);
	
    MatSetFromOptions(matrix);
	
	MatSetUp(matrix);
}

petscpp::Matrix::~Matrix() {
    PetscErrorCode err = MatDestroy(&matrix);
    petscpp::handleError(err);
}

void petscpp::Matrix::set(PetscInt row, PetscInt col, PetscScalar value) {
    PetscErrorCode err = MatSetValue(matrix, row, col, value, INSERT_VALUES);
    petscpp::handleError(err);
}

PetscInt petscpp::Matrix::getColumns() {
	PetscInt columns;
	MatGetSize(matrix, nullptr, &columns);
	return columns;
}

PetscInt petscpp::Matrix::getRows() {
	PetscInt rows;
	MatGetSize(matrix, &rows, nullptr);
	return rows;
}

void petscpp::Matrix::assemble() {
    assemblyBegin();
    assemblyEnd();
}

void petscpp::Matrix::assemblyBegin() {
    PetscErrorCode err = MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
    petscpp::handleError(err);
}

void petscpp::Matrix::assemblyEnd() {
    PetscErrorCode err = MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
    petscpp::handleError(err);
}

Mat petscpp::Matrix::getPetscMatrix() {
    return matrix;
}

void petscpp::Matrix::print() {
    MatView(matrix, PETSC_VIEWER_STDOUT_WORLD);
}

void petscpp::Matrix::setNullSpace() {
	MatNullSpace nullspace;
	
	PetscErrorCode err = MatNullSpaceCreate(petscpp::worldComm(), PETSC_TRUE, 0, 0, &nullspace);
	petscpp::handleError(err);
	
	MatSetNullSpace(matrix, nullspace);
	
	MatNullSpaceDestroy(&nullspace);
}
