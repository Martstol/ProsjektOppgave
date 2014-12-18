// This file is created by students of the HPC-Lab at the Norwegian University of Science and Technology (NTNU)
// and is part of the HPC-Lab Snow Simulator application distributed under the GPL license. 
// Copyright (c) 2006-2013 High Performance Lab at the Norwegian University of Science and Technology (NTNU)
// Department of Computer and Information Science (IDI). All rights reserved.
// See the file README.md for more information.

#include <stdexcept>

#include <petscsys.h>

#include "petscpp.h"
#include "petscppvector.h"
#include "petscpperror.h"

petscpp::Vector::Vector(PetscInt n) {
    PetscErrorCode err = VecCreate(petscpp::worldComm(), &vector);
    petscpp::handleError(err);

    err = VecSetSizes(vector, PETSC_DECIDE, n);
    petscpp::handleError(err);

    err = VecSetFromOptions(vector);
    petscpp::handleError(err);
}

petscpp::Vector::~Vector() {
    PetscErrorCode err = VecDestroy(&vector);
    petscpp::handleError(err);
}

void petscpp::Vector::set(PetscInt i, PetscScalar value) {
    PetscErrorCode err = VecSetValue(vector, i, value, INSERT_VALUES);
    petscpp::handleError(err);
}

void petscpp::Vector::fill(PetscScalar value) {
    PetscErrorCode err = VecSet(vector, value);
    petscpp::handleError(err);
}

void petscpp::Vector::assemble() {
    assemblyBegin();
    assemblyEnd();
}

void petscpp::Vector::assemblyBegin() {
    PetscErrorCode err = VecAssemblyBegin(vector);
    petscpp::handleError(err);
}

void petscpp::Vector::assemblyEnd() {
    PetscErrorCode err = VecAssemblyEnd(vector);
    petscpp::handleError(err);
}

PetscInt petscpp::Vector::size() {
	PetscInt size;
	VecGetSize(vector, &size);
	return size;
}

Vec petscpp::Vector::getPetscVector() {
    return vector;
}

PetscReal petscpp::Vector::getNorm(NormType type) {
    PetscReal val;
    VecNorm(vector, type, &val);
    return val;
}

void petscpp::Vector::print() {
    VecView(vector, PETSC_VIEWER_STDOUT_WORLD);
}

PetscScalar * petscpp::Vector::getArray() {
    PetscScalar * array;
	
    PetscErrorCode err = VecGetArray(vector, &array);
    petscpp::handleError(err);
	
    return array;
}

void petscpp::Vector::restoreArray(PetscScalar ** array) {
    PetscErrorCode err = VecRestoreArray(vector, array);
    petscpp::handleError(err);
}