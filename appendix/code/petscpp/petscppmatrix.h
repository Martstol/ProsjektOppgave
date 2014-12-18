// This file is created by students of the HPC-Lab at the Norwegian University of Science and Technology (NTNU)
// and is part of the HPC-Lab Snow Simulator application distributed under the GPL license. 
// Copyright (c) 2006-2013 High Performance Lab at the Norwegian University of Science and Technology (NTNU)
// Department of Computer and Information Science (IDI). All rights reserved.
// See the file README.md for more information.

#ifndef PETSC_CPP_MATRIX_H
#define PETSC_CPP_MATRIX_H

#include <petscmat.h>

#include "petscppdm.h"

namespace petscpp {
    class Matrix {
    private:
        Mat matrix;
    public:
        Matrix(PetscInt rows, PetscInt columns);
        ~Matrix();
		
        Matrix & operator=(const Matrix&) = delete;
        Matrix(const Matrix&) = delete;
		
        void set(PetscInt row, PetscInt col, PetscScalar value);
		PetscInt getRows();
		PetscInt getColumns();
		
        void assemble();
        void assemblyBegin();
        void assemblyEnd();
		
		void setNullSpace();
		
        Mat getPetscMatrix();
        void print();
    };
}

#endif
