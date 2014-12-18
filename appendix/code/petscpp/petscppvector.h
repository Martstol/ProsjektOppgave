// This file is created by students of the HPC-Lab at the Norwegian University of Science and Technology (NTNU)
// and is part of the HPC-Lab Snow Simulator application distributed under the GPL license. 
// Copyright (c) 2006-2013 High Performance Lab at the Norwegian University of Science and Technology (NTNU)
// Department of Computer and Information Science (IDI). All rights reserved.
// See the file README.md for more information.

#ifndef PETSC_CPP_VECTOR_H
#define PETSC_CPP_VECTOR_H

#include <petsc.h>
#include <petscvec.h>

namespace petscpp {

    namespace norm {
        NormType const inf = NORM_INFINITY;
        NormType const norm1 = NORM_1;
        NormType const norm2 = NORM_2;
    }

    class Vector {
    private:
        Vec vector;

    public:
        Vector(PetscInt n);
        ~Vector();
		
        Vector & operator=(const Vector&) = delete;
        Vector(const Vector&) = delete;
		
        Vector( Vector&& mv ) : vector( mv.vector ) { mv.vector = nullptr; }
        Vector & operator=(Vector&& mv) {vector = mv.vector; mv.vector = nullptr;}
        
        void set(PetscInt i, PetscScalar value);
        void fill(PetscScalar value);
		PetscInt size();
        PetscReal getNorm(NormType type);
		
        void assemble();
        void assemblyBegin();
        void assemblyEnd();
		
        Vec getPetscVector();
		
        void print();
		
        PetscScalar * getArray();
        void restoreArray(PetscScalar ** array);
    };

}

#endif
