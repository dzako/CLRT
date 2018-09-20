/* 
 *Lagrangian Reachtube Computations
 *Authors: Md Ariful Islam
 *Contact: mdarifui@andrew.cmu.edu
 */

/* Standard C++ library*/
#include <unistd.h>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <iostream>
#include <vector>
#include <functional>
#include <numeric>
#include <set>
#include <iterator>
/* Eigen Library*/
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>

/* CAPD Library*/
#include "capd/capdlib.h"
#include <capd/vectalg/Vector_inline.h>
#include <capd/matrixAlgorithms/floatMatrixAlgorithms.hpp>

/* Library defined by myself */
#ifndef INTEIG_H
#define INTEIG_H
    #include "computeIntEig.h"
#endif

/*Define namespace*/
using namespace Eigen;
using namespace capd;
using namespace std;

// compute Jordan-Wielandt Matrix
// (0, A^T, A, 0)
IMatrix computeJW(IMatrix iA){
	int dim = iA.numberOfRows();
	IMatrix iB(2*dim, 2*dim);
	// compute 0's in the diagonal blocks: (0,0) and (1,1) blocks
	for(int i=0; i<dim; i++){
		for(int j=0; j<dim; j++){
			iB[i][j] = 0.0;
			iB[dim+i][dim+j] = 0;

		}
	}

	// write (0,1)-block: A^T
	IMatrix iAtrans = transpose(iA);
	for(int i=0; i<dim; i++){
		for(int j=0; j<dim; j++){
			iB[i][dim+j] = iAtrans[i][j];
		}
	}

	// write (1,0)-block: A
	for(int i=0; i<dim; i++){
		for(int j=0; j<dim; j++){
			iB[dim+i][j] = iA[i][j];
		}
	}

	return iB;

}

double computeSVD(IMatrix iA){
	int dim = iA.numberOfRows();

	IMatrix iB = computeJW(iA);
	double svd = computeIntEig(iB);
	
	return svd;

}
