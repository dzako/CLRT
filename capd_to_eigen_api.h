/* 
 *Lagrangian Reachtube Computations
 *Authors: Md Ariful Islam and Jacek Cyranks
 *Contact: ariful.islam@ttu.edu
 */

/* Standard C++ Library */
#include <unistd.h>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
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

/*Define namespace*/
using namespace Eigen;
using namespace capd;
using namespace std;

////////////// BEGIN: API between CAPD and EIGEN ////////////////////////////////
/* CAPD to Eigen */
MatrixXd capd_to_eigen(DMatrix capd_mat, int dim) {
    MatrixXd eig_mat(dim, dim);
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            eig_mat(i,j) = capd_mat[i][j];
        }
    }
    return eig_mat;
}

/* Eigen to CAPD */
// MatrixXd -> IMatrix
DMatrix eigen_to_capd(MatrixXd eig_mat, int dim) {
    DMatrix capd_mat(dim, dim);
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            capd_mat[i][j] = eig_mat(i,j);
        }
    }
    return capd_mat;
}
// VectorXd -> IVector
IVector eigen_to_capd(VectorXd eig_vec, int dim){
    IVector capd_vec(dim);
    for(int i=0; i<dim; i++){
        capd_vec[i] = eig_vec(i);
    }
    return capd_vec;
}
