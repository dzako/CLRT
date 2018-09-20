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


/*function computes double change of coordinates matrix C (inverse)*/
DMatrix computeNorm(const DMatrix& F, const int DIM){
    DVector rVal(DIM), iVal(DIM);
    DMatrix rVec(DIM, DIM), iVec(DIM, DIM);
    
    /* Compute Eigenvalue and Eigenvector in CAPD */
    
    capd::alglib::computeEigenvaluesAndEigenvectors(F, rVal, iVal, rVec, iVec);
    if(_DEBUG_) {
        cout << "rVec = [" ;
        for(int i=0; i<DIM; i++){
            for(int j=0; j<DIM; j++){
                cout << rVec[i][j] << " ";
            }
            if(i < DIM-1)
                cout << ";";
        }
        cout << "]" << endl;
        cout << "iVec = [";
        for(int i=0; i<DIM; i++){
            for(int j=0; j<DIM; j++){
                cout << iVec[i][j] << " ";
            }
            if(i < DIM-1)
                cout << ";";
        }
        cout << "]" << endl;
    }

    DMatrix r(DIM, DIM);
    int i, j;
    for(i = 0; i < DIM - 1; i += 2){
        for(j = 0; j < DIM; j++) {
            if(iVec[i][j] != 0 || iVec[i + 1][j] != 0) { //block with eigenvectors in form v+iw, v-iw is detected
            //instead of block [v+iw,v-iw] we take block [v,w]
            r[i][j] = rVec[i][j]; //v_1
            r[i + 1][j] = rVec[i + 1][j]; //v_2
            r[i][j + 1] = iVec[i][j]; //w_1
            r[i + 1][j + 1] = iVec[i + 1][j]; //w_2
            j++;
          } else { //we rewrite real eigenvectors
            r[i][j] = rVec[i][j];
            r[i + 1][j] = rVec[i + 1][j];
          }
        }
    }
    //in case of odd dimension we check last row
    if(DIM % 2 == 1){
        for(j = 0; j < DIM; j++) {
          if(iVec[i][j] != 0){
            r[i][j + 1] = iVec[i][j];
            r[i][j] = rVec[i][j];
            j++;
          } else {
            r[i][j]=rVec[i][j];
          }
        }
    }
    return capd::matrixAlgorithms::gaussInverseMatrix(r);
}

/*function computes double change of coordinates matrix C (inverse)*/
DMatrix computeNorm_eigen(const DMatrix& F, const int DIM){
    
    /* Compute Eigenvalue and Eigenvector in CAPD */
    MatrixXd A = capd_to_eigen(F, DIM);
    EigenSolver<MatrixXd> es(A);

    DMatrix r(DIM, DIM);
    int i, j;
    for(i = 0; i < DIM - 1; i += 2){
        VectorXd iv1 = es.eigenvectors().col(i).imag();
        VectorXd iv2 = es.eigenvectors().col(i+1).imag();

        VectorXd rv1 = es.eigenvectors().col(i).real();
        VectorXd rv2 = es.eigenvectors().col(i+1).real();

        for(j = 0; j < DIM; j++) {

            if(iv1[j] != 0 || iv2[j] != 0) { //block with eigenvectors in form v+iw, v-iw is detected
            //instead of block [v+iw,v-iw] we take block [v,w]
            r[i][j] = rv1[j]; //v_1
            r[i + 1][j] = rv2[j]; //v_2
            r[i][j + 1] = iv1[j]; //w_1
            r[i + 1][j + 1] = iv2[j]; //w_2
            j++;
          } else { //we rewrite real eigenvectors
            r[i][j] = rv1[j];
            r[i + 1][j] = rv2[j];
          }
        }
    }
    //in case of odd dimension we check last row
    if(DIM % 2 == 1){
        VectorXd iv1 = es.eigenvectors().col(i).imag();
        VectorXd rv1 = es.eigenvectors().col(i).real();
        for(j = 0; j < DIM; j++) {
          if(iv1[j] != 0){
            r[i][j + 1] = iv1[j];
            r[i][j] = rv1[j];
            j++;
          } else {
            r[i][j]=rv1[j];
          }
        }
    }
    return capd::matrixAlgorithms::gaussInverseMatrix(r);
}


DMatrix chol(DMatrix dA, int DIM) {
    MatrixXd A(DIM,DIM);
    for(int i=0; i < DIM; i++){
        for(int j=0; j<DIM; j++){
            A(i,j) = dA[i][j]; 
        }
    }
    
    LLT<MatrixXd> lltOfA(A); // compute the Cholesky decomposition of A
    MatrixXd L = lltOfA.matrixL(); // retrieve factor L  in the decomposition

    DMatrix dC(DIM,DIM);
    for(int i=0; i < DIM; i++){
        for(int j=0; j<DIM; j++){
            dC[i][j] = L(i,j); 
        }
    }

    return transpose(dC);
}

DMatrix sqrt(DMatrix dF, int DIM) {
    MatrixXd F(DIM,DIM);
    //cout << "dF: " << dF << endl;
    for(int i=0; i < DIM; i++){
        for(int j=0; j<DIM; j++){
            F(i,j) = dF[i][j]; 
        }
    }
    MatrixXd A = F.transpose()*F;
    //cout << "Here is a random positive-definite matrix, A:" << endl << A << endl << endl;
    SelfAdjointEigenSolver<MatrixXd> es(A);
    MatrixXd sqrtA = es.operatorSqrt();
    //cout << "The square root of A is: " << endl << sqrtA << endl;
    //cout << "If we square this, we get: " << endl << sqrtA*sqrtA << endl;

    DMatrix dsqrtC(DIM,DIM);
    for(int i=0; i < DIM; i++){
        for(int j=0; j<DIM; j++){
            dsqrtC[i][j] = sqrtA(i,j); 
        }
    }

    return dsqrtC;
}