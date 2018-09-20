/* 
 *Lagrangian Reachtube Computations
 *Authors: Md Ariful Islam and Jacek Cyranks
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

/* std::vector<double> -> DVector*/
DVector stdVec_to_dVec(std::vector<double> std_vec, int begin, int end){

    int dim = end - begin;
    DVector d_vec(dim);
    for(int i=begin; i<end; i++)
        d_vec[i] = std_vec[i];
    return d_vec;    
}

/* DVector -> std::vector<double>*/
std::vector<double> dVec_to_stdVec(DVector d_vec){
    int dim = d_vec.dimension();
    std::vector<double> std_vec(dim);
    for(int i=0; i<dim; i++)
        std_vec[i] = d_vec[i];
    return std_vec;    
}

/* CAPD: DMatrix -> IMatrix */
IMatrix d2iMatrix(DMatrix &m){
    const int R = m.numberOfRows();
    const int C = m.numberOfColumns();
    
    IMatrix r(R,C);
    for(int i=0; i < R; i++){
        for(int j=0; j < C; j++){
            r[i][j] = m[i][j];
        }
    }
    return r;
}

DMatrix capd_lMatrix(IMatrix iA){
    int R = iA.numberOfRows();
    int C = iA.numberOfColumns();
    DMatrix r(R,C);
    for(int i=0; i<R; i++){
        for(int j=0; j<C; j++)
            r[i][j] = leftBound(iA[i][j]);
    }
    return r;
}

DMatrix capd_rMatrix(IMatrix iA){
    int R = iA.numberOfRows();
    int C = iA.numberOfColumns();
    DMatrix r(R,C);
    for(int i=0; i<R; i++){
        for(int j=0; j<C; j++)
            r[i][j] = rightBound(iA[i][j]);
    }
    return r;
}

/* Matlab's mag(Y): Interval Vector version*/
IVector mag(IVector ivec, int dim){
    IVector res(dim);
    //int n = ivec.numOfRows();
    for(int i=0; i < dim; i++){
        double lb = abs(leftBound(ivec[i]));
        double rb = abs(rightBound(ivec[i]));
        double maxVal = std::max(lb, rb);   
        res[i] = maxVal;
    }
    return res;
}

DMatrix mag(IMatrix &m){
    int R = m.numberOfRows();
    int C = m.numberOfColumns();
    DMatrix r(R,C);
    for(int i=0; i<R; i++){
        for(int j=0; j<C; j++){
            r[i][j] = std::max(abs(leftBound(m[i][j])),abs(rightBound(m[i][j])));
        }
    }
    return r;
}

/* Matlab's mag(Y): Interval version*/
Interval mag(Interval iV){
    Interval res;
    double lb = abs(leftBound(iV));
    double rb = abs(rightBound(iV));
    double maxVal = std::max(lb, rb);   
    res = maxVal;
    return res;
}
/* Matlab's all(in0(A,B)): checking inclusion of iA in iB componentwise*/
bool  all_in(IVector iA, IVector iB, int dim){
    for(int i=0; i< dim; i++){
        // res = ( b.inf < a.inf ) & ( a.sup < b.sup );
        if(!(leftBound(iB[i]) < leftBound(iA[i]) && rightBound(iA[i]) < rightBound(iB[i])))
            return false;
    }
    return true;
}

/*Matlab's B = flipud(A) method -- non-rigorous: flip upside down or reverse in case of vector */
IVector flipud(IVector vecA, int dim){
    IVector vecB(dim);
    for(int i = 0; i<dim; i++){
        vecB[i] = vecA[dim-1-i];
    }
   return vecB;
}

DMatrix takeMidM(const IMatrix& m){
    const int R = m.numberOfRows();
    const int C = m.numberOfColumns();
    
    DMatrix r(R,C);
    for(int i=0; i < R; i++){
        for(int j=0; j < C; j++){
            r[i][j] = rightBound(mid(m[i][j]));
        }
    }
    return r;
}

DVector takeMidIV(const IVector& m){
    const int dim = m.dimension();
    
    DVector r(dim);
    for(int i=0; i < dim; i++){
        r[i] = 0.5*(leftBound(m[i])+rightBound(m[i]) );
    }
    return r;
}

DMatrix takeRadM(const IMatrix& m){
    const int R = m.numberOfRows();
    const int C = m.numberOfColumns();
    
    DMatrix r(R,C);
    for(int i=0; i < R; i++){
        for(int j=0; j < C; j++){
            r[i][j] = rightBound( 0.5* ( rightBound(m[i][j]) - leftBound(m[i][j]) ) );
        }
    }
    return r;
}

IMatrix getIMatrix(DMatrix lm, DMatrix rm){
    int R = rm.numberOfRows();
    int C = rm.numberOfColumns();

    IMatrix im(R,C);
    for(int i=0; i<R; i++){
        for(int j=0; j<C; j++)
            im[i][j] = interval(lm[i][j], rm[i][j]);
    }
    return im;
}

IMatrix midRad(DMatrix dAc, DMatrix dAd){
    const int R = dAc.numberOfRows();
    const int C = dAc.numberOfColumns();
    IMatrix r(R,C);
    for(int i=0; i < R; i++){
        for(int j=0; j < C; j++){
            r[i][j] = interval(dAc[i][j]-dAd[i][j], dAc[i][j]+dAd[i][j]);
        }
    }
    return r;
}

void printIMatrix(IMatrix& m){

    int R = m.numberOfRows();
    int C = m.numberOfColumns();
    //cout << endl;
    for(int i=0; i<R; i++){
        cout << "\t";
        for(int j=0; j<C; j++){
            cout << leftBound(m[i][j]) << " ";
        }
        cout << endl;
    }
    cout << endl;
    
    for(int i=0; i<R; i++){
        cout << "\t";
        for(int j=0; j<C; j++){
            cout << rightBound(m[i][j]) << " ";
        }
        cout << endl;
    }
    cout << endl;

}

void printDMatrix(DMatrix& m){

    int R = m.numberOfRows();
    int C = m.numberOfColumns();
    //cout << endl;
    for(int i=0; i<R; i++){
        cout << "\t";
        for(int j=0; j<C; j++){
            cout << m[i][j] << " ";
        }
        cout << endl;
    }  

}

void printIVector(IVector& m){
    int dim = m.dimension();
    //cout << endl;
    cout << "\t";
    for(int i=0; i<dim; i++){
        cout << leftBound(m[i]) << " ";
    }  
    cout << endl;
    cout << "\t";
    for(int i=0; i<dim; i++){
        cout << rightBound(m[i]) << " ";
    }  
    cout << endl;

}

void printDVector(DVector& m){
    int dim = m.dimension();
    //cout << endl;
    cout << "\t";
    for(int i=0; i<dim; i++){
        cout << m[i] << " ";
    }  
    cout << endl;
}

void writeIMatrix(IMatrix m, std::ofstream &fout){
    int R = m.numberOfRows();
    int C = m.numberOfColumns();
    //leftbound
    for(int i=0; i<R; i++){
        for(int j=0; j<C; j++){
            fout << leftBound(m[i][j]) << " ";
        }
    }

    //rightbound
    for(int i=0; i<R; i++){
        for(int j=0; j<C; j++){
            fout << rightBound(m[i][j]) << " ";
        }
    }
}

void writeDMatrix(DMatrix m, std::ofstream &fout){
    int R = m.numberOfRows();
    int C = m.numberOfColumns();
    //leftbound
    for(int i=0; i<R; i++){
        for(int j=0; j<C; j++){
            fout << m[i][j] << " ";
        }
    }
}

void writeIVector(IVector& m, std::ofstream &fout){
    int dim = m.dimension();
    for(int i=0; i<dim; i++){
        fout << leftBound(m[i]) << " ";
    }  
    for(int i=0; i<dim; i++){
        fout << rightBound(m[i]) << " ";
    }  
    cout << endl;

}

IMatrix deleteRowCol(IMatrix& m, int idx){
    const int R = m.numberOfRows();
    const int C = m.numberOfColumns();
    
    //cout << endl << "Inside deleteRowCol: " << endl;
    //cout << "iA: " << endl;
    //printIMatrix(m);
    IMatrix r(R-1,C-1);
    //r = m;
    int cCount, rCount;
    rCount = 0;
    for(int i=0; i < R; i++){
        if(i == idx){
            continue;
        }
        cCount = 0;
        for(int j=0; j < C; j++){
            if(j == idx){
                continue;
            }
            r[rCount][cCount] = m[i][j];
            cCount++;
        }
        rCount++;
    }
    //cout << "Idx: " << idx << endl;
    //cout << "new Mat: " << endl;
    //printIMatrix(r);
    return r;
        
}

IMatrix getRowsCols(IMatrix& m, vector<int> idxVec){
    const int R = m.numberOfRows();
    const int C = m.numberOfColumns();
    
    //cout << "Inside getRowsCols: " << endl;
    // cout << "Index vector: " << endl;
    // for(int i: idxVec){
    //     cout << i << " ";
    // }
    //cout << endl;
    int nDim = idxVec.size();
    IMatrix r(nDim,nDim);
    //r = m;
    //cout << "nDIM: " << nDim << endl;
    int cCount, rCount;
    rCount = 0;
    for(int i: idxVec){
        cCount = 0;
        for(int j: idxVec){
            r[rCount][cCount] = m[i][j];
            cCount++;
        }
        rCount++;
    }
    // cout << "new Mat: " << endl;
    // printIMatrix(r);
    return r;
        
}

DMatrix deleteRowCol(DMatrix& m, int idx){
    const int R = m.numberOfRows();
    const int C = m.numberOfColumns();
    
    DMatrix r(R-1,C-1);
    //r = m;
    int cCount, rCount;
    rCount = 0;
    for(int i=0; i < R; i++){
        if(i == idx){
            continue;
        }
        cCount = 0;
        for(int j=0; j < C; j++){
            if(j == idx){
                continue;
            }
            r[rCount][cCount] = m[i][j];
            cCount++;
        }
        rCount++;
    }
    return r;
        
}