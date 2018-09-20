/* 
 *Lagrangian Reachtube Computations
 *Authors: Md Ariful Islam and Jacek Cyranka
 *Contact: ariful.islam@ttu.edu
 */

/* Standard C++ library*/
#include <unistd.h>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <functional>
#include <numeric>
#include <set>
#include <iterator>
#include <stdio.h>

/* Eigen Library*/
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

/* CAPD Library*/
#include "capd/capdlib.h"
#include <capd/vectalg/Vector_inline.h>
#include <capd/matrixAlgorithms/floatMatrixAlgorithms.hpp>
#include "capd/intervals/lib.h"
/*GiNaC Library*/
#include <ginac/ginac.h>

/* NLOPT Library*/
#include<nlopt.hpp>

/* Library defined by myself */
#ifndef INTEIG_H
#define INTEIG_H
    #include "computeIntEig.h"
#endif

#ifndef NORM_H
#define NORM_H
    #include "computeNorm.h"
#endif

#ifndef SVD_H
#define SVD_H
    #include "computeSVD.h"
#endif

/* GLOBAL VARIABLES */
double TILDE_FACTOR = 35;
int DREAL_MAX_ITER = 1000;
double DREAL_FACTOR = 1.05;
double DREAL_DELTA = 1e-6;
double DUBINS_V = 1;
/*Define namespace*/
using namespace Eigen;
using namespace capd;
using namespace std;
using namespace GiNaC;


/************************* Symbolic Computation ******************************/
// define symbolic vector: ex type vector
std::vector<ex> sym_vec(string varName, int dim){

    std::vector<ex> symVec(dim);
    for(int i =0; i < dim; i++){ // naming start from 1, ex: x1, x2, ..., xdim
        symbol x(varName + std::to_string(i+1));
        ex poly = pow(x,1); // x^1 = x
        symVec[i] = poly;
    }
    return symVec;
}

// Matrix-vector multiplication: outVec = M*xVec
vector<ex> sym_mult(DMatrix &M, std::vector<ex> xVec){

    int C = M.numberOfColumns();
    int R = M.numberOfRows();
    int dim = xVec.size();

    assert(dim == C);

    std::vector<ex> outVec(dim);

    for(int i=0; i<R; i++){
        ex poly;
        for(int j=0; j<C; j++){
            poly += M[i][j]*xVec[j];
        }
        outVec[i] = poly;
    }

    return outVec;
}

// Inner product of two vector: xVec'*yVec
ex sym_inner(std::vector<ex> xVec, std::vector<ex> yVec){

    int dim_x = xVec.size();
    int dim_y = yVec.size();
    assert(dim_x == dim_y);
    ex out;
    for(int i=0; i<dim_x; i++){
        out += xVec[i]*yVec[i];
    }

    return out;
}

// Return symbolic expression: xVec'*M*xVec - rad^2
// M and rad: numeric
// xVec: symbolic
ex sym_ellipse(std::vector<ex> xVec, DMatrix M, double rad){

    std::vector<ex> mxVec = sym_mult(M, xVec); // M*xVec
    ex ellipse = sym_inner(mxVec, xVec) - rad*rad; // xVec'*(M*xVec);

    return ellipse;
}
/************************* End of Symbolic Computation ****************************/


/******************* Definition of Model's dynamic in symbolic expression *******************/

/* Brusselator dynamics: f(x)=[f1(x); f2(x)], where x = [x1;x2] */
std::vector<ex> sym_bruss(std::vector<ex> xvec, ex tvar){ // tvar is used for time-varying system only
    int dim = xvec.size(); // 2 here
    assert(dim == 2);
    std::vector<ex> fdyn(dim);

    fdyn[0] = 1 + xvec[0]*xvec[0]*xvec[1] - 2.5*xvec[0]; //1+x*x*y-2.5*x
    fdyn[1] = 1.5*xvec[0] - xvec[0]*xvec[0]*xvec[1]; // 1.5*x-x*x*y
    return fdyn;

}

std::vector<ex> sym_ms(std::vector<ex> xvec, ex tvar){ // tvar is used for time-varying system only
    int dim = xvec.size(); // 2 here
    assert(dim == 2);
    ex x1 = xvec[0];
    ex x2 = xvec[1];

    std::vector<ex> fdyn(dim);
 
    fdyn[0] = x2*x1*x1* (1 - x1)/0.3 - x1/6; //x2*x1^2* (1 - x1)/0.3 - x1/6
    fdyn[1] = (0.5*(1 + (1-exp(-100*x1+10))/(1-exp(-100*x1+10))))*x2/150 + 
              (1 - (0.5*(1 + (1-exp(-100*x1+10))/(1-exp(-100*x1+10)))))*(1 - x2)/20;

    return fdyn;

}

std::vector<ex> sym_vdp(std::vector<ex> xvec, ex tvar){ // tvar is used for time-varying system only
    int dim = xvec.size(); // 2 here
    assert(dim == 2);
    ex x1 = xvec[0];
    ex x2 = xvec[1];

    std::vector<ex> fdyn(dim);
    //vectorString = "var:x,y;fun: -y, (x^2-1)*y+x;";
    fdyn[0] = -x2;
    fdyn[1] = (x1*x1 - 1)*x2 + x1;

    return fdyn;

}

/* Robot Arm dynamics: f(x)=[f1(x); f2(x); f3(x); f4(x)], where x = [x1;x2;x3;x4] */
std::vector<ex> sym_robot(std::vector<ex> xvec, ex tvar){ // tvar is used for time-varying system only
    
    //cout << "Inside sym_robot function ..." << endl;
    int dim = xvec.size(); 
    assert(dim == 4); // 4-dimensional systems
    std::vector<ex> fdyn(dim);

    /* Define Parameters*/
    int m = 1;
    int l = 3;
    int kp1 = 2;
    int kp2 = 1;
    int kd1 = 2;
    int kd2 = 1;
    int u1 = 2;
    int u2 = 1;

    /* State variables */
    ex x1 = xvec[0];
    ex x2 = xvec[1];
    ex x3 = xvec[2];
    ex x4 = xvec[3];

    /* Define dynamics */
    fdyn[0] = x3; //1+x*x*y-2.5*x
    fdyn[1] = x4; // 1.5*x-x*x*y
    fdyn[2] = (-2*m*x2*x3*x4-kp1*x1-kd1*x3)/(m*x2*x2+l/3)+(kp1*kp1)/(m*x2*x2+l/3);
    fdyn[3] = x2*x3*x3-kp2*x2/m-kd2*x4/m+kp2*kp2/m;

    return fdyn;
}

/* Twelve dynamics: dimension = 12 */
std::vector<ex> sym_twelve(std::vector<ex> xvec, ex tvar){ // tvar is used for time-varying system only
    
    //cout << "Inside sym_robot function ..." << endl;
    int dim = xvec.size(); 
    assert(dim == 12); // 12-dimensional systems
    std::vector<ex> fdyn(dim);

    /* State variables */
    ex x1 = xvec[0];
    ex x2 = xvec[1];
    ex x3 = xvec[2];
    ex x4 = xvec[3];
    ex x5 = xvec[4];
    ex x6 = xvec[5];
    ex x7 = xvec[6];
    ex x8 = xvec[7];
    ex x9 = xvec[8];
    ex x10 = xvec[9];
    ex x11 = xvec[10];
    ex x12 = xvec[11];
 
    /* Define dynamics */
    fdyn[0] = -1.47*x1 - 0.574*x2 - 0.583*x3 + 0.942*x4 + 0.714*x5 \
            +1.35*x6 + 0.861*x7 + 1.51*x9 + 1.10*x11 + 1.14*x12 + x2*x6*x9;
    fdyn[1] = -0.76*x2 + 1.54*x3 - 1.63*x4 - 0.819*x5 - 0.925*x6 - 1.06*x7 - \
            0.535*x10 + 0.681*x11 + x7*x10*x11; 
    fdyn[2] = 1.77*x1 - 0.983*x2 - 0.766*x3 - 0.930*x6 - 1.13*x9 - 0.998*x10 \
            + 1.84*x11 + 0.772*x12 + x3*x3*x3;
    fdyn[3] = -1.37*x1 + 1.31*x2 + 0.790*x3 - 0.572*x4 - 1.07*x5 - 0.783*x6 - \
            0.938*x8 + 2.03*x9 - 0.857*x11 + x4*x8*x11;
    fdyn[4] = -1.04*x1 + 0.755*x2 + 1.12*x4 - 0.538*x5 - 0.563*x8 + 1.40*x9 - \
            1.32*x11 + x1*x2*x5;
    fdyn[5] = 1.26*x2 + 1.11*x4 - 0.885*x6 - 0.935*x8 + 1.05*x9 - 1.30*x11 + x2*x4*x11;
    fdyn[6] = - 0.845*x1 - 1.67*x7 + 1.72*x9 + 0.839*x10 - 1.25*x11 - 0.600*x12 + x1*x7*x11;
    fdyn[7] = 0.583*x1 + 1.20*x4 + 0.816*x5 + 0.599*x6 - 0.735*x9 + x1*x1*x1;
    fdyn[8] = - 1.25*x1 + 1.04*x2 + 0.990*x3 - 1.44*x4 - 1.44*x5 - 1.29*x6 - 1.20*x9 \
            - 1.16*x10 - 3.15*x11 - 0.998*x12 + x1*x6*x10;
    fdyn[9] = 1.43*x2 + 0.835*x3 + 0.545*x7 - 0.890*x8 - 0.973*x10 + x2*x3*x8;
    fdyn[10] = - 2.29*x3 + 1.41*x4 + 1.43*x5 + 0.512*x6 + 2.00*x7 - 0.801*x8 + \
            2.46*x9 - 0.659*x10 - 0.833*x11 + 1.29*x12 + x6*x6*x6;
    fdyn[11] = - 0.850*x1 - 0.919*x3 + 0.515*x7 + 1.20*x9 - 1.25*x11 + x1*x9*x11;

    return fdyn;
}

/* Tweleve dynamics: dimension = 7 */
std::vector<ex> sym_biology(std::vector<ex> xvec, ex tvar){ // tvar is used for time-varying system only
    
    //cout << "Inside sym_robot function ..." << endl;
    int dim = xvec.size(); 
    assert(dim == 7); // 7-dimensional systems
    std::vector<ex> fdyn(dim);

    /* State variables */
    ex x1 = xvec[0];
    ex x2 = xvec[1];
    ex x3 = xvec[2];
    ex x4 = xvec[3];
    ex x5 = xvec[4];
    ex x6 = xvec[5];
    ex x7 = xvec[6];

    /* Define dynamics */
    fdyn[0] = -0.4*x1 + 50*x3*x4;
    fdyn[1] = 0.4*x1 - x2; 
    fdyn[2] = x2 - 50*x3*x4;
    fdyn[3] = 50*x5*x6 - 50*x3*x4;
    fdyn[4] = -50*x5*x6 + 50*x3*x4;
    fdyn[5] = 0.5*x7 - 50*x5*x6;
    fdyn[6] = -0.5*x7 + 50*x5*x6;

    return fdyn;
}

/* Dubins Car dynamics: f(x)=[f1(x); f2(x); f3(x)], where x = [x1;x2;x3] */
std::vector<ex> sym_dubins(std::vector<ex> xvec, ex tvar){ // tvar is used for time-varying system only
    
    //cout << "Inside sym_robot function ..." << endl;
    int dim = xvec.size(); // 2 here
    assert(dim == 3); // 3-dimensional systems
    std::vector<ex> fdyn(dim);

    /* Define Parameters*/
    double v = 1;
    
    /* State variables */
    ex x1 = xvec[0];
    ex x2 = xvec[1];
    ex x3 = xvec[2];
    //ex x4 = xvec[3];

    /* Define dynamics */
    fdyn[0] = v*cos(x3); //1+x*x*y-2.5*x
    fdyn[1] = v*sin(x3); // 1.5*x-x*x*y
    fdyn[2] = x1*sin(tvar);
    //fdyn[3] = 0;

    return fdyn;
}

/* Dubins Car dynamics: f(x)=[f1(x); f2(x); f3(x)], where x = [x1;x2;x3] */
std::vector<ex> sym_dubins_2(std::vector<ex> xvec, ex tvar){ // tvar is used for time-varying system only
    
    //cout << "Inside sym_robot function ..." << endl;
    int dim = xvec.size(); // 2 here
    assert(dim == 3); // 3-dimensional systems
    std::vector<ex> fdyn(dim);

    /* Define Parameters*/
    double v = 1;
    
    /* State variables */
    ex x1 = xvec[0];
    ex x2 = xvec[1];
    ex x3 = xvec[2];
    //ex x4 = xvec[3];

    /* Define dynamics */
    fdyn[0] = v*cos(x3); //1+x*x*y-2.5*x
    fdyn[1] = v*sin(x3); // 1.5*x-x*x*y
    fdyn[2] = x1*x2*sin(tvar);
    //fdyn[3] = 0;

    return fdyn;
}

std::vector<ex> sym_dubins_ti(std::vector<ex> xvec, ex tvar){ // tvar is used for time-varying system only
    
    //cout << "Inside sym_robot function ..." << endl;
    int dim = xvec.size(); // 2 here
    assert(dim == 3); // 3-dimensional systems
    std::vector<ex> fdyn(dim);

    /* Define Parameters*/
    double v = 1;
    
    /* State variables */
    ex x1 = xvec[0];
    ex x2 = xvec[1];
    ex x3 = xvec[2];
    //ex x4 = xvec[3];

    /* Define dynamics */
    fdyn[0] = v*cos(x3); //1+x*x*y-2.5*x
    fdyn[1] = v*sin(x3); // 1.5*x-x*x*y
    fdyn[2] = 1;
    //fdyn[3] = 0;

    return fdyn;
}
/******************* End of Definition of Model's dynamic in symbolic expression *******************/

/********************* Local Optimization  **************************/
struct constraint_data {
    double rad;
    DVector cx;
    DMatrix M;
    constraint_data(DVector c, double r, DMatrix m){
        rad = r;
        cx = c; 
        M = m;
    }
    void print(){
        cout << "cx: " << endl;
        printDVector(cx);
        cout << endl;
        cout << "rad: " << rad << endl;
        cout << "M: " << endl;
        printDMatrix(M);
        cout << endl;
    }
};


double cons_func(const std::vector<double> &x, std::vector<double> &grad, void *data){   
    int dim = x.size();
    constraint_data* consData = reinterpret_cast<constraint_data*>(data); 
    DVector cx = consData->cx;
    double rad = consData-> rad;
    DMatrix M = consData->M;
    DVector d_x = stdVec_to_dVec(x, 0, dim-1); // from index 0 to dim-1: time is excluded 
    double val = (d_x - cx)*(M*(d_x - cx));
    return val;
}

// objective function for brusselator
double obj_bruss(const std::vector <double> &x, std::vector<double> &grad, void *data){
    DMatrix* pM = reinterpret_cast<DMatrix*>(data); 
    DMatrix M = *pM;
    int dim = x.size();
    
    // compute fdyn for brusselator system
    DVector fdyn(dim-1); // remove time dimension: x[2]
    fdyn[0] = 1 + x[0]*x[0]*x[1] - 2.5*x[0]; // 1+x1*x1*x2 - 2.5*x1
    fdyn[1] = 1.5*x[0] - x[0]*x[0]*x[1];  // 1.5*x1 - x1*x1*x2
    double val = fdyn*(M*fdyn);
    return val;

}

// objective function for MS model
double obj_ms(const std::vector <double> &x, std::vector<double> &grad, void *data){
    DMatrix* pM = reinterpret_cast<DMatrix*>(data); 
    DMatrix M = *pM;
    int dim = x.size();
    
    double x1 = x[0];
    double x2 = x[1];

    // compute fdyn for brusselator system
    DVector fdyn(dim-1); // remove time dimension: x[2]

    fdyn[0] = x2*x1*x1* (1 - x1)/0.3 - x1/6; //x2*x1^2* (1 - x1)/0.3 - x1/6
    fdyn[1] = (0.5*(1 + (1-exp(-100*x1+10))/(1-exp(-100*x1+10))))*x2/150 + 
              (1 - (0.5*(1 + (1-exp(-100*x1+10))/(1-exp(-100*x1+10)))))*(1 - x2)/20;
    double val = fdyn*(M*fdyn);
    return val;

}

// objective function for MS model
double obj_vdp(const std::vector <double> &x, std::vector<double> &grad, void *data){
    DMatrix* pM = reinterpret_cast<DMatrix*>(data); 
    DMatrix M = *pM;
    int dim = x.size();
    
    double x1 = x[0];
    double x2 = x[1];

    // compute fdyn for brusselator system
    DVector fdyn(dim-1); // remove time dimension: x[2]

    fdyn[0] = -x2;
    fdyn[1] = (x1*x1 - 1)*x2 + x1;
    double val = fdyn*(M*fdyn);
    return val;

}

double obj_robot(const std::vector<double> &xvec, std::vector<double> &grad, void *data){
    DMatrix* pM = reinterpret_cast<DMatrix*>(data); 
    DMatrix M = *pM;
    int dim = xvec.size();

    // compute fdyn for brusselator system
    DVector fdyn(dim-1); // remove time dimension: x[4]

    /* Define Parameters*/
    int m = 1;
    int l = 3;
    int kp1 = 2;
    int kp2 = 1;
    int kd1 = 2;
    int kd2 = 1;
    int u1 = 2;
    int u2 = 1;

    /* State variables */
    double x1 = xvec[0];
    double x2 = xvec[1];
    double x3 = xvec[2];
    double x4 = xvec[3];

    /* Define dynamics */
    fdyn[0] = x3; //1+x*x*y-2.5*x
    fdyn[1] = x4; // 1.5*x-x*x*y
    fdyn[2] = (-2*m*x2*x3*x4-kp1*x1-kd1*x3)/(m*x2*x2+l/3)+(kp1*kp1)/(m*x2*x2+l/3);
    fdyn[3] = x2*x3*x3-kp2*x2/m-kd2*x4/m+kp2*kp2/m;

    double val = fdyn*(M*fdyn);
    //cout << "Local obj fun: fval = " << val << endl;
    return val;
}

double obj_dubins(const std::vector<double> &xvec, std::vector<double> &grad, void *data){
    DMatrix* pM = reinterpret_cast<DMatrix*>(data); 
    DMatrix M = *pM;
    int dim = xvec.size();

    // compute fdyn for brusselator system
    DVector fdyn(dim-1); // remove time dimension: x[4]

    /* Define Parameters*/
    double v = 1;

    /* State variables */
    double x1 = xvec[0];
    double x2 = xvec[1];
    double x3 = xvec[2];
    double t = xvec[3];
    /* Define dynamics */
    fdyn[0] = v*std::cos(x3); //1+x*x*y-2.5*x
    fdyn[1] = v*std::sin(x3); // 1.5*x-x*x*y
    fdyn[2] = x1*std::sin(t);

    double val = fdyn*(M*fdyn);
    return val;
}

double obj_dubins_2(const std::vector<double> &xvec, std::vector<double> &grad, void *data){
    DMatrix* pM = reinterpret_cast<DMatrix*>(data); 
    DMatrix M = *pM;
    int dim = xvec.size();

    // compute fdyn for brusselator system
    DVector fdyn(dim-1); // remove time dimension: x[4]

    /* Define Parameters*/
    double v = 1;

    /* State variables */
    double x1 = xvec[0];
    double x2 = xvec[1];
    double x3 = xvec[2];
    double t = xvec[3];
    /* Define dynamics */
    fdyn[0] = v*std::cos(x3); //1+x*x*y-2.5*x
    fdyn[1] = v*std::sin(x3); // 1.5*x-x*x*y
    fdyn[2] = x1*x2*std::sin(t);

    double val = fdyn*(M*fdyn);
    return val;
}

double obj_dubins_ti(const std::vector<double> &xvec, std::vector<double> &grad, void *data){
    DMatrix* pM = reinterpret_cast<DMatrix*>(data); 
    DMatrix M = *pM;
    int dim = xvec.size();

    // compute fdyn for brusselator system
    DVector fdyn(dim-1); // remove time dimension: x[4]

    /* Define Parameters*/
    double v = 1;

    /* State variables */
    double x1 = xvec[0];
    double x2 = xvec[1];
    double x3 = xvec[2];
    double t = xvec[3];
    /* Define dynamics */
    fdyn[0] = v*std::cos(x3); //1+x*x*y-2.5*x
    fdyn[1] = v*std::sin(x3); // 1.5*x-x*x*y
    fdyn[2] = 1;

    double val = fdyn*(M*fdyn);
    return val;
}

double obj_biology(const std::vector<double> &xvec, std::vector<double> &grad, void *data){
    DMatrix* pM = reinterpret_cast<DMatrix*>(data); 
    DMatrix M = *pM;
    int dim = xvec.size();

    // compute fdyn for brusselator system
    DVector fdyn(dim-1); // remove time dimension: x[7]

    /* State variables */
    double x1 = xvec[0];
    double x2 = xvec[1];
    double x3 = xvec[2];
    double x4 = xvec[3];
    double x5 = xvec[4];
    double x6 = xvec[5];
    double x7 = xvec[6];

    double t = xvec[7];
    /* Define dynamics */
    fdyn[0] = -0.4*x1 + 50*x3*x4;
    fdyn[1] = 0.4*x1 - x2; 
    fdyn[2] = x2 - 50*x3*x4;
    fdyn[3] = 50*x5*x6 - 50*x3*x4;
    fdyn[4] = -50*x5*x6 + 50*x3*x4;
    fdyn[5] = 0.5*x7 - 50*x5*x6;
    fdyn[6] = -0.5*x7 + 50*x5*x6;

    double val = fdyn*(M*fdyn);
    return val;
}

double obj_twelve(const std::vector<double> &xvec, std::vector<double> &grad, void *data){
    DMatrix* pM = reinterpret_cast<DMatrix*>(data); 
    DMatrix M = *pM;
    int dim = xvec.size();

    // compute fdyn for brusselator system
    DVector fdyn(dim-1); // remove time dimension: x[13]

    /* State variables */
    double x1 = xvec[0];
    double x2 = xvec[1];
    double x3 = xvec[2];
    double x4 = xvec[3];
    double x5 = xvec[4];
    double x6 = xvec[5];
    double x7 = xvec[6];
    double x8 = xvec[7];
    double x9 = xvec[8];
    double x10 = xvec[9];
    double x11 = xvec[10];
    double x12 = xvec[11];

    double t = xvec[12];
    /* Define dynamics */
    fdyn[0] = -1.47*x1 - 0.574*x2 - 0.583*x3 + 0.942*x4 + 0.714*x5 \
            +1.35*x6 + 0.861*x7 + 1.51*x9 + 1.10*x11 + 1.14*x12 + x2*x6*x9;
    fdyn[1] = -0.76*x2 + 1.54*x3 - 1.63*x4 - 0.819*x5 - 0.925*x6 - 1.06*x7 - \
            0.535*x10 + 0.681*x11 + x7*x10*x11; 
    fdyn[2] = 1.77*x1 - 0.983*x2 - 0.766*x3 - 0.930*x6 - 1.13*x9 - 0.998*x10 \
            + 1.84*x11 + 0.772*x12 + x3*x3*x3;
    fdyn[3] = -1.37*x1 + 1.31*x2 + 0.790*x3 - 0.572*x4 - 1.07*x5 - 0.783*x6 - \
            0.938*x8 + 2.03*x9 - 0.857*x11 + x4*x8*x11;
    fdyn[4] = -1.04*x1 + 0.755*x2 + 1.12*x4 - 0.538*x5 - 0.563*x8 + 1.40*x9 - \
            1.32*x11 + x1*x2*x5;
    fdyn[5] = 1.26*x2 + 1.11*x4 - 0.885*x6 - 0.935*x8 + 1.05*x9 - 1.30*x11 + x2*x4*x11;
    fdyn[6] = - 0.845*x1 - 1.67*x7 + 1.72*x9 + 0.839*x10 - 1.25*x11 - 0.600*x12 + x1*x7*x11;
    fdyn[7] = 0.583*x1 + 1.20*x4 + 0.816*x5 + 0.599*x6 - 0.735*x9 + x1*x1*x1;
    fdyn[8] = - 1.25*x1 + 1.04*x2 + 0.990*x3 - 1.44*x4 - 1.44*x5 - 1.29*x6 - 1.20*x9 \
            - 1.16*x10 - 3.15*x11 - 0.998*x12 + x1*x6*x10;
    fdyn[9] = 1.43*x2 + 0.835*x3 + 0.545*x7 - 0.890*x8 - 0.973*x10 + x2*x3*x8;
    fdyn[10] = - 2.29*x3 + 1.41*x4 + 1.43*x5 + 0.512*x6 + 2.00*x7 - 0.801*x8 + \
            2.46*x9 - 0.659*x10 - 0.833*x11 + 1.29*x12 + x6*x6*x6;
    fdyn[11] = - 0.850*x1 - 0.919*x3 + 0.515*x7 + 1.20*x9 - 1.25*x11 + x1*x9*x11;

    double val = fdyn*(M*fdyn);
    return val;
}

/* Local optimization using NLOPT */

double localOptimization(DVector cx, double delta, DMatrix M, double dt, IVector box,
    double (*fpObj)(const std::vector<double> &, std::vector<double> &, void *)){

    double fval = 0;
    
    // Need to convert (cx: DVector) to (cx: vector<double>) for NLOPT library
    int dim = cx.dimension();
    std::vector<double> cVec(dim);
    for(int i=0; i<dim; i++)
        cVec[i] = cx[i];

    /* Setting the problem */
    //cout << "before initiation..." << endl; 
    // create the optimization object
    nlopt::opt opt(nlopt::LN_COBYLA, dim+1); // time is (dim+1)-th dimension 
    
    // setting up lower bound
    vector<double> lbVec(dim+1);
    for(int i=0; i<dim; i++)
        lbVec[i] = leftBound(box[i]);
    lbVec[dim] = 0; // adding lower bound for time: 0         
    opt.set_lower_bounds(lbVec);
    
    // setting upper bound
    vector<double> ubVec(dim+1);
    for(int i=0; i<dim; i++)
        ubVec[i] = rightBound(box[i]);
    ubVec[dim] = dt; // adding upper bound for time: dt         
    opt.set_upper_bounds(ubVec);
    
    // setting objective function
    //cout << "before setting obj..." << endl;
    opt.set_max_objective(fpObj, &M); // additional data is M
    
    // adding inequality constraint: x'*M*x - delta^2 < 0

    constraint_data consData(cx, delta, M);
    opt.add_inequality_constraint(cons_func, &consData, 1e-8);

    // set tolerance 
    opt.set_xtol_rel(1e-4);
    
    // initial guess
    vector<double> x0(dim+1);
    for(int i=0; i<dim; i++)
        x0[i] = 0.5*(lbVec[i]+ubVec[i]);
    x0[dim] = 0; // time is initially 0
    
    //cout << "before call optimization..." << endl;
    nlopt::result res = opt.optimize(x0, fval);
    //cout << "After call optimization..." << endl;
    return fval;
}
/********************* End of Local Optimization ***************************/

/********** Validated global optimization using dReal **********************************/

void gen_script(ex fcon1, double fval, ex fcon2, std::vector<ex> xvec, IVector box, ex tvar, double dt){
    std:: ofstream fsrcipt;
    fsrcipt.open("script.dr", std::ofstream::out);
    int dim = box.dimension();

    /* Define variables */
    fsrcipt << "var: " << endl;
    for(int i=0; i < dim; i++){
        fsrcipt << "[ " << leftBound(box[i]) << ", " << rightBound(box[i]) << "] " << xvec[i] << ";" << endl;
    }
    // adding time variable
    fsrcipt << "[0, " << dt << "] " << tvar << ";" << endl;
    // define pi variable -- symbolic expression may contain it
    fsrcipt << "[3.141592653589793, 3.141592653589794] pi;" << endl;

    /* Define constraints */
    fsrcipt << endl << "ctr: " << endl;
    // first constraint: fdyn'*M*fdyn > fval;
    fsrcipt << fcon1 << " >= " << fval << ";" << endl;
    // first constraint: fdyn'*M*fdyn > fval;
    fsrcipt << fcon2 << " <= 0;" << endl;

    fsrcipt.close();


}

double dRealMax(double rad, double fval0, DVector cx, DMatrix M, IVector box, double dt,
    std::vector<ex> (*fpDyn)(std::vector<ex>, ex), int max_iter){
    
    int dim = cx.dimension();
    // create symbolic xVec(dim)
    std::vector<ex> xvec = sym_vec("x", dim);
    symbol tsym("t");
    ex tvar = pow(tsym,1);

    // Constructing symbolic Constraint: fdyn'*M*fdyn  
    std::vector<ex> fdyn = fpDyn(xvec, tvar);
    std::vector<ex> mxvec = sym_mult(M, fdyn); // M*fdyn
    ex sym_fcon1;
    sym_fcon1 = sym_inner(mxvec, fdyn); // fdyn'*(M*fdyn);
    
    // Constructing equation: ((x-cx)'*M*(x-cx)-r^2), x is symbolic vector 
    std::vector<ex> new_xvec(dim); // new_xvec = xvec - cx
    for(int i=0; i < dim; i++){
        new_xvec[i] = xvec[i] - cx[i];   
    }
    mxvec = sym_mult(M, new_xvec); // M*(xvec-cx)
    ex sym_fcon2 = sym_inner(mxvec, new_xvec) - rad*rad;
    
    /* running dReal in a loop */
    double fval = fval0;
    int iter;
    FILE *fp;
    for(iter = 0; iter <= max_iter; iter++){
        // generate dReal script 
        //cout << "dReal loop: fval = " << fval << endl;
        gen_script(sym_fcon1, fval, sym_fcon2, xvec, box, tvar, dt);
        // run the dreal command in a separate process and save the output in fout
        string command = "dReal_new script.dr --precision " + std::to_string(DREAL_DELTA);
        FILE *fp = popen(command.c_str(), "r");

        char sys_out[1024];
        if(fp != NULL){
            fgets(sys_out, 1024, fp); // get the command output in sys_out
        } else {
            cout << "something is wrong here! " << endl;
        }
        if(sys_out[0] == 'u'){
            pclose(fp);
            return fval;
        }
        
        fval = fval*DREAL_FACTOR;
        pclose(fp);
        //system("rm script.dr script.out");
    } 
    
    cout << "Warning: dReal could not find maximum with provided maximum iteration!!!";
    exit(1);
    //return fval;    
}
/********** End of Validated global optimization using dReal **********************************/

/***************************** Continuous Flowpipe Computation *************************/
void  bloatBox(IVector &box, double dtilde){
    int dim = box.dimension();
    for(int i=0; i<dim; i++){
        box[i] = box[i]+interval(-dtilde, dtilde);
    }
}

double computeCset(DVector cx, double delta0, DMatrix M, double dt, IVector box,
    double (*fpObj)(const std::vector<double> &, std::vector<double> &, void *),
    std::vector<ex> (*fpDyn)(std::vector<ex>, ex)){

    //cout << "Inside computeCset..." << endl;
    int dim = cx.dimension();
    IVector cSet(dim);
    // dtilde: dtilde_factor = f(delta0) = 10*delta0 + 0.1 ... some linear function of delta0
    //double dtilde = TILDE_FACTOR*delta0; 
    double dtilde_factor = TILDE_FACTOR*delta0 + 0.15;
    double dtilde = delta0*dtilde_factor;
    double Delta = delta0 + dtilde; // bloat by dtilde
    cout << "Initial Delta: " << delta0 << endl;
    // local optimization as a initial guess
    double fval0, fval;
    fval0 = localOptimization(cx, Delta, M, dt, box, fpObj);
    // validate this maximum value using dReal
    fval = dRealMax(Delta, fval0, cx, M, box, dt, fpDyn, DREAL_MAX_ITER);
    cout << "fval0 = " << fval0 << "\t" << "fval = " << fval << endl;
    interval i_fval = std::abs(fval);
    double  fmax = rightBound(dt*sqrt(i_fval)); // sqrt needs to be rigorous
    while(Delta < (delta0 + fmax)){
        cout << "Delta: " << Delta << endl;
        dtilde = delta0*dtilde_factor*0.5; // dtilde decrease by half
        Delta = Delta + dtilde;  
        bloatBox(box, dtilde); // bloating box by dtilde amount
        fval = dRealMax(Delta, fval, cx, M, box, dt, fpDyn, DREAL_MAX_ITER);
        i_fval = std::abs(fval);
        fmax = rightBound(dt*sqrt(i_fval));
    }
    cout << "Final Delta: " << Delta << endl;
    return Delta;
}


/*********************************     MAIN ALGORITHM  ************************/
typedef capd::dynsys::Solver<IMap> OdeSolver;

int main(int argc, char *argv[]){
        
/******************** Tool options ***********************************/
        const bool EFF_EIG = true;
        double EPS = 0.05;
        
        double ist_th = 0.01;
        double norm_th = 0.001;
        

        const int MAXITER = 10;
        bool alw_change = false;
        int counter = 0; // count the number of integration step
        int step = 1; // Norm change after this number of steps
        // When to change Norm 
        
/******************** Debuging option *********************************/
        _DEBUG = true; // first-level debug output
        _DEBUG_ = false; // second-level debug output

/********* Processing command line argument: Norm and Grad *************/
        if(argc < 3){
            cout << "usage:\t./mainCLRT <model_name> <Time_horizon> <ist_th> <norm_th> <delta> "<< endl;
            cout << "       \tModel name example: bruss, fvdp, robot etc." << endl;
            cout << "       \tTIME_horizon > 0 : Computation horizon: " << endl;
            cout << "       \tOptional: Threshold for IST in (0,1) with default value 0.01" << endl; 
            cout << "       \tOptional: Threshold for norm change (0, 0.1) with default value 0.001" << endl;
            cout << "       \tExample run: ./mainCLRT bruss 5 0.02 0.001 1e-5"<< endl;
            return EXIT_FAILURE;
        }

        string model = argv[1]; 
        double TIMEHORIZON = atof(argv[2]);

       //  IST Threshold 
       if(argc > 3){
            ist_th = atof(argv[3]);
       } 

       // Norm Changing threshold
       if(argc > 4){
            norm_th = atof(argv[4]);
       }
       // dreal delta
       if(argc > 5){
            DREAL_DELTA = atof(argv[5]);
       }
       
       if(_DEBUG){
            cout << "First level debug output is chosen..." << endl;
            cout << "Model: " << model << endl;
            cout << "Time Horizon: " << TIMEHORIZON << endl;
            cout << "IST Threshold: " << ist_th << endl;
            cout << "Norm Change Threshold: " << norm_th << endl;
            cout << endl << endl;
       }

/******** Creating output files *******************************************************************/
        std::ostringstream file_suffix;
        file_suffix << model << "_" << TIMEHORIZON;
        file_suffix <<"_" << ist_th;
        file_suffix <<"_" << norm_th;
        file_suffix << ".txt";
        ofstream fout_dSet, fout_cSet;
        fout_dSet.open("outputs/dset_" + file_suffix.str(), std::ofstream::out);
        fout_cSet.open("outputs/cset_" + file_suffix.str(), std::ofstream::out);
        /************** Setting precision of output files ************************************/
        fout_dSet.precision(16);
        fout_cSet.precision(16);
        cout.precision(16);

/******************** Model Input *******************************************/
        int DIM;
        double initR, timeStep;
        
        // creating function pointer for objective function
        double (*fpObj)(const std::vector<double> &, std::vector<double> &, void *);
        

        // creating function pointer for symbolic model equation
        std::vector<ex> (*fpDyn)(std::vector<ex>, ex);

        string vectorString;
        if(model == "vdp") {
            fpObj = &obj_vdp;
            fpDyn = &sym_vdp;
            TILDE_FACTOR = 35;
            DIM = 2;
            initR = 0.01;    // initial radius
            timeStep = 0.005; // time-step in CAPD solver
            //TIMEHORIZON = 2; // integration duration
            
            //vectorField("time:t;var:x,y;fun: y, (1-x^2)*y-x+1.2*sin(0.6283*t);");
            vectorString = "var:x,y;fun: -y, (x^2-1)*y+x;";
            //x0[0] = -1; x0[1] = -1;
        } else if(model == "fivdp") {
            DIM = 2;
            initR = 0.01;    // initial radius
            timeStep = 0.01; // time-step in CAPD solver
            EPS = 0.02;
            //TIMEHORIZON = 2; // integration duration
            //vectorField("time:t;var:x,y;fun: y, (1-x^2)*y-x+1.2*sin(0.6283*t);");
            vectorString = "time:t;var:x,y;fun: -y, (x^2-1)*y+x+1.2*sin(0.6283*t);";
            //x0[0] = -1; x0[1] = -1;
        } else if (model == "bruss") {
            fpObj = &obj_bruss;
            fpDyn = &sym_bruss;
            DIM = 2;
            initR = 0.01;    // initial radius
            EPS = 0.03;
            timeStep = 0.01; // time-step in CAPD solver
            //TIMEHORIZON = 5; // integration duration
            vectorString = "var:x,y;fun: 1+x*x*y-2.5*x, 1.5*x-x*x*y;";
            cout << "Running brusselator model: var:x,y;fun: 1+x*x*y-2.5*x, 1.5*x-x*x*y" << endl;
            cout << "Time Horizon = " << TIMEHORIZON << endl;
            cout << "Initial Time Step = " << timeStep << endl; 
            cout << endl;

        } else if (model == "dubins"){
            fpObj = &obj_dubins;
            fpDyn = &sym_dubins;
            DREAL_FACTOR = 1.001; // special case
            TILDE_FACTOR = 5;

            DIM = 3;
            initR = 0.01;
            timeStep = 0.01;
            EPS = 0.01;
            //TIMEHORIZON = 0.03;
            double v = 1; // velocity
            DUBINS_V = v;
            vectorString = "time: t;var:x,y,th;fun: "+ to_string(v)+"*cos(th),"+ to_string(v) +"*sin(th), x*sin(t);";
            cout << "Running Dubins car model: " << vectorString << endl;
            cout << "Time Horizon = " << TIMEHORIZON << endl;
            cout << "Initial Time Step = " << rightBound(timeStep) << endl; 
            cout << endl;
            
        } else if (model == "dubins_2"){
            fpObj = &obj_dubins_2;
            fpDyn = &sym_dubins_2;
            DREAL_FACTOR = 1.002; // special case
            TILDE_FACTOR = 15;

            DIM = 3;
            initR = 0.01;
            timeStep = 0.01;
            //EPS = 0.01;
            //TIMEHORIZON = 0.03;
            double v = 1; // velocity
            DUBINS_V = v;
            vectorString = "time: t;var:x,y,th;fun: "+ to_string(v)+"*cos(th),"+ to_string(v) +"*sin(th), x*y*sin(t);";
            cout << "Running Dubins car model: " << vectorString << endl;
            cout << "Time Horizon = " << TIMEHORIZON << endl;
            cout << "Initial Time Step = " << rightBound(timeStep) << endl; 
            cout << endl;
            
        } 
        else if (model == "dubins_ti"){
            fpObj = &obj_dubins_ti;
            fpDyn = &sym_dubins_ti;
            DREAL_FACTOR = 1.0001; // special case
            TILDE_FACTOR = 2;

            DIM = 3;
            initR = 0.01;
            timeStep = 0.01;
            EPS = 0.01;
            //TIMEHORIZON = 0.03;
            double v = 1; // velocity
            DUBINS_V = v;
            vectorString = "var:x,y,th;fun: "+ to_string(v)+"*cos(th),"+ to_string(v) +"*sin(th), 1;";
            cout << "Running Dubins car model: " << vectorString << endl;
            cout << "Time Horizon = " << TIMEHORIZON << endl;
            cout << "Initial Time Step = " << rightBound(timeStep) << endl; 
            cout << endl;
            
        } 
        else if (model == "robot"){
            fpObj = &obj_robot;
            fpDyn = &sym_robot;
            DIM = 4;
            initR = 0.005;
            timeStep = 0.01;
            EPS = 0.01;
            //TIMEHORIZON = 1;
            vectorString = "par:m,l,kp1,kp2,kd1,kd2;var:x1,x2,x3,x4; \
            fun:x3,x4,(-2*m*x2*x3*x4-kp1*x1-kd1*x3)/(m*x2*x2+l/3)+(kp1*kp1)/(m*x2*x2+l/3), \
            x2*x3*x3-kp2*x2/m-kd2*x4/m+kp2*kp2/m;";
    
            cout << "Running Robotic Arm Model: " << vectorString << endl;
            cout << "Time Horizon = " << TIMEHORIZON << endl;
            cout << "Initial Time Step = " << rightBound(timeStep) << endl; 
            cout << endl;
        } else if (model == "ms") {
            fpObj = &obj_ms;
            fpDyn = &sym_ms;
            DIM = 2;
            initR = 0.001;
            timeStep = 0.005;
            //TIMEHORIZON = 0.03;
            vectorString = "var:x1,x2;fun: x2*x1^2* (1 - x1)/0.3 - x1/6, (0.5*(1 + (1-exp(-100*x1+10))/(1-exp(-100*x1+10))))*x2/150 + (1 - (0.5*(1 + (1-exp(-100*x1+10))/(1-exp(-100*x1+10)))))*(1 - x2)/20;";
            cout << "MS model: " << vectorString << endl;
            cout << "Time Horizon = " << TIMEHORIZON << endl;
            cout << "Initial Time Step = " << rightBound(timeStep) << endl; 
            cout << endl;     

        } else if (model == "biology") {
            fpObj = &obj_biology;
            fpDyn = &sym_biology;
            DIM = 7;
            initR = 0.0001;
            timeStep = 0.005;
            //EPS = 0.05;
            //TIMEHORIZON = 0.03;
            vectorString = "var:x1,x2,x3,x4,x5,x6,x7;\
                fun: -0.4*x1 + 50*x3*x4,\
                0.4*x1 - x2,\
                x2 - 50*x3*x4,\
                50*x5*x6 - 50*x3*x4,\
               -50*x5*x6 + 50*x3*x4,\
                0.5*x7 - 50*x5*x6,\
               -0.5*x7 + 50*x5*x6;";
            cout << "Biology model: " << vectorString << endl;
            cout << "Time Horizon = " << TIMEHORIZON << endl;
            cout << "Initial Time Step = " << rightBound(timeStep) << endl; 
            cout << endl;     

        } else if(model=="twelve"){
            fpObj = &obj_twelve;
            fpDyn = &sym_twelve;
            DIM = 12;
            initR = 0.0001;
            timeStep = 0.005;
            //EPS = 0.05;
            TILDE_FACTOR = 50;
            vectorString = "var:x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12;\
            fun: -1.47*x1 - 0.574*x2 - 0.583*x3 + 0.942*x4 + 0.714*x5 \
            +1.35*x6 + 0.861*x7 + 1.51*x9 + 1.10*x11 + 1.14*x12 + x2*x6*x9,\
            -0.76*x2 + 1.54*x3 - 1.63*x4 - 0.819*x5 - 0.925*x6 - 1.06*x7 - \
            0.535*x10 + 0.681*x11 + x7*x10*x11,\
            1.77*x1 - 0.983*x2 - 0.766*x3 - 0.930*x6 - 1.13*x9 - 0.998*x10 \
            + 1.84*x11 + 0.772*x12 + x3^3,\
            -1.37*x1 + 1.31*x2 + 0.790*x3 - 0.572*x4 - 1.07*x5 - 0.783*x6 - \
            0.938*x8 + 2.03*x9 - 0.857*x11 + x4*x8*x11,\
            -1.04*x1 + 0.755*x2 + 1.12*x4 - 0.538*x5 - 0.563*x8 + 1.40*x9 - \
            1.32*x11 + x1*x2*x5,\
            1.26*x2 + 1.11*x4 - 0.885*x6 - 0.935*x8 + 1.05*x9 - 1.30*x11 + x2*x4*x11,\
            - 0.845*x1 - 1.67*x7 + 1.72*x9 + 0.839*x10 - 1.25*x11 - 0.600*x12 + x1*x7*x11,\
            0.583*x1 + 1.20*x4 + 0.816*x5 + 0.599*x6 - 0.735*x9 + x1^3,\
            - 1.25*x1 + 1.04*x2 + 0.990*x3 - 1.44*x4 - 1.44*x5 - 1.29*x6 - 1.20*x9 \
            - 1.16*x10 - 3.15*x11 - 0.998*x12 + x1*x6*x10,\
            1.43*x2 + 0.835*x3 + 0.545*x7 - 0.890*x8 - 0.973*x10 + x2*x3*x8,\
            - 2.29*x3 + 1.41*x4 + 1.43*x5 + 0.512*x6 + 2.00*x7 - 0.801*x8 + \
            2.46*x9 - 0.659*x10 - 0.833*x11 + 1.29*x12 + x6^3,\
            - 0.850*x1 - 0.919*x3 + 0.515*x7 + 1.20*x9 - 1.25*x11 + x1*x9*x11;";
        } else {// default is also bruss
            DIM = 2;
            initR = 0.01;    // initial radius
            timeStep = 0.005; // time-step in CAPD solver
            TIMEHORIZON = 1; // integration duration
            
            vectorString = "var:x,y;fun: y, 0;";
            //x0[0] = 1;  x0[1] = 1;

        }
        
/********************* Setting up initial region *****************************************************************/
        try{
            
            IVector x0(DIM);
            if(model == "bruss"){
                x0[0] = 1;  x0[1] = 1;
            } else if (model == "fvdp") {
                x0[0] = -1;  x0[1] = -1;
            } else if (model == "fivdp"){
                x0[0] = -1; x0[1] = -1;
            } 
            else if(model == "dubins") {
                x0[0] = 0; x0[1] = 0; x0[2] = 0.7854;  
            } else if(model == "dubins_2") {
                x0[0] = 0; x0[1] = 0; x0[2] = 0.7854;  
            } 
            else if(model == "dubins_ti") {
                x0[0] = 0; x0[1] = 0; x0[2] = 0.7854;  
            } 
            else if (model == "robot") {
                x0[0] = 1.505;  x0[1] = 1.505; x0[2] = 0.005; x0[3] = 0.005;
            } else if (model == "ms") {
                x0[0] = 0.8;    x0[1] = 0.5;
            } else if (model == "vdp") {
                x0[0] = -1;    x0[1] = -1;
            } else if(model == "biology"){
                for(int i=0; i<DIM; i++)
                    x0[i] = 0.1;
            } else if(model=="twelve"){
                for(int i=0; i<DIM; i++){
                    x0[i] = 0.01;
                }
            }              
            else {
                x0[0] = 1;  x0[1] = 1;
            }
            /*intial radius*/
            interval currentR = initR;
            
            /* initial set */
            IVector currentSet(DIM);
            IVector currentSet_new(DIM);
            
            /* Initialize currentSet */
            for(int i=0; i < DIM; i++) {
                currentSet[i] = x0[i] + Interval(-initR, initR);
            }
            //cout << "you are before CAPD settings" << endl;
/************* Settings for CAPD solver *******************************************/
            IMap vectorField(vectorString);

            // Setting parameter values
            if(model == "robot") {
                vectorField.setParameter("m",interval(1)); 
                vectorField.setParameter("l", interval(3));
                vectorField.setParameter("kp1", interval(2)); 
                vectorField.setParameter("kp2", interval(1)); 
                vectorField.setParameter("kd1", interval(2)); 
                vectorField.setParameter("kd2", interval(1));
            }

            OdeSolver solver(vectorField,10); // solver object
            solver.setAbsoluteTolerance(1e-10);
            solver.setRelativeTolerance(1e-10);
            solver.setStep(timeStep);

/********************* Intermediate varaible declarations ****************************************/
            /* variable related to center; center is computed using CAPD here*/
            IVector currentCenter(DIM), CcurrentCenter(DIM), newCcurrentCenter(DIM), currentCcurrentCenter(DIM);  
            IMatrix igrad(DIM,DIM) , icgm(DIM,DIM), imcgm(DIM,DIM), imcgm_cur(DIM,DIM), imcgm_new(DIM,DIM),centergrad(DIM,DIM), currentC(DIM, DIM), currentCi(DIM, DIM), newC(DIM, DIM), newCi(DIM, DIM);
            DMatrix dCi(DIM, DIM), dC(DIM,DIM);
            IMatrix dgrad(DIM,DIM), sq_dgrad(DIM,DIM), epsMat(DIM, DIM), omegaMat(DIM, DIM);
            double lambdaMax, lambdaMaxM, lambdaMaxM_cur, lambdaMaxM_new;
            interval eig_dgrad;
            double norm_dgrad; 
            
            interval stretch2m, stretch2, stretch2m_cur, stretch2m_new;
            
            /* currentSet representation based on GRAD computation option*/
            C1Rect2Set cfullSet(currentSet);
            C1Rect2Set temp_cfullSet(currentSet);

            currentCenter = x0;
            /* C1 representaion of center_set*/
            C1Rect2Set centerPt(currentCenter);
            
/************** Initializing coordiante system with Identity matrix **********************/
            
            currentC = IMatrix::Identity(DIM); // C interval matrix
            currentCi = IMatrix::Identity(DIM); // C inverse interval matrix
            newC = IMatrix::Identity(DIM); // C interval matrix
            newCi = IMatrix::Identity(DIM); // C inverse interval matrix

/****************** Dumping initial discrete set and norm in the output file *********************/
           // initial discrete set
            fout_dSet << "0.0 "; // initial time = 0;
            for(int i=0; i< DIM; i++){
                fout_dSet << leftBound(currentCenter[i]) <<" " << rightBound(currentCenter[i])<<" "; 
            }
            fout_dSet << initR << " ";

            // initial Norm for discrete set (Euclidean)
            for(int i=0; i< DIM; i++){
                for(int j=0; j < DIM; j++){
                    fout_dSet << rightBound(currentC[i][j]) <<" "; 
                } 
            }
            fout_dSet << endl;
/************************ Main integration loop *******************************************/
            int nsteps = 1;
            /* Integration parameter*/
            double adapt_dt = timeStep; // initial value of adaptive time step
            double totalTime = 0.0; // increases by T in each iteration
            double temp_totalTime = 0;
            bool isfstep = false; // to handle boundary condition in cSet computation
            while(totalTime < TIMEHORIZON){
                if(nsteps % 2 == 1) { // change timeStep in odd steps
                    adapt_dt = timeStep; // initialize with user provided value
                    int iter = 0; // IST iteration
                    /********************** IST Loop for  to Compute adaptive dt **************************/
                    do { 
                        iter++;
                        temp_totalTime = totalTime + adapt_dt;
                        temp_cfullSet = C1Rect2Set(currentSet);
                        solver.setStep(adapt_dt); 
                        temp_cfullSet.setCurrentTime(temp_totalTime);
                        temp_cfullSet.move(solver); // integrate cfullSet using C1 algorithm
                        igrad = (IMatrix)temp_cfullSet;

                        // compute Displacement Gradient Tensor: [F] - I
                        dgrad = igrad - IMatrix::Identity(DIM); // need to define identity matrix 
                        sq_dgrad = transpose(dgrad)*dgrad; 
                        eig_dgrad = computeMaxEigPD(sq_dgrad);
                        
                        // compute norm of Displacement Gradient Tensor
                        norm_dgrad = rightBound(sqrt(eig_dgrad));
                        if(_DEBUG){
                            if(MAXITER == iter)
                            {
                                std::cout << "Warning: Max iteration reached!";
                            }
                        }
                        // test whether dt needs to be reduced
                       if(norm_dgrad > ist_th)
                        {
                            adapt_dt = 0.5*adapt_dt;
                        }

                    } while(iter <= MAXITER && norm_dgrad > ist_th); // end of IST Loop
                    /********************** End of IST Loop *************************************/
                }
                // final step
                if(totalTime + adapt_dt > TIMEHORIZON){
                    adapt_dt = TIMEHORIZON - totalTime;
                    isfstep = true;
                }
            
                std::cout << endl << "################ current time = [" << totalTime << "," << totalTime + adapt_dt << "] ################### \n";
                cout << "current dt: " << adapt_dt << endl;
                if(_DEBUG){ // output M0
                    cout << "Current Center: " << endl;
                    printIVector(currentCenter);
                    cout << "Current Radius: " << rightBound(currentR) << endl;
                    cout << endl;
                    cout << "Current Norm M0: " << endl;
                    IMatrix iM0 = transpose(currentC)*currentC;
                    DMatrix M0 = takeMidM(iM0);
                    printDMatrix(M0);
                    cout << endl;
                }

                // update totaltime with proper dt
                totalTime = totalTime + adapt_dt; 
        
                /********** Compute New Center at t+dt using CAPD ***************************/
                solver.setStep(adapt_dt); // update solver timestep with new dt
                centerPt.setCurrentTime(totalTime); // for time-varying system
                centerPt.move(solver); // integrate one-step
                centergrad = (IMatrix)centerPt; // gradient of center at t+dt, used for optimal norm computaion
                
                /******* Convert |B|_M0 into |B|_2 for CAPD Integration *******/
                IVector currentCBds(DIM);
                currentCcurrentCenter = currentC * currentCenter;  
                for(int i=0; i < DIM; i++){
                    currentCBds[i] = currentCcurrentCenter[i] + interval(-rightBound(currentR) , rightBound(currentR) );
                }
                currentSet = currentCi * currentCBds; 
                
                /****** Computing stretching factor from CG matrix *********/   
                // step 1: Gradient computation using C1 integration of |B|_2 in CAPD
                //cfullSet = C1Rect2Set(currentSet);
                cfullSet.setCurrentTime(totalTime);
                cfullSet.move(solver); 
                igrad = (IMatrix)cfullSet;
                
                // step 2: CG matrix in current norm
                IMatrix U(DIM, DIM);
                U = currentC * igrad * currentCi;
                imcgm_cur = transpose(U)*U;

                // step 3: Compute lambda max of CG
                lambdaMaxM_cur = computeIntEig(imcgm_cur);
                // step 4: Compute stretching factor in current norm 
                stretch2m_cur = sqrt(lambdaMaxM_cur);
                
                // step-5: Compute New Norm
                
                DMatrix F = takeMidM(igrad);
                if(_DEBUG){
                    cout << "Computation of new Norm:" << endl;
                    cout << "F_mid: " << endl;
                    printDMatrix(F); 
                    cout << endl;
                }

                /* Compute new norm: M1 = C1^T*C1, here, C1 = newC and C0 = currentC */
                dC = computeNorm(F, DIM);
                IMatrix iCi(dCi), iC(dC);   //interval versions
                dCi = capd::matrixAlgorithms::gaussInverseMatrix(dC);
                iCi = capd::matrixAlgorithms::gaussInverseMatrix(iC);
            
                //set newC and newCi
                newC = iC;
                newCi = iCi;
                if(_DEBUG){
                    cout << "New Norm M1: " << endl;
                    IMatrix iM1 = transpose(newC)*newC;
                    DMatrix M1 = takeMidM(iM1);
                    printDMatrix(M1);
                    cout << endl;
                }

                //CG matrix in M1 (new) norm
                U = (newC * igrad) * newCi;
                imcgm_new = transpose(U)*U;
                
                // CG matrix in M0-M1 norm
                U = (newC * igrad) * currentCi;
                imcgm = transpose(U) * U;

                // compute eigenvalues
                lambdaMaxM = computeIntEig(imcgm);
                lambdaMaxM_new = computeIntEig(imcgm_new);

                // compute stretching factor
                stretch2m_new = sqrt(lambdaMaxM_new); 
                stretch2m = sqrt(lambdaMaxM);                        
                if(_DEBUG) {
                    cout << "CG0: ";
                    printIMatrix(imcgm_cur);
                    cout << endl;

                    cout << "CG1: ";
                    printIMatrix(imcgm_new);
                    cout << endl;

                    cout << "CG: ";
                    printIMatrix(imcgm);
                    cout << endl;

                    cout << "M0 Stretching factor: " << rightBound(stretch2m_cur) << endl;
                    cout << "M1 Stretching factor: " << rightBound(stretch2m_new) << endl;
                    cout << "M0/M1 Stretching factor: " << rightBound(stretch2m) << endl;
                }

                // checking norm changing condition
                double th = rightBound(stretch2m_cur) - rightBound(stretch2m_new);
                if( th >= norm_th ) {
                    cout << "*** Changing to new norm ***" << endl << endl;
                    currentC = newC;
                    currentCi = newCi;
                } else {
                    cout << "*** Keeping current norm ***" << endl;
                    stretch2m = stretch2m_cur;
                }

                // step 6: Compute the new radius
                cout << "Used SF: " << rightBound(stretch2m) << endl << endl;
                currentR = initR * stretch2m; 
                //currentR = currentR*stretch2m;   
                /* Reinitialize center and set for next step t1*/
                currentCenter = (IVector)centerPt; // new center at t+dt
                
                /***** write discrete set in output file *****/
                fout_dSet << totalTime << " ";
                for(int i=0; i< DIM; i++){
                    fout_dSet << leftBound(currentCenter[i]) <<" " << rightBound(currentCenter[i])<<" "; 
                }
                fout_dSet << rightBound(currentR) << " ";
                // Norm
                IMatrix normMat(DIM,DIM);
                normMat = transpose(currentC)*currentC; 
                for(int i=0; i<DIM; i++){
                    for(int j=0; j<DIM; j++){
                        fout_dSet << rightBound(normMat[i][j]) << " ";
                    }
                } 
                fout_dSet << endl;

                if(_DEBUG){
                    cout << "New Center: " << endl;
                    printIVector(currentCenter);
                    cout << "New Radius: " << rightBound(currentR) << endl;
                    cout << endl;
                    cout << "New Norm M: " << endl;
                    IMatrix iM = transpose(currentC)*currentC;
                    DMatrix M = takeMidM(iM);
                    printDMatrix(M);
                    cout << endl;
                }   
                if((nsteps % 2) == 1 ){ // Continuous parts in every other steps
                    //if(_DEBUG){
                    cout << "Begining of Continuous parts: " << endl;
                    DVector cx = takeMidIV(currentCenter);
                    double delta = rightBound(currentR);
                    DMatrix M = takeMidM(normMat);
                    double dt = adapt_dt;
                    IVector box = (IVector)cfullSet;
                    double contRad = computeCset(cx, delta, M, dt, box,fpObj, fpDyn);

                    /* dumping result to file */
                    // time interval 
                    fout_cSet << totalTime - adapt_dt << " " << totalTime + adapt_dt << " ";
                    
                    // center cx
                    for(int i=0; i < DIM; i++)
                        fout_cSet << cx[i] << " ";

                    // radius
                    fout_cSet << contRad << " ";

                    // M Matrix
                    for(int i=0; i < DIM; i++){
                        for(int j=0; j < DIM; j++){
                            fout_cSet << M[i][j] << " ";
                        }
                    }
                    fout_cSet << endl;
                }
                nsteps++;
            } // END of main while loop
            
        } catch(exception& e) {
            cout << "\n\nException caught!\n" << e.what() << endl << endl;
        }
        
        /* Closing output files */
        fout_dSet.close();
        fout_cSet.close();
        return EXIT_SUCCESS;
    } // END
