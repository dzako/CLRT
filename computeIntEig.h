/* 
 *Lagrangian Reachtube Computations
 *Authors: Md Ariful Islam and Jacek Cyranks
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

/* Library defined by myself*/
#include "capd_api.h"
#include "capd_to_eigen_api.h"

/*Define namespace*/
using namespace Eigen;
using namespace capd;
using namespace std;

/* Define some GLOBAL VARIABLES*/
#define REAL_MIN 2.225073858507201e-308
#define MMAX 100
#define INF 1.7976931348623157E+308

bool _DEBUG; // only important 
bool _DEBUG_; // many more

////////////////// WRAPPER for STD::IOTA //////////////////////////////////////
template < class T > struct  IotaWrapper
{
    typedef T type;
    typedef std::function<type(const type&)> IncrFunction;

    type value;
    IncrFunction incrFunction;
    IotaWrapper () = delete;
    IotaWrapper(const type& n, const IncrFunction incrFunction): value(n), incrFunction(incrFunction){};
    operator type() { return value; }
    IotaWrapper& operator++ () { value = incrFunction(value); return *this;} 
};

////////////////// Begin of sorting with index map /////////////////////////////////////
/* Comparison struct used by sort method*/
template <class T> struct index_cmp 
{
  index_cmp(const T arr) : arr(arr) {}
  bool operator()(const int a, const int b) const
  { 
    return arr[a] < arr[b];
  }
  const T arr;
};

/* This implementation is O(n), but also uses O(n) extra memory*/
template< class T >
void reorder(
  std::vector<T> & unordered, 
  std::vector<int> const & index_map, 
  std::vector<T> & ordered)
{
  // copy for the reorder according to index_map, because unsorted may also be
  // sorted
  std::vector<T> copy = unordered;
  ordered.resize(index_map.size());
  for(int i = 0; i<index_map.size();i++)
  {
    ordered[i] = copy[index_map[i]];
  }
}

/* Sorting that gives Index
Same as Matlab's: [B, I] = sort(A)
*/
template <class T>
void sort(
  std::vector<T> & unsorted,
  std::vector<T> & sorted,
  std::vector<int> & index_map)
{
  // Original unsorted index map
  index_map.resize(unsorted.size());
  for(size_t i=0;i<unsorted.size();i++)
  {
    index_map[i] = i;
  }
  // Sort the index map, using unsorted for comparison
  sort(
    index_map.begin(), 
    index_map.end(), 
    index_cmp<std::vector<T>& >(unsorted));

  sorted.resize(unsorted.size());
  reorder(unsorted,index_map,sorted);
}

////////////// BEGIN: UTIL FUNCTIONS ////////////////////////////////

/*Matlab's [eigVal, eigVec] = eig(A) method -- non-rigorous*/
void getEigen(MatrixXd A, VectorXd& eigVal, MatrixXd& eigVec){
    
   SelfAdjointEigenSolver<MatrixXd> eigensolver(A);
   if (eigensolver.info() != Success) abort();
   eigVal = eigensolver.eigenvalues();
   eigVec = eigensolver.eigenvectors();
}


/* Matlan's: C = setdiff(A,B)*/
vector<int> setdiff(vector<int> first, vector<int> second, int dim){
    //int first[] = {5,10,15,20,25};
    //int second[] = {50,40,30,20,10};
    //assert(dim >= first.size() && dim >= second.size());
    std::vector<int> v(dim);                      
    std::vector<int>::iterator it;

    std::sort (first.begin(),first.end());     //  5 10 15 20 25
    std::sort (second.begin(),second.end());   // 10 20 30 40 50

    it=std::set_difference (first.begin(), first.end(), second.begin(), second.end(), v.begin());
                                               //  5 15 25  0  0  0  0  0  0  0
    v.resize(it-v.begin());                      //  5 15 25
    assert(v.size() <= dim);

    // std::cout << "The difference has " << (v.size()) << " elements:\n";
    // for (it=v.begin(); it!=v.end(); ++it)
    //     std::cout << ' ' << *it;
    // std::cout << '\n';

  return v;
}

/* Gershgorin method*/
double gershgorinEig(const IMatrix& m){
    interval eig1 = capd::matrixAlgorithms::Gershgorin< IMatrix, true >::maxEigenValueOfSymMatrix(m);
    return rightBound(eig1);
}

////////////// BEGIN: Rump's Method for point Matrix -- Rigorous ////////////////////
void verifyeig(MatrixXd A, double lam, VectorXd rvec, Interval &iLam, IVector &iVec, int dim){
// consider rvecs contains single eigenvector corresponding to eigen value lam
    //cout << "Inside verifyeig: " << endl;
    //cout << "A: " << A << endl;
    /* sorting the element in rvecs */
    int rows = dim;
    int cols = 1;
    //int cols = rvecs.cols();
    //std::vector<double> row_sum(cols);
    VectorXd row_sum = rvec;
    //cout << endl << endl;
    //cout << "rvec: " << endl << rvec << endl;
    //cout << "row sum: " << endl << row_sum << endl;

    // copy row_sum to std::vector (vec_sum)
    vector<double> vec_sum;
    vec_sum.resize(row_sum.size());
    VectorXd::Map(&vec_sum[0], row_sum.size()) = row_sum;
    
    // sorting [N,I] = sort(sum(abs(xs),2));
    std::vector<int> sort_idx; // I: std::vector
    std::vector<double> sorted_sum; // N
    sort(vec_sum,sorted_sum,sort_idx);
    // converting std::vec to Eigen::Vector
    VectorXi idx_vec = Eigen::Map<Eigen::VectorXi,Eigen::Unaligned>(sort_idx.data(), sort_idx.size());
    
    // u = I(1:n-k)
    VectorXi u_idx = idx_vec.head(rows-cols);
    // v = I(n-k+1:n)
    VectorXi v_idx = idx_vec.tail(cols);
    MatrixXd midA = A; // A is point matrix, so midA is same as A

    // R = midA - lambda*speye(n);
    MatrixXd R = midA - lam*MatrixXd::Identity(dim, dim);
    //cout << "R: " << R << endl;
    //cout << "u_idx: " << endl << u_idx << endl;
    //cout << "v_idx: " << endl << v_idx << endl;
    
    //R(:,v) = -xs;
    R.col(v_idx(0)) = -rvec;
    //cout << "R: " << R << endl;

    //y = R\(midA*xs-lambda*xs);
    VectorXd y = R.colPivHouseholderQr().solve(midA*rvec-lam*rvec);
    //cout << "y: " << y << endl;

    //xs(u,:) = xs(u,:) - y(u,:);
    for(int i=0; i< rows-cols; i++){
        rvec(u_idx(i)) = rvec(u_idx(i)) - y(u_idx(i));
    }
    //cout << "rvec: " << rvec << endl;

    //lambda = lambda - sum(diag(y(v,:)))/k;
    /******** From here, we assume k = 1 */
    lam = lam - y(v_idx(0));
    //R = midA - lambda*speye(n);
    R = midA - lam*MatrixXd::Identity(dim, dim);
    //R(:,v) = -xs;
    R.col(v_idx(0)) = -rvec;
    R = R.inverse();
    //cout << "R_inv: " << endl << R << endl;
    
    //C = A - intval(lambda)*speye(n);
    iLam = lam;
    // convert A to iA
    DMatrix dA = eigen_to_capd(A, dim);
    IMatrix iA = d2iMatrix(dA);
    // for(int i=0; i<dim; i++){
    //  for(int j=0; j<dim; j++){
    //      iA[i][j] = A(i,j);
    //  }
    // }
    IMatrix iC = iA - iLam*IMatrix::Identity(dim);
    //cout << "C Matrix: " << iC << endl;

    //Z = - R * ( C * xs ); rowsXcols matrix
    // convert R to IMatrix
    DMatrix dR = eigen_to_capd(R, dim);
    IMatrix iR = d2iMatrix(dR);
    // for(int i=0; i<dim; i++){
    //  for(int j=0; j<dim; j++){
    //      iR[i][j] = R(i,j);
    //  }
    // }
    // rvecs to IMatrix 
    IVector iRvec = eigen_to_capd(rvec, dim);
    // for(int i=0; i<rows; i++){
    //  iRvec[i] = rvec(i);
    // }
    IVector iZ = - iR * (iC * iRvec);
    //cout << "Z: " << iZ << endl;

    // C(:,v) = -xs;
    for(int i=0; i < rows; i++){
        iC[i][v_idx[0]] = -iRvec[i];
    }

    //  C = speye(n) - R * C;
     iC = IMatrix::Identity(rows) - iR * iC; 
    
    //  Y = Z;
    IVector iY = iZ; 
    
    //Eps = 0.1*mag(Y)*hull(-1,1) + midrad(0,realmin);
    IVector iEps = mag(iY, dim)*interval(-0.1,0.1) + interval(-REAL_MIN, REAL_MIN);
    //cout << "Eps: " << iEps << endl;
    
    bool ready = false; 
    int m = 0;
    IVector iX(dim);
    IVector iXX(dim);
    std::vector<bool> isSubset;
    while(!ready && m < MMAX){
        m++;
        iX = iY + iEps;
        iXX = iX;
        iXX[v_idx[0]] = 0.0;
    //  // Y = Z + C*X + R*(XX*X(v,:));
        iY = iZ + iC*iX + iR*(iXX*iX[v_idx[0]]);
        //ready = all(all(in0(Y,X)));
         ready = all_in(iY,iX, dim);
        
     }

     if(ready){ // inclusion found
        //cout << "inclusion found in Rump's method" << endl;
        //cout << "m = " << m << endl;
         Interval rad = mag(iY[v_idx[0]]); // Eigenvalue Correction
         iLam = iLam + interval(-rightBound(rad), rightBound(rad));
         //Y(v,:) = 0;
         iY[v_idx[0]] = 0;
         iVec = iRvec + iY;
         //cout << "Interval Eigenvalue: " << iLam << endl;
        //cout << "Interval eigenvector: " << iVec << endl;

     } else {
        // use CAPD version
        cout << "Verified eigenvalue inclusion not found using Rump's method!!" << endl;
        cout << "Using Gershgorin Method ..." << endl;
        
        DMatrix dA = eigen_to_capd(A, dim); 
        IMatrix iA(dA);
        double lam = gershgorinEig(iA);
        //cout << "lam: " << lam << endl;
        iLam = lam;
        //exit(1);
     }
 
}
////////////// UTILITY FUNCTION for INTERLACING METHOD////////////////////////////////

//Rohn's unverified method
double rohnunverified (MatrixXd Ac, MatrixXd Ad, int dim){
    // non-rigourous eigenvalues for Ac
    VectorXd evalAc; 
    MatrixXd evecAc;
    getEigen(Ac,evalAc, evecAc); 

    // Non-rigorous Eigenvalue for Ad
    VectorXd evalAd; 
    MatrixXd evecAd;
    getEigen(Ad,evalAd, evecAd); // non-rigourous eigenvalue and eigenvector

    double spectral = evalAd[dim-1];

    return evalAc[dim-1]+spectral;

}

double unverifiedUpper(IMatrix iA, int dim){
    DMatrix magA = mag(iA);
    MatrixXd A = capd_to_eigen(magA, dim);

    // non-rigourous eigenvalue and eigenvector
    VectorXd evalA; 
    MatrixXd evecA;
    getEigen(A,evalA, evecA); 
    double eigAbs = evalA[dim-1];

    // Rohn's unverified method
    DMatrix dAc = takeMidM(iA);
    DMatrix dAd = takeRadM(iA);
    MatrixXd Ac = capd_to_eigen(dAc, dim);
    MatrixXd Ad = capd_to_eigen(dAd, dim);
    double maxEig = rohnunverified(Ac, Ad, dim);
    return std::min(maxEig, eigAbs);

}

// Courant Fischer Theorem
// λ_m(iA) ≤ λ_m(mag(iA)).
double cfEigMax(IMatrix iA, int dim) {
    double lamMax;
    
    //cout << "inside cfEigMax " << endl << "iA: " << endl;
    //printIMatrix(iA);
    DMatrix magA = mag(iA);
    MatrixXd A = capd_to_eigen(magA, dim);

    // non-rigourous eigenvalue and eigenvector
    VectorXd evalA; 
    MatrixXd evecA;
    getEigen(A,evalA, evecA); 

    Interval i_evalA;
    IVector i_evecA(dim);
    if(evalA[dim-1] != 0)
        verifyeig(A, evalA(dim-1), evecA.col(dim-1), i_evalA, i_evecA, dim);
    else
        i_evalA = 0.0;

    lamMax = rightBound(i_evalA);    
    //cout << "Lam: " << lamMax << endl;
    return lamMax;
}

////////////// BEGIN: Rohn's Method for Interval Matrix -- Rigorous ////////////////////
Interval rohnEig(MatrixXd Ac, MatrixXd Ad, int dim){
    // compute eigenvalue of Ad
    VectorXd evalAc; 
    MatrixXd evecAc;
    getEigen(Ac,evalAc, evecAc); // non-rigourous eigenvalue and eigenvector
    // Maximum eigVal = 
    Interval i_evalAc; // need to compute spectral radius, so all eigenvalue is need
    IVector i_evecAc(dim);

    // only maximum rigorous value is computed: max is in (dim-1)-th index
    if(evalAc[dim-1] != 0)
        verifyeig(Ac, evalAc(dim-1), evecAc.col(dim-1), i_evalAc, i_evecAc, dim);
    else
        i_evalAc = 0.0;
    
    if(_DEBUG_){
        cout << endl << "Ac: " << endl << Ac << endl;
        cout << "Max Eigenvalue of Ac: " << endl;
        cout << "[" << leftBound(i_evalAc) << ", " << rightBound(i_evalAc) << "]" << endl;
    }
    
    // compute eigenvalue of Ad
    VectorXd evalAd; 
    MatrixXd evecAd;
    getEigen(Ad,evalAd, evecAd); // non-rigourous eigenvalue and eigenvector
    // Maximum eigVal = 
    Interval i_evalAd; // need to compute spectral radius, so all eigenvalue is need
    IVector i_evecAd(dim);
    if(evalAd[dim-1] != 0)
        verifyeig(Ad, evalAd(dim-1), evecAd.col(dim-1), i_evalAd, i_evecAd, dim);
    else
        i_evalAd = 0.0;

    if(_DEBUG_) {
        cout << endl << "Ad: " << endl << Ad << endl;
        cout << "Max Eigenvalue of Ad: " << endl;
        for(int i=0; i<dim; i++){
            cout << "[" << leftBound(i_evalAd) << ", " << rightBound(i_evalAd) << "]" << endl;
        }
    }

    // compute spectral of i_evalAd_vec
    // spectral = max([abs(inf(iEigsAd)); abs(sup(iEigsAd))]);
    interval spectral = i_evalAd; // for symmetric matrix it is upper bound of max Eig
    //Apply the theorem
    //ienc = infsup( inf(iEigsAc) - spectral, sup(iEigsAc) + spectral );
    Interval iLam = interval(leftBound(i_evalAc)-rightBound(spectral),rightBound(i_evalAc)+rightBound(spectral));

    if(_DEBUG_) {
        cout << endl << "iA: " << endl << Ad << endl;
        cout << "Max Eigenvalue of iA: " << endl;
        cout << "[" << leftBound(iLam) << ", " << rightBound(iLam) << "]" << endl;
    }
    return iLam;

}

double rohnEigMax(IMatrix iA, int dim) {
    double lamMax;
    //cout << "inside rohnEigMax " << endl;
    DMatrix dAc = takeMidM(iA);
    DMatrix dAd = takeRadM(iA);
    MatrixXd Ac = capd_to_eigen(dAc, dim);
    MatrixXd Ad = capd_to_eigen(dAd, dim);

    if(_DEBUG_){
        cout << endl << "Inside rohnEig: " << endl;
        cout << "Ac: " << endl << Ac << endl;
        cout << "Ad: " << endl << Ad << endl;
    }

    Interval iLam = rohnEig(Ac, Ad, dim);
    lamMax = rightBound(iLam);    
    return lamMax;
}

double combineRohnCfEigMax(IMatrix iA){
    int dim = iA.numberOfRows();
    double lam1 = rohnEigMax(iA, dim);
    double lam2 = cfEigMax(iA, dim);

    return lam1 < lam2? lam1: lam2;
}

///////////////// BEGIN Hladik's Method /////////////////////////////////
int selectIndex(IMatrix iA){
    int dim = iA.numberOfRows();

    if(dim == 1)
        return 0;

    double best = INF;
    IMatrix iB(dim-1,dim-1);
    int bestIdx;
    for(int i=0; i<dim; i++){
        iB = deleteRowCol(iA, i);
        double lam = unverifiedUpper(iB, dim-1);
        if(lam < best){
            bestIdx = i;
            best = lam;
        }
    }
    return bestIdx;
}

//// DIRECT METHOD ///////////////
double directInterlacing(IMatrix iA){
    int dim = iA.numberOfRows();
    IMatrix iB(iA);
    DVector eigsUp(dim);

    //cout << "iA: " << endl;
    //printIMatrix(iA);
    // eigsUp: contain eigVal in descending order
    for(int k=0; k < dim; k++){
        //cout << "k: " << k << endl;
        //cout << "current iB: " << endl;
        //printIMatrix(iB);
        //cout << "*********************" << endl;
        eigsUp[k] = combineRohnCfEigMax(iB); 
        //cout << "max Eig: " << eigsUp[k] << endl;
        int idx = selectIndex(iB);
        //cout << "Best IDX: " << idx << endl << endl;
        iB = deleteRowCol(iB, idx);
        //cout << "deleting idx from iB: " << endl;
        //printIMatrix(iB);
        //cout << endl << endl;
    }
    //eigsUp[dim-1] = combineRohnCfEigMax(iB);
    //cout << "Forward Direction: " << endl;
    //cout << "Max Eig: " << eigsUp[0] << endl << endl;
    vector<int> idxVec(dim);
    IotaWrapper<int> fW(0, [](const int& n){ return n+1;}); // [0, 1, 2, ..., n-1]
    std::iota(idxVec.begin(), idxVec.end(), fW);

    vector<int> newIdxVec;
    //cout << endl << "reverse steps" << endl;
    
    for(int k=0; k<dim; k++){
        //cout << "K = " << k << endl;
        //cout << "***************************" << endl;
        vector<int> remind = setdiff(idxVec, newIdxVec, dim);
        //cout << "remind in hladikDirectEigMax" << endl;
        // cout << "size of remind: " << remind.size() << endl;
        // for(int i: remind){
        //     cout << i << " ";
        // }
        
        IMatrix itemp = getRowsCols(iA, remind);
        //cout << "item: ";
        //printIMatrix(itemp);
        int bestIdx = selectIndex(itemp);
        //cout << "best IDx: " << remind[bestIdx] << endl;
        newIdxVec.push_back(remind[bestIdx]);
        // cout << "current Idx Set: ";
        // for(int i: newIdxVec)
        //     cout << i << " ";
        // //cout << endl << "iA: " << endl;
        //printIMatrix(iA);
        iB = getRowsCols(iA, newIdxVec);
        //printIMatrix(iB);
        
        double newUp = combineRohnCfEigMax(iB);
        //cout << "Max Eig of iB: " << newUp << endl;
        eigsUp[dim-1-k] = std::min(eigsUp[dim-1-k],newUp);
        //cout << endl;

    }
    //cout << "Max Eig of iB: " << eigsUp[0] << endl;

    return eigsUp[0]; // double check whether max is in 0- or (dim-1)-th index
}


double hladikDirectEigMax(IMatrix iA, int dim){
    //int dim = iA.numberOfRows();
    double directEigMax;
    directEigMax = directInterlacing(iA);

    return directEigMax; // send Max
}

//// INDIRECT METHOD /////////////// -- Not complete yet
double indirectInterlacing(IMatrix iA, int dim){ 
    //int dim = iA.numberOfRows();

    DMatrix dAc = takeMidM(iA);
    DMatrix dAd = takeRadM(iA);
    MatrixXd Ac = capd_to_eigen(dAc, dim);
    MatrixXd Ad = capd_to_eigen(dAd, dim); 
    // Compute bound for Ac
    VectorXd evalAc; 
    MatrixXd evecAc;
    getEigen(Ac,evalAc, evecAc); // non-rigourous eigenvalue and eigenvector
    // Maximum eigVal = 
    IVector i_evecAc(dim); // need to compute spectral radius, so all eigenvalue is need
    Interval i_evalAc; // need to compute spectral radius, so all eigenvalue is need
    //i = dim-1;
    //for(int i=dim-1; i>=0; i--){
    if(evalAc[dim-1] != 0)
        verifyeig(Ac, evalAc(dim-1), evecAc.col(dim-1), i_evalAc, i_evecAc, dim);
    else
        i_evalAc = 0.0;
    
    //}

    // compute bound for [-Ad, Ad]
    IMatrix iAd = getIMatrix(-dAd, dAd);
    //cout << "iAd: ";
    //printIMatrix(iAd);
    double evalAd = combineRohnCfEigMax(iAd);

    return rightBound(i_evalAc) + evalAd; // double check whether max is in 0- or (dim-1)-th index
}


double hladikIndirectEigMax(IMatrix iA, int dim){
    //int dim = iA.numberOfRows();
    double indirectEigMax;
    indirectEigMax = indirectInterlacing(iA, dim);

    return indirectEigMax; // send Max
}


// Combine All Methods
double effEigMax(IMatrix iA, int dim){
    double best = INF;

    double lam = rohnEigMax(iA, dim);
    if(lam < best)
        best = lam;
    
    lam = cfEigMax(iA, dim);
    if(lam < best)
        best = lam;

    
    // lam = hladikDirectEigMax(iA, dim);
    // if(lam < best)
    //     best = lam;

    // lam = hladikIndirectEigMax(iA, dim);
    // if(lam < best)
    //     best = lam;
    return best;
}

double hladikDiagMax(IMatrix iA){
    int dim = iA.numberOfRows();
    IMatrix iB(iA);
    for(int i=0; i < dim; i++){
        double diagMax = rightBound(iA[i][i]);
        iB[i][i] = interval(diagMax, diagMax);
    }
    double lam = effEigMax(iA, dim);
    return lam;
}

// CombinedApproach
double computeIntEig(IMatrix iA){

    int dim = iA.numberOfRows();
    double best = INF;
    double lam = effEigMax(iA, dim);
    if(lam < best)
        best = lam;
    lam = hladikDiagMax(iA);
    if(lam < best)
        best = lam;
    return best;

}

///////// CAPD Baed EIG computation //////////////
// Eigenvalue for PD Matrix: only for Dim = 2 at the momemt
void computeEigvPD(const IMatrix& m, interval& eig1, interval& eig2){
    interval delta2 = (m[0][0] + m[1][1]) * (m[0][0] + m[1][1]) - 4 * (m[0][0] * m[1][1] - m[0][1] * m[1][0]);
    interval stretch2;
    
    if(delta2 >= 0.){
        
    }else{
        //we may set delta2 leftBound to 0., as we asume about the set of matrices (Cauchy-Green stress tensor)
        //that they are symmetric.
        delta2.setLeftBound(0.);
    }
    eig1 = ( m[0][0] + m[1][1] - sqrt(delta2) ) / 2. ;
    if(!(eig1 >= 0.))
        eig1.setLeftBound(0.);
    eig2 = ( m[0][0] + m[1][1] + sqrt(delta2) ) / 2. ;
    if(!(eig2 >= 0.))
        eig2.setLeftBound(0.);
}


// Max Eig of Interval Matrix in CAPD: Gershgorin Method
double computeMaxEigPD(const IMatrix& m){
    interval eig1 = capd::matrixAlgorithms::Gershgorin< IMatrix, true >::maxEigenValueOfSymMatrix(m);
    if(!(eig1 >= 0.))
        eig1.setLeftBound(0.);
    return rightBound(eig1);
}
