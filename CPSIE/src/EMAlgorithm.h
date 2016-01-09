#ifndef EMALGORITHM_H
#define EMALGORITHM_H

#include <iostream>
#include <vector>
#include <thread>
#include <ctime>
#include <cstdlib>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "KmerHash.h"

using namespace std;

extern "C" void slsqp_(int* M, int* MEQ, int* LA, int* N,
                        double* X, double* XL, double* XU, 
                        double* F, double* C, double* G, double* A,
                        double* ACC, int* ITER, int* MODE,
                        double* W, int* L_W, int* JW, int* L_JW);

class EMAlgorithm {
    public:
        const double EPS = 1e-20;
        int K;
        int readLength;

        long long NW;
        int NG;
        vector <int> NE;
        vector <int> NX;
        vector <int> NA;

        vector <Eigen::MatrixXd> L;
        Eigen::MatrixXd W;
        vector <Eigen::SparseMatrix <double> > Tau;
        vector <Eigen::MatrixXd> A;

        Eigen::MatrixXd Z;
        vector <Eigen::MatrixXd> X;
        vector <Eigen::MatrixXd> Mu;

        vector <vector <double> > PSI;

        EMAlgorithm();
        EMAlgorithm(KmerHash&);
        ~EMAlgorithm();

        void initialIndice();
        void initialCoefficients(KmerHash&);
        void initialConstraints();
        Eigen::MatrixXd initialX(int);
        void initialVariables();

        void eStep();
        void mStep();
        void optimizeQ(int);
        double QFunction(Eigen::MatrixXd, int);
        Eigen::MatrixXd QGradient(Eigen::MatrixXd, int);
        Eigen::MatrixXd QConstraints(Eigen::MatrixXd, int);
        Eigen::MatrixXd QConstraintsNormal(Eigen::MatrixXd, int);
        void work(int);

        void computePSI();
};

#endif
