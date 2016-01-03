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

extern "C" {
    void slsqp(const int* M, const int* MEQ);
}

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

        EMAlgorithm();
        EMAlgorithm(const KmerHash&);
        ~EMAlgorithm();

        void initialIndice();
        void initialCoefficients(const KmerHash&);
        void initialConstraints();

        void initialVariables();
        void eStep();
        void mStep(int);
        void optimizeQ(int, int);
        double QFunction(Eigen::MatrixXd, int);
        Eigen::MatrixXd QDerivate(Eigen::MatrixXd, int);
        void work(int);

        void computePSI();
};

#endif
