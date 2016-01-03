#include "EMAlgorithm.h"

EMAlgorithm::EMAlgorithm() {}

EMAlgorithm::EMAlgorithm(const KmerHash& kmerHasher) {
    K = kmerHasher.K;
    readLength = kmerHasher.readLength;
    NW = kmerHasher.NW;
    NG = kmerHasher.NG;
    NE.clear();
    for(int g = 0; g < NG; g ++) {
        NE.push_back(kmerHasher.NE[g]);
    }
    
    initialIndice();
    initialCoefficients(kmerHasher);
    initialConstraints();
}

EMAlgorithm::~EMAlgorithm() {}

void EMAlgorithm::initialIndice() {
    NX.clear();
    for(int g = 0; g < NG; g ++) {
        NX.push_back(NE[g] * (NE[g] + 1) / 2);
    }
}

void EMAlgorithm::initialCoefficients(const KmerHash& kmerHasher) {
    L.clear();
    for(int g = 0; g < NG; g ++) {
        L.push_back(Eigen::MatrixXd(1, NX[g]));
        int col = 0;
        for(int e = 0; e < NE[g]; e ++) {
            int st = kmerHasher.geneBoundary[g][e].first;
            int ed = kmerHasher.geneBoundary[g][e].second + 1;
            L[g](0, col) = ed - st - readLength + 1;
            ++ col;
        }
        for(int ei = 0; ei < NE[g]; ei ++) {
            for(int ej = ei + 1; ej < NE[g]; ej ++) {
                L[g](0, col) = readLength - 1;
                ++ col;
            }
        }
    }
    cerr << "Length initialized" << endl;
    
    Tau.clear();
    for(int g = 0; g < NG; g ++) {
        Tau.push_back(Eigen::SparseMatrix<double>(NW, NX[g]));
    }
    W.resize(NW, 1);
    int row = 0;
    for(auto& it: kmerHasher.kmerTable) {
        W(row, 0) = it.first;
        for(auto& x: it.second) {
            vector <string> subid;
            boost::split(subid, x.first, boost::is_any_of(","));
            int g = boost::lexical_cast<int> (subid[0]);
            int ei = boost::lexical_cast<int> (subid[1]);
            int ej = boost::lexical_cast<int> (subid[2]);
            int col = 0;
            if(ei == ej) {
                col = ei;
            } else {
                col = (2 * NE[g] - ei) * (ei + 1) / 2 + ej - ei - 1;
            }
            Tau[g].insert(row, col) = x.second;
        }
        ++ row;
    }
    cout << "Tau initalization finished!" << endl;
}

void EMAlgorithm::initialConstraints() {
    NA.clear();
    for(int g = 0; g < NG; g ++) {
        if(NE[g] > 1) {
            NA.push_back(3 * NE[g] - 4);
        } else {
            NA.push_back(NE[g]);
        }
    }

    A.clear();
    for(int g = 0; g < NG; g ++) {
        A.push_back(Eigen::MatrixXd::Zero(NA[g], NX[g]));

        if(NA[g] == 1) {
            A[g](0, 0) = 1;
            break;
        }

        int row = 0;
        
        int NRightJunc = NE[g] - 1;
        int l = NE[g];
        int r = l + NE[g] - 1;
        for(; row < NRightJunc ;row ++) {
            A[g](row, row) = 1;
            for(int i = l; i < r; i ++) {
                A[g](row, i) = -1;
            }
            l = r;
            r = r + NE[g] - (row + 2);
        }

        int NLeftJunc = NRightJunc + NE[g] - 1;
        for(; row < NLeftJunc ; row ++) {
            A[g](row, row - NRightJunc + 1) = 1;
            l = NE[g];
            r = row - NRightJunc;
            for(int i = 1; r >= 0; i ++, r --) {
                A[g](row, l + r) = -1;
                l += NE[g] - i;
            }
        }

        int NPsi = NLeftJunc + NE[g] - 2;
        for(; row < NPsi ; row ++) {
            for(int i = 0; i < NX[g]; i ++) {
                if(i >= NE[g]) {
                    A[g](row, i) = -1;
                } else if(i != row - NLeftJunc + 1) {
                    A[g](row, i) = 1;
                } else {
                    A[g](row, i) = 0;
                }
            }
        }

        for(int i = 0 ; i < NA[g] ; i ++) {
            A[g].row(i) = A[g].row(i).cwiseQuotient(L[g]);
        }
    }
    cout << "Constraints initialization finished!" << endl;
}

void EMAlgorithm::initialVariables() {
    Z.resize(NG, 1);
    Z.setRandom(NG, 1);
    double tot = 0.0;
    for(int g = 0; g < NG; g ++) {
        Z(g, 0) += 1;
        tot += Z(g, 0);
    }
    for(int g = 0; g < NG; g ++) {
        Z(g, 0) /= tot;
    }

    for(int g = 0; g < NG; g ++) {
        X.push_back(Eigen::MatrixXd::Random(NX[g], 1));
        tot = 0.0;
        for(int i = 0; i < NX[g]; i ++) {
            X[g](i, 0) += 1;
            if(i < NE[g]) {
                X[g](i, 0) *= 100;
            }
            tot += X[g](i, 0);
        }
        for(int i = 0; i < NX[g]; i ++) {
            X[g](i, 0) /= tot;
        }
    }
}

void EMAlgorithm::eStep() {
    Mu.clear();
    Eigen::MatrixXd tot(NW, 1);
    Eigen::MatrixXd temp(NW, 1);
    for(int g = 0; g < NG; g ++) {
        Mu.push_back(Eigen::SparseMatrix<double>(NW, 1));
        temp = Tau[g] * X[g] * Z(g, 0);
        tot += temp;
    }
    for(int g = 0; g < NG; g ++) {
        temp = temp.cwiseQuotient(tot);
    }
    for(int g = 0; g < NG; g ++) {
        temp = temp.cwiseProduct(W);
        Mu[g] = temp.sparseView();
    }
}

void EMAlgorithm::mStep(int T) {
    double tot = 0.0;
    for(int g = 0; g < NG; g ++) {
        Z(g, 0) = Mu[g].sum();
        tot += Z(g, 0);
    }
    for(int g = 0; g < NG; g ++) {
        Z(g, 0) /= tot;
    }

    for(int g = 0; g < NG; g ++) {
        optimizeQ(g, T);
    }
}

void EMAlgorithm::optimizeQ(int T, int g) {
    
}

double EMAlgorithm::QFunction(Eigen::MatrixXd X, int g) {
    Eigen::MatrixXd temp = Tau[g] * X;
    temp = temp.array().log();
    for(int w = 0; w < NW; w ++) {
        if(isinf(-temp(w, 0))) {
            temp(w, 0) = 0;
        }
    }
    cout << temp << endl;
    cout << Mu[g].rows() << ' ' << Mu[g].cols() << endl;
    cout << temp.rows() << ' ' << temp.cols() << endl;

    return - (Mu[g].transpose() * temp)(0, 0);
}

Eigen::MatrixXd EMAlgorithm::QDerivate(Eigen::MatrixXd X, int g) {
    Eigen::MatrixXd temp = Tau[g] * X;
    temp = Mu[g].cwiseQuotient(temp);
    Eigen::MatrixXd jac = Tau[g].transpose() * temp;
    jac /= jac.array().sum();
    return - jac;
}

void EMAlgorithm::work(int T) {
    initialVariables();

    //Eigen::MatrixXd prevZ(Z.data());
    T = 1;
    for(int proc = 0; proc < T; proc ++) {
        cout << "\n\n++++++++++" << proc << " iteration processed..." << endl;
        eStep();
        mStep(T);
        cout << QFunction(X[0], 0) << endl;
        //cout << "here?" << endl;
        //cout << QDerivate(X[0], 0) << endl;
    }
    computePSI();
}

void EMAlgorithm::computePSI() {

}
