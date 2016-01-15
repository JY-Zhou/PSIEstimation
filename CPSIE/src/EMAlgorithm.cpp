#include "EMAlgorithm.h"
#define EIGEN_DONT_PARALLELIZE

EMAlgorithm::EMAlgorithm() {}

EMAlgorithm::EMAlgorithm(KmerHash& kmerHasher) {
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

void EMAlgorithm::initialCoefficients(KmerHash& kmerHasher) {
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
    cout << "Length initialized" << endl;
    
    Tau.clear();
    for(int g = 0; g < NG; g ++) {
        Tau.push_back(Eigen::SparseMatrix<double>(NW, NX[g]));
    }
    W.resize(NW, 1);
    int row = 0;
    for(auto it: kmerHasher.kmerTable) {
        W(row, 0) = kmerHasher.kmerCount[it.first];
        for(auto x: it.second) {
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

Eigen::MatrixXd EMAlgorithm::initialX(int g) {
    Eigen::MatrixXd ret(NX[g], 1);
    //ret.setRandom(NX[g], 1);
    //int tot = 0;
    //for(int i = 0; i < NX[g]; i ++) {
    //    ret(i, 0) += 1;
    //    if(i < NE[g]) {
    //        ret(i, 0) *= 10;
    //    }
    //    tot += ret(i, 0);
    //}
    //for(int i = 0; i < NX[g]; i ++) {
    //    ret(i, 0) /= tot;
    //}

    //while((A[g] * ret).minCoeff() < -EPS) {
    //    ret.setRandom(NX[g], 1);
    //    tot = 0;
    //    for(int i = 0; i < NX[g]; i ++) {
    //        ret(i, 0) += 1;
    //        if(i < NE[g]) {
    //            ret(i, 0) *= 10;
    //        }
    //        tot += ret(i, 0);
    //    }
    //    for(int i = 0; i < NX[g]; i ++) {
    //        ret(i, 0) /= tot;
    //    }
    //}
    for(int e = 0; e < NE[g]; e ++) {
        ret(e, 0) = 0.8 / NE[g];
    }
    for(int e = NE[g]; e < NX[g]; e ++) {
        ret(e, 0) = 0.2 / (NX[g] - NE[g]);
    }

    return ret;
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
        X.push_back(initialX(g));
    }
}

void EMAlgorithm::eStep() {
    Mu.clear();
    Eigen::MatrixXd tot(NW, 1);
    tot.setZero();
    Eigen::MatrixXd temp(NW, 1);
    for(int g = 0; g < NG; g ++) {
        Mu.push_back(Eigen::SparseMatrix<double>(NW, 1));
        temp = Tau[g] * X[g] * Z(g, 0);
        tot += temp;
    }
    for(int g = 0; g < NG; g ++) {
        temp = Tau[g] * X[g] * Z(g, 0);
        temp = temp.cwiseQuotient(tot);
        temp = temp.cwiseProduct(W);
        for(int w = 0; w < NW; w ++) {
            if(std::isnan(temp(w, 0))) {
                temp(w, 0) = 0;
            }
        }
        Mu[g] = temp.sparseView();
    }

    //for(int i = 0; i < Mu[0].rows() ; i ++) {
    //    cout << Mu[0](i, 0) << endl;
    //    if(i % 50 == 0) {
    //        getchar();
    //    }
    //}
    //cout << "Look at me !!!! " << Mu[0].array().sum() << endl;
    //cout << W.array().sum() << endl;
    //cout << Mu[0].size() << ' ' << Mu[0].nonZeros() << endl;
}

void EMAlgorithm::mStep() {
    double tot = 0.0;
    for(int g = 0; g < NG; g ++) {
        Z(g, 0) = Mu[g].sum();
        tot += Z(g, 0);
    }
    for(int g = 0; g < NG; g ++) {
        Z(g, 0) /= tot;
    }

    //#pragma omp parallel for schedule(dynamic)
    for(int g = 0; g < NG; g ++) {
        optimizeQ(g);
    }
    //thread task[8];
    //for(int i = 0 ; i < 8 ; i ++) {
    //    task[i] = thread([i] {cout <<">>>>Thread " << i << " run.." << endl;});
    //}
    //int g = 0; 
    //while(g < NG) {
    //    task[g % 8].join();
    //    task[g % 8] = thread([this, g] {optimizeQ(g);} );
    //    g ++;
    //}
    //for(int i = 0 ; i < 8 ; i ++) {
    //    task[i].join();
    //}
}

void EMAlgorithm::optimizeQ(int g) {
    //cout << "\nGene-" << g << " started ...." << endl;
    int f_M = NA[g] + 1;
    int f_MEQ = 1;
    int f_LA = f_M;
    int f_N = NX[g];
    Eigen::MatrixXd f_X = initialX(g);
    //Eigen::MatrixXd f_X = X[g];
    Eigen::MatrixXd f_XL;
    Eigen::MatrixXd f_XU;
    f_XL.setZero(f_N, 1);
    f_XU.setOnes(f_N, 1);
    double f_F = QFunction(f_X, g);
    Eigen::MatrixXd f_C = QConstraints(f_X, g);
    Eigen::MatrixXd f_G = QGradient(f_X, g);
    Eigen::MatrixXd f_A = QConstraintsNormal(f_X, g);
    double f_ACC = 1e-9;
    int f_ITER = 100;
    int f_MODE = 0;
    int f_N1 = f_N + 1;
    int f_MINEQ = f_M - f_MEQ + 2 * f_N1;
    int f_L_W = (3 * f_N1 + f_M) * (f_N1 + 1)
              + (f_N1 - f_MEQ + 1) * (f_MINEQ + 2) + 2 * f_MINEQ
              + (f_N1 + f_MINEQ) * (f_N1 - f_MEQ) + 2 * f_MEQ + f_N1
              + f_N1 * f_N / 2 + 2 * f_M + 3 * f_N + 3 * f_N1 + 1;
    Eigen::MatrixXd f_W(f_L_W, 1);
    f_W.setZero(f_L_W, 1);
    int f_L_JW = f_MINEQ;
    Eigen::MatrixXi f_JW(f_L_JW, 1);
    f_JW.setZero(f_L_JW, 1);
    int t = 0;
    
    //ostringstream log;

    //log << "Now at gene " << g << endl;

    while(true) {
        t ++;
        //log << g << ": " <<  "-Monitor::::"<< t << ": " <<  f_X.sum() << endl;
        //log << g << ": " <<  "==" << endl;
        //log << g << ": " <<  "MODE " << f_MODE << endl;
        if(f_MODE == 0 || f_MODE == 1) {
            f_F = QFunction(f_X, g);
            f_C = QConstraints(f_X, g);
        }

        if(f_MODE == 0 || f_MODE == -1) {
            f_G = QGradient(f_X, g);
            f_A = QConstraintsNormal(f_X, g);
        }

        //log << g << ": " <<  "X==" << endl;
        //log << g << ": " <<  f_X << endl;
        //log << g << ": " <<  "F==" << endl;
        //log << g << ": " <<  f_F << endl;
        //log << g << ": " <<  "C==" << endl;
        //log << g << ": " <<  f_C << endl;
        //log << g << ": " <<  "G==" << endl;
        //log << g << ": " <<  f_G << endl;
        //log << g << ": " <<  "A==" << endl;
        //log << g << ": " <<  f_A << endl;

        slsqp_(&f_M, &f_MEQ, &f_LA, &f_N,
              f_X.data(), f_XL.data(), f_XU.data(),
              &f_F, f_C.data(), f_G.data(), f_A.data(),
              &f_ACC, &f_ITER, &f_MODE,
              f_W.data(), &f_L_W, f_JW.data(), &f_L_JW);

        if(f_MODE != -1 and f_MODE != 1) {
            //log << g << ": " <<  f_MODE << endl;
            break;
        }
    }
    //cout << "iter time " << t << endl;
    //cout << "x sumation " << f_X.sum() << endl;
    //
    //cout << "Optimal is " << f_F << endl;
    //cout << "Exit Mode is " << f_MODE << endl;
    //if(std::isnan(f_F) or std::isnan(f_X.sum())) {
    //    cerr << log.str() << endl;
    //} else {
    //    log.str("");
    //}
    X[g] = f_X;

    return;
}

double EMAlgorithm::QFunction(Eigen::MatrixXd X, int g) {
    ostringstream log;
    Eigen::MatrixXd temp = Tau[g] * X;
    temp = temp.array().log();
    int i = 0;
    for(int w = 0; w < NW; w ++) {
        if(std::isinf(-temp(w, 0)) || std::isinf(temp(w, 0))) {
            i ++;
            temp(w, 0) = 0;
        }
        if(std::isnan(temp(w, 0))) {
            temp(w, 0) = 0;
        }
    }
    //cout << "The number of inf is: " << i << endl;
    //cout << Mu[g].rows() << ' ' << Mu[g].cols() << endl;
    //cout << temp.rows() << ' ' << temp.cols() << endl;

    double ret = - (Mu[g].transpose() * temp)(0, 0);
    log << "Mu[g].sum() = " << Mu[g].sum() << endl;
    log << "temp.sum() = " << temp.sum() << endl;
    if(std::isnan(ret)) {
        cerr << log.str() << endl;
    }
    return ret;
}

Eigen::MatrixXd EMAlgorithm::QGradient(Eigen::MatrixXd X, int g) {
    //cout << "++++++++++++++calGrad ++++++++" << endl;
    //cout << "X\n" << endl;
    //cout << "[";
    for(int i = 0; i < NX[g]; i ++) {
        //cout << X(i, 0) << (i == NX[g] - 1 ? "]\n" : ", ");
    }
    //cout << endl;
    Eigen::MatrixXd temp = Tau[g] * X;
    temp = Mu[g].cwiseQuotient(temp);
    for(int w = 0 ; w < NW; w ++) {
        if(std::isnan(temp(w, 0))) {
            temp(w, 0) = 0;
            ////cout << "temp NAN at " << w << endl;
        }
        if(std::isinf(temp(w, 0)) || std::isinf(-temp(w, 0))) {
            //cout << "temp INF at " << w << endl;
            double s = 0;
            for(int i = 0; i < NX[g]; i ++) {
                //cout << Tau[g].coeffRef(w, i) << " * " << X(i, 0) << " = " << Tau[g].coeffRef(w, i) * X(i, 0) << endl;
                s += Tau[g].coeffRef(w, i) * X(i, 0);
            }
            //cout << "total " << s << endl;
            temp(w, 0) = 0;
        }
    }
    Eigen::MatrixXd jac = Tau[g].transpose() * temp;

    //cout << "jac before " << jac << endl;
    //cout << "jac summation" << jac.sum() << endl;
    jac /= jac.sum();
    //cout << "jac " << jac << endl;
    ////cout << "jac size " << jac.rows() << ' ' << jac.cols() << endl;
    ////cout << "Summation is: " <<  jac.array().sum() << endl;
    Eigen::MatrixXd ret(NX[g] + 1, 1);
    ret.setZero(NX[g] + 1, 1);
    ret.block(0, 0, NX[g], 1) = -jac;
    return ret;
}

Eigen::MatrixXd EMAlgorithm::QConstraints(Eigen::MatrixXd X, int g) {
    Eigen::MatrixXd cons(NA[g] + 1, 1);
    cons(0, 0) = X.sum() - 1;
    cons.block(1, 0, NA[g], 1) = A[g] * X;
    //for(int i = 0 ; i < NA[g] + 1; i ++ ) {
    //    if(!(cons(i, 0) > -EPS)) {
    //        cout << "FALSE" << endl;
    //        cout << cons(i, 0) << endl;
    //    }
    //}
    return cons;
}

Eigen::MatrixXd EMAlgorithm::QConstraintsNormal(Eigen::MatrixXd X, int g) {
    Eigen::MatrixXd normals(NA[g] + 1, NX[g] + 1);
    normals.setZero(NA[g] + 1, NX[g] + 1);
    normals.block(0, 0, 1, NX[g]) = Eigen::MatrixXd::Ones(1, NX[g]);
    normals.block(1, 0, NA[g], NX[g]) = A[g];
    //int cnt = 0;
    //for(int i = 0; i < NA[g] + 1; i ++) {
    //    for(int j = 0; j < NX[g] + 1; j ++) {
    //        normals(i, j) = cnt ++;
    //    }
    //}
    //return normals.transpose();
    return normals;
}

void EMAlgorithm::work(int T) {
    initialVariables();
    
    for(int g = 0; g < NG; g ++) {
        Z(g, 0) = 1.0 / NG;
        for(int e = 0; e < NE[g]; e ++) {
            X[g](e, 0) = 0.8 / NE[g];
        }
        for(int e = NE[g]; e < NX[g]; e ++) {
            X[g](e, 0) = 0.2 / (NX[g] - NE[g]);
        }
    }

    Eigen::MatrixXd preZ = Z;
    int proc = 0;
    //while(true) {
    while(proc < 2) {
        cout << "\n\n++++++++++" << proc ++ << " iteration processed..." << endl;
        //cout << Z << endl;
        cout << "E-step...." << endl;
        eStep();
        cout << "M-step...." << endl;
        mStep();
        //cout << "==================Debug==============" << endl;
        //cout << X[0] << endl << endl;
        //cout << A[0] * X[0] << endl << endl;
        //cout << X[0].array().sum() << endl;

        //cout << "Tau[0] " << Eigen::MatrixXd(Tau[0]).array().sum() << endl;
        //cout << "mu" << endl;
        //cout << Mu[0].rows() << " " << Mu[0].cols() << endl;
        //cout << Mu[0].array().sum() << endl;
        //int s = 0;
        //for(int i = 0; i < Mu[0].rows(); i ++) {
        //    if(Mu[0](i,0) != 0) s ++;
        //}
        //cout << s << endl;

        //cout << "==================QFunc==============" << endl;
        //cout << "QFunc(0) = " << QFunction(X[0], 0) << endl;
        //cout << "==================QGrad==============" << endl;
        //cout << "QGrad(0) = [\n" << QDerivate(X[0], 0) << "\n]" << endl;
        //cout << "==================Debug==============" << endl;
        //cout << preZ << endl;
        //cout << "====" << endl;
        //
        //cout << Z << endl;
        if(((preZ - Z).array() < 1e-5).all()) {
            break;
        } else {
            preZ = Z;
        }
        //cout << Z << endl;
    }
    computePSI();
}

void EMAlgorithm::computePSI() {
    PSI.clear();
    for(int g = 0; g < NG; g ++) {
        Eigen::MatrixXd tempPSI = X[g].cwiseQuotient(L[g].transpose());
        double sumExon = EPS;
        double sumJunction = EPS;
        //cout << tempPSI << endl;
        //cout << X[g] << endl;
        //cout << "Summation is " << X[g].sum() << endl;
        //cout << "================" << endl;
        //cout << QFunction(X[g], g) << endl;
        //cout << "================" << endl;
        //cout << QGradient(X[g], g) << endl;
        //cout << "================" << endl;
        //cout << QConstraints(X[g], g) << endl;
        //cout << "================" << endl;
        //cout << QConstraintsNormal(X[g], g) << endl;
        //cout << "================" << endl;
        //Eigen::MatrixXd a(2,3);
        //a << 1, 2, 3,
        //    4, 5, 6;
        //cout << a << endl;
        //double* pnt = a.data();
        //for(int i = 0 ; i < 9 ; i++) {
        //    cout << *pnt << " ";
        //    pnt ++;
        //}
        //cout << endl;
        //getchar();
        for(int i = 0; i < NX[g]; i ++) {
            if(i < NE[g]) {
                sumExon += tempPSI(i, 0);
            } else {
                sumJunction += tempPSI(i, 0);
            }
        }
        tempPSI /= (sumExon - sumJunction);
        //cout << sumExon << endl;
        //cout << sumJunction << endl;
        //getchar();
        vector <double> temp;
        for(int e = 0; e < NE[g]; e ++) {
            temp.push_back(tempPSI(e, 0));
        }
        PSI.push_back(temp);
    }
    for(int g = 0; g < NG; g ++) {
        cout << "[ ";
        for(int e = 0; e < NE[g]; e ++) {
            cout << PSI[g][e] << ", ";
        }
        cout << "]" << endl;
    }
}
