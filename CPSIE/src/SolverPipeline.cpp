#include <iostream>
#include <algorithm>
#include "KmerHash.h"
#include "EMAlgorithm.h"

using namespace std;

int main() {
    cout << "Compile passed!" << endl;
    string readPath = "../input/reads.fq";
    string genomePath = "../input/genome.fa";
    string exonBoundaryPath = "../input/exonBoundary.bed";
    int K = 15;
    int readLength = 75;
    KmerHash kmerHasher(K, readLength, genomePath, exonBoundaryPath, readPath);
    cout << kmerHasher.kmerCount.size() << " kmers" << endl;
    //cout << kmerHasher.kmers.size() << endl;
    //for(int i = 0; i < kmerHasher.kmers.size(); i++) {
    //    long long id = kmerHasher.kmers[i];
    //    string mer = "";
    //    for(int j = 0; j < K; j++) {
    //        int c = id & 3;
    //        if(c == 0) mer = "A" + mer;
    //        if(c == 1) mer = "C" + mer;
    //        if(c == 2) mer = "G" + mer;
    //        if(c == 3) mer = "T" + mer;
    //        id >>= 2;
    //    }
    //    id = kmerHasher.kmers[i];
    //    cout << mer << endl;
    //    cout << "#" << i << "\t";
    //    cout << kmerHasher.kmerCount[id] << "\t";
    //    vector <string> temp;
    //    for(auto& x: kmerHasher.kmerTable[id]) {
    //        temp.push_back(x.first);
    //    }
    //    sort(temp.begin(), temp.end());
    //    for(int j = 0; j < temp.size(); j++) {
    //        //cout << x.first << ":" << (int)(x.second*100000) << "\t";
    //        cout << temp[j] << ":" << (int)(kmerHasher.kmerTable[id][temp[j]] * 100000) << "\t";
    //    }
    //    cout << endl;
    //}
    EMAlgorithm solver(kmerHasher);
    solver.work(10);
    cout << "Run passed!" << endl;
    return 0;
}
