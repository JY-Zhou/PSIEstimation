#ifndef KMERHASH_H
#define KMERHASH_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <ctime>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;

class KmerHash {
    public:
        int K;
        int readLength;
        long long mask;
        long long NW;
        int NG;
        vector <int> NE;
        //vector <long long> kmers;

        unordered_map <long long, long long> kmerCount;
        unordered_map <long long, unordered_map<string, double> > kmerTable;
        vector <vector<pair<int, int> > > geneBoundary;

        KmerHash();
        KmerHash(int, string, string, string);
        ~KmerHash();
        
        void readReads(string);
        void readReadsFasta(string);
        void readReadsFastq(string);
        void readGenome(string);
        void readFromGtf(string);
        void readFromBed(string);
        void buildKmerTable(string);
        void mergeKmerTable();

        inline int parseBP(char);
        inline char antiBP(char);
        inline void readKmer(string);
        inline string getAntiSense(string);
        inline double contribution(int, int, int, int, int);
};

#endif
