#include <iostream>
#include <algorithm>
#include <ctime>
#include <cstring>
#include <boost/lexical_cast.hpp>
#include "KmerHash.h"
#include "EMAlgorithm.h"

using namespace std;

void work(string genomePath, string exonBoundaryPath, string readPath, int K, string outputPath) {
    time_t st, ed;
    time(&st);
    KmerHash kmerHasher(K, genomePath, exonBoundaryPath, readPath);
    time(&ed);
    cout << "Elasped time " << difftime(ed, st) << "s. " << endl;

    time(&st);
    EMAlgorithm solver(kmerHasher); 
    solver.work(10, outputPath);
    time(&ed);
    cout << "Elasped time " << difftime(ed, st) << "s. " << endl;
    cout << "=== Finished! ===" << endl;
}

void showHelp() {
    cout << "Allowed options" << endl;
    cout << endl;
    cout << "\t-g <FASTA file>\tThe reference genome." << endl;
    cout << endl;
    cout << "\t-a <BED/GTF file>\tThe annotation of exons." << endl;
    cout << endl;
    cout << "\t-r <FASTA/FASTQ file>\tThe read sequences." << endl;
    cout << endl;
    cout << "\t-k <Integer>\tThe length of kmer when computing. (Default: 25)" << endl;
    cout << endl;
    cout << "\t-o <JSON file>\tThe result of the exon inclusion ratio." << endl;
    cout << endl;
}

int main(int argc, char** argv) {
    //ios_base::sync_with_stdio(false);
    //cin.tie(0);

    if(argc != 11) {
        showHelp();
        return 0;
    }

    string genomePath = "";
    string exonBoundaryPath = "";
    string readPath = "";
    int K = 25;
    string outputPath = "";

    for(int i = 1; i < argc; i ++) {
        if(strcmp(argv[i], "-g") == 0) {
            genomePath = argv[++i];
        } else if(strcmp(argv[i], "-a") == 0) {
            exonBoundaryPath = argv[++i];
        } else if(strcmp(argv[i], "-r") == 0) {
            readPath = argv[++i];
        } else if(strcmp(argv[i], "-k") == 0) {
            K = boost::lexical_cast<int>(argv[++i]);
        } else if(strcmp(argv[i], "-o") == 0) {
            outputPath = argv[++i];
        }
    }

    long st = clock();
    time_t n_st, n_ed;
    time(&n_st);
    work(genomePath, exonBoundaryPath, readPath, K, outputPath);
    time(&n_ed);
    cout << "CPU Time elapsed: " << (clock() - st) / CLOCKS_PER_SEC << "s. "<< endl;

    cout << "Natural Time elapsed: " << difftime(n_ed, n_st) << "s. " << endl;

    return 0;
}
