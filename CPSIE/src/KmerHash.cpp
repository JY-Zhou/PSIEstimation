#include "KmerHash.h"

KmerHash::KmerHash() {}

KmerHash::KmerHash(int K, string genomePath, string exonBoundaryPath, string readPath) {
    this -> K = K;
    this -> readLength = -1;
    mask = 0;
    for(int i = 0; i < K; i++) {
        mask = (mask << 2) + 3;
    }
    NW = 0;
    NG = 0;
    NE.clear();
    kmerCount.clear();
    kmerTable.clear();
    geneBoundary.clear();

    readReads(readPath);
    readGenome(exonBoundaryPath);
    buildKmerTable(genomePath);
    mergeKmerTable();
}

KmerHash::~KmerHash() {}

inline int KmerHash::parseBP(char c) {
    if(c == 'A' or c == 'a') return 0;
    if(c == 'C' or c == 'c') return 1;
    if(c == 'G' or c == 'g') return 2;
    if(c == 'T' or c == 't') return 3;
    return -1;
}

inline char KmerHash::antiBP(char c) {
    if(c == 'A' or c == 'a') return 'T';
    if(c == 'C' or c == 'c') return 'G';
    if(c == 'G' or c == 'g') return 'C';
    if(c == 'T' or c == 't') return 'A';
    return 'N';
}

inline string KmerHash::getAntiSense(string reads) {
    //return "";
    //string ret = "";
    //for(int i = reads.length() - 1; i >= 0; i --) {
    //    char c = reads[i];
    //    if(c == 'A' or c == 'a') ret += 'T';
    //    if(c == 'C' or c == 'c') ret += 'G';
    //    if(c == 'G' or c == 'g') ret += 'C';
    //    if(c == 'T' or c == 't') ret += 'A';
    //}
    //return ret;
    reverse(reads.begin(), reads.end());
    for(int i = 0; i < reads.length(); i ++) {
        reads[i] = antiBP(reads[i]);
    }
    return reads;
}

inline void KmerHash::readKmer(string reads) {
    long long id = 0;
    int st = 0;
    while(st < K - 1) {
        id = (id << 2) + parseBP(reads[st]);
        st ++;
    }

    while(st < reads.length()) {
        id = ((id << 2) + parseBP(reads[st])) & mask;
        if(kmerCount.count(id) > 0) {
            kmerCount[id] = kmerCount[id] + 1;
        } else {
            kmerCount[id] = 1;
            kmerTable[id] = unordered_map<string, double> ();
        }
        st ++;
    }
}

void KmerHash::readReads(string readPath) {
    vector <string> fileInfo;
    boost::split(fileInfo, readPath, boost::is_any_of("."));
    string ext = fileInfo.back();
    boost::algorithm::to_lower(ext);
    if(ext == "fa" || ext == "fasta") {
        cout << ".FASTA reads file format detected" << endl;
        readReadsFasta(readPath);
    } else if(ext == "fq" || ext == "fastq") {
        cout << ".FASTQ reads file format detected" << endl;
        readReadsFastq(readPath);
    }
}

void KmerHash::readReadsFastq(string readPath) {
    ifstream readsFile(readPath.c_str(), ios::in);
    int proc = 0;
    string line;
    cout << "Processing reads..." << endl;
    while(getline(readsFile, line)) {
        proc ++;
        if(proc % 100000 == 0)
            cout << ">" << flush;

        string identifier = "";
        string reads = "";
        string strand = "";
        string score = "";

        if(line[0] == '@') {
            identifier = boost::trim_copy(line);
            getline(readsFile, line);
            reads = boost::trim_copy(line);
            getline(readsFile, line);
            strand = boost::trim_copy(line);
            getline(readsFile, line);
            score = boost::trim_copy(line);
        }

        if(readLength == -1) {
            readLength = reads.length();
        }

        readKmer(reads);
        readKmer(getAntiSense(reads));
    }
    cout << endl;
    readsFile.close();
    cout << "Total " << proc << " reads." << endl;
    cout << kmerCount.size() << " kmers before merging"<< endl;
    cout << "End of read reads." << endl;
}

void KmerHash::readReadsFasta(string readPath) {
    ifstream readsFile(readPath.c_str(), ios::in);
    int proc = 0;
    string line;
    cout << "Processing reads... (FASTA Format)" << endl;
    while(getline(readsFile, line)) {
        proc ++;
        if(proc % 100000 == 0)
            cout << ">" << flush;

        string identifier = "";
        string reads = "";

        if(line[0] == '>') {
            identifier = boost::trim_copy(line);
            getline(readsFile, line);
            reads = boost::trim_copy(line);
        }

        if(readLength == -1) {
            readLength = reads.length();
        }

        readKmer(reads);
        readKmer(getAntiSense(reads));
    }
    cout << endl;
    readsFile.close();
    cout << "Total " << proc << " reads." << endl;
    cout << kmerCount.size() << " kmers before merging"<< endl;
    cout << "End of read reads." << endl;

}

inline double KmerHash::contribution(int st, int ed, int l, int r, int len) {
    int cil = min(len - readLength + 1, readLength - K + 1);
    int ret = min(l - st + 1, ed - r + 1);
    return min(ret, cil);
}

void KmerHash::readGenome(string exonBoundaryPath) {
    vector <string> fileInfo;
    boost::split(fileInfo, exonBoundaryPath, boost::is_any_of("."));
    string ext = fileInfo.back();
    boost::algorithm::to_lower(ext);
    if(ext == "bed") {
        cout << ".BED annotation file format detect." << endl;
        readFromBed(exonBoundaryPath);
    } else if(ext == "gtf") {
        cout << ".GTF annotation file format detect." << endl;
        readFromGtf(exonBoundaryPath);
    }
}

void KmerHash::readFromGtf(string exonBoundaryPath) {

}

void KmerHash::readFromBed(string exonBoundaryPath) {
    const int COL_CHROMST = 1;
    const int COL_CHROMED = 2;
    const int COL_STRAND = 5;
    const int COL_EXONNUM = 9;
    const int COL_EXONST = 11;
    const int COL_EXONLEN = 10;
    cout << "Read genome ..." << endl;

    ifstream exonBoundaryFile(exonBoundaryPath.c_str(), ios::in);
    string transcript;
    while(getline(exonBoundaryFile, transcript)) {
        transcript = boost::trim_copy(transcript);
        vector <string> subcol;
        vector <string> subst;
        vector <string> sublen;
        boost::split(subcol, transcript, boost::is_any_of("\t"));
        boost::split(subst, subcol[COL_EXONST], boost::is_any_of(","));
        boost::split(sublen, subcol[COL_EXONLEN], boost::is_any_of(","));

        vector <pair <int, int> > exons;
        int exonnum = boost::lexical_cast<int> (subcol[COL_EXONNUM]);
        int chromst = boost::lexical_cast<int> (subcol[COL_CHROMST]);
        for(int i = 0; i < exonnum; i++) {
            int exonst = chromst + boost::lexical_cast<int>(subst[i]);
            int exoned = exonst + boost::lexical_cast<int>(sublen[i]) - 1;
            exons.push_back(make_pair(exonst, exoned));
        }
        geneBoundary.push_back(exons);
    }
    exonBoundaryFile.close();

    NG = geneBoundary.size();
    for(int g = 0; g < NG; g++) {
        NE.push_back(geneBoundary[g].size());
    }
    cout << NG << " genes" << endl;
    //for(int g = 0; g < NG; g++) {
    //    for(int e = 0; e < NE[g]; e++) {
    //        cout << geneBoundary[g][e].first << "\t" << geneBoundary[g][e].second << endl;
    //    }
    //}
}

void KmerHash::buildKmerTable(string genomePath) {
    ifstream genomeFile(genomePath.c_str(), ios::in);
    string geneSeq = "";
    string line;
    while(getline(genomeFile, line)) {
        if(line[0] != '>') {
            geneSeq += boost::trim_copy(line);
        }
    }
    genomeFile.close();

    //int gp = 0;
    //long long gid = 0;
    //for(; gp < K - 1; gp ++) {
    //    gid = (gid << 2) + parseBP(geneSeq[gp]);
    //}
    //for(; gp < geneSeq.size(); gp ++) {
    //    gid = ((gid << 2) + parseBP(geneSeq[gp])) & mask;
    //    kmers.push_back(gid);
    //}

    for(int g = 0; g < NG; g++) {
        //cout << "Gene-" << g << " processed..." << endl;

        for(int e = 0; e < NE[g]; e++) {
            ostringstream oss;
            oss << g << "," << e << "," << e;
            string exonid = oss.str();
            int exonst = geneBoundary[g][e].first;
            int exoned = geneBoundary[g][e].second + 1;
            
            int p = 0;
            long long id = 0;
            double tot = (exoned - exonst - readLength + 1) * (readLength - K + 1);
            
            //cout << "Exon start site: " << exonst << endl;
            //for(int i = exonst ; i < exoned; i++) {
            //    cout << geneSeq[i];
            //} cout << endl;

            for(; p < K - 1; p ++) {
                id = (id << 2) + parseBP(geneSeq[exonst + p]);
            }

            for(; exonst + p < exoned; p ++) {
                id = ((id << 2) + parseBP(geneSeq[exonst + p])) & mask;
                double contr = contribution(0, exoned - exonst, p - K + 1, p + 1, exoned - exonst) / tot;
                if(kmerTable.count(id) > 0) {
                    if(kmerTable[id].count(exonid) > 0) {
                        kmerTable[id][exonid] = kmerTable[id][exonid] + contr;
                    } else {
                        kmerTable[id][exonid] = contr;
                    }
                }
            }
        }
        for(int ei = 0; ei < NE[g]; ei++) {
            for(int ej = ei + 1; ej < NE[g]; ej++) {
                ostringstream oss;
                oss << g << "," << ei << "," << ej;
                string exonid = oss.str();
                int exonedi = geneBoundary[g][ei].second + 1;
                int exonstj = geneBoundary[g][ej].first;

                int p = 0;
                long long id = 0;
                double tot = (readLength - 1) * (readLength - K + 1);
                for(; p < K - 1; p ++) {
                    id = (id << 2) + parseBP(geneSeq[exonedi - readLength + 1 + p]);
                }
                for(; p < readLength - 1; p ++) {
                    id = ((id << 2) + parseBP(geneSeq[exonedi - readLength + 1 + p])) & mask;
                    double contr = contribution(0, 2*readLength - 2, p - K + 1, p + 1, 2*readLength - 2) / tot;
                    if(kmerTable.count(id) > 0) {
                        if(kmerTable[id].count(exonid) > 0) {
                            kmerTable[id][exonid] = kmerTable[id][exonid] + contr;
                        } else {
                            kmerTable[id][exonid] = contr;
                        }
                    }
                }
                for(p = 0; p < readLength - 1; p ++) {
                    id = ((id << 2) + parseBP(geneSeq[exonstj + p])) & mask;
                    double contr = contribution(0, 2*readLength - 2, readLength + p - K, readLength + p, 2*readLength - 2) / tot;
                    if(kmerTable.count(id) > 0) {
                        if(kmerTable[id].count(exonid) > 0) {
                            kmerTable[id][exonid] = kmerTable[id][exonid] + contr;
                        } else {
                            kmerTable[id][exonid] = contr;
                        }
                    }
                }
            }
        }
    }

    cout << "End of building kmer table" << endl;
}

void KmerHash::mergeKmerTable() {
    unordered_map<string, long long> temp;
    auto it = kmerTable.begin();
    while(it != kmerTable.end()) {
        long long id = it -> first;
        ostringstream oss;
        vector <string> keys;
        for(auto x: it -> second) {
            keys.push_back(x.first);
        }
        sort(keys.begin(), keys.end());
        for(auto x: keys) {
            oss << x << (int)(it -> second[x] * 10000000);
        }
        string value = oss.str();
        if(temp.count(value) > 0) {
            long long presentid = temp[value];
            kmerCount[presentid] += kmerCount[id];
            for(auto x: keys) {
                kmerTable[presentid][x] += it -> second[x];
            }
            kmerCount.erase(id);
            it = kmerTable.erase(it);
        } else {
            temp[value] = id;
            it ++;
        }
    }
    NW = temp.size();
}
