#ifndef CLUSTER_H
#define CLUSTER_H

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <sstream>
#include <iomanip>
#include <set>
#include <cmath>
#include <fstream>

using namespace std;

class Cluster
{
public:
    double value;
    set<unsigned> items;

    Cluster();
    Cluster(unsigned);
    Cluster(const set<unsigned> &);
    Cluster(const vector<unsigned> &);

    bool operator == (const Cluster &) const;

    static void getClusterInfo(const Cluster &, const Cluster &, int&, int&, int&);
    static double exCF(const Cluster &, const Cluster &, const vector<vector<unsigned>> &);
    static string getBestClusterDotCode(const vector<string> &, const vector<string> &, const vector<vector<unsigned>> &, const vector<Cluster> &);
    static void generateRSF(ofstream&, const vector<Cluster> &, const vector<string> &);
};

#endif // CLUSTER_H
