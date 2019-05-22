#include "cluster.h"

Cluster::Cluster(){
    this->value = 0;
}

Cluster::Cluster(unsigned i){
    this->items.insert(i);
    this->value = 0;
}


Cluster::Cluster(const set<unsigned> &c){
    this->items = c;
}

Cluster::Cluster(const vector<unsigned> &c){
    for(auto cl : c)
        this->items.insert(cl);
}

bool Cluster::operator == (const Cluster &r) const{
    return r.items == this->items;
}

void Cluster::generateRSF(ofstream &output, const vector<Cluster> &clusters, const vector<string> &names){
    for(unsigned i = 0; i < clusters.size(); i++)
        for(auto j : clusters[i].items)
            output << "contain " << i << " " << names[j] << endl;
}

string Cluster::getBestClusterDotCode(const vector<string> &names, const vector<string> &oenames, const vector<vector<unsigned>> &oasymatrix, const vector<Cluster> &clusters) {
    string res = "digraph bestClustering\n{\n";
    stringstream temp;
    vector<int> cls(names.size());
    int matrix[clusters.size()][clusters.size()];
    for(unsigned i = 0; i < clusters.size(); i++)
        for(unsigned j = 0; j < clusters.size(); j++)
            matrix[i][j] = 0;
    for(unsigned i = 0; i < clusters.size(); i++) {
        temp << "subgraph cluster" << i << "{\nlabel=\"Cluster " << i + 1 << "\";" << endl;
        temp << "style = bold;" << endl;
        for(auto j : clusters[i].items) {
            cls[j] = i;
            temp << j << " [label=\"";
            temp << names[j];
            temp << "\"];" << endl;
        }
        temp << "}" << endl;
    }
    if(names[0] != oenames[0]){
        for(unsigned i = 0; i < clusters.size(); i++){
            temp << "subgraph cluster" << clusters.size() + i << "{\nlabel=\"Cluster " << i + 1 << "\";" << endl;
            temp << "style = bold;" << endl;
            for(auto j : clusters[i].items){
                temp << "e" << j << " [label=\"";
                temp << oenames[j];
                temp << "\"];" << endl;
            }
            temp << "}" << endl;
        }
        for(unsigned i = 0; i < oasymatrix.size(); i++)
            for(unsigned j = 0; j < oasymatrix[i].size(); j++)
                if(oasymatrix[i][j] > 0){
                    temp << "e" << i << " -> e" << j;
                    temp << " [label=\"" << oasymatrix[i][j] << "\"]";
                    temp << ";" << endl;
                }
    }
    for(unsigned i = 0; i < clusters.size(); i++){
        temp << "C" << i << " [label=\"";
        temp << i + 1;
        temp << "\"];" << endl;
    }
    for(unsigned i = 0; i < oasymatrix.size(); i++)
        for(unsigned j = 0; j < oasymatrix[i].size(); j++)
            if(oasymatrix[i][j] > 0){
                matrix[cls[i]][cls[j]] += oasymatrix[i][j];
                temp << i << " -> " << j;
                temp << " [label=\"" << oasymatrix[i][j] << "\"]";
                temp << ";" << endl;
            }
    for(unsigned i = 0; i < clusters.size(); i++)
        for(unsigned j = 0; j < clusters.size(); j++)
            if(matrix[i][j] > 0){
                temp << "C" << i << " -> C" << j;
                temp << " [label=\"" << matrix[i][j] << "\"]";
                temp << ";" << endl;
            }
    auto ttttt = temp.str();
    res += temp.str();
    res += "}";
    return res;
}
double Cluster::exCF(const Cluster &cl, const Cluster &clp, const vector<vector<unsigned>> &matrix){
    double inner = 0, ios, in, outer = 0, out, io = 0;
    double cf = 0;
    if(cl.items.size() == 0)
        return 0;
    for(set<unsigned>::iterator i = cl.items.begin(); i != cl.items.end(); i++){
        in = out = ios = 0;
        for(unsigned j = 0; j < matrix.size(); j++)
            if(matrix[*i][j] > 0){
                if(cl.items.find(j) != cl.items.end())
                    in += matrix[*i][j];
                else if(clp.items.find(j) != clp.items.end())
                    ios += matrix[*i][j];
                else
                    out += matrix[*i][j];
            }
        inner += in;
        outer += out;
        io += ios;
    }
    if(outer > io)
        return -1;
    if(inner + io + outer != 0)
        cf = (inner / (inner + 2 * io + outer));
    return cf;
}
