#include <iostream>
#include <vector>
#include <unordered_map>
#include <set>
#include <fstream>
#include <time.h>
#include <random>
#include <chrono>
#include <algorithm>
#include <queue>
#include <cmath>

#include "chromosome.h"
#include "prg.h"
#include "cluster.h"

using namespace std;

struct GAParams{
    unsigned populationSize;
    unsigned iterationCount;
    double pC;
    double pM;
    vector<vector<unsigned>> matrix;
} params;


PRG randomGenerator;

void Mutate(Chromosome &ch)
{
    unsigned r1 = randomGenerator.nextInt() % ch.data.size(), r2;
    do{
        r2 = randomGenerator.nextInt() % ch.data.size();
    }while(ch.data[r2] == ch.data[r1]);
    swap(ch.data[r1], ch.data[r2]);
}

void Crossover(Chromosome &ch1, Chromosome &ch2){
    unsigned size = ch1.data.size();
    unsigned i =randomGenerator.nextInt() % size;
    if(ch1.data[i] == ch2.data[i]){
        unsigned j = i;
        do{
            i = (i + 1) % size;
        }while(ch1.data[i] == ch2.data[i] && i != j);
        if(i == j)
            return;
    }
    unsigned c = ch1.data[i];
    while(true){
        unsigned temp = ch1.data[i];
        ch1.data[i] = ch2.data[i];
        ch2.data[i] = temp;
        if(ch1.data[i] == c)
            break;
        unsigned j = (i + 1) % size;
        while(ch1.data[i] != ch1.data[j])
            j = (j + 1) % size;
        i = j;
    }
}

void getTreeDetails(const Chromosome &ch, unordered_map<unsigned, Cluster> &TDHCclrs, unordered_map<unsigned, vector<unsigned>> &TDHCchilds){
    TDHCchilds.clear();
    TDHCclrs.clear();
    vector<unsigned> d;
    unsigned n = (ch.data.size() + 3) / 2;
    for(unsigned i = 0; i < n; i++)
    {
        d.push_back(1);
        TDHCclrs[i].items.insert(i);
    }
    for(unsigned i = 0; i < n - 2; i++)
        d.push_back(3);
    d.push_back(2);
    for(unsigned i = 0; i < ch.data.size(); i++){
        unsigned ind;
        for(ind = 0; ind < d.size(); ind++)
            if(d[ind] == 1)
                break;
        TDHCchilds[ch.data[i]].push_back(ind);
        for(set<unsigned>::iterator j = TDHCclrs[ind].items.begin(); j != TDHCclrs[ind].items.end(); j++)
            TDHCclrs[ch.data[i]].items.insert(*j);
        d[ind]--;
        d[ch.data[i]]--;
    }
    unsigned in1, in2;
    for(in1 = 0; in1 < d.size(); in1++)
        if(d[in1] == 1)
            break;
    for(in2 = in1 + 1; in2 < d.size(); in2++)
        if(d[in2] == 1)
            break;
    TDHCchilds[in2].push_back(in1);
    for(set<unsigned>::iterator j = TDHCclrs[in1].items.begin(); j != TDHCclrs[in1].items.end(); j++)
        TDHCclrs[in2].items.insert(*j);
}


void Evaluate(Chromosome &ch){
    unordered_map<unsigned, Cluster> TDHCclrs;
    unordered_map<unsigned, vector<unsigned>> TDHCchilds;
    ch.fitness = 0;
    getTreeDetails(ch, TDHCclrs, TDHCchilds);
    queue<unsigned> q;
    unsigned idx = TDHCclrs.size() - 1, idx1, idx2;
    TDHCclrs[idx].value = 1;
    q.push(idx);
    do{
        idx = q.front();
        q.pop();
        if(TDHCchilds[idx].size() != 2){
            ch.fitness += TDHCclrs[idx].value;
            continue;
        }
        idx1 = TDHCchilds[idx][0];
        idx2 = TDHCchilds[idx][1];
        TDHCclrs[idx1].value = Cluster::exCF(TDHCclrs[idx1], TDHCclrs[idx], params.matrix);
        TDHCclrs[idx2].value = Cluster::exCF(TDHCclrs[idx2], TDHCclrs[idx], params.matrix);
        if(TDHCclrs[idx1].value + TDHCclrs[idx2].value >= TDHCclrs[idx].value){
            q.push(idx1);
            q.push(idx2);
        }
        else{
            ch.fitness += TDHCclrs[idx].value;
        }
    }while(!q.empty());
    ch.fitness--;
}

Chromosome InitPop(vector<Chromosome> &pop){
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    int index = 0;
    double fit = -1;
    unsigned n = params.matrix.size();
    vector<unsigned> temp;
    pop.resize(params.populationSize);
    for(unsigned i = n; i < 2 * n - 2; i++){
        temp.push_back(i);
        temp.push_back(i);
    }
    temp.push_back(2 * n - 2);
    for(unsigned i = 0; i < params.populationSize; i++){
        shuffle(temp.begin(), temp.end(), std::default_random_engine(seed));
        pop[i].data = temp;
        Evaluate(pop[i]);
        if(pop[i].fitness > fit){
            fit = pop[i].fitness;
            index = i;
        }
    }
    return pop[index];
}

vector<Chromosome> RouletteSelectPop(vector<Chromosome> &source) {
    vector<Chromosome> dest;
    unsigned size = source.size();
    double p[size];
    double sum = 0;
    dest.resize(size);
    for(unsigned i = 0; i < size; i++){
        sum += source[i].fitness;
        p[i] = sum;
    }
    for(unsigned i = 0; i < size; i++){
        double r = sum * randomGenerator.nextDouble();
        for(unsigned j = 0; j < size; j++){
            if(r <= p[j]){
                dest[i] = source[j];
                break;
            }
        }
    }
    return dest;
}

Chromosome RunAlgorithm() {
    vector<Chromosome> pop;
    Chromosome best = InitPop(pop);
    vector<double> bests;
    for(unsigned i = 0; i < params.iterationCount; i++) {
        Chromosome b;
        pop = RouletteSelectPop(pop);
        for(unsigned it = 0; it < params.populationSize; it++){
            Chromosome &ch1 = pop[it++];
            Chromosome &ch2 = pop[it];
            if(randomGenerator.nextDouble() < params.pC)
                Crossover(ch1, ch2);
            if(randomGenerator.nextDouble() < params.pM)
                Mutate(ch1);
            if(randomGenerator.nextDouble() < params.pM)
                Mutate(ch2);
            Evaluate(ch1);
            if(ch1.fitness > b.fitness)
                b = ch1;
            Evaluate(ch2);
            if(ch2.fitness > b.fitness)
                b = ch2;
        }
        if(b.fitness > best.fitness)
            best = b;
        bests.push_back(best.fitness);
        pop.push_back(best);
        cout << i + 1 << ": " << best.fitness << endl;
    }
    return best;
}

vector<Cluster> getCluster(const Chromosome &ch){
    vector<Cluster> cl;
    unordered_map<unsigned, Cluster> TDHCclrs;
    unordered_map<unsigned, vector<unsigned>> TDHCchilds;
    getTreeDetails(ch, TDHCclrs, TDHCchilds);
    queue<unsigned> q;
    unsigned idx = TDHCclrs.size() - 1, idx1, idx2;
    TDHCclrs[idx].value = 1;
    idx1 = TDHCchilds[idx][0];
    idx2 = TDHCchilds[idx][1];
    TDHCclrs[idx1].value = Cluster::exCF(TDHCclrs[idx1], TDHCclrs[idx], params.matrix);
    TDHCclrs[idx2].value = Cluster::exCF(TDHCclrs[idx2], TDHCclrs[idx], params.matrix);
    if((TDHCclrs[idx1].value + TDHCclrs[idx2].value) < TDHCclrs[idx].value){
        cl.push_back(TDHCclrs[idx]);
        return cl;
    }
    q.push(idx1);
    q.push(idx2);
    do{
        idx = q.front();
        q.pop();
        if(TDHCchilds[idx].size() != 2){
            cl.push_back(TDHCclrs[idx]);
            continue;
        }
        idx1 = TDHCchilds[idx][0];
        idx2 = TDHCchilds[idx][1];
        TDHCclrs[idx1].value = Cluster::exCF(TDHCclrs[idx1], TDHCclrs[idx], params.matrix);
        TDHCclrs[idx2].value = Cluster::exCF(TDHCclrs[idx2], TDHCclrs[idx], params.matrix);
        if(TDHCclrs[idx1].value + TDHCclrs[idx2].value >= TDHCclrs[idx].value){
            q.push(idx1);
            q.push(idx2);
        }
        else{
            cl.push_back(TDHCclrs[idx]);
        }
    }while(!q.empty());
    return cl;
}

string getTotalDotCode(Chromosome ch, vector<Cluster> best, const vector<string> &names){
    unordered_map<unsigned, Cluster> TDHCclrs;
    unordered_map<unsigned, vector<unsigned>> TDHCchilds;
    getTreeDetails(ch, TDHCclrs, TDHCchilds);
    string res = "graph treeResult\n{\n";
    set<unsigned> tcls;
    stringstream temp;
    queue<unsigned> q;
    set<unsigned>::iterator it;
    for(auto from : TDHCchilds)
        for(auto to : from.second){
            temp << "t" << from.first << " -- t" << to << ";" << endl;
            TDHCclrs[to].value = Cluster::exCF(TDHCclrs[to], TDHCclrs[from.first], params.matrix);
        }
    for(auto cl : TDHCclrs){
        tcls.clear();
        temp << "t" << cl.first << " [label=\"";
        for(it = cl.second.items.begin(); it != cl.second.items.end(); it++)
            tcls.insert(*it);
        it = tcls.begin();
        temp << names[*it];
        for(it++; it != tcls.end(); it++)
            temp << ", " << names[*it];
        temp << "\"";
        if(find(best.begin(), best.end(), cl.second) != best.end())
            temp << " , color=lightblue, fontcolor=black, style=filled";
        temp << "];" << endl;
    }
    for(unsigned i = 0; i < names.size(); i++){
        temp << "c" << i << " [label=\"" << names[i] << "\"];" << endl;
        for(unsigned j = i + 1; j < names.size(); j++)
            if(params.matrix[i][j] == 1)
                temp << "c" << i << " -- c" << j << ";" << endl;
    }
    res += temp.str();
    res += "}";
    return res;
}

int main(int argc, char *argv[])
{
    if (argc < 2)
        return 1;
    string filename(argv[1]);
    ifstream file(filename);
    int ind = 0;
    unordered_map<string, int> classes;
    vector<pair<string, string>> strs;
    vector<vector<unsigned>> matrix;
    vector<string> names;
    set<string> tempNames;
    string from, to;
    while(file >> from >> to){
        tempNames.insert(from);
        tempNames.insert(to);
        strs.push_back(pair<string, string>(from, to));
    }
    file.close();
    for(auto name : tempNames){
        classes[name] = ind++;
        names.push_back(name);
    }
    int n = classes.size();
    if(n == 0)
        throw "Invalid file content.";
    for(int i = 0; i < n; i++){
        matrix.push_back(vector<unsigned>(n));
        fill(matrix[i].begin(), matrix[i].end(), 0);
    }
    for(unsigned i = 0; i < strs.size(); i++){
        matrix[classes[strs[i].first]][classes[strs[i].second]]++;
        matrix[classes[strs[i].second]][classes[strs[i].first]]++;
    }
    params.matrix = matrix;
    params.populationSize = 300 * matrix.size();
    params.iterationCount = 20 * matrix.size();
    params.pM = 0.04 * (log(params.populationSize) / log(2));
    if(params.populationSize <= 1000)
        params.pC = 0.7;
    else if(params.populationSize >= 10000)
        params.pC = 0.9;
    else
        params.pC = 0.7 + (params.populationSize - 1000) / 45000;
    Chromosome result = RunAlgorithm();
    auto c = getCluster(result);
    ofstream output(filename + ".dot");
    output << getTotalDotCode(result, c, names);
    return 0;
}
