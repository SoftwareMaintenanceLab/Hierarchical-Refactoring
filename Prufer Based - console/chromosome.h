#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include<vector>

using namespace std;

class Chromosome
{
public:
    vector<unsigned> data;
    double fitness;

    Chromosome(){
        this->fitness = -1;
    }

    Chromosome(const Chromosome &ch){
        this->fitness = ch.fitness;
        this->data = ch.data;
    }

    Chromosome(const vector<unsigned> &d){
        this->data = d;
    }
};
#endif // CHROMOSOME_H
