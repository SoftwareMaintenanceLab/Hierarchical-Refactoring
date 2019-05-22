#ifndef PRG_H
#define PRG_H

#include <random>
#include <ctime>

class PRG
{
    std::mt19937 _gen;
    unsigned generate();

public:
    PRG(){
        this->_gen.seed(std::time(NULL));
    }

    unsigned nextInt();
    double nextDouble();
    bool nextBoolean();
};

#endif // PRG_H
