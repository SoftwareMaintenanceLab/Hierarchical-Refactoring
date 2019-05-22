#include "prg.h"

unsigned PRG::generate(){
    return this->_gen();
}

unsigned PRG::nextInt(){
    return this->generate();
}

double PRG::nextDouble(){
    double v = this->generate() / (this->_gen.max() * 1.0);
    return v;
}

bool PRG::nextBoolean(){
    return this->generate() % 2;
}
