//
//  Implicit.hpp
//  Assignment4
//
//  Created by Sagar Dolas on 16/01/16.
//  Copyright © 2016 Sagar Dolas. All rights reserved.
//

#ifndef Implicit_hpp
#define Implicit_hpp

#include <stdio.h>
#include "ParabolicSolver.h"
#include <numeric>

class Implicit  : public ParabolicSolver{
    
public:
    
    std::vector<double> d_;
    std::vector<double> r_;
    std::vector<double> z_;
    
    Implicit(int argc, char *argv[]);
    Implicit();
    ~Implicit() ;
    void allocMemory() override;
    virtual void solver() override;
    //void computeCgIter();
};

#endif /* Implicit_hpp */
