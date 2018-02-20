//
//  explicit.hpp
//  Assignment4
//
//  Created by Sagar Dolas on 16/01/16.
//  Copyright © 2016 Sagar Dolas. All rights reserved.
//

#ifndef Explicit_hpp
#define Explicit_hpp

#include <stdio.h>
#include <iostream>
#include "ParabolicSolver.h"

class Explicit : public ParabolicSolver
{
    
    public:
    
    Explicit(int, char* []);
    ~Explicit();
    void allocMemory() override;
    virtual void solver() override;

};

#endif /* explicit_hpp */
