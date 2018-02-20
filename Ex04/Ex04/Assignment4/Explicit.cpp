//
//  explicit.cpp
//  Assignment4
//
//  Created by Sagar Dolas on 16/01/16.
//  Copyright © 2016 Sagar Dolas. All rights reserved.
//

#include "Explicit.hpp"

Explicit::Explicit(int argc, char *argv[]) :ParabolicSolver(argc,argv)
{
    std::cout << "Explicit Solver Construction started\n";
}

//Implicit::Implicit(){}

Explicit::~Explicit()
{
    std::cout << "Destructor called for Explicit case\n";
    MPI_Finalize();
}

// Allocate memory
void Explicit::allocMemory()
{
    
    if(this->cart_rank ==0)
    {
        for(int i=0; i<numGridPointsX_*numGridPointsY_; ++i)
        {
            ureal_.push_back(0.);
        }
    }    
    
    if (isBoundaryProcessor()) {
        for (auto i = 0; i<sizeproc_ + numGridPointsX_ ; ++i) {
            u_.push_back(0);
            unew_.push_back(0);
        }
    }
    else{
        for (auto i = 0; i<sizeproc_ + 2* numGridPointsX_ ; ++i) {
            u_.push_back(0);
            unew_.push_back(0);
        }
    }
}

void Explicit::solver(){
    
    communicateGhostCells(this->u_);
                          
    // Explicit Euler solver
    for (auto it = this->index_.begin(); it!= this->index_.end(); ++it) {
        this->unew_[*it] =  this->u_[*it] + ((tau_ * kappa_) * laplacianOperator((*it),this->u_));
    }
    this->u_.swap(this->unew_);
}
    



