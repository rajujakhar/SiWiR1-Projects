//
//  Implicit.cpp
//  Assignment4
//
//  Created by Sagar Dolas on 16/01/16.
//  Copyright © 2016 Sagar Dolas. All rights reserved.
//

#include "Implicit.hpp"

Implicit::Implicit(int argc, char *argv[]) :ParabolicSolver(argc,argv){
     //std::cout << "In the Constructor of Implicit\n";
}

Implicit::~Implicit(){
   // std::cout << "Destructor called for Implicit\n";
    MPI_Finalize();
}

// Allocate memory
void Implicit::allocMemory(){
   
   // Need to uncomment this
   // if(this->cart_rank ==0)
   // {
        for(int i=0; i<numGridPointsX_*numGridPointsY_; ++i)
        {
            this->ureal_.push_back(0.);
        }
  //  }
    
    
    if (isBoundaryProcessor()) {
        for (auto i = 0; i< this->sizeproc_ + numGridPointsX_ ; ++i) {
            u_.push_back(0.0);
            d_.push_back(0.0);
            r_.push_back(0.0);
            z_.push_back(0.0);
            
        }
    }
    else{
        for (auto i = 0; i< this->sizeproc_ + 2* numGridPointsX_ ; ++i) {
            u_.push_back(0.0);
            d_.push_back(0.0);
            r_.push_back(0.0);
            z_.push_back(0.0);
            
        }
    }
}

/*void Implicit::computeCgIter(){
    
    
}*/

void Implicit::solver(){
    
    isSolutionConverged =0;
    double u_right,u_left;
    cgCounter_ =0;
    scalarProductrr_  =0.0;
    communicateGhostCells(u_);
   
    for(auto it=this->index_.begin(); it!=this->index_.end();++it){
       u_right =tau_*(1-alpha_)*kappa_*laplacianOperator(*it,this->u_) + this->u_[*it];
       u_left = u_[*it] - tau_*alpha_*kappa_*laplacianOperator(*it,this->u_);
       this->r_[*it] = u_right - u_left;
    }
    

    for(int i=0; i<r_.size(); ++i)
    {
        scalarProductrr_ += r_[i]*r_[i];
    }
    
     MPI_Allreduce(&scalarProductrr_, &deltazero_, 1, MPI_DOUBLE, MPI_SUM, new_comm);
    
    if(discreteL2Norm(deltazero_) < tolerance_)
    {
        isSolutionConverged=1;
    }
    
    this->d_=this->r_;
    
    
    while ((cgCounter_ <num_iter) && (this->isSolutionConverged==0)) {
        //computeCgIter();
        
        
        scalarProductrr_ = 0.0;
        scalarProductDZ_ = 0.0;
        
        communicateGhostCells(this->d_);
        
        for(auto it=this->index_.begin(); it!= this->index_.end();++it){
            this->z_[*it] = this->d_[*it] - tau_*alpha_*kappa_*laplacianOperator(*it,this->d_);
        }
       
        for(int i=0; i<d_.size(); ++i)
        {
            scalarProductDZ_ += d_[i]*z_[i];
        }
        
        
        MPI_Allreduce(&scalarProductDZ_, &globalScalarProduct_, 1, MPI_DOUBLE, MPI_SUM, new_comm);
        
        this->alphacg_ = deltazero_ / globalScalarProduct_;
        
        for (auto it = this->index_.begin(); it!= this->index_.end(); ++it) {
            
            u_[*it] += this->alphacg_ * this->d_[*it];
            r_[*it] -= this->alphacg_ * this->z_[*it];
        }
        
        for(int i=0; i<r_.size(); ++i)
        {
            scalarProductrr_ += r_[i]*r_[i];
        }
        
        
        MPI_Allreduce(&scalarProductrr_, &deltaone_, 1, MPI_DOUBLE,MPI_SUM, new_comm);
        
        if(discreteL2Norm(deltaone_) < tolerance_) {
            this->isSolutionConverged = 1;
            break;
        }
        
        betacg_ = deltaone_ / deltazero_;
        
        for (auto it = this->index_.begin(); it!= this->index_.end(); ++it) {
            this->d_[*it] = this->r_[*it] + betacg_ * this->d_[*it];
        }
        
        deltazero_ = deltaone_;
        
        ++cgCounter_;
        MPI_Barrier(new_comm);
    }
   
    if(cart_rank ==0)
    {
   // std::cout << "Iter taken to converge: " << countGlobal_  << "  is " << cgCounter_ << std::endl;
    std::cout << "The CG iterations are "<<cgCounter_<<" and the Discrete L2 norm: at global time " << countGlobal_ << "is := " << discreteL2Norm(deltazero_) << std::endl;
    }
    
    
}
