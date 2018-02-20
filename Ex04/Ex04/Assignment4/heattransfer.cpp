//
//  main.cpp
//  Assignment4
//
//  Created by Sagar Dolas on 05/01/16.
//  Copyright © 2016 Sagar Dolas. All rights reserved.
//

#include "ParabolicSolver.h"
#include "Implicit.hpp"
#include "Explicit.hpp"
#include <iostream>
#include <cmath>
#include <mpi.h>

template <typename T>
void runSimulation(int &argc, char *argv[])
{
   T sim(argc, argv);
   sim.initValues();
   sim.createCartCord();
   sim.domainDecompose();
   sim.distributeMemory();
   sim.allocMemory();
   sim.init_u_();
   sim.correctIndexCalculator();
   sim.starttime();
   sim.callSolver();
   sim.endtime();
};

int main(int argc,char * argv[]) {
    
    
    if (argc==10) {
        std::cout << "You have chosen for default boundary condition, to change please enter d for donut like boundary condition at the end of command line arguments " << std::endl;
        if (std::stod(argv[8]) == 0) {
            runSimulation<Explicit>(argc,argv);
        }
        else
            runSimulation<Implicit>(argc,argv);
    }
    else if(argc ==11){
          std::cout<<"You have chosen for donut like bounary condition"<<std::endl;
          if (std::stod(argv[8]) == 0) {
             runSimulation<Explicit>(argc,argv);
          }
          else
             runSimulation<Implicit>(argc,argv);
    }
    else{
        std::cout<<"The number of arguments do not match minimum requirements for program to run , pleas enter the 10 number of arguments , Program exiting "<<std::endl;
        exit(1);
    }


    return 0;
}
