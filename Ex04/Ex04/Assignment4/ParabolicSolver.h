
#ifndef Parabolic_h
#define Parabolic_h
#pragma once

#include "global.hpp"
#include <vector>
#include <mpi.h>
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include "Timer.h"
#include <locale>
#include <string>
#include <algorithm>


class ParabolicSolver {
    
public:
    
    double alpha_;
    int numGridPointsX_;
    int numGridPointsY_;
    int rank_;
    int size_;
    double hx_;
    double hy_;
    double hxsqvinv_;
    double hysqvinv_;
    double tau_;
    size_t numTimeSteps_;
    double kappa_;
    size_t sizeproc_;
    std::vector<size_t> domain_;
    std::vector<double> index_;
    size_t vtkSpacing;
    siwir::Timer timer;
    std::string boundary;
    double radius1;
    double radius2;
    
    // Simulation parameters
    
    std::vector<double> u_;
    std::vector<double> unew_;
    std::vector<double> ureal_;
    double timetaken;
    
    // Cartesian coordinate variables
    
    MPI_Comm new_comm;
    MPI_Status status_;
    int ndims;
    int reorder;
    int up;
    int down;
    int cart_rank;
    std::vector<int> coords;
    std::vector <int> dim;
    std::vector <int> periods;
    
    // cg variables
    
    double scalarProductDZ_;
    double globalScalarProduct_;
    double scalarProductrr_;
    double alphacg_;
    double deltazero_;
    double deltaone_;
    double tolerance_;
    int isSolutionConverged;
    double betacg_;
    int num_iter;
    int countGlobal_;
    int cgCounter_;
    std::ofstream file_;

    // Initial torus condition
    bool isTorusInit;
    // Member Functions
    
    ParabolicSolver(int argc, char *argv[]);
    virtual ~ParabolicSolver()=0;

    void correctIndexCalculator();
    void domainDecompose();
    void createCartCord();
    void initValues();
    void distributeMemory();
    virtual void allocMemory()=0;
    void init_u_();
    size_t getNumRows();
    void communicateGhostCells(std::vector<double> &ghostCellVec);
    double laplacianOperator(const size_t &index, const std::vector<double> &vec );
    virtual void solver()=0;
    bool isBoundaryProcessor();
    void displayData(const std::vector<double> &data);
    double discreteL2Norm(const double value);
    void getGlobal_u_();
    void writeSolution();
    void callSolver();
    void starttime();
    void endtime();

    //VTK Helper function
    void writeVTU_file();
    void writePVTU_file();
    void writePVD_file();

    void torusInitCondition();
    inline bool isInsideconcentricCircles(const double r1, const double r2,const double x,const double y);
    
};

#endif /* EllipticSolver_h */
