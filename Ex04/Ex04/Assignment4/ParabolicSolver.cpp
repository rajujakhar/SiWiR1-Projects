//
//  ParabolicSolver.cpp
//  Assignment4
//
//  Created by Raju Ram on 17/01/16.
//  Copyright © 2016  Raju Ram. All rights reserved.
//


#include "ParabolicSolver.h"

ParabolicSolver::ParabolicSolver(int argc, char *argv[]) : numGridPointsX_(std::stoi(argv[1]) +1),numGridPointsY_(std::stoi(argv[2]) + 1),tau_(std::stod(argv[6])),numTimeSteps_(std::stoi(argv[5])),alpha_(std::stod(argv[8])),tolerance_(std::stod(argv[4])),num_iter(std::stoi(argv[3])),kappa_(std::stod(argv[7])),vtkSpacing(std::stoi(argv[9])){
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size( MPI_COMM_WORLD , &size_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    
    for(int i=0;i<5;++i)
        domain_.push_back(0);
    
    hx_ = 1./(numGridPointsX_-1);
    hy_ = 1./(numGridPointsY_-1);
    hxsqvinv_ = 1.0 / (hx_ * hx_ );
    hysqvinv_ = 1.0 / (hy_ * hy_);
    isSolutionConverged =0;
    countGlobal_ = 0;
    
    this->sizeproc_ =0;
    this->radius1 = 0.2;
    this->radius2 = 0.4;

    if(argc ==11){
        boundary = argv[10];
        if(this->rank_ == 0)
            std::cout<<"The inner radius of donut like shape is:="<<radius1<<" and the outer radius is :="<<radius2<<" You can change it in constructor of parabolicsolver"<<std::endl;
    }


    std::transform(boundary.begin(), boundary.end(),boundary.begin(), ::toupper);

    if(boundary =="D")
        this->isTorusInit = true;
    else
        this->isTorusInit = false;

    if(this->rank_ ==0){
        std::cout << "alpha is:=" << alpha_ << std::endl;
        std::cout << "tau is :="<<tau_<<std::endl;
        std::cout<< "kappa is:="<<kappa_<<std::endl;
    }
}

ParabolicSolver::~ParabolicSolver(){

}

size_t ParabolicSolver::getNumRows(){
    
    size_t rowsize = u_.size() / (numGridPointsX_);
    
    // If size_1 is 1 we are unnecessarily adding one row
    if(size_ ==1) rowsize -=1;
    return rowsize;
    
}

void ParabolicSolver::correctIndexCalculator(){
    
    size_t numGridY = getNumRows();
    
    for(auto i=1; i<numGridY - 1; ++i)
    {
        for(auto j=1 ;j<numGridPointsX_-1;++j)
        {
            index_.push_back(i*numGridPointsX_+j);
        }
    }
    
}

void ParabolicSolver::initValues(){
    ndims = 2;
    dim.push_back(1);
    dim.push_back(size_);
    
    for(int i=0;i <ndims; i++){
    periods.push_back(0);
    coords.push_back(0);
    }
    
    reorder = 1;
    up=0;
    down =0;
}

void ParabolicSolver::distributeMemory(){
    if(cart_rank == 0)
    {
        sizeproc_ = domain_[1] * numGridPointsX_;
        if (domain_[2]==0)
        {
            for (int i = 1; i < size_; ++i)
            {
                MPI_Send(&sizeproc_,1,MPI_INT,i,0,MPI_COMM_WORLD);
            }
        }
        else
        {
            size_t sizeproc_2_= domain_[3] * numGridPointsX_;
            for (int i = 1; i < size_-1; ++i)
            {
                MPI_Send(&sizeproc_,1,MPI_INT,i,0,MPI_COMM_WORLD);
            }
            MPI_Send(&sizeproc_2_,1,MPI_INT,size_-1,0,MPI_COMM_WORLD);
        }
        
         //MPI_Barrier(new_comm);
    }
    
   
    
    else
    {
        MPI_Recv(&sizeproc_,1,MPI_INT,0,0,MPI_COMM_WORLD,&status_);
      //  std::cout<<"I am Rank "<<rank_<<" and I have received "<<sizeproc_<<"memeory for computation"<<std::endl;
    }
    
    int size = (int)domain_.size();
    MPI_Bcast(domain_.data(),size, MPI_INT, 0, new_comm);

}

// Initialize u values
void ParabolicSolver::init_u_(){
    double x,y;
    
    size_t begin = numGridPointsX_*domain_[1]*cart_rank;
    size_t end = begin + sizeproc_;
    double r1 =radius1;
    double r2 =radius2;
   
    for(size_t i=begin; i<end; ++i){
        
        x=hx_*(i%numGridPointsX_);
        y=hy_*(i/numGridPointsX_);
        
        // The following if checks if we need to include torus/donut like init condition
        if(this->isTorusInit) {
            if (isInsideconcentricCircles(r1, r2, x, y)){
                if(cart_rank == 0){
                    u_[i-begin]=1;
                }
                else{
                    u_[i-begin+numGridPointsX_]=1;
                }
            }
            else{
                
                 if(cart_rank == 0){
                    u_[i-begin]=0;
                 }
                 else{
                    u_[i-begin+numGridPointsX_]=0;
                 }
               }
        }
        else {
              if(cart_rank == 0){
               u_[i-begin]=sin(PI*x)*sin(PI*y);
               }
               else{
               u_[i-begin+numGridPointsX_]=sin(PI*x)*sin(PI*y);
               }
        }
       
    }
    
    // Initialize the global u value to the initial given solution u(x,y,t=0)
    if(cart_rank ==0)
    {
        for(int i=0; i<numGridPointsX_*numGridPointsY_; ++i)
        {
            x=hx_*(i%numGridPointsX_);
            y=hy_*(i/numGridPointsX_);
       
            if(isTorusInit) {
                if(isInsideconcentricCircles(r1, r2, x, y))
                    this->ureal_[i] = 1;
                else
                    this->ureal_[i] = 0;
                 
            }
            else
            this->ureal_[i] = sin(PI*x)*sin(PI*y);
        }
    }
}

// Create Cartesian Coordinate
void ParabolicSolver::createCartCord(){
    new_comm = MPI_COMM_NULL;
    MPI_Cart_create(MPI_COMM_WORLD,ndims,dim.data(),periods.data(),reorder,&new_comm);
    MPI_Comm_rank(new_comm,&cart_rank);
    MPI_Cart_coords(new_comm,cart_rank,2,coords.data());
    std::cout << "For cart_rank: " <<  cart_rank <<" cords[0] " <<  coords[0] << " cords[1] " << coords[1] << std::endl;
    MPI_Cart_shift(new_comm,1,1,&down,&up);
}

// Communicate the u value to nearby processors
void ParabolicSolver::communicateGhostCells(std::vector<double> &ghostCellVec){
    //Send up
    if(this->cart_rank ==0){
        MPI_Send(ghostCellVec.data() + sizeproc_-numGridPointsX_, numGridPointsX_,MPI_DOUBLE, up,0,new_comm);
    }
    else if( this->up >0){
        MPI_Send(ghostCellVec.data() + sizeproc_, numGridPointsX_,MPI_DOUBLE, up,0,new_comm);
    }
    
    // Receive from down
    if(this->down >= 0){
        MPI_Recv(ghostCellVec.data(), numGridPointsX_, MPI_DOUBLE, down,0, new_comm, &status_);
    }
    
    MPI_Barrier(new_comm);
    
    // Send down
    if(this->down >=0){
        MPI_Send(ghostCellVec.data()+numGridPointsX_, numGridPointsX_, MPI_DOUBLE,down,1,new_comm);
    }
    
    // Receive from up
    if(this->cart_rank ==0){
        MPI_Recv(ghostCellVec.data()+sizeproc_, numGridPointsX_, MPI_DOUBLE, up,1, new_comm, &status_);
    }
    else if(this->up>0){
        MPI_Recv(ghostCellVec.data()+sizeproc_+numGridPointsX_, numGridPointsX_, MPI_DOUBLE, up,1, new_comm, &status_);
    }
    
    
    
}

// Laplace Operator
double ParabolicSolver::laplacianOperator(const size_t &index, const std::vector<double> &vec ){
    
    return (( ((vec[index+1] - (2 * vec[index]) + vec[index-1]) * hxsqvinv_) +  ((vec[index + numGridPointsX_] - (2 * vec[index]) + vec[index - numGridPointsX_]) * hysqvinv_)));
}

void ParabolicSolver::domainDecompose(){
    
    if(cart_rank ==0)
    {
        size_t dataprocydirect = 0;
        if(numGridPointsY_ %size_ ==0)
        {
            dataprocydirect = numGridPointsY_/ size_;
            domain_[0] = size_;
            domain_[1] = dataprocydirect;
            domain_[2] = 0;
            domain_[3] = 0;
            domain_[4] = 1;
            std::cout<<"Domain is Decomposed into equal subdomains "<<std::endl;
        }
        else
        {
            if(numGridPointsY_ %(size_-1)!=0 && (!((numGridPointsY_ -((numGridPointsY_ /(size_-1))*(size_-1))) == 1 )))
            {
                dataprocydirect = numGridPointsY_ / (size_-1);
                domain_[0] = size_-1;
                domain_[1] = dataprocydirect;
                domain_[2] = 1;
                domain_[3] = numGridPointsY_- (domain_[0] * domain_[1]);
                domain_[4] = 2;
                std::cout<<"Domain is  Decomposed into irregular domains 2nd form "<<std::endl;
            }
            else
            {
                dataprocydirect = numGridPointsY_ / size_;
                domain_[0] = size_-1;
                domain_[1] = dataprocydirect;
                domain_[2] = 1;
                domain_[3] = numGridPointsY_ - (domain_[0]*domain_[1]);
                domain_[4] = 3;
                std::cout<<"Domain is Decomposed into irregular subdomains 3rd form"<<std::endl;
            }
            
        }
        
    }
}

bool ParabolicSolver::isBoundaryProcessor(){
    
    if (this->up <0 || this->down <0) {
        return true;
    }
    else
        return false;
}

void ParabolicSolver::displayData(const std::vector<double> &data){
    
    size_t numy = data.size() / numGridPointsX_ ;
    for (size_t i =0 ; i < numy; ++i) {
        std::cout << "cart_rank: " << cart_rank << " : ";
        for (size_t j =0 ; j <numGridPointsX_; ++j) {
            std::cout<<data[i*numGridPointsX_ + j]<<"\t";
        }
        std::cout<<std::endl;
    }
    
}

double ParabolicSolver::discreteL2Norm(const double value){
    
    size_t interiorpoints = (numGridPointsX_ -2) * (numGridPointsY_ -2);
    return sqrt(value /interiorpoints );
}

void ParabolicSolver::getGlobal_u_(){
    
    if(cart_rank ==0){
        for(int i=0; i< sizeproc_; ++i)
        {
         this->ureal_[i]= this->u_[i];
        }
    
        if(size_ !=1) {
            
        for(int i=1;i<size_-1; ++i)
        MPI_Recv(this->ureal_.data()+(i*sizeproc_), (int)sizeproc_, MPI_DOUBLE, i, 3,new_comm, &status_);
        
        // If domain is decomposed in first form
        if(domain_[4] == 1)
            MPI_Recv(this->ureal_.data()+(size_-1)*sizeproc_, (int)sizeproc_, MPI_DOUBLE, size_-1, 3,new_comm, &status_);
        else
        MPI_Recv(this->ureal_.data()+(size_-1)*sizeproc_, (int)(domain_[3]*numGridPointsX_), MPI_DOUBLE, size_-1, 3,new_comm, &status_);
        }
        
    }
    
    else{
        if(size_ !=1)
        MPI_Send(this->u_.data()+numGridPointsX_, (int)sizeproc_, MPI_DOUBLE, 0, 3, new_comm);
    }
    
   if(size_ !=1)
   MPI_Barrier(new_comm);
   
    if(cart_rank == 0) {
        writeSolution();
    }
    
}

void ParabolicSolver::writeVTU_file() {
    std::string str("vtk/time_step_");
    str.append(std::to_string(countGlobal_));
    str.append("_");
    str.append(std::to_string(cart_rank));
    str.append(".vtu");
    file_.open(str);
    
    int begin,end;
    double x,y;
    begin = (int) (cart_rank*domain_[1]*numGridPointsX_);
    end = begin + (int)sizeproc_;
    std::string q = "'";
    std::string s= "<?xml version=\"1.0\"?> \n <VTKFile type=\"UnstructuredGrid\"  version=\"0.1\"> \n <UnstructuredGrid> \n <Piece NumberOfPoints=\"";
    //s=  s +q +q ;
    s.append(std::to_string(sizeproc_));
   // s=  s +q +q ;
    s = s + "\" NumberOfCells=\"";
    //s=  s +q +q ;
    s= s+ std::to_string(sizeproc_);
    //s=  s +q +q ;
    s = s + "\">\n <Points>\n";
    s = s + "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    
    file_ << s ;
    
    for (int i = begin; i < end; ++i)
    {
        x = hx_ * (i % numGridPointsX_);
        y = hy_ * (i / numGridPointsY_);
        file_<<x<<"  "<<y<<"  "<<0.0<<std::endl;
    }
    
    s = "</DataArray> \n </Points> \n <Cells> \n <DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\">\n";
    
    file_ <<s ;
    
    for(int i=0; i<sizeproc_; ++i)
    {
       file_ << std::to_string(i) << "\t";
    }
    
    s = "\n </DataArray> \n <DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\">\n";
    
    file_ <<s;
    
    for(int i=0; i<sizeproc_; ++i)
    {
        file_ << std::to_string(i+1) << "\t";
    }
    
    s= "\n</DataArray> \n <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    
    file_ <<s;
    
    for(int i=0; i<sizeproc_; ++i)
    {
       file_ << "1" << "\t";
    }
    
    s= "\n</DataArray> \n </Cells> <PointData> \n <DataArray type=\"Float64\" Name=\"Temperature\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    
    file_ << s;
    
    for (int i = 0; i < sizeproc_; ++i)
    {
       // x = hx_ * (i % numGridPointsX_);
       // y = hy_ * (i / numGridPointsY_);
        if(cart_rank ==0)
        file_<< u_[i]<< "\t";
        else
        file_<< u_[i+numGridPointsX_]<< "\t";
    }
    
    s ="\n</DataArray> \n </PointData> \n </Piece> \n </UnstructuredGrid> \n </VTKFile>\n";
    
    file_ << s;
    
    
    file_.close();
    
}

void ParabolicSolver::writePVTU_file()
{
    std::string str("vtk/time_step_");
    str.append(std::to_string(countGlobal_));
    str.append(".pvtu");
    file_.open(str);
    
    str = "<?xml version=\"1.0\"?> \n <VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"> \n <PUnstructuredGrid> \n <PPoints> \n <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"/> \n </PPoints> \n <PPointData> \n <PDataArray type=\"Float64\" Name=\"Temperature\" NumberOfComponents=\"1\" format=\"ascii\"/> \n </PPointData> \n ";
    file_ << str;
    
    // Append .vtu files for each processor at given time step
    for(int i=0; i< size_; ++i)
    {
        file_ << "<Piece Source=\"time_step_" << std::to_string(this-> countGlobal_) << "_" << std::to_string(i) << ".vtu\"/> \n";
    }
    
    str = "\n </PUnstructuredGrid> \n </VTKFile>" ;
    file_ << str << "\n";
    file_ .close();

}

void ParabolicSolver::writePVD_file()
{
    std::string str("vtk/solution.pvd");
    file_.open(str);
    
    str = " <?xml version=\"1.0\"?> \n <VTKFile type=\"Collection\" version=\"0.1\"> \n <Collection> \n ";
    file_ << str <<std::endl;
    
    for(int i=0; i<this->numTimeSteps_; i=i+int(this->vtkSpacing)) {
       file_ << "<DataSet timestep=\"0\" file=\"time_step_"<<i<<".pvtu\"/> \n";
     }
    
    str =  " </Collection> \n </VTKFile> \n";
    file_ << str << std::endl;

}

 void ParabolicSolver::writeSolution(){
     double x,y;
     std::string str("sol/time_step_");
     str.append(std::to_string(countGlobal_));
     str.append(".txt");
     file_.open(str);
     const int numGridPoints = (this->numGridPointsX_)*(this->numGridPointsY_);
     
    for (auto i = 0; i < numGridPoints; ++i)
    {
        x = hx_ * (i % numGridPointsX_);
        y = hy_ * (i / numGridPointsY_);
        file_<<x<<"  "<<y<<"  "<<this->ureal_[i]<<std::endl;
    }
    file_.close();
   //std::cout << "The solution has been written in solution.txt" << std::endl;
    
}

void ParabolicSolver::callSolver(){
   
    // This time the count Global is zero
    writePVD_file();
    writePVTU_file();
    getGlobal_u_();
    
    while(this->countGlobal_ < this->numTimeSteps_)
    {
     this->solver();
     this->countGlobal_++;
      
      if((this->countGlobal_ % this->vtkSpacing) == 0)
        {
            writeVTU_file();
        }
        MPI_Barrier(new_comm);
        
        if((cart_rank == 0) && (this->countGlobal_ % this->vtkSpacing) == 0)
        {
            writePVTU_file();
        }
        
        // to check gnuplot. delete this for actual solution
        if((this->countGlobal_ % this->vtkSpacing) == 0)
            getGlobal_u_();
    }
    
    if(cart_rank ==0)
    writePVD_file();
    
}

inline bool ParabolicSolver::isInsideconcentricCircles(const double r1, const double r2,const double x,const double y){
    
    if (((pow(x-0.5, 2)) + (pow(y-0.5, 2)) <= (pow(r2, 2))) &&
        (((pow(x-0.5, 2)) + (pow(y-0.5, 2)) >= (pow(r1, 2)))  ))
    {
        return true;
    }
    else
        return false;
    
}

void ParabolicSolver::starttime(){
    timer.reset();
}

void ParabolicSolver::endtime(){

    this->timetaken= timer.elapsed();

    if(this->cart_rank == 0)
        std::cout<<"The time taken to complete simulation is :="<<this->timetaken<<std::endl;

}


