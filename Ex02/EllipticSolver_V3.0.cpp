#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>
#include <omp.h>
#include <fstream>

void calIndexes(const int numGridX, const int numGridY, std::vector<int> &redIndex,std::vector<int> &blackIndex)
{
 bool toggle = false;
 int numTotGrids = numGridX*numGridY;

 if((numGridX % 2) == 1)
 {
   for(int i=0;i < numTotGrids; i+=2)
   {
     redIndex.push_back(i);
     if(i!= (numTotGrids - 1))
     blackIndex.push_back(i+1);
   }

 }
 else
 {
   for(size_t i=0; i<numGridY; ++i)
   {
     toggle = !toggle;

     for(size_t j=0;j<numGridX;j+=2)
       {
         if(toggle)
           {
             redIndex.push_back(i*numGridX +j);
             blackIndex.push_back(1+i*numGridX +j);
           }
         else
           {
             redIndex.push_back(1+ i*numGridX +j);
             blackIndex.push_back(i*numGridX +j);
           }
       }
   }
 }

}

void Solve_Serial(std::vector<double> u , const std::vector<double> f, const int numgridpoints_X, const int numiter){

  const double PI = 4.0 * atan(1.0);
  std::cout << "Starting Serial Solver " << std::endl;
  const size_t numgridpoints = u.size();
  const int numgridpoints_x = numgridpoints_X;
  const int numgridpoints_y = numgridpoints/numgridpoints_x;
  const double hx  = 2.0 / (numgridpoints_x-1);
  const double hy = 1.0 / (numgridpoints_y-1);
  const double hxsq = hx * hx ;
  const double hysq = hy * hy ;
  const double K = 2 * PI;
  long double resi_dual;
  std::vector<double> Resi_D(numgridpoints_x * numgridpoints_y ,0.0);
  int size = Resi_D.size();
  long double resi_temp =1;
  int counter=0;
  const int total_points = numgridpoints_x * numgridpoints_y;

  int numiter_s = 0 ;
  double tolerance = pow(10, -3);
  while (numiter_s < numiter){
    resi_dual = 0.0;

    for (size_t i = numgridpoints_x +1  ; i < numgridpoints_x * (numgridpoints_y - 1)-1 ; i++){
      if (!( ((i+1)% numgridpoints_x == 0) || (i% numgridpoints_x ==0))){
        Resi_D[i]= f[i] - ((u[i] * ((2/hxsq) + (2/(hysq)) + (K* K))) - (((1/hxsq)*(u[i-1] + u[i+1])) + (1/hysq)*( u[i - numgridpoints_x] + u[i + numgridpoints_x])));
        resi_dual += Resi_D[i] * Resi_D[i];
      }
    }
    for (size_t i = numgridpoints_x +1  ; i < numgridpoints_x * (numgridpoints_y - 1)-1 ; i++) {
      if (!( ((i+1)% numgridpoints_x == 0) || (i% numgridpoints_x ==0))){
        u[i] = (1 /((2/hxsq) + (2/(hysq)) + (K* K))) * ( f[i] + (1/hxsq)*(u[i-1] + u[i+1]) + (1/hysq)*( u[i - numgridpoints_x] + u[i + numgridpoints_x]));
      }
    }
    //std::cout << "The Residual for Iteration number ---   "<<numiter_s+1<< "   is   "<<sqrt((resi_dual/total_points))<< std::endl;
    resi_temp = resi_dual;
  ++numiter_s;
  }
  std::ofstream datafile;
  datafile.open("Solution_serial.txt");
  for (size_t i = 0; i < numgridpoints_y; ++i) {
    for (size_t j = 0; j < numgridpoints_x; ++j) {

      datafile << hx * (j% numgridpoints_x) <<"       "<< hy * (i% numgridpoints_y)<<"      "<< u[i*numgridpoints_x + j]<< std::endl ;
    }
  }
  datafile.close();



}

void Solve_parallel(std::vector<double> u , const std::vector<double> f, const int numgridpoints_X,const std::vector<int>& redvec, const std::vector<int>& blackvec , const int numiter){

  const double PI = 4.0 * atan(1.0);
  std::cout << "Starting parallel Solver " << std::endl;
  const size_t numgridpoints = u.size();
  const int numgridpoints_x = numgridpoints_X;
  const int numgridpoints_y = numgridpoints/numgridpoints_x;
  const double hx  = 2.0 / (numgridpoints_x-1);
  const double hy = 1.0 / (numgridpoints_y-1);
  const double hxsq = hx * hx ;
  const double hysq = hy * hy ;
  const double K = 2 * PI;
  long double resi_dual;
  std::vector<double> Resi_D(numgridpoints_x * numgridpoints_y ,0.0);
  int size = Resi_D.size();
  long double resi_temp =1;
  int counter=0;
  int total_points = numgridpoints_x * numgridpoints_y;
  int num_threads;

  int numiter_p = 0 ;
  double tolerance = pow(10, -3);
  while (numiter_p < numiter){

    resi_dual = 0.0;
    /*
    for (size_t i = 0; i < redvec.size(); i++) {
      std::cout << "The element of redvec are " << redvec[i]<<std::endl;
    }
    */
    //#pragma omp parallel for schedule (static)
    for (size_t i = numgridpoints_x +1  ; i < numgridpoints_x * (numgridpoints_y - 1)-1 ; i++){
      if (!( ((i+1)% numgridpoints_x == 0) || (i% numgridpoints_x ==0))){
        Resi_D[i]= f[i] - ((u[i] * ((2/hxsq) + (2/(hysq)) + (K* K))) - (((1/hxsq)*(u[i-1] + u[i+1])) + (1/hysq)*( u[i - numgridpoints_x] + u[i + numgridpoints_x])));
        resi_dual += Resi_D[i] * Resi_D[i];
        //std::cout << "checking" <<resi_dual<< std::endl;
      }
    }
    // updating red values
    #pragma omp parallel for schedule (static)
    for (int i = 0; i < redvec.size(); ++i){
      if (!(( ((redvec[i]+1)% numgridpoints_x == 0) || (redvec[i]% numgridpoints_x ==0)) || (((redvec[i]/numgridpoints_x == 0)|| (redvec[i]/numgridpoints_x == (numgridpoints_y-1)))))){
        u[redvec[i]] = (1 /((2/hxsq) + (2/(hysq)) + (K* K))) * ( f[redvec[i]] + (1/hxsq)*(u[redvec[i]-1] + u[redvec[i]+1]) + (1/hysq)*( u[redvec[i] - numgridpoints_x] + u[redvec[i] + numgridpoints_x]));
        //std::cout << "The value at u is "<<u[i]<<"  at "<<redvec[i] << std::endl;
      }
    }
     //updating the black values
    #pragma omp parallel for schedule (static)
    for (size_t i = 0; i < blackvec.size(); ++i) {
      if (!(( ((blackvec[i]+1)% numgridpoints_x == 0) || (blackvec[i]% numgridpoints_x ==0)) || (((blackvec[i]/numgridpoints_x == 0)|| (blackvec[i]/numgridpoints_x == (numgridpoints_y-1)))))) {
        u[blackvec[i]] = (1 /((2/hxsq) + (2/(hysq)) + (K* K))) * ( f[blackvec[i]] + (1/hxsq)*(u[blackvec[i]-1] + u[blackvec[i]+1]) + (1/hysq)*( u[blackvec[i] - numgridpoints_x] + u[blackvec[i] + numgridpoints_x]));
      }

    }

    std::cout << "The Residual for Iteration number ---   "<<numiter_p+1<< "   is   "<<sqrt((resi_dual/total_points))<< std::endl;
    resi_temp = resi_dual;
    //if (resi_temp < pow(10,-11))
      //std::cout << "The Solution is converged" << std::endl;
  ++numiter_p;
  }
  std::ofstream datafile_p;
  datafile_p.open("Solution_parallel.txt");
  for (size_t i = 0; i < numgridpoints_y; ++i) {
    for (size_t j = 0; j < numgridpoints_x; ++j) {

      datafile_p << hx * (j% numgridpoints_x) <<"       "<< hy * (i% numgridpoints_y)<<"      "<< u[i*numgridpoints_x + j]<< std::endl ;
    }
  }
  datafile_p.close();

}

void StartInit(int const numgridpoints_x,int const numgridpoints_y, const int numiter ) {
  /* initilization based on the numgridpoints_x and numgridpoints_y */
  // numgridpoints_x = numgrid_interval

  const double PI = 4.0 * atan(1.0);
  const double hx  = 2.0 / (numgridpoints_x-1);
  const double hy = 1.0 / (numgridpoints_y-1);
  std::cout << "Initilization Starts......." << std::endl;
  std::vector<double> u(numgridpoints_x * numgridpoints_y, 0.0);
  std::vector<double> f(numgridpoints_x * numgridpoints_y, 0.0);

  // initilization of upper Boundary Values
  std::cout << "Initilizating Upper Boundary Values" << std::endl;
  #pragma omp parallel for schedule (static)
  for (int i = (numgridpoints_y-1) * numgridpoints_x; i < (numgridpoints_x * numgridpoints_y) ; ++i ) {
      u[i] = sin(2 * PI * (hx *(i - ((numgridpoints_y - 1.0) * numgridpoints_x))) * sinh(2.0 * PI));
      std::cout << "the values  u_i for boundary " <<u[i]<<std::endl;
  }

 
 std::ofstream datafile1;
  datafile1.open("u_boundary.txt");
  
    for (size_t i = 0; i < numgridpoints_x; ++i) {

      datafile1 << i <<"    "<< u[i + (numgridpoints_y-1) * numgridpoints_x ]<< std::endl;
    
  }
  datafile1.close();


  // evauating the function
  double x ,y ;
  std::cout << "Initilizating the Domain with Function values" << std::endl;
  #pragma omp parallel for schedule (static)
  for (int  i = 0; i < numgridpoints_x * numgridpoints_y; ++i) {
    x = hx * (i % numgridpoints_x);
    y = hy * (i / numgridpoints_x);
    f[i] = 4 * PI *PI * sin(2 * PI * x) * sinh(2* PI * y );
    //std::cout << "tHE VALUE f_i" <<f[i]<< std::endl;

  }
  // Founding the Red and Black array
  std::vector<int> redIndArr;
  std::vector<int> blackIndArr;
  calIndexes(numgridpoints_x,numgridpoints_y,redIndArr, blackIndArr);
  //std::cout << "Red Index Array\n";
  //for(int i= 0 ; i < redIndArr.size(); ++i)
    //std::cout << "Red "

  std::cout << "Now initilization ends......Starting the solver " << std::endl;
  std::chrono::time_point<std::chrono::system_clock> Start, End;



  Start = std::chrono::system_clock::now();
  Solve_Serial(u,f, numgridpoints_x, numiter);
  End = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_sec = End - Start;
  std::cout << "The time taken for Serial computation is " <<elapsed_sec.count()<< std::endl;




  /*Start = std::chrono::system_clock::now();
  Solve_parallel(u,f, numgridpoints_x,redIndArr, blackIndArr, numiter);
  End = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_sec_1 = End - Start;
  std::cout << "The time taken for parallel computation is " <<elapsed_sec_1.count()<< std::endl;*/




}

int main(int argc, char const *argv[]) {

  std::string s1 = argv[1],s2= argv[2],s3 = argv[3];
  int numgridpoints_x = std::stoi(s1) ;
  int numgridpoints_y = std::stoi(s2);
  int numiter = std::stoi(s3);
  StartInit(numgridpoints_x, numgridpoints_y,numiter);
  std::cout << "Input Data has been taken...........Initilization in progess " << std::endl;

  return 0 ;

}
