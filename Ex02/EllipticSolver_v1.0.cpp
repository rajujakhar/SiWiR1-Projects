#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>
#include <omp.h>
#include <string>
#include <fstream>

void Solve_Serial(std::vector<double> u , const std::vector<double> f, const int numgridpoints_X, const int numiter_){

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
  int numiter = numiter_;
  double tolerance = pow(10, -3);
  
  while (counter< numiter)
  {
    resi_dual = 0.0;

    for (size_t i = numgridpoints_x +1  ; i < numgridpoints_x * (numgridpoints_y - 1)-1 ; i++)
    {
      if (!( ((i+1)% numgridpoints_x == 0) || (i% numgridpoints_x ==0))){
        Resi_D[i]= f[i] - ((u[i] * ((2/hxsq) + (2/(hysq)) + (K* K))) - (((1/hxsq)*(u[i-1] + u[i+1])) + (1/hysq)*( u[i - numgridpoints_x] + u[i + numgridpoints_x])));
        resi_dual += Resi_D[i] * Resi_D[i];
        ++counter;
      }
    }

    for (size_t i = numgridpoints_x +1  ; i < numgridpoints_x * (numgridpoints_y - 1)-1 ; i++) 
    {
      if (!( ((i+1)% numgridpoints_x == 0) || (i% numgridpoints_x ==0))){
        u[i] = (1 /((2/hxsq) + (2/(hysq)) + (K* K))) * ( f[i] + (1/hxsq)*(u[i-1] + u[i+1]) + (1/hysq)*( u[i - numgridpoints_x] + u[i + numgridpoints_x]));
      }
    }

    // Residual Function
    //std::cout << "The counter is " <<counter<< std::endl;
    //std::cout << "The Residual for Iteration number ---   "<<counter+1<< "   is   "<<sqrt((resi_dual/counter))<< std::endl;
    resi_temp = resi_dual;
    //if (resi_temp < pow(10,-11))
      //std::cout << "The Solution is converged" << std::endl;
 

  ++counter;
  }

  // This following code writes the result of u in  the file
  std::ofstream mystream;
  mystream.open("uVal.txt"); 
  for(int i= 0; i< numgridpoints_y; ++i )
  {
    for(int j=0; j<numgridpoints_x ; ++j)
    {
      mystream << hx*(j%numgridpoints_x)  << "           " << hy*(i%numgridpoints_y) << "          " << u[i*numgridpoints_x + j]  << "\n";
    }
  }
  mystream.close();

}

void Solve_parallel(std::vector<double> u , const std::vector<double> f, const int numgridpoints_X, int const numiter_){

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
  int numiter = numiter_;
  std::cout<<numiter<<std::endl;
  double tolerance = pow(10, -3);
  while (counter < numiter){

    resi_dual = 0.0;
    //#pragma omp parallel for schedule (static)
    for (size_t i = numgridpoints_x +1  ; i < numgridpoints_x * (numgridpoints_y - 1)-1 ; i++){
      if (!( ((i+1)% numgridpoints_x == 0) || (i% numgridpoints_x ==0))){
        Resi_D[i]= f[i] - ((u[i] * ((2/hxsq) + (2/(hysq)) + (K* K))) - (((1/hxsq)*(u[i-1] + u[i+1])) + (1/hysq)*( u[i - numgridpoints_x] + u[i + numgridpoints_x])));
        resi_dual += Resi_D[i] * Resi_D[i];
       // std::cout << "checking" <<resi_dual<< std::endl;

      }
    }
    // updating red values
    #pragma omp parallel for schedule (static)
    for (size_t i = numgridpoints_x +1  ; i < numgridpoints_x * (numgridpoints_y - 1)-1 ; i+=2) {
      if (!( ((i+1)% numgridpoints_x == 0) || (i% numgridpoints_x ==0))){
        u[i] = (1 /((2/hxsq) + (2/(hysq)) + (K* K))) * ( f[i] + (1/hxsq)*(u[i-1] + u[i+1]) + (1/hysq)*( u[i - numgridpoints_x] + u[i + numgridpoints_x]));
      }
    }
    // updating black values
    #pragma omp parallel for schedule (static)
    for (size_t i = numgridpoints_x +2  ; i < numgridpoints_x * (numgridpoints_y - 1)-1 ; i+=2) {
      if (!( ((i+1)% numgridpoints_x == 0) || (i% numgridpoints_x ==0))){
        u[i] = (1 /((2/hxsq) + (2/(hysq)) + (K* K))) * ( f[i] + (1/hxsq)*(u[i-1] + u[i+1]) + (1/hysq)*( u[i - numgridpoints_x] + u[i + numgridpoints_x]));
      }
    }
    // Residual Function
    //std::cout << "The counter is " <<counter<< std::endl;
    std::cout << "The Residual for Iteration number ---   "<<counter+1<< "   is   "<<sqrt((resi_dual/total_points))<< std::endl;
    resi_temp = resi_dual;
    //if (resi_temp < pow(10,-11))
      std::cout << "The Solution is converged" << std::endl;
  ++counter;
  }


}

void StartInit(int const numgridpoints_x,int const numgridpoints_y , int const numiter) {
   //initilization based on the numgridpoints_x and numgridpoints_y 
  // numgridpoints_x = numgrid_interval

  const double PI = 4.0 * atan(1.0);
  const double hx  = 2.0 / (numgridpoints_x-1);
  const double hy = 1.0 / (numgridpoints_y-1);
  int const numiter_ = numiter;
  std::cout << "Initilization Starts......." << std::endl;
  std::vector<double> u(numgridpoints_x * numgridpoints_y, 0.0);
  std::vector<double> f(numgridpoints_x * numgridpoints_y);

  // initilization of upper Boundary Values
  std::cout << "Initilizating Upper Boundary Values" << std::endl;
  #pragma omp parallel for schedule (static)
  for (int i = (numgridpoints_y-1) * numgridpoints_x; i < (numgridpoints_x * numgridpoints_y) ; ++i ) {
      u[i] = sin(2 * PI * (hx *(i - ((numgridpoints_y - 1.0) * numgridpoints_x))) * sinh(2.0 * PI));
      //std::cout << "the values  u_i " <<u[i]<<std::endl;
  }
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

  std::cout << "Now initilization ends......Starting the solver " << std::endl;
  std::chrono::time_point<std::chrono::system_clock> Start, End;



  Start = std::chrono::system_clock::now();
  Solve_Serial(u,f, numgridpoints_x, numiter_);
  End = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_sec = End - Start;
  std::cout << "The time taken for Serial computation is " <<elapsed_sec.count()<< std::endl;



 /* Start = std::chrono::system_clock::now();
  Solve_parallel(u,f, numgridpoints_x, numiter_);
  End = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_sec_1 = End - Start;
  std::cout << "The time taken for parallel computation is " <<elapsed_sec_1.count()<< std::endl;*/ 



}
// implementing the exact program 

int main(int argc,const char *argv[])
{
  std::string s1 = argv[1],s2= argv[2],s3 = argv[3];
  int numgridpoints_x = std::stoi(s1);
  int numgridpoints_y = std::stoi(s2);
  int numiter = std::stoi(s3);

 
  std::cout << "Input Data has been taken...........Initilization in progess\n numgridpoints_x:  " <<numgridpoints_x << "  numgridpoints_y : "<<numgridpoints_y<<"  numiter : "<<numiter<< std::endl;
  StartInit(numgridpoints_x, numgridpoints_y ,numiter);

  return 0 ;


  // include time
}
