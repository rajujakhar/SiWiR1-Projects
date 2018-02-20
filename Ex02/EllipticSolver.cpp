#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>

void Solve_Serial(std::vector<double>& u , const std::vector<double>& f, const int numgridpoints_X){

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

  int numiter = 0 ;
  double tol = pow(10,-5);
  while (resi_temp > tol){
    for (size_t i = numgridpoints_x +1  ; i < numgridpoints_x * (numgridpoints_y - 1)-1 ; i++) {
      if (!( ((i+1)% numgridpoints_x == 0) || (i% numgridpoints_x ==0))){
        u[i] = (1 /((2/hxsq) + (2/(hysq)) + (K* K))) * ( f[i] + (1/hxsq)*(u[i-1] + u[i+1]) + (1/hysq)*( u[i - numgridpoints_x] + u[i + numgridpoints_x]));
      }
    }
    resi_dual = 0.0;
    for (size_t i = numgridpoints_x +1  ; i < numgridpoints_x * (numgridpoints_y - 1)-1 ; i++){
      if (!( ((i+1)% numgridpoints_x == 0) || (i% numgridpoints_x ==0))){
        Resi_D[i]= f[i] - ((u[i] * ((2/hxsq) + (2/(hysq)) + (K* K))) - (((1/hxsq)*(u[i-1] + u[i+1])) + (1/hysq)*( u[i - numgridpoints_x] + u[i + numgridpoints_x])));
        resi_dual += Resi_D[i] * Resi_D[i];
        }
    }
    // Residual Function
    std::cout << "The Residual for Iteration number ---   "<<numiter+1<< "   is   "<<sqrt(resi_dual)<< std::endl;
    resi_temp = resi_dual;
    if (resi_temp < pow(10,-11))
      std::cout << "The Solution is converged" << std::endl;
  ++numiter;
  }

}

void StartSolver(int const numgridpoints_x,int const numgridpoints_y ) {
  /* initilization based on the numgridpoints_x and numgridpoints_y */
  // numgridpoints_x = numgrid_interval

  const double PI = 4.0 * atan(1.0);
  const double hx  = 2.0 / (numgridpoints_x-1);
  const double hy = 1.0 / (numgridpoints_y-1);
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
  Solve_Serial(u,f, numgridpoints_x);
  End = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_sec = End - Start;
  std::cout << "The time taken for computation is " <<elapsed_sec.count()<< std::endl;

}

int main()
{
  // include time
  int numgridpoints_x;
  int numgridpoints_y;
  std::cout << "Enter the Division you want in X and Y direction as follows " << std::endl;
  std::cout << "Enter the division in X direction" << std::endl;
  std::cin >> numgridpoints_x;
  std::cout << "Enter the division in Y direction" << std::endl;
  std::cin >> numgridpoints_y;
  StartSolver(numgridpoints_x, numgridpoints_y );
  //std::cout << "Input Data has been taken...........Initilization in progess " << std::endl;

  return 0 ;


  // include time
}
