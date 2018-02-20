//*************************************************************************************
//******* 2D POISSON's EQUATION USING MPI PARALLELISATION****************************//
//*************************************************************************************
/*
Explanation .....How we are solving the particular problem 
*/

//***************************
//*** INCLUDES ************//
//***************************

#include <iostream>
#include <cmath>
#include <mpi.h>

//***************************************************
//****** GLOBAL VARAIBLES *************************//
//*************************************************** 

const double PI = 4.0 * atan(1.0);
const double con_sin = sinh(2*PI); // Has to be changed afterwards
const double tol = 0.0001;

//**************************************************
//******* FUNCTION DECLARATION *******************//
//**************************************************
void domaindecompose(const int &numx,const int &numy,const int &size,int domain[]);

//*****************
//** MAIN()******//
//*****************
int main(int argc, char *argv[])
{
	int size(0);
	int rank(0);

	// *************************************
	// *** MPI Initialisation **************
	//**************************************
	MPI_Init(&argc,&argv);
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  	// Checking whether number of parameters passed during main call are as per the requirement or not and if 
  	// they are not , program is terminated and this is checked by processor 0 

  	if (rank==0)
	{
  		if (argc!=4) {
   		std::cout << "Aborting the Program , Please check the number of the input arguments" << std::endl;
    	exit(1);
  	}
  		else
    	std::cout << "OK .....Parameters Taken successfully" << std::endl;
	}

	//********************************************************************
	//****** PARAMETER INITIALISATION **********************************//
	//********************************************************************
	
	std::string s1 = argv[1],s2= argv[2],s3 = argv[3];
  	int numx_interval = std::stoi(s1) ;
  	int numy_interval = std::stoi(s2);
  	int num_iter = std::stoi(s3);
  	int numgridpoints_x = numx_interval+1;
  	int numgridpoints_y = numy_interval+1;
  	const double hx  = 2.0 / (numx_interval);
  	const double hy = 1.0 / (numy_interval);
  	const double hxsq = hx * hx ;
  	const double hysq = hy * hy ;
  	const double hxsqinv = 1.0/ hxsq;
  	const double hysqinv = 1.0/ hysq;
  	const double K = 2 * PI;
  	const double constant = (2*hxsqinv) + (2*hysqinv) + (K*K);
  	double local_residual= 0.0;
  	int sizeproc=0;
  	int domain[4] ={0};

  	//***************************************
  	//** Domain Decompose *****************//
  	//***************************************
  	if (rank==0)
  	{
		  domaindecompose(numgridpoints_x,numgridpoints_y,size,domain);

  	}
    
  	//***************************************
  	//******Cartesian Cordinate System ****** 
  	//***************************************

    MPI_Comm new_comm(MPI_COMM_NULL);
    int ndims = 2;
    int dim[ndims] = {1,size};
    int periods[ndims] = {0};
    int reorder = 1;
    int up,down = 0;
    int coords[2] ={0,0};
    int cart_rank(0);

    MPI_Cart_create(MPI_COMM_WORLD,ndims,dim,periods,reorder,&new_comm);
    MPI_Comm_rank(new_comm,&cart_rank);
    MPI_Cart_coords(new_comm,cart_rank,2,coords);
    MPI_Cart_shift(new_comm,1,1,&down,&up);
    std::cout<<"this is my rank"<<cart_rank<<std::endl;
    std::cout<<"This is my neighbours"<<up<<"  "<<down<<std::endl;
    
  	//********************************************
  	// Allocating memory in each processor *****//
  	//********************************************
  	// Memory in each processor allocated according to domain array 
    
    // Communicate 
    if(cart_rank == 0)
    {
      if (domain[2]==0)
      {
        sizeproc = domain[1] * numgridpoints_x;
        for (int i = 1; i < size; ++i)
        {
          MPI_Send(&sizeproc,1,MPI_INT,i,0,MPI_COMM_WORLD);
        }
      }
      else
      {
        sizeproc= domain[1] * numgridpoints_x;
        int sizeproc_2= domain[3] * numgridpoints_x;
        for (int i = 1; i < size-1; ++i)
        {
          MPI_Send(&sizeproc,1,MPI_INT,i,0,MPI_COMM_WORLD);
        }
          MPI_Send(&sizeproc_2,1,MPI_INT,size-1,0,MPI_COMM_WORLD);
      }
    }
	  else
	  {
  		MPI_Status status;
  		MPI_Recv(&sizeproc,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
  		std::cout<<"I am Rank "<<rank<<" and I have received "<<sizeproc<<"memeory for computation"<<std::endl;
  	}

    // Memory allocation for each processor according to cartesian topology
    if(down<0 || up <0)
    {
      //std::cout<<"I am rank "<<cart_rank<<std::endl;
      //std::cout<<"I am allocating memeory    "<<sizeproc+numgridpoints_x<<" for computation"<<std::endl;
      // Allocating for u
      double u[sizeproc+numgridpoints_x]={0};
      // Allocating for f
      double f[sizeproc]={0};
      //Allocating for r
      double r[sizeproc]={0};
      //Allocating for d
      double d[sizeproc+(1*numgridpoints_x)]={0};
    

    }
    else
    {
      //std::cout<<"I am rank "<<cart_rank<<std::endl;
      //std::cout<<"I am allocating memory    "<<sizeproc+(2*numgridpoints_x)<<"for computation"<<std::endl;
      //Allocating for u
      double u[sizeproc+(2*numgridpoints_x)]={0};
      //Allocating for f
      double f[sizeproc]={0};
      //Allocating for r
      double r[sizeproc]={0};
      //Allocating for d
      double d[sizeproc+(2*numgridpoints_x)]={0};
    }

    //***********************************************
    //**** Boundary conditions and Initalisation***//
    //***********************************************


    

  	MPI_Finalize();

return 0;
}

void domaindecompose(const int &numx,const int &numy,const int &size,int domain[])
{
	std::cout<<"Domain Decmposition starts for 1D array of processor"<<std::endl;	
	int dataprocydirect = 0;
	int numprocxdirect = 1;
	if(numy%size ==0)
	{
		dataprocydirect = numy / size;
		domain[0] = size;
		domain[1] = dataprocydirect;
		domain[2] = 0;
		domain[3] = 0;
		std::cout<<"Domain Decomposed "<<std::endl;
	}
  else
  {
    if(numy%(size-1)!=0)
    {
      dataprocydirect = numy / (size-1);
      domain[0] = size-1;
      domain[1] = dataprocydirect;
      domain[2] = 1;
      domain[3] = numy- (domain[0] * domain[1]);
      std::cout<<"Domain Decomposed "<<std::endl;
    }
    else
    {
      dataprocydirect = numy / size;
      domain[0] = size-1;
      domain[1] = dataprocydirect;
      domain[2] = 1;
      domain[3] = numy - (domain[0]*domain[1]);
      std::cout <<"Domain Decomposed\n";
    }

  }
  
	
}



