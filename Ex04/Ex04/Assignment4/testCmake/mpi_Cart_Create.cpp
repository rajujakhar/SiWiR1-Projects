#include <iostream>
#include <mpi.h>
#include <vector>

int main(int argc, char *argv[])
{
	int size,rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	//std::cout << "No of processors: " << size;

std::cout << "Initiating Values \n";
    
    auto ndims = 2;
    std::vector<int> dim,periods,coords;
    dim.push_back(1);
    dim.push_back(size);
    
    for(int i=0;i <ndims; i++)
    {
    periods.push_back(0);
    coords.push_back(0);
    }
    
    int reorder = 1;
    int up=0;
    int down = 0;

    MPI_Comm new_comm(MPI_COMM_NULL);
    MPI_Cart_create(MPI_COMM_WORLD,ndims,dim.data(),periods.data(), reorder, &new_comm);
    //MPI_Cart_create(MPI_COMM_WORLD,ndims,&dim[0],&periods[0],reorder,&new_comm);
    
    std::cout << "Cart_Create Successfull\n";

	MPI_Finalize();

}
