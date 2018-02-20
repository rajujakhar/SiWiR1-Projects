// Include Header files 
#include <iostream> 
#include <fstream>
#include <string>
#include <vector>
#include <assert.h>
#include "Timer.h"

// Using Likwid 

#ifdef USE_LIKWID
extern "C"{ 
#include <likwid.h>
}
#endif

using namespace std;

// Main function passing argument
int main(int argc, char *argv[])

{
    
    //Declaring the vector
    
    if(argc != 4)
    cout<<"Please pass the executable, A matrix, B matrix and C matrix properly\n";

    vector<double> v_A;
    vector<double> v_B;
    vector<double> v_C;
    
    
    //Declaring the variables 
    double value_A;
    double value_B;
    
    
    // Declaring object of Fstream header file 
    ifstream file_A;
    ifstream file_B;
    ofstream file_C;
    

    // Opening The files 
    file_A.open(argv[1]);
    file_B.open(argv[2]);
    file_C.open(argv[3]);


    // Reading till untill end of file has reached 
    while(file_A>>value_A)
    {   
	    v_A.push_back(value_A);
    }
    // Reading the second file 
    
    while (file_B>> value_B)
    {
       v_B.push_back(value_B);
    }

    // Closing the files 
    file_A.close();
    file_B.close();
    
    // Fetching the number of rows and columns 

    int numRow_A= v_A[0];
    int numCol_A = v_A[1];
    int numCol_B = v_B[1];
    int numRow_B = v_B[0];
    int numRow_C = numRow_A;
    int numCol_C = numCol_B ;  
    int size_C = numRow_C*numCol_C;
    
    // Pushing thw rows and columns

    v_C.push_back(numRow_C);
    v_C.push_back(numCol_C);

    assert(numCol_A == numRow_B);

    siwir:: Timer timer;
    
    // Naive Matrix Matrix multiplication  O(n3)

    double timeTaken = 0.0;
    double temp = 0.0;
    int row = 0;
    int col = 0;
#ifdef USE_LIKWID
	likwid_markerInit();
	likwid_markerStartRegion("vector");
#endif
    
    for(int i=0; i< size_C ; ++i)
        {
            //cout<<"i is: "<<i <<endl;
            temp = 0.0;
            row = i/numCol_B;
            col = i % numCol_B; 

            for(int j=0 ; j < numCol_A ; ++j)
            {

                temp +=  v_A[2+ row*numCol_A + j] * v_B[2 + col+ j*numCol_B];
            }
          
            v_C.push_back(temp);
            
	   }

        timeTaken = timer.elapsed();
#ifdef USE_LIKWID
	
	likwid_markerStopRegion("vector");
	likwid_markerClose();
#endif
        
	file_C<< v_C[0]<<" "<<v_C[1]<<endl;

      cout<<"\nTime taken in matrix multiplication is:       " <<  timeTaken << endl;
        
	// Writing into the output file   
	  for (int i=2;i<(size_C + 2);i++)
        	file_C <<v_C[i]<<endl;
        	
        	
        	
        file_C.close();

    return 0;

}

