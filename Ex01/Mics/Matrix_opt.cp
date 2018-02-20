// Include Header files 
#include <iostream> 
#include <fstream>
#include <string>
#include <vector>
#include <assert.h>
#include "Timer.h"

#ifdef USE_LIKWID
extern "C"{ 
#include <likwid.h>
}
#endif

using namespace std;


// Main function passing argument
int main(int argc, char *argv[])

{
      //Declaring the variables 

    if(argc != 4)
    cout<<"Please pass the executable, A matrix, B matrix and C matrix properly\n";  

    double value_A;
    double value_B;
    
    
    vector<double> v_A;
    vector<double> v_B;
    
    
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
        //file>> arr[n]; // putting it in array 
        v_A.push_back(value_A);
        //++n;

    }
    // Reading the second file 
    while (file_B>> value_B)
    {
        //Reading the contents of the second file 

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
    int size_B = numCol_A*numCol_B;
    
    //Declaring the vector
    
    if(argc != 4)
    cout<<"There is some error";
    
    vector<double> v_C(size_C + 2);
    vector<double> v_B_T(size_B + 2);

    v_C[0]=numRow_C;
    v_C[1] = numCol_C;
    
    
    for(int i=0;i<numCol_A;++i)
	{
	for(int j=0;j<numCol_B;++j)
		{
		v_B_T[2+j*numCol_A + i ]= v_B[2+i*numCol_B+j];
		}
	}

    assert(numCol_A == numRow_B);

    siwir:: Timer timer;
    double timeTaken;
    
    // Block Matrix multiplication  O(n^(2.78))
   
#ifdef USE_LIKWID
	likwid_markerInit();
	likwid_markerStartRegion("vector");
#endif
    int s = 16;
    double temp=0.0;

  // Blocking and Copy optimization technique . 
	

	 for(int jj=0;jj<numCol_B;jj+= s)
 	{
        	for(int kk=0;kk<numCol_A;kk+= s)
        	{
                	for(int i=0;i<numRow_A;i++)
                	         for(int j = jj; j<((jj+s)>numCol_B?numCol_B:(jj+s)); j++) 
                 	       {
                        	        temp = 0.0;
                                        for(int k = kk; k<((kk+s)>numCol_A?numCol_A:(kk+s)); k++)
                        	        {
							  
							  	    
						//temp += v_A[2+i*numCol_A + k]*v_B[2+k*numCol_B+ j];  //  B is read in row major format  
						temp += v_A[2+i*numCol_A + k]*v_B_T[2+j*numCol_A+ k];    // B is read in column major format                          
					        }
                        	       
                        	        v_C[2+i*numCol_C+ j]+= temp ;
                         }
             }
      }
	
	
        timeTaken = timer.elapsed();
#ifdef USE_LIKWID
	
	likwid_markerStopRegion("vector");
	likwid_markerClose();
#endif
	// Writing into output files 

	file_C<< v_C[0]<<" "<<v_C[1]<<endl;

      cout<<"\nTime taken in matrix multiplication is:       " <<  timeTaken << endl;
        for (int i=2;i<(size_C + 2);i++)
        	file_C <<v_C[i]<<endl;
        	
        	
     file_C.close();
     
     return 0;

}

