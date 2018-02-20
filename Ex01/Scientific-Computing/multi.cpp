#include <iostream>
#include <ctime>
#include <cmath>

int main()
{
long size = 4096*4*4*4;
double a[size]; 
double b[size]; 
double c[size];

for(long i=0;i<size;i++)
	{
	a[i] = i%3;
	b[i] = i%2;	
	c[i] = 0.0;
	}
int N = 64*2*2*2;
int s = 16*2*2*2; 
double temp = 0.0;
time_t start , end ,diff;
start = time (NULL); 

/*for(int jj=0;jj<N;jj+= s){
        for(int kk=0;kk<N;kk+= s){
                for(int i=0;i<N;i++){
                        for(int j = jj; j<((jj+s)>N?N:(jj+s)); j++){
                                temp = 0;
                                for(int k = kk; k<((kk+s)>N?N:(kk+s)); k++){
                                        //temp += a[i][k]*b[k][j];
					temp += a[i*N + k]*b[k*N+ j];	                                
					}
                                c[i*N+ j] += temp;
                        }
                }
        }
}*/

        

end = time (NULL);
diff = difftime(end, start);
double dif = diff*1000;
std::cout<<dif<<std::endl;

for(int i=0;i<N;++i){
	for(int j=0;j<N;++j){
		//std::cout<<c[i*N +j]<<"\t";
		}
	//std::cout<<std::endl;
	}

return 0;


}
