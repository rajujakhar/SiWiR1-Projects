#include<iostream>
#include<vector>

//This function stores red and black index for red black Guass Siedel Method
 void calIndexes(const int numGridX, const int numGridY, std::vector<int> &redIndex,std::vector<int> &blackIndex)
{
 bool toggle = false;
 
 // We have to exclude boundary indexed in the red and black array
   for(size_t i=1; i<numGridY - 1; ++i)
   {
   		toggle = !toggle;
  	 	for(size_t j=1 ;j<numGridX-1;j+=2)
   		{
   			if(toggle)
           {
             redIndex.push_back(i*numGridX +j);
             if((numGridX % 2) == 1)
              {
                if (!((j+2) % numGridX == 0)) blackIndex.push_back(1+i*numGridX +j);
              } 
              else  blackIndex.push_back(1+i*numGridX +j);
           }
         else
           {
           	 blackIndex.push_back(i*numGridX +j);
             if((numGridX % 2) == 1)
              {
              if (!((j+2) % numGridX == 0)) redIndex.push_back(1+ i*numGridX +j);
              }
              else redIndex.push_back(1+ i*numGridX +j);
           }

   		} 
   }

}


int main()
{

	int numGridX,numGridY;
	std::cout<<"Please Enter Grid Point in X and Y direction respectively \n";
	std::cin>>numGridX >> numGridY ;

	std::vector<int> redInd;
	std::vector<int> blackInd;
	calIndexes(numGridX,numGridY,redInd, blackInd);

	std::cout << "Red Indexs are : \n";
	for(size_t i =0; i<redInd.size() ; ++i)
	{
	std::cout << redInd[i] << "  " ;
	}

	std::cout << "\nBlack Indexs are : \n";
	for(size_t i =0; i<blackInd.size() ; ++i)
	{
	std::cout << blackInd[i] << "  " ;
	}

    std::cout<< std::endl; 

	
	return 0;
}