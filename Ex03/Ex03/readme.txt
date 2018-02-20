// instruction to run the program

compile : type make
run : mpirun -np 4 ./cg 100 100 1000 0.0001
