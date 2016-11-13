# Simulation and Scientific Computing 1 
#Assignment 1
#Group: Mayank Patel, Vinayak Gholap

To run the code follow the steps. 
1. Open the terminal in the code directory. 
2. First type "make clean" (without the inverted commas) to remove all the existing object files (.o).
3. Then to compile all the files type "make" (without the inverted commas). This will compile all the .h and .cpp files. This is done by the Makefile. Appropriate flags are given in Makefile as required.
4. We now will have an executable file "matmult". 
5. Here to execute this file you will need to pass three arguments, type "./likwidScript". Here you can also give other input files instead of 2048*2048-1 and 2048*2048-2 inside the likwidScript file.
6. Now you will have the runtimes for different implimentations in seconds.
7. An output file named C.out will be generated. This file contains the number of Rows and number of Columns in initial line and data elements below. You can compare your solver solution with reference solution using- vimdiff C.out Ref.out. Here Ref.out should be the reference file you need to check the solution with the solver result.
8. A Performance report is present named "Performance.pdf".
9. We have also documented the code using Doxygen tool. To generate the html files and the latex files, type in terminal doxygen. This will generate necessary files automatically.

