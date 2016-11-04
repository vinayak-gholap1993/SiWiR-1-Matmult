
#include "IncludeFiles.hpp"
#include "Matrix.hpp"



int main(int argc, char *argv[])
{
try
{
  std::cout << BOLD(FBLU("\n ----------------- Program starts from here ------------- \n"));
    
  std::string A_inputFile=argv[1];
  std::string B_inputFile=argv[2];
  std::string C_outputFile=argv[3];
  
 siwir::Timer time; 
  
Matrix A,B,D;	// object of matrix class
Matrix call; 
 
  std::cout<<"no of parameter is... "<<argc<<std::endl;
  A.readFile(A_inputFile);
  B.readFile(B_inputFile);

  Matrix C(A.getNumberOfRow(),B.getNumberOfColumn());	//object create
  call.navierImplimentation(A,B,C);
  //time.reset();   //timer start
 // C = A * B;
  //real runtime = time.elapsed();  //runtime
  /*  
  time.reset();   //timer start
  call.optiImlimentation(A,B,C);
  real runtime1 = time.elapsed();  //runtime
  */
  //std::cout << BOLD(FGRN("\n runtime of navier is... ")) << runtime;
  //std::cout << BOLD(FGRN("\n runtime of opti. is... ")) << runtime1;
  
  
  C.writeFile(C_outputFile);
}
  
catch (const std::exception& e)
 {
   std::cout << BOLD(FRED(" a standard exception was caught, with message "))<< e.what()<<"\n";
 }
  
  std::cout << BOLD(FBLU("\n \n ------------------ end program ------------------"));
  std::cout<<"\n";

 return 0;  
}

