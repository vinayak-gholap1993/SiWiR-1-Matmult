
#include "IncludeFiles.hpp"
#include "Matrix.hpp"

#ifdef USE_LIKWID
extern "C" {
#include <likwid.h>
}
#endif


int main(int argc, char *argv[])
{
try
{
    std::cout << BOLD(FBLU("\n ----------------- Program starts from here ------------- \n"));
    
  std::string A_inputFile=argv[1];
  std::string B_inputFile=argv[2];
  std::string C_outputFile=argv[3];
  
  
Matrix A,B;	// object of matrix class
Matrix call; 
 std::cout<<" \n "<<B.getNumberOfColumn();
 
 std::cout<<"no of parameter is... "<<argc<<std::endl;
  A.readFile(A_inputFile);
  B.readFile(B_inputFile);

  Matrix C(A.getNumberOfRow(),B.getNumberOfColumn());	//object create
 
#ifdef USE_LIKWID
   likwid_markerInit();
   likwid_markerStartRegion( "vec" );
#endif   
 
   call.navierImplimentation(A,B,C);
 
#ifdef USE_LIKWID
   likwid_markerStopRegion( "vec" );
    likwid_markerClose();
#endif   
   
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

