
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
  
Matrix A,B;	// object of matrix class
Matrix call; 
 
  std::cout<<"no of parameter is... "<<argc<<std::endl;
  A.readFile(A_inputFile);
  B.readFile(B_inputFile);

  Matrix C(A.getNumberOfRow(),B.getNumberOfColumn());	//object create

/*
  time.reset();   //timer start
  call.NaiveImplimentation(A,B,C);
  real runtime = time.elapsed();  //runtime
  std::cout << BOLD(FGRN("\n runtime of naive is... ")) << runtime;
  */
  /*
  time.reset();   //timer start
  call.OptiImplimentation(A,B,C);
  real runtime1 = time.elapsed();  //runtime
  std::cout << BOLD(FGRN("\n runtime of opti. is... ")) << runtime1;
  */
  

  time.reset();   //timer start
  call.BlockImplimentation(A,B,C);
  real runtime2 = time.elapsed();  //runtime
  std::cout << BOLD(FGRN("\n runtime of opti. is... ")) << runtime2;

  
  
 // std::cout<<C;
  

  
 // 
  //std::cout << BOLD(FGRN("\n runtime of blas is... ")) << runtimeBlas;
  
  
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

/*
#include "IncludeFiles.hpp"
#include "Matrix.hpp"


extern "C" {
#include <cblas.h>
}

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

//------------------------------------- Start likwid -------  
#ifdef USE_LIKWID
   likwid_markerInit();
#endif
//---------------------------------------- navierImplimentation ----------------------------------------

#ifdef USE_LIKWID
   likwid_markerStartRegion( "Naive" );
#endif  
   
   
 siwir::Timer time; 
  
Matrix A,B,D;	// object of matrix class
Matrix call; 
 
  std::cout<<"no of parameter is... "<<argc<<std::endl;
  A.readFile(A_inputFile);
  B.readFile(B_inputFile);

  Matrix C(A.getNumberOfRow(),B.getNumberOfColumn());	//object create

  const int lda = A._columns;
  const int ldb = 1024;
  const int ldc = 931;

  time.reset();   //timer start
  call.NaiveImplimentation(A,B,C);
  real runtime = time.elapsed();  //runtime
    
  
#ifdef USE_LIKWID
   likwid_markerStopRegion( "Naive" );
#endif
//---------------------------------------------- optiImlimentation --------------------------------------------- 

#ifdef USE_LIKWID
   likwid_markerStartRegion( "Opti" );
#endif
   
   
  time.reset();   //timer start
  call.OptiImplimentation(A,B,C);
  real runtime1 = time.elapsed();  //runtime

  
#ifdef USE_LIKWID
   likwid_markerStopRegion( "Opti" );
#endif

//---------------------------------------------- optiImlimentation --------------------------------------------- 

#ifdef USE_LIKWID
   likwid_markerStartRegion( "Block" );
#endif
   
   
  time.reset();   //timer start
  call.BlockImplimentation(A,B,C);
  real runtimeB = time.elapsed();  //runtime

  
#ifdef USE_LIKWID
   likwid_markerStopRegion( "Block" );
#endif
   
   
//------------------------------------------------------- Blas ------------------------------------------------------------   
 
#ifdef USE_LIKWID
   likwid_markerStartRegion( "Blas" );
#endif
   
   
  time.reset();   //timer start
  cblas_dgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, A._rows, A._columns, B._columns, 1.0, A._vecObject.data(), lda, B._vecObject.data(), ldb, 0.0, C._vecObject.data(), ldc );
  real runtimeBlas = time.elapsed();  //runtime

  
#ifdef USE_LIKWID
   likwid_markerStopRegion( "Blas" );
#endif
  
   
  std::cout << BOLD(FGRN("\n runtime of naive is... ")) << runtime;
  std::cout << BOLD(FGRN("\n runtime of opti. is... ")) << runtime1;
 std::cout << BOLD(FGRN("\n runtime of opti. is... ")) << runtimeB;
  std::cout << BOLD(FGRN("\n runtime of blas is... ")) << runtimeBlas;
  
  
  C.writeFile(C_outputFile);
}
  
catch (const std::exception& e)
 {
   std::cout << BOLD(FRED(" a standard exception was caught, with message "))<< e.what()<<"\n";
 }
  
  std::cout << BOLD(FBLU("\n \n ------------------ end program ------------------"));
  std::cout<<"\n";

  
//--------------------------------- close likwid -------------------------  
#ifdef USE_LIKWID
  // likwid_markerStopRegion( "likwid" );
   likwid_markerClose();
#endif   
     
//--------------------------------------------------------------------  
  
 return 0;  
}
*/
