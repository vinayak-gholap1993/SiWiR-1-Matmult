

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

//---------------------------------------------------------------Strassen Implimentation--------------------------------------------------------------------------



/**
 * multiplication of two matrix with blocking and unrolling
 * @param ::vector&A object of vector
 * @param ::vector&B object of vector
 * @param ::vector&C object of vector
 * @param loopIter blocking size
 * **/
void matrixMulti(::vector& A,::vector& B,::vector& C,const u_int& loopIter)
{
  {
   u_int j , istart, iend, jstart, jend,kstart, kend, ii, jj, kk ,iblock , jblock ,kblock;
   u_int ia = 30 , ib = 300 , ic = 10;
   u_int x = loopIter , y = loopIter ,z = loopIter;

   iblock = std::min(ia,x); jblock = std::min(ib,y); kblock = std::min(ic,z);

    for (ii = 0; ii < loopIter; ii += iblock){       //// blocking for i loop with blocksize of "iblock"

            istart = ii; iend = std::min(loopIter,ii+iblock);

            for ( kk = 0; kk < loopIter; kk += kblock) {     ////// blocking for k loop with blocksize of "kblock"

                kstart = kk; kend = std::min(loopIter,kk+kblock);

                for ( jj = 0; jj < loopIter; jj += jblock){       ////// blocking for j loop with blocksize of "jblock"

                                jstart = jj; jend = std::min(loopIter,jj+jblock);

                        for (u_int i = istart; i < iend ; ++i) {

                            for (u_int k = kstart; k < kend; ++k) {

                             //We are doing unrolling 8 time so check if it is less than 7 than go to remaining loop
                                for ( j = jstart; j > ((jend-7) && j < jend); j+=8) {

                                    // Loop unrolling with stride
                                       C[ i * loopIter + j]  += A[ i*loopIter + k] * B[ k*loopIter + j];
                                       C[ i * loopIter + (j+1)]  += A[ i*loopIter + k] * B[ k*loopIter + (j+1)];
                                       C[ i * loopIter + (j+2)]  += A[ i*loopIter + k] * B[ k*loopIter + (j+2)];
                                       C[ i * loopIter + (j+3)]  += A[ i*loopIter + k] * B[ k*loopIter + (j+3)];

                                       C[ i * loopIter + (j+4)]  += A[ i*loopIter + k] * B[ k*loopIter + (j+4)];
                                       C[ i * loopIter + (j+5)]  += A[ i*loopIter + k] * B[ k*loopIter + (j+5)];
                                       C[ i * loopIter + (j+6)]  += A[ i*loopIter + k] * B[ k*loopIter + (j+6)];
                                       C[ i * loopIter + (j+7)]  += A[ i*loopIter + k] * B[ k*loopIter + (j+7)];

                            } //j
                            //std::cout << "\n remaining loop " << (jend - unrollJEND);
                            for(u_int ju = j ; ju < jend ; ju+=4)
                            {	      C[ i * loopIter + ju]  += A[ i*loopIter + k] * B[ k*loopIter + ju];
                                      C[ i * loopIter + (ju+1)]  += A[ i*loopIter + k] * B[ k*loopIter + (ju+1)];
                                      C[ i * loopIter + (ju+2)]  += A[ i*loopIter + k] * B[ k*loopIter + (ju+2)];
                                      C[ i * loopIter + (ju+3)]  += A[ i*loopIter + k] * B[ k*loopIter + (ju+3)];
                            } //new j
                        } //k
                    } //i
                } //jj
            } //kk
     } //ii

}
}

/**
 * Addition of two matrix
 * @param ::vector&A object of vector
 * @param ::vector&B object of vector
 * @param ::vector&C object of vector
 * @param loopIter blocking size
 * **/
void sum(::vector& A, ::vector& B, ::vector& C,const u_int& loopIter)
{
  //std::cout<<"\n I am in sum function... ";
   for(u_int i = 0 ;i < loopIter ; ++i)
   {
      for(u_int j = 0 ;j < loopIter ; j+=8)
      {
        //if (loopIter==2) std::cout<<"\n "<<A[i * loopIter + j] <<" "<< B[i * loopIter + j];

        C[i * loopIter + j] = A[i * loopIter + j] + B[i * loopIter + j];
        C[i * loopIter + (j+1)] = A[i * loopIter + (j+1)] - B[i * loopIter + (j+1)];
        C[i * loopIter + (j+2)] = A[i * loopIter + (j+2)] - B[i * loopIter + (j+2)];
        C[i * loopIter + (j+3)] = A[i * loopIter + (j+3)] - B[i * loopIter + (j+3)];

        C[i * loopIter + (j+4)] = A[i * loopIter + (j+4)] - B[i * loopIter + (j+4)];
        C[i * loopIter + (j+5)] = A[i * loopIter + (j+5)] - B[i * loopIter + (j+5)];
        C[i * loopIter + (j+6)] = A[i * loopIter + (j+6)] - B[i * loopIter + (j+6)];
        C[i * loopIter + (j+7)] = A[i * loopIter + (j+7)] - B[i * loopIter + (j+7)];
        //std::cout<<"\n "<<A[i * loopIter + j] <<" "<< B[i * loopIter + j];
      } //j
   }//i
}

/**
 * Substraction of two matrix
 * @param ::vector&A object of vector
 * @param ::vector&B object of vector
 * @param ::vector&C object of vector
 * @param loopIter blocking size
 * **/
void sub(::vector& A, ::vector& B, ::vector& C,const u_int& loopIter)
{
  //std::cout<<"\n I am in sub function... ";
   for(u_int i = 0 ;i < loopIter ; ++i)
   {
      for(u_int j = 0 ;j < loopIter ; j+=8)
      {
        //if (loopIter==2) std::cout<<"\n "<<A[i * loopIter + j] <<" "<< B[i * loopIter + j];

        C[i * loopIter + j] = A[i * loopIter + j] - B[i * loopIter + j];
        C[i * loopIter + (j+1)] = A[i * loopIter + (j+1)] - B[i * loopIter + (j+1)];
        C[i * loopIter + (j+2)] = A[i * loopIter + (j+2)] - B[i * loopIter + (j+2)];
        C[i * loopIter + (j+3)] = A[i * loopIter + (j+3)] - B[i * loopIter + (j+3)];

        C[i * loopIter + (j+4)] = A[i * loopIter + (j+4)] - B[i * loopIter + (j+4)];
        C[i * loopIter + (j+5)] = A[i * loopIter + (j+5)] - B[i * loopIter + (j+5)];
        C[i * loopIter + (j+6)] = A[i * loopIter + (j+6)] - B[i * loopIter + (j+6)];
        C[i * loopIter + (j+7)] = A[i * loopIter + (j+7)] - B[i * loopIter + (j+7)];

      } //j
   }//i
}

//be care fully pass vector not matrix object
void strassen( ::vector& A , ::vector& B , ::vector& C ,const u_int& rowA ,const u_int& colA ,const u_int& colB )
{
  u_int block = rowA; 	//rowA = 1024
 u_int x = colA , y = colB;
  u_int unblock = x , size;  //unblock index is fix as no of col and row are fix and A , B index fix so...
  unblock = y;
  if( block <= 64 )
  {
    //std::cout<<"\n I am ready for multiplication "<<std::endl;
    matrixMulti(A,B,C,block);
  }
  else
  {

    unblock = block; 	//to access vectors of elements
    block /= 2;		//to run loop
    size = block * block;

    //std::cout<<"\n I am not ready for multiplication...\n ";
    //divide matrix into 4 sub matrixes
    ::vector a11, a12, a21, a22,b11, b12, b21, b22 , c11, c12, c21, c22,addAresult, addBresult , M1 , M2 , M3 , M4 , M5 , M6 , M7;

    a11.resize(size,0.0),a12.resize(size,0.0),a21.resize(size,0.0),a22.resize(size,0.0),
    b11.resize(size,0.0),b12.resize(size,0.0),b21.resize(size,0.0),b22.resize(size,0.0),addAresult.resize(size,0.0),addBresult.resize(size,0.0),
    c11.resize(size,0.0),c12.resize(size,0.0),c21.resize(size,0.0),c22.resize(size,0.0),
    M1.resize(size,0.0),M2.resize(size,0.0),M3.resize(size,0.0),M4.resize(size,0.0),M5.resize(size,0.0),M6.resize(size,0.0),M7.resize(size,0.0);

                            //fill-up vector with data
                          for(u_int i = 0 ; i < block ; ++i)
                            for(u_int j = 0; j < block ; ++j)
                            {
                              a11[ i * block + j] = A[ i * unblock + j];
                              a12[ i * block + j] = A[ (i * unblock) + (j + block)];
                              a21[ i * block + j] = A[ (i+block) * unblock + j];
                              a22[ i * block + j] = A[ (i+block) * unblock + (j + block)];

                              b11[ i * block + j] = B[ i * unblock + j];
                              b12[ i * block + j] = B[ i * unblock + (j + block)];
                              b21[ i * block + j] = B[ (i+block) * unblock + j];
                              b22[ i * block + j] = B[ (i+block) * unblock + (j + block)];

                            }

                           sum(a11,a22,addAresult,block);
                           sum(b11,b22,addBresult,block);
                           strassen(addAresult,addBresult,M1,block,block,block);

                          sum(a21,a22,addAresult,block);	// a21 + a22
                          strassen(addAresult,b11,M2,block,block,block);	//p2 = (a21 + a22) * (b11)


                          sub(b12, b22, addBresult, block); // b12 - b22
                          strassen(a11, addBresult, M3, block,block,block); // p3 = (a11) * (b12 - b22)

                          sub(b21, b11, addBresult, block); // b21 - b11
                          strassen(a22, addBresult, M4, block,block,block); // p4 = (a22) * (b21 - b11)

                          sum(a11, a12, addAresult, block); // a11 + a12
                          strassen(addAresult, b22, M5, block,block,block); // p5 = (a11+a12) * (b22)

                          sub(a21, a11, addAresult, block); // a21 - a11
                          sum(b11, b12, addBresult, block); // b11 + b12
                          strassen(addAresult, addBresult, M6, block,block,block); // p6 = (a21-a11) * (b11+b12)

                          sub(a12, a22, addAresult, block); // a12 - a22
                          sum(b21, b22, addBresult, block); // b21 + b22
                          strassen(addAresult, addBresult, M7, block,block,block); // p7 = (a12-a22) * (b21+b22)

                          //finally add all partial results
                          sum(M3, M5, c12, block); // c12 = p3 + p5
                          sum(M2, M4, c21, block); // c21 = p2 + p4

                          sum(M1, M4, addAresult, block); // p1 + p4
                          sum(addAresult, M7, addBresult, block); // p1 + p4 + p7
                          sub(addBresult, M5, c11, block); // c11 = p1 + p4 - p5 + p7

                          sum(M1, M3, addAresult, block); // p1 + p3
                          sum(addAresult, M6, addBresult, block); // p1 + p3 + p6
                          sub(addBresult, M2, c22, block); // c22 = p1 + p3 - p2 + p6


                          for(u_int i = 0 ; i < block ; ++i)
                            for(u_int j = 0; j < block ; ++j)
                            {
                               C[ i * unblock + j] = c11[ i * block + j];
                               C[ (i * unblock) + (j + block)] = c12[ i * block + j];
                               C[ (i+block) * unblock + j] = c21[ i * block + j];
                               C[ (i+block) * unblock + (j + block)] = c22[ i * block + j];
                            }

  } //else

} //function
//---------------------------------------------------------------------------------------------------------------------------------------------------------------



int main(int argc, char *argv[])
{

  
try
{
  std::cout << BOLD(FBLU("\n ----------------- Program starts from here ------------- \n"));
    
   
  siwir::Timer time; 
  
  Matrix A,B;	// object of matrix class
  Matrix call; 
 
  if( argc != 4 ) 	std::cerr << BOLD(FRED(" Enter Appropriate Number of Arguments "))<<std::endl;
  else
  {
  //std::cout<<"no of parameter is... "<<argc<<std::endl;
  A.readFile(argv[1]);
  B.readFile(argv[2]);

  
  Matrix C(A.getNumberOfRow(),B.getNumberOfColumn());	//object create
  Matrix C_U(A.getNumberOfRow(),B.getNumberOfColumn());
  Matrix C_B(A.getNumberOfRow(),B.getNumberOfColumn());
  Matrix C_U_B(A.getNumberOfRow(),B.getNumberOfColumn());
  

  Matrix dummy(A.getNumberOfRow(),B.getNumberOfColumn());

      //heating up processor
    for(u_int i = 0; i<2 ; ++i)
        C.BlockImplimentation(A,B,dummy);

    //------------------------------------- Start likwid -------
  #ifdef USE_LIKWID
     likwid_markerInit();
  #endif

  
//---------------------------------------- NaiveImplimentation ----------------------------------------

#ifdef USE_LIKWID
   likwid_markerStartRegion( "Naive" );
#endif

  time.reset();   //timer start
  call.NaiveImplimentation(A,B,C);
  //real runtimeN = time.elapsed();  //runtime
  //std::cout << BOLD(FGRN("\n runtime for Naive ...		")) << runtimeN;

#ifdef USE_LIKWID
   likwid_markerStopRegion( "Naive" );
#endif

  //---------------------------------------------- UnrollingImplimentation ---------------------------------------------
#ifdef USE_LIKWID
   likwid_markerStartRegion( "Unrolling" );
#endif
   
  time.reset();   //timer start
  call.UnrollingImplimentation(A,B,C_U);
  //real runtimeU = time.elapsed();  //runtime
  //std::cout << BOLD(FGRN("\n runtime for Unroll ...	        ")) << runtimeU;

#ifdef USE_LIKWID
   likwid_markerStopRegion( "Unrolling" );
#endif
   //---------------------------------------------- BlockImplimentation ---------------------------------------------
#ifdef USE_LIKWID
   likwid_markerStartRegion( "Block" );
#endif
 
  time.reset();   //timer start
  call.BlockImplimentation(A,B,C_B);
  //real runtimeB = time.elapsed();  //runtime
  //std::cout << BOLD(FGRN("\n runtime for Blocking ...	")) << runtimeB;

#ifdef USE_LIKWID
   likwid_markerStopRegion( "Block" );
#endif

//---------------------------------------------- BlockingAndUnrollingImplimentation ---------------------------------------------

#ifdef USE_LIKWID
   likwid_markerStartRegion( "Blocking&Unrolling" );
#endif

     time.reset();   //timer start
     call.BlockingAndUnrollingImplimentation(A,B,C_U_B);
     //real runtimeBU = time.elapsed();  //runtime
     //std::cout << BOLD(FGRN("\n runtime for Block And Unroll...")) << runtimeBU;

#ifdef USE_LIKWID
   likwid_markerStopRegion( "Blocking&Unrolling" );
#endif

//---------------------------------------- StrassenImplimentation ----------------------------------------
  


   if((A._columns == A._rows) && (B._columns == B._rows))
   {

#ifdef USE_LIKWID
   likwid_markerStartRegion( "Strassen" );
#endif

      Matrix C_S(A.getNumberOfRow(),B.getNumberOfColumn());

        time.reset();   //timer start
        strassen( A._vecObject , B._vecObject , C_S._vecObject , A._rows , A._columns , B._columns);
        //real runtimeS = time.elapsed();  //runtime
        //std::cout << BOLD(FGRN("\n runtime for Strassen ...	")) << runtimeS;

#ifdef USE_LIKWID
   likwid_markerStopRegion( "Strassen" );
#endif


//------------------------------------------------------- Blas ------------------------------------------------------------   
 
#ifdef USE_LIKWID
   likwid_markerStartRegion( "Blas" );
#endif

  time.reset();   //timer start
  cblas_dgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, A._rows, A._columns, B._columns, 1.0, A._vecObject.data(), A._columns, B._vecObject.data(), B._columns, 0.0, C._vecObject.data(), C._columns );
  //real runtimeBlas = time.elapsed();  //runtime
  //std::cout << BOLD(FGRN("\n runtime for BLAS ...		")) << runtimeBlas;

#ifdef USE_LIKWID
   likwid_markerStopRegion( "Blas" );
#endif

   } //if strassen
  
   else
   {
                std::cout << BOLD(FMAG(" Strassen and BLAS needs Square matrix ... "))<<std::endl;
   }

   C_U_B.writeFile(argv[3]);	//output data file
    
  } //else
} //try
  
catch (const std::exception& e)
 {
   std::cerr << BOLD(FRED(" a standard exception was caught, with message "))<< e.what()<<std::endl;
 }
  
  std::cout << BOLD(FBLU("\n \n ------------------ end program ------------------"));
  std::cout<<"\n";

#ifdef USE_LIKWID
   likwid_markerClose();
#endif

  
 return 0;  
}





