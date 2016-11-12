
#include "Matrix.hpp"

Matrix::Matrix(const uint16_t& numRow, const uint16_t& numCol)
{
  //std::cout<<"\n Matrix configuration starts... "; 
  
  this->_rows = numRow;		//assign no of _rows
  this->_columns = numCol;	//assign no of _columns
  
 // std::cout<<"\n Number of rows is... "<<_rows<<"\t Number of columns is..."<<_columns; 
  /**
   * void resize (size_type n, const value_type& val) means new element n are initialized as copies of val
   ***/
  this->_vecObject.resize(_rows * _columns , 0.0);

 // std::cout<<"\n Matrix Configured...\n";
}

// loop interchange gives better result compare to NaiveImplimentation

void Matrix::NaiveImplimentation( const Matrix& A,const Matrix& B,Matrix& C)
{
  assert((A._columns == B._rows)|| " Matrix A columns and B rows do not match -NaiveImplimentation ");	///assert


  //std::cout<<" I am inside NaiveImplimentation "<<std::endl;

  for( u_int i = 0; i <	A._rows ; ++i)	//Matrix A
    for( u_int k = 0; k< B._rows ; k++)   //Matrix B
    {
     for( u_int j = 0; j < B._columns ; ++j)	//Matrix B
      {
        C._vecObject[ i * B._columns + j] += A(i,k) * B(k,j);
      } //j
    }  //k
}

// Loop unrolling Optimization
// j loop unrolled with unrolling factor =8

void Matrix::UnrollingImplimentation(const Matrix& A, const Matrix& B, Matrix& C)
{
  assert((A._columns == B._rows) || " Matrix A columns and B rows do not match -navierImplimentation ");	///assert
  register real temp1 = 0.0,temp2 = 0.0 , temp3 = 0.0,temp4 = 0.0,temp5 = 0.0 , temp6 = 0.0,temp7 = 0.0 ,temp8 = 0.0;


  for( u_int i = 0; i <	A._rows ; ++i)	//Matrix A

      for( u_int k = 0; k< B._rows ; ++k)	//Matrix C
      {
          for( u_int j = 0; j < B._columns ; j+=8)	//Matrix B with stride for loop unrolling
          {
              temp1 += A(i,k) * B(k,j);
              temp2 += A(i,k) * B(k,j+1);
              temp3 += A(i,k) * B(k,j+2);
              temp4 += A(i,k) * B(k,j+3);
              temp5 += A(i,k) * B(k,j+4);
              temp6 += A(i,k) * B(k,j+5);
              temp7 += A(i,k) * B(k,j+6);
              temp8 += A(i,k) * B(k,j+7);
              C._vecObject[ i * B._columns + j] = temp1 + temp2 + temp3 + temp4 + temp5 + temp6 + temp7 + temp8;
              } //j
              temp1 = 0.0 , temp2 = 0.0 ,temp3 = 0.0,temp4 = 0.0,temp5 = 0.0 , temp6 = 0.0,temp7 = 0.0 ,temp8 = 0.0;

            }  //k
        }

// Blocking Optimization
//Blocking the matrix with some block factors for each matrix

void Matrix::BlockImplimentation(const Matrix& A, const Matrix& B, Matrix& C)
{
   u_int j , istart, iend, jstart, jend,kstart, kend, ii, jj, kk ,iblock , jblock ,kblock;
   u_int ia = 30, ib = 300 , ic = 10;
   u_int x = A._rows , y = B._rows ,z = B._columns;

   iblock = std::min(ia,x); jblock = std::min(ib,y); kblock = std::min(ic,z);

    for (ii = 0; ii < A._rows; ii += iblock){       //// blocking for i loop with blocksize of "iblock"

            istart = ii; iend = std::min(A._rows,ii+iblock);

            for ( kk = 0; kk < B._rows; kk += kblock) {     ////// blocking for k loop with blocksize of "kblock"

                kstart = kk; kend = std::min(B._rows,kk+kblock);

                for ( jj = 0; jj < B._columns; jj += jblock){       ////// blocking for j loop with blocksize of "jblock"

                                jstart = jj; jend = std::min(B._columns,jj+jblock);

                        for (u_int i = istart; i < iend ; ++i) {

                            for (u_int k = kstart; k < kend; ++k) {

                                for ( j = jstart; j < jend; ++j) {

                                                          C._vecObject[ i * B._columns + j]  += A(i,k) * B(k,j);

                            } //j
                        } //k
                    } //i
                } //jj
            } //kk
     } //ii

}


//Blocking and Unrolling Optimization
//Blocking the matrix and unrolling the j loop for efficient utilization of Matrix B Cache data
// Also included the Unrolled Remainder j loop

void Matrix::BlockingAndUnrollingImplimentation(const Matrix& A, const Matrix& B, Matrix& C)
{
   u_int j , istart, iend, jstart, jend,kstart, kend, ii, jj, kk ,iblock , jblock ,kblock;
   u_int ia = 30 , ib = 300 , ic = 10;
   u_int x = A._rows , y = B._rows ,z = B._columns;
   
   iblock = std::min(ia,x); jblock = std::min(ib,y); kblock = std::min(ic,z);
   
    for (ii = 0; ii < A._rows; ii += iblock){       //// blocking for i loop with blocksize of "iblock"
 
            istart = ii; iend = std::min(A._rows,ii+iblock);
 
            for ( kk = 0; kk < B._rows; kk += kblock) {     ////// blocking for k loop with blocksize of "kblock"
 
                kstart = kk; kend = std::min(B._rows,kk+kblock);
 
                for ( jj = 0; jj < B._columns; jj += jblock){       ////// blocking for j loop with blocksize of "jblock"
 
                                jstart = jj; jend = std::min(B._columns,jj+jblock);
 
                        for (u_int i = istart; i < iend ; ++i) {
                            u_int temp = i * B._columns;
                            for (u_int k = kstart; k < kend; ++k) { 
			       
			     //We are doing unrolling 8 time so check if it is less than 7 than go to remaining loop
                                for ( j = jstart; j> (jend-7); j+=8) {

                                    // Loop unrolling with stride
                                                          C._vecObject[ temp + j]  += A(i,k) * B(k,j);
                                                          C._vecObject[ temp + (j+1)] += A(i,k) * B(k,j+1);
                                                          C._vecObject[ temp + (j+2)]  += A(i,k) * B(k,j+2);
                                                          C._vecObject[ temp + (j+3)] += A(i,k) * B(k,j+3);
                                                          C._vecObject[ temp + (j+4)]  += A(i,k) * B(k,j+4);
                                                          C._vecObject[ temp + (j+5)] += A(i,k) * B(k,j+5);
                                                          C._vecObject[ temp + (j+6)]  += A(i,k) * B(k,j+6);
                                                          C._vecObject[ temp + (j+7)] += A(i,k) * B(k,j+7);



                            } //j
                            //std::cout << "\n remaining loop " << (jend - unrollJEND);
			    for(u_int ju = j ; ju < jend ; ju+=4)
			    {
                                C._vecObject[ temp + ju]  += A(i,k) * B(k,ju);
                                C._vecObject[ temp + (ju+1)]  += A(i,k) * B(k,ju+1);
                                C._vecObject[ temp + (ju+2)]  += A(i,k) * B(k,ju+2);
                                C._vecObject[ temp + (ju+3)]  += A(i,k) * B(k,ju+3);
			    } //new j 
                        } //k
                    } //i
                } //jj
            } //kk
     } //ii
 
}


//read data from file
bool Matrix::readFile(const std::string& fileName)
{
  std::ifstream file(fileName);
  register u_int temp = 0;
  
  
  if(file.is_open())
  {
    file>>_rows;	//get rows
    file>>_columns;	//get columns
    
    //newColumns = _columns + pad;
    
    this->_vecObject.resize(_rows * _columns,0.0);	//allocate
    
    for(u_int i = 0;i < _rows ;++i)
    {
      temp = i * _columns;
      for(u_int j = 0;j < _columns ; ++j)
      {
	file >> this->_vecObject[ temp + j ]; 
      }
    }


    return true;
  }//if
  else
  {
   std::cerr << BOLD(FRED("\n Unable to open file  ")) <<std::endl;
    return false;
  }
  file.close();
}

// write data to file
bool Matrix::writeFile(const std::string& fileName)
{
  std::ofstream file(fileName);		//object of ofstream
  //register u_int newColumns = 0 , temp = 0;

  if(file.is_open())
  {   
    file << this->_rows << " " << _columns<<"\n";
    
    for(auto iter = this->_vecObject.begin();  iter<this->_vecObject.end();++iter)
      file << *iter <<"\n";
    
   //std::cout<<"\n file writeing done...";
    return true;
  }//if
  else
  {
    std::cerr << BOLD(FRED("\n Unable to open file  ")) <<std::endl;
    return false;
  }
  file.close(); 
}

