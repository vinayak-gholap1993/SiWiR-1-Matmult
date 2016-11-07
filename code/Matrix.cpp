
#include "Matrix.hpp"

Matrix::Matrix(const uint16_t& numRow, const uint16_t& numCol)
{
  std::cout<<"\n Matrix configuration starts... "; 
  
  this->_rows = numRow;		//assign no of _rows
  this->_columns = numCol;	//assign no of _columns
  
  std::cout<<"\n Number of rows is... "<<_rows<<"\t Number of columns is..."<<_columns; 
  /**
   * void resize (size_type n, const value_type& val) means new element n are initialized as copies of val
   ***/
  this->_vecObject.resize(_rows * _columns , 0.0);
 /* 
   for(u_int row=0; row<numRow ;++row)
	for(u_int col=0 ;col<numCol ;++col)
	  this->_vecObject[row * numCol + col] = row*numCol + col;
  */
  std::cout<<"\n Matrix Configured...\n";
}
/*
Matrix Matrix::operator*(const Matrix& B)
{
  assert(this->getNumberOfColumn() == B.getNumberOfRow());	///assert
  Matrix C(this->_rows,B._columns);		//temporary object of Matrix
  register real temp = 0.0;
 
  for( u_int i = 0; i <	this->_rows ; ++i)	//Matrix A 
    for( u_int j = 0; j < B._columns ; ++j)	//Matrix B
    {
      for( u_int k = 0; k< B._rows ; k++)	//Matrix C
      {	
	temp += (this->operator()(i,k)) * B(k,j);
      } //k
      C._vecObject[ i * B._columns + j] = temp;
      temp = 0.0;
    }  //j
 return C;    
}
*/
/*
//simple code that do not give good performance
void Matrix::NaiveImplimentation(const Matrix& A, const Matrix& B, Matrix& C)
{
  assert((A._columns == B._rows) || " Matrix A columns and B rows do not match -NaiveImplimentation ");	///assert
  register real temp = 0.0;
  
  //C._vecObject.resize(A._rows * B._columns,0.0);
  //std::cout<<" I am inside NaiveImplimentation "<<std::endl;
  
  for( u_int i = 0; i <	A._rows ; ++i)	//Matrix A
  {
    for( u_int j = 0; j < B._columns ; ++j)	//Matrix B
    {
      for( u_int k = 0; k< B._rows ; k++)	//Matrix C
      {	
	temp += A(i,k) * B(k,j);
      } //k
     
      C._vecObject[ i * B._columns + j] = temp;
      temp = 0.0;
    }  //j
  } //i
}
*/



void Matrix::BlockImplimentation(const Matrix& A, const Matrix& B, Matrix& C)
{
   u_int istart, iend, jstart, jend,kstart, kend, ii, jj, kk;
 
    for (ii = 0; ii < A._rows; ii += iblock){       //// blocking for i loop with blocksize of "iblock"
 
            istart = ii; iend = std::min(A._rows,ii+iblock);
 
            for ( kk = 0; kk < B._rows; kk += kblock) {     ////// blocking for k loop with blocksize of "kblock"
 
                kstart = kk; kend = std::min(B._rows,kk+kblock);
 
                for ( jj = 0; jj < B._columns; jj += jblock){       ////// blocking for j loop with blocksize of "jblock"
 
                                jstart = jj; jend = std::min(B._columns,jj+jblock);
 
 
                        for (u_int i = istart; i < iend ; i+=strideforBlocking) {
 
                            for (u_int k = kstart; k < kend; k++) {
 
                                for (u_int j = jstart; j < jend; j++) {

                                    // Loop unrolling with stride

                                                          C._vecObject[ i * B._columns + j]  += A(i,k) * B(k,j);
                                                          //C._vecObject[ (i+1) * B._columns + j] += A(i+1,k) * B(k,j);



                            } //j

                        } //k
                    } //i
                } //jj
            } //kk
     } //ii
 
}


/*
//simple code that do not give good performance
void Matrix::NaiveImplimentation(const Matrix& A, const Matrix& B, Matrix& C)
{
  assert((A._columns == B._rows) || " Matrix A columns and B rows do not match -navierImplimentation ");	///assert
  u_int istart, iend, jstart, jend,kstart, kend, ii, jj, kk;
  
  
  for( ii = 0; ii <A._rows; ii+=iblock)
    {                                       					// blocking for i loop with blocksize of "block"
       istart = ii; iend = std::min(A._rows,ii + iblock); 		// start value and end value of block in i loop
 
      for( jj = 0; jj <B._rows; jj+=jblock)
      {                                     					// blocking for j loop with blocksize of "block"
        jstart = jj; jend = std::min(B._columns,jj + jblock); 	// start value and end value of block in j loop
	  
	for( kk = 0; kk <B._rows; kk+=kblock)
	  {                                     					// blocking for j loop with blocksize of "block"
             kstart = kk; kend = std::min(B._rows,kk + kblock); 	// start value and end value of block in j loop
	   
	    for( u_int i = istart; i <= iend ; ++i)	//Matrix A
	      {
              for( u_int j = jstart; j < jend ; ++j)	//Matrix B
                {
                for( u_int k = kstart; k< kend ; k++)	//Matrix C
                  {
		    C._vecObject[ i * B._columns + j] += A(i,k) * B(k,j);
                  } //j
              }  //k
          } //i
	  } //kk
      } // jj
  } //ii
}
*/

void Matrix::NaiveImplimentation(const Matrix& A, const Matrix& B, Matrix& C)
{
  assert((A._columns == B._rows) || " Matrix A columns and B rows do not match -navierImplimentation ");	///assert
  register real temp1 = 0.0, temp2 = 0.0, temp3 = 0.0, temp4 = 0.0, temp5 = 0.0, temp6 = 0.0, temp7 = 0.0, temp8 = 0.0;
 
  //std::cout<<" I am inside NaiveImplimentation "<<std::endl;
  
  for( u_int i = 0; i <	A._rows ; ++i)	//Matrix A 
    for( u_int j = 0; j < B._columns ; j+=strideforNaive)	//Matrix B with stride for loop unrolling
    {
      for( u_int k = 0; k< B._rows ; k++)	//Matrix C
      {	
          //loop unrolling on j loop

        temp1 += A(i,k) * B(k,j);
        temp2 += A(i,k) * B(k,j+1);
        temp3 += A(i,k) * B(k,j+2);
        temp4 += A(i,k) * B(k,j+3);
        temp5 += A(i,k) * B(k,j+4);
        temp6 += A(i,k) * B(k,j+5);
        temp7 += A(i,k) * B(k,j+6);
        temp8 += A(i,k) * B(k,j+7);
      } //k
      C._vecObject[ i * B._columns + j] = temp1;
      C._vecObject[ i * B._columns + j+1] = temp2;
      C._vecObject[ i * B._columns + j+2] = temp3;
      C._vecObject[ i * B._columns + j+3] = temp4;
      C._vecObject[ i * B._columns + j+4] = temp5;
      C._vecObject[ i * B._columns + j+5] = temp6;
      C._vecObject[ i * B._columns + j+6] = temp7;
      C._vecObject[ i * B._columns + j+7] = temp8;
      temp1 = 0.0, temp2 =0, temp3 =0, temp4 = 0.0, temp5 = 0.0, temp6 = 0.0, temp7 = 0.0, temp8 = 0.0;
    }  //j 
}

// loop interchange gives better result compare to navierImplimentation
void Matrix::OptiImplimentation( const Matrix& A,const Matrix& B,Matrix& C)
{
  assert((A._columns == B._rows)|| " Matrix A columns and B rows do not match -optiImplimentation ");	///assert
//  register real temp = 0.0;
 
  for( u_int i = 0; i <	A._rows ; ++i)	//Matrix A 
    for( u_int k = 0; k< B._rows ; k++)   //Matrix B
    {
     for( u_int j = 0; j < B._columns ; ++j)	//Matrix B
      {	
	C._vecObject[ i * B._columns + j] += A(i,k) * B(k,j);
      } //j
      //C._vecObject[ i * B._columns + j] = temp;
      //temp = 0.0;
    }  //k  
}



//read data from file
bool Matrix::readFile(const std::string& fileName)
{
  std::ifstream file(fileName);
  
  if(file.is_open())
  {
    file>>_rows;	//get rows
    file>>_columns;	//get columns
    this->_vecObject.resize(_rows*_columns,0.0);	//allocate
    
    
    for(auto iter = this->_vecObject.begin();  iter<this->_vecObject.end();++iter)
      file >> *iter;
    
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
  //u_int temp = 0.0;
  if(file.is_open())
  {   
    file << this->_rows << " " << this->_columns;
    /*
      for(u_int i = 0; i< this->_rows ;++i)
    {
      temp = i * this->_rows;		//in order to reduce pressure on complier 
      for(u_int j=0 ;j< this->_columns ;++j)
      {
	file <<"\n" << this->_vecObject[temp + j];		//reading file may take some time during that time we initialized vector.
      }//j
    }//i*/
    
    for(auto iter = this->_vecObject.begin();  iter<this->_vecObject.end();++iter)
      file <<"\n" << *iter;
    
    std::cout<<"\n file writeing done...";
    return true;
  }//if
  else
  {
    std::cerr << BOLD(FRED("\n Unable to open file  ")) <<std::endl;
    return false;
  }
  file.close(); 
}

