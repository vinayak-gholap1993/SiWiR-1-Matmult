
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

//simple code that do not give good performance
void Matrix::NaiveImplimentation(const Matrix& A, const Matrix& B, Matrix& C)
{
  assert((A._columns == B._rows) || " Matrix A columns and B rows do not match -navierImplimentation ");	///assert
  register real temp = 0.0;
  //uint16_t block = 8000;
  
  for(u_int ii = 0; ii <A._rows; ii+=block)
    {                                       					// blocking for i loop with blocksize of "block"
      u_int istart = ii, iend = std::min(ii + block -1,A._rows); 		// start value and end value of block in i loop

      for(u_int jj = 0; jj <B._columns; jj+=block)
      {                                     					// blocking for j loop with blocksize of "block"

          u_int jstart = jj, jend = std::min(jj + block -1,B._columns); 	// start value and end value of block in j loop

          for( u_int i = istart; i < iend ; ++i)	//Matrix A
          {

              for( u_int j = jstart; j < jend ; ++j)	//Matrix B
              {
                  for( u_int k = 0; k< B._rows ; k++)	//Matrix C
                  {

                      temp += A(i,k) * B(k,j);
                  } //k
                  C._vecObject[ i * B._columns + j] = temp;
                  temp = 0.0;
              }  //j
          } //i
      } // jj
  } //ii
}

/*
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
*/

void Matrix::OptiImplimentation(const Matrix& A, const Matrix& B, Matrix& C)
{

  assert((A._columns == B._rows) || " Matrix A columns and B rows do not match -navierImplimentation ");	///assert
  register real temp = 0.0;
  //uint16_t block = 8000;
  
  for(u_int ii = 0; ii <A._rows; ii+=block)
    {                                       					// blocking for i loop with blocksize of "block"
      u_int istart = ii, iend = std::min(ii + block -1,A._rows); 		// start value and end value of block in i loop

      for(u_int kk = 0; kk <B._columns; kk+=block)
      {                                     					// blocking for j loop with blocksize of "block"

          u_int kstart = kk, kend = std::min(kk + block -1,B._columns); 	// start value and end value of block in j loop

          for( u_int i = istart; i < iend ; ++i)	//Matrix A
          {
	      temp = i * B._columns;
              for( u_int k = kstart; k < kend ; ++k)	//Matrix B
              {
                  for( u_int j = 0; j< B._rows ; j++)	//Matrix C
                  {

                      C._vecObject[ temp + j] += A(i,k) * B(k,j);
                  } //k
 
              }  //j
          } //i
      } // jj
  } //ii
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
  
  if(file.is_open())
  {   
    file << this->_rows << " " << this->_columns;
    
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

