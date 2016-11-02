
#include "Matrix.hpp"

Matrix::Matrix(const u_int& numRow, const u_int& numCol)
{
  std::cout<<"\n Matrix configuration starts... "; 
  
  this->_rows = numRow;		//assign no of _rows
  this->_columns = numCol;	//assign no of _columns
  
  std::cout<<"\n Number of rows is... "<<_rows<<"\t Number of columns is..."<<_columns; 
  /**
   * void resize (size_type n, const value_type& val) means new element n are initialized as copies of val
   ***/
  this->_vecObject.resize(_rows * _columns , 0.0);
  
   for(u_int row=0; row<numRow ;++row)
	for(u_int col=0 ;col<numCol ;++col)
	  this->_vecObject[row * numCol + col] = row*numCol + col;
  
  std::cout<<"\n Matrix Configured...\n";
}

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

//simple code that do not give good performance
void Matrix::navierImplimentation(const Matrix& A, const Matrix& B, Matrix& C)
{
  assert((A._columns == B._rows) || " Matrix A columns and B rows do not match -navierImplimentation ");	///assert
  register real temp = 0.0;
  
  //C._vecObject.resize(A._rows * B._columns,0.0);
  //std::cout<<" I am inside navierImplimentation "<<std::endl;
  
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
/*
// loop interchange gives better result compare to navierImplimentation
void optiImlimentation( const Matrix& A,const Matrix& B,Matrix& C)
{
  assert((A._columns == B._rows)|| " Matrix A columns and B rows do not match -optiImplimentation ");	///assert
  register real temp = 0.0;
 
  for( u_int i = 0; i <	A._rows ; ++i)	//Matrix A 
    for( u_int k = 0; k< B._rows ; k++)   //Matrix B
    {
     for( u_int j = 0; j < B._columns ; ++j)	//Matrix B
      {	
	temp += A(i,k) * B(k,j);
      } //j
      C._vecObject[ i * B._columns + j] = temp;
      temp = 0.0;
    }  //k  
}
*/
//read data from file
bool Matrix::readFile(const std::string& fileName)
{
  std::ifstream file(fileName);
  
  if(file.is_open())
  {
    file>>_rows;	//get rows
    file>>_columns;	//get columns
    register int temp = 0;
    this->_vecObject.resize(_rows*_columns,0.0);	//allocate
    
    for(u_int i = 0; i< _rows ;++i)
    {
      temp = i * _columns;		//in order to reduce pressure on complier
      for(u_int j=0 ;j< _columns ;++j)
      {
	file>>this->_vecObject[temp + j];		//reading file may take some time during that time we initialized vector.
      }//j
    }//i
    std::cout<<"\n file reading done...";
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
  register int temp = 0;
  std::ofstream file(fileName);		//object of ofstream
  
  if(file.is_open())
  {   
    file << this->_rows << "\t" << this->_columns;
    
    for(u_int i = 0; i< this->_rows ;++i)
    {
      temp = i * this->_rows;		//in order to reduce pressure on complier 
      for(u_int j=0 ;j< this->_columns ;++j)
      {
	file <<"\n" << this->_vecObject[temp + j];		//reading file may take some time during that time we initialized vector.
      }//j
    }//i

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

/*
bool Matrix::readFile(const std::string& fileName)
{
  u_int dummyRow,dummyCol,j=0;
  register int temp = 0;
  std::cout<<" file name in readFile function "<< filename <<std::endl;
  std::ifstream file("A.in");
  
  if(file.is_open())
  {
    file>>dummyRow;
    file>>dummyCol;
    
    this->_rows = dummyRow;this->_columns = dummyCol;
    
    this->_vecObject.resize(dummyRow*dummyCol,0.0);	//allocate
    
    for(u_int i = 0; i< dummyRow ;++i)
    {
      //remove bubble in multiply pipeline
      temp = i * dummyCol;		//in order to reduce pressure on complier
      //this->_vecObject[temp + j] = 0.0;
      //this->_vecObject[temp + j + 1] = 0.0;
      
      for(j=0 ;j< dummyCol ;++j)
      {
	//this->_vecObject[temp + j + 2] = 0.0;
	//this->_vecObject[temp + j + 3] = 0.0;
	file>>this->_vecObject[temp + j];		//reading file may take some time during that time we initialized vector.
      }//j
      j = 0;
    }//i
    std::cout<<"file reading done...";
    return true;
  }//if
  else
  {
   std::cerr << BOLD(FRED("\n Unable to open file  ")) <<std::endl;
    return false;
  }
  file.close();
}

bool Matrix::writeFile(const std::string& fileName)
{
  register int temp = 0;
  std::cout<<" file name in writeFile function "<< filename <<std::endl;
  std::ofstream file("C.in");		//object of ofstream
  
  if(file.is_open())
  {   
    std::cout<<"\n\n"<< this->_rows << "\t" << this->_columns;
    file << this->_rows << "\t" << this->_columns;
    
    for(u_int i = 0; i< this->_rows ;++i)
    {
      //remove bubble in multiply pipeline
      temp = i * this->_rows;		//in order to reduce pressure on complier
      //this->_vecObject[temp + j] = 0.0;
      //this->_vecObject[temp + j + 1] = 0.0;
      
      for(u_int j=0 ;j< this->_columns ;++j)
      {
	//this->_vecObject[temp + j + 2] = 0.0;
	//this->_vecObject[temp + j + 3] = 0.0;
	file <<"\n" << this->_vecObject[temp + j];		//reading file may take some time during that time we initialized vector.
      }//j
    }//i

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


*/