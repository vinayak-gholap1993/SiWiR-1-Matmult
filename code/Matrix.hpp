
/**
 * All matrix operation done in this folder
 * 
 **/
#pragma once

#include "IncludeFiles.hpp"

class Matrix
{  
//------------------------------- friend ---------------------------------------  

  /**
     *  Print matrix
     *	@param object of ostream(standard class) with reference.
     *	@param object of Matrix class with reference.
     * 	@return object of ostream(standard class) with reference.
     * 
    **/

  friend std::ostream& operator<<(std::ostream& out,const Matrix& obj)
 {  
  for(uint16_t row= 0 ; row<obj._rows ;++row)
    {
    for(uint16_t col= 0 ;col<obj._columns  ;++col)
      {
	out<< obj._vecObject[row * obj._columns + col]<<"\t"; 
      }
    std::cout<<std::endl;
    }
return out;
 }

//------------------------------- Public functions of class ------------------------

public:
  /**
   * Default Constructor
   * **/
  Matrix(){};
    /**
     * Cosntructor for assign dimension of matrix and all element set to 0;
     * @param numRow Number of row in matrix
     * @param numCol Number of column in matrix
     * **/
  Matrix(const uint16_t& numRow, const uint16_t& numCol);
    /**
     * Default Destructor 
     * **/
  ~Matrix(){};
    
 //--------------------------------------- Some Basic Functions ------------------------------

    /**
     * 	Get number of rows
     * **/
  inline uint16_t getNumberOfRow(void)const
    {
      //std::cout<< " number of rows "<<this->_rows;
      return this->_rows;
    }
    
    /**
     * 	Get number of columns
     * **/
  inline u_int getNumberOfColumn(void)const
    {
      //std::cout<< " number of rows "<<this->_columns;
      return this->_columns;
    }
    /**
     * Check matrix is square or not
     * **/
  inline bool isSquareMatrix(void)const
    {
      return (this->_rows == this->_columns);
    }

 //----------------------  Printing Functions using call operator overloading --------------------------------------------
 /**
  * This function need to change instead of mutliply in Printing time think in for loop
  * @param rowIndexStart index of row from where we want start Printing
  * @param rowIndexend index of row from where we want end Printing
  * @param colIndexStart index of column from where we want start Printing
  * @param colIndexend index of column from where we want end Printing
  * 
  * **/
  void operator()(const u_int &rowIndexStart, const u_int &rowIndexend, const u_int &colIndexStart, const u_int &colIndexend)
    {
      std::cout<<"\n";
      for(u_int row=rowIndexStart; row<=rowIndexend ;++row)
      {
	for(u_int col=colIndexStart ;col<=colIndexend ;++col)
	{
	  std::cout<<this->_vecObject[row * this->_columns + col]<<"\t";
	}
       std::cout<<"\n";
      }
    }

    /**
     * Map from 2D to 1D return index of 1D vector
     * **/
  inline real operator()(const u_int &rowIndex,const u_int &colIndex)
    {
      assert(rowIndex < _rows && colIndex < _columns);			///assert
      return this->_vecObject[rowIndex * (this->_columns) + colIndex];
    }
    
  inline real operator()(const u_int &rowIndex,const u_int &colIndex)const
    {
       assert(rowIndex < _rows && colIndex < _columns);			//assert
      return this->_vecObject[rowIndex * (this->_columns) + colIndex];
    }
   
//-----------------------------  Perform Arithmetic Operation ----------------------------------------
 /**
  *  Matrix multiplication 
  * @param *this assume *this is pointer that point to matrix A.
  * @param B matrix B.
  * @return Matrix object where result of matrix A and matrix B stores.
  * **/
  Matrix operator*(const Matrix& B);
  
  /**
   * Simple matrix - matrix multiplication and padding
   * @param A  constant of object of matrix class with reference, const becuase it does't change ,we just read data
   * @param B  constant of object of matrix class with reference, const becuase it does't change ,we just read data
   * @param C  object of matrix class with reference as we don't want to copy whole matrix
   * **/
  
  void NaiveImplimentation(const Matrix& A, const Matrix& B, Matrix& C);
  
  /**
   * Basic Matrix Matrix Implementation without any optimization
   * @param A  constant of object of matrix class with reference, const becuase it does't change ,we just read data
   * @param B  constant of object of matrix class with reference, const becuase it does't change ,we just read data
   * @param C  object of matrix class with reference as we don't want to copy whole matrix
   * **/
  
  void UnrollingImplimentation(const Matrix& A,const Matrix& B,Matrix& C);
  
  /**
   * Unrolling j loop
   * @param A  constant of object of matrix class with reference, const becuase it does't change ,we just read data
   * @param B  constant of object of matrix class with reference, const becuase it does't change ,we just read data
   * @param C  object of matrix class with reference as we don't want to copy whole matrix
   * **/
  void BlockImplimentation(const Matrix& A, const Matrix& B, Matrix& C);
  void BlockingAndUnrollingImplimentation(const Matrix& A, const Matrix& B, Matrix& C);
  
  
// --------------------------------- Read and write data from file -----------------------------------
  /**
   * Read data from file named fileName
   * @param fileName name of file from which we read our data
   * @return success read data from file than true else false. 
   * **/
  bool readFile( const std::string& );

    /**
   * Read data from file named fileName
   * @param fileName name of file from which we read our data
   * @return success read data from file than true else false. 
   * **/
  bool writeFile( const std::string& );
    
//--------------------------------------- private function ---------------------------------    
  

  std::vector< real > _vecObject;  
  u_int _rows;
  u_int _columns;
    
};  //class Matrix



  
  
  
  
  
 
