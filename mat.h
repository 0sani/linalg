#include <stdio>


template<typename T>
class Matrix {
private:
	T* contents;
	const int rows;
	const int cols;
public:
	// basic matrix constructor	
	Matrix(int rows, int cols);

	// vector constructor
	explicit Matrix(int rows);

	Matrix(const Matrix&); // copy constructor
	Matrix(Matrix&&); // move constructor

	Matrix& operator=(const Matrix&); // copy assignment
	Matrix& operator=(X&&); // move assignment

	~Matrix();	
}
