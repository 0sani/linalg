#include <iostream>


template<typename T>
class Matrix {
private:
	T* contents;
	const int rows;
	const int cols;
public:
	
	Matrix(int rows, int cols);


	// vector constructor
	explicit Matrix(int rows);

	
};
