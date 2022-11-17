#include <iostream>


template<typename T, size_t rows, size_t cols>
class Matrix {
private:
	T* _contents;
	constexpr static size_t _rows = rows;
	constexpr static size_t _cols = cols;
	constexpr static size_t _size = rows*cols;
public:
	Matrix() {
		_contents = new T[_size];
	}

	T& operator[](size_t index) {
		assert(index < _size);
		return _contents[index];
	}

	T& operator()(size_t row, size_t col) {
		assert(row < _rows && col < _cols);
		return _contents[row*_rows + col];
	}

	int operator*=(int scalar) {
		for (size_t i = 0; i < _size; i++)
			_contents[i] *= scalar;
		return scalar;
	}

};
