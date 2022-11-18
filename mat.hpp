#include <iostream>
#include <cassert>

template<typename T, size_t rows, size_t cols>
class Matrix {
private:
    constexpr static size_t _size = rows * cols;
	T _contents[cols][rows];
public:
	T& operator[](size_t index) {
		assert(index < _size);
		return _contents[index];
	}

	T& operator()(size_t x, size_t y) {
		assert(x < rows && y < cols);
		return _contents[x*rows + y];
	}


	Matrix& operator*=(T scalar) {
		for (size_t i = 0; i < _size; i++)
			_contents[i] *= scalar;
		return *this;
	}


	Matrix& operator+=(Matrix<T, rows, cols> other) {
		for (size_t i = 0; i < _size; i++)
			_contents[i] += other[i];
		return *this;
	}

	Matrix& operator-=(Matrix<T, rows, cols> other) {
		for (size_t i = 0; i < _size; i++)
			_contents[i] -= other[i];
		return *this;
	}

	Matrix& operator/=(T scalar) {
		for (size_t i = 0; i < _size; i++)
			_contents[i] /= scalar;
		return *this;
	}

	Matrix operator+(Matrix<T, rows, cols> other) {
		Matrix<T, rows, cols> res;
		for (size_t i = 0; i < _size; i++) {
			res[i] = _contents[i] + other[i];
		}
		return res;
	}

	Matrix operator-(Matrix<T, rows, cols> other) {
		Matrix<T, rows, cols> res;
		for (size_t i = 0; i < _size; i++) {
			res[i] = _contents[i] - other[i];
		}
		return res;
	}

	template<size_t C>
	Matrix operator*(Matrix<T, cols, C> other) {
		Matrix<T, rows, C> res;
		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < C; j++) {
				T sum = 0;
				for (size_t k = 0; k < rows; k++){

					sum += _contents[i*rows+k] * other(k, j);
				}
				res(i,j) = sum;
			}
			
		}
		return res;
	}

};
