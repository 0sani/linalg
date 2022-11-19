#include <iostream>
#include <cassert>
#include <array>

template<typename T, size_t rows, size_t cols>
class Matrix {
private:
    constexpr static size_t _size = rows * cols;
	T _contents[_size]{0};
public:
	Matrix() = default;

	//Diagonal Matrix Constructor
	Matrix(const std::array<T, std::min(rows, cols)>& arr) { 
		for (size_t i = 0; i < std::min(rows, cols); ++i)
			_contents[(i * cols) + i] = arr[i];
	}

	explicit Matrix(const std::array<T, _size>& data) {
		for (size_t i = 0; i < _size; i++) {
			_contents[i] = data[i];
		}
	}

	T& operator[](size_t index) {
		assert(index < _size);
		return _contents[index];
	}

	T& operator()(size_t row, size_t col) {
		assert(row < rows && col < cols);
		return _contents[row*cols + col];
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
	Matrix<T, rows, C> operator*(Matrix<T, cols, C> other) {
		Matrix<T, rows, C> res;
		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < C; j++) {
				T sum = 0;
				for (size_t k = 0; k < rows; k++){

					sum += _contents[i*cols+k] * other(k, j);
				}
				res(i,j) = sum;
			}
		}
		return res;
	}

	friend std::ostream &operator<<(std::ostream &output, Matrix<T, rows, cols> mat) {
	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < cols; j++) {
			output << mat(i,j) << " ";
		}
		output << std::endl;
	}
	return output;
	}
};
