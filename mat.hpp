#include <iostream>
#include <cassert>
#include <cstring>
#include <array>

template<typename T, size_t rows, size_t cols>
class Matrix {
private:
    constexpr static size_t _size = rows * cols;
	T _contents[_size]{};
public:
    Matrix() = default;

	Matrix(const std::array<T, _size>& arr) { memcpy(_contents, arr.data(), _size * sizeof(T)); }
	
	// Diagonal Matrix Constructor
	template <size_t dim, std::enable_if_t<dim == std::min(rows, cols), bool> = true>  
	Matrix(const std::array<T, dim>& arr) {
		for (size_t i = 0; i < dim; ++i) {
			if constexpr (rows == cols)
			{
				_contents[(i * dim) + i] = arr[i];
			}
			else
			{
				_contents[i * (cols + 1)] = arr[i];
			}
		}
	}

	T& operator[](size_t index) {
		assert(index < _size);
		return _contents[index];
	}

	const T& operator[](size_t index) const {
		assert(index < _size);
		return _contents[index];
	}

	T& operator()(size_t y, size_t x) {
		assert(x < cols && y < rows);
		return _contents[y*cols + x];
	}

	const T& operator()(size_t y, size_t x) const {
		assert(x < cols && y < rows);
		return _contents[y*cols + x];
	}

	Matrix& operator*=(T scalar) {
		for (size_t i = 0; i < _size; ++i)
			_contents[i] *= scalar;
		return *this;
	}


	Matrix& operator+=(Matrix<T, rows, cols> other) {
		for (size_t i = 0; i < _size; ++i)
			_contents[i] += other[i];
		return *this;
	}

	Matrix& operator-=(Matrix<T, rows, cols> other) {
		for (size_t i = 0; i < _size; ++i)
			_contents[i] -= other[i];
		return *this;
	}

	Matrix& operator/=(T scalar) {
		for (size_t i = 0; i < _size; ++i)
			_contents[i] /= scalar;
		return *this;
	}

	Matrix operator+(Matrix<T, rows, cols> other) {
		Matrix<T, rows, cols> res;
		for (size_t i = 0; i < _size; ++i) {
			res[i] = _contents[i] + other[i];
		}
		return res;
	}

	Matrix operator-(Matrix<T, rows, cols> other) {
		Matrix<T, rows, cols> res;
		for (size_t i = 0; i < _size; ++i) {
			res[i] = _contents[i] - other[i];
		}
		return res;
	}

	template<size_t C>
	Matrix operator*(Matrix<T, cols, C> other) {
		Matrix<T, rows, C> res;
		for (size_t i = 0; i < rows; ++i) {
			for (size_t j = 0; j < C; ++j) {
				T sum = 0;
				for (size_t k = 0; k < rows; ++k){

					sum += _contents[i*rows+k] * other(k, j);
				}
				res(i,j) = sum;
			}
			
		}
		return res;
	}

	friend std::ostream &operator<<(std::ostream &os, const Matrix<T, rows, cols>& mat) {
		for (size_t i = 0; i < rows; ++i) {
			for (size_t j = 0; j < cols; ++j) {
				os << mat(i, j) << ' ';
			}
			os << std::endl;
		}
		return os;
	}
};
