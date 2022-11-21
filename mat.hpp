#pragma once

#include <iostream>
#include <cassert>
#include <array>
#include <cmath>
#include "basis.hpp"

template<typename T, size_t rows, size_t cols>
class Matrix {
private:
    constexpr static size_t _size = rows * cols;
	size_t _pivots = 0;
	T _contents[_size]{0};
public:
	Matrix() = default;

	//Diagonal Matrix Constructor
	explicit Matrix(const std::array<T, std::min(rows, cols)>& arr) { 
		for (size_t i = 0; i < std::min(rows, cols); ++i)
			_contents[(i * cols) + i] = arr[i];
	}

	explicit Matrix(const std::array<T, _size>& data) {
		for (size_t i = 0; i < _size; i++) {
			_contents[i] = data[i];
		}
	}

	explicit Matrix(const Basis<T, rows>& basis) {
		size_t size = basis.size();
		for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < size; ++j) {
                _contents[i*cols+j] = basis[j][i];
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

	T& operator()(size_t row, size_t col) {
		assert(row < rows && col < cols);
		return _contents[row*cols + col];
	}

	const T& operator()(size_t row, size_t col) const{
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
				for (size_t k = 0; k < cols; k++){
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
			output << mat(i,j) << "\t";
		}
		output << std::endl;
	}
	return output;
	}

	Matrix row_interchange(size_t a, size_t b) {
		assert(a < rows && b < rows);
		if (a == b) return *this;
		T temp;
		for (size_t i = 0; i < cols; i++) {
			temp = _contents[a*cols + i];
			_contents[a*cols + i] = _contents[b*cols + i];
			_contents[b*cols + i] = temp;
		}
		return *this;
	}

	Matrix row_scaling(size_t row, T scalar) {
		assert(row < rows);
		for (size_t i = 0; i < cols; ++i) {
			_contents[row*cols + i] *= scalar;
		}
		return *this;
	}

	Matrix row_replacement(size_t a, size_t b, T scalar) {
		assert(a < rows && b < rows && a != b);
		for (size_t i = 0; i < cols; ++i) {
			_contents[a*cols+i] -= _contents[b*cols+i] * scalar;
		}
		return *this;
	}


	Matrix RREF() {
		// Using Guassian Elimination so not 100% numerically stable
		// Partial pivoting used though less likely to be unstable
		// If something breaks in 6 months then it might be here

		Matrix<T, rows, cols> res = *this;
		int p_row = 0;
		int p_col = 0;

		// Forward phase
		while (p_row < rows && p_col < cols) {			
			// find k-th pivot
			// argmax
			size_t i_max = 0;
			T highest = 0;
			for (size_t i = p_row; i < rows; ++i) {
				if (std::abs(res(i, p_col)) >= highest) {
					highest = std::abs(res(i, p_col));
					i_max = i;
				}
			}

			if (res(i_max, p_col) == 0) {
				++p_col;
			} else {
				res.row_interchange(p_row, i_max);

				for (size_t i = p_row + 1; i < rows; ++i) {
					T f = res(i, p_col) / res(p_row, p_col);
					res.row_replacement(i, p_row, f);
				}
				++p_row;
				++p_col;
			}
		}

		// Backward phase
		p_row = rows - 1;
		p_col = cols - 1;

		while (p_row >= 0 && p_col >= 0) {
			// find pivot col
			for (size_t i = 0; i < cols; ++i) {
				if (res(p_row, i)) {
					p_col = i;
					++res._pivots;
					break;
				}
			}

			// scale pivot row so leading num is 1
			T f = (res(p_row, p_col)) ? 1 / res(p_row, p_col) : 0;

			// create zeros above pivot col
			res.row_scaling(p_row, f);
			for (int i = p_row-1; i >= 0; --i) {
				res.row_replacement(i, p_row, res(i, p_col));
			}
			p_row--;
			p_col--;
		}
		
		return res;
	}

	size_t rank() {
		if (_pivots == 0) {
			return RREF()._pivots;
		}
		return _pivots;
	}
	
	size_t nullity() {
		return cols - rank();
	}
};
