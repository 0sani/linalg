#pragma once

#include <iostream>
#include <cassert>
#include <array>
#include <cmath>

template<typename T>
class Matrix {
private:
	const size_t _rows;
	const size_t _cols;
    const size_t _size = _rows * _cols;
	size_t _pivots = 0;
	T* _contents;

	// not sure where to put it but for floating point errors
	double epsilon = 1e-10;
public:

	Matrix(size_t rows, size_t cols)
		: _rows {rows}, _cols{cols} {
		_contents = new T[_size];
	}

	Matrix(size_t rows, size_t cols, const std::vector<T>& arr)
		: _rows {rows}, _cols{cols} {
		
		assert(arr.size() == _size || arr.size() == std::min(_rows, _cols));
		_contents = new T[_size];
		
		if (arr.size() == _size) {
			for (size_t i = 0; i < _size; ++i)
				_contents[i] = arr[i];
		} else {
			for (size_t i = 0; i < arr.size(); ++i)
				_contents[i*(_cols + 1)] = arr[i];
		}
	}

	Matrix operator=(const Matrix& arr) {
		assert(arr._rows == _rows && arr._cols == _cols);
		for (size_t i = 0; i < _size; i++) {
			_contents[i] = arr[i];
		}
		return *this;
	}

	bool operator==(const Matrix& arr) {
		assert(arr._rows == _rows && arr._cols == _cols);
		for (size_t i = 0; i < _size; i++) {
			if (std::abs(_contents[i] - arr[i]) > epsilon) {
				return false;
			}
		}
		return true;
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
		assert(row < _rows && col < _cols);
		return _contents[row*_cols + col];
	}

	const T& operator()(size_t row, size_t col) const{
		assert(row < _rows && col < _cols);
		return _contents[row*_cols + col];
	}

	Matrix& operator*=(T scalar) {
		for (size_t i = 0; i < _size; ++i)
			_contents[i] *= scalar;
		return *this;
	}

	Matrix& operator/=(T scalar) {
		for (size_t i = 0; i < _size; ++i)
			_contents[i] /= scalar;
		return *this;
	}

	Matrix& operator+=(Matrix<T> other) {
		assert(other._rows == _rows && other._cols == _cols);
		for (size_t i = 0; i < _size; ++i)
			_contents[i] += other[i];
		return *this;
	}

	Matrix& operator-=(Matrix<T> other) {
		assert(other._rows == _rows && other._cols == _cols);
		for (size_t i = 0; i < _size; ++i)
			_contents[i] -= other[i];
		return *this;
	}

	Matrix operator+(Matrix<T> other) {
		Matrix res(_rows, _cols);
		for (size_t i = 0; i < _size; ++i) {
			res[i] = _contents[i] + other[i];
		}
		return res;
	}

	Matrix operator-(Matrix<T> other) {
		Matrix res(_rows, _cols);
		for (size_t i = 0; i < _size; ++i) {
			res[i] = _contents[i] - other[i];
		}
		return res;
	}

	Matrix operator*(Matrix<T> other) {
		assert(_cols == other._rows);
		Matrix res(_rows, other._cols);
		for (size_t i = 0; i < _rows; ++i) {
			for (size_t j = 0; j < other._cols; j++) {
				T sum = 0;
				for (size_t k = 0; k < _cols; k++){
					sum += _contents[i*_cols+k] * other(k, j);
				}
				res(i,j) = sum;
			}
		}
		return res;
	}

	friend std::ostream &operator<<(std::ostream &output, Matrix mat) {
	for (size_t i = 0; i < mat._rows; ++i) {
		for (size_t j = 0; j < mat._cols; j++) {
			output << mat(i,j) << "\t";
		}
		output << std::endl;
	}
	return output;
	}

	Matrix row_interchange(size_t a, size_t b) {
		assert(a < _rows && b < _rows);
		if (a == b) return *this;
		T temp;
		for (size_t i = 0; i <_cols; ++i) {
			temp = _contents[a*_cols + i];
			_contents[a*_cols + i] = _contents[b*_cols + i];
			_contents[b*_cols + i] = temp;
		}
		return *this;
	}

	Matrix row_scaling(size_t row, T scalar) {
		assert(row <_rows);
		for (size_t i = 0; i <_cols; ++i) {
			_contents[row*_cols + i] *= scalar;
		}
		return *this;
	}

	Matrix row_replacement(size_t a, size_t b, T scalar) {
		assert(a <_rows && b <_rows && a != b);
		for (size_t i = 0; i <_cols; ++i) {
			_contents[a*_cols+i] -= _contents[b*_cols+i] * scalar;
		}
		return *this;
	}


	Matrix RREF() {
		// Using Guassian Elimination so not 100% numerically stable
		// Partial pivoting used though less likely to be unstable
		// If something breaks in 6 months then it might be here

		// Update two days later:
		// It is numerically unstable :D

		Matrix res(_rows, _cols);
		res = *this;
		int p_row = 0;
		int p_col = 0;

		// Forward phase
		while (p_row <_rows && p_col <_cols) {			
			// find k-th pivot
			// argmax
			size_t i_max = 0;
			T highest = 0;
			for (size_t i = p_row; i <_rows; ++i) {
				if (std::abs(res(i, p_col)) >= highest) {
					highest = std::abs(res(i, p_col));
					i_max = i;
				}
			}

			if (res(i_max, p_col) == 0) {
				++p_col;
			} else {
				res.row_interchange(p_row, i_max);

				for (size_t i = p_row + 1; i <_rows; ++i) {
					T f = res(i, p_col) / res(p_row, p_col);
					res.row_replacement(i, p_row, f);
				}
				++p_row;
				++p_col;
			}
		}

		// zero any potential floating point errors
		for (size_t i = 0; i < res._size; ++i) {
			if (res[i] < epsilon) res[i] = 0;
		}

		// Backward phase
		p_row =_rows - 1;
		p_col =_cols - 1;

		while (p_row >= 0 && p_col >= 0) {
			// find pivot col
			for (size_t i = 0; i <_cols; ++i) {
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
		return _cols - rank();
	}


};
