#pragma once

#include <iostream>
#include <cassert>
#include <array>
#include <vector>
#include <cmath>
#include <limits>

template<typename T>
class Matrix {
private:
	T* _contents;

	// not sure where to put it but for floating point errors
	T epsilon = std::numeric_limits<T>::epsilon();
public:
	const size_t rows;
	const size_t cols;
    const size_t size = rows * cols;
	size_t _pivots = 0;
	std::vector<size_t> pivot_cols;

	Matrix(size_t rows, size_t cols)
		: rows {rows}, cols{cols} {
		_contents = new T[size];
	}

	Matrix(size_t rows, size_t cols, const std::vector<T>& arr)
		: rows {rows}, cols{cols} {
		
		assert(arr.size() == size || arr.size() == std::min(rows, cols));
		_contents = new T[size];
		
		if (arr.size() == size) {
			for (size_t i = 0; i < size; ++i)
				_contents[i] = arr[i];
		} else {
			for (size_t i = 0; i < size; ++i) {
				if (i % cols == i / rows) {
					_contents[i] = arr[i%cols];
				} else {
					_contents[i] = 0;
				}
			}
		}
	}

	T *begin() {
		return &_contents[0];
	}

	T *end() {
		return &_contents[size] + 1;
	}

	Matrix operator=(const Matrix& arr) {
		assert(arr.rows == rows && arr.cols == cols);
		for (size_t i = 0; i < size; i++) {
			_contents[i] = arr[i];
		}
		return *this;
	}

	bool operator==(const Matrix& arr) {
		assert(arr.rows == rows && arr.cols == cols);
		for (size_t i = 0; i < size; i++) {
			// compensate for floating point errors
			if (std::fabs(_contents[i] - arr[i]) > epsilon * std::fabs(_contents[i] + arr[i]) * 2) {
				return false;
			}
		}
		return true;
	}

	bool operator!=(const Matrix& arr) {
		return !(*this==arr);
	}

	T& operator[](size_t index) {
		assert(index < size);
		return _contents[index];
	}

	const T& operator[](size_t index) const {
		assert(index < size);
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
		for (size_t i = 0; i < size; ++i)
			_contents[i] *= scalar;
		return *this;
	}

	Matrix& operator/=(T scalar) {
		for (size_t i = 0; i < size; ++i)
			_contents[i] /= scalar;
		return *this;
	}

	Matrix& operator+=(Matrix<T> other) {
		assert(other.rows == rows && other.cols == cols);
		for (size_t i = 0; i < size; ++i)
			_contents[i] += other[i];
		return *this;
	}

	Matrix& operator-=(Matrix<T> other) {
		assert(other.rows == rows && other.cols == cols);
		for (size_t i = 0; i < size; ++i)
			_contents[i] -= other[i];
		return *this;
	}

	Matrix operator+(Matrix<T> other) {
		Matrix res(rows, cols);
		for (size_t i = 0; i < size; ++i) {
			res[i] = _contents[i] + other[i];
		}
		return res;
	}

	Matrix operator-(Matrix<T> other) {
		Matrix res(rows, cols);
		for (size_t i = 0; i < size; ++i) {
			res[i] = _contents[i] - other[i];
		}
		return res;
	}

	Matrix operator*(Matrix<T> other) {
		assert(cols == other.rows);
		Matrix res(rows, other.cols);
		for (size_t i = 0; i < rows; ++i) {
			for (size_t j = 0; j < other.cols; j++) {
				T sum = 0;
				for (size_t k = 0; k < cols; k++){
					sum += _contents[i*cols+k] * other(k, j);
				}
				res(i,j) = sum;
			}
		}
		return res;
	}

	void round() {
		for (size_t i = 0; i < size; ++i) {
			if (std::fabs(_contents[i]) < 1e-12) {
				_contents[i] = 0;
			}
		}
	}

	friend std::ostream &operator<<(std::ostream &output, Matrix mat) {
	for (size_t i = 0; i < mat.rows; ++i) {
		for (size_t j = 0; j < mat.cols; j++) {
			output << mat(i,j) << "\t";
		}
		output << std::endl;
	}
	return output;
	}

	T determinant() {
		if (is_lower_triangular() || is_upper_triangular()) {
			// can do this since triangular matrices must be square
			T res = 1;
			for (size_t i = 0; i < rows; ++i) {
				res *= _contents[i*(cols+1)];
			} 
			return res;
		}
		// TODO Add proper determinant function
		return -1;
	}

	Matrix row_interchange(size_t a, size_t b) {
		assert(a < rows && b < rows);
		if (a == b) return *this;
		T temp;
		for (size_t i = 0; i <cols; ++i) {
			temp = _contents[a*cols + i];
			_contents[a*cols + i] = _contents[b*cols + i];
			_contents[b*cols + i] = temp;
		}
		return *this;
	}

	Matrix row_scaling(size_t row, T scalar) {
		assert(row <rows);
		for (size_t i = 0; i <cols; ++i) {
			_contents[row*cols + i] *= scalar;
		}
		return *this;
	}

	Matrix row_replacement(size_t a, size_t b, T scalar) {
		assert(a <rows && b <rows && a != b);
		for (size_t i = 0; i <cols; ++i) {
			_contents[a*cols+i] -= _contents[b*cols+i] * scalar;
		}
		return *this;
	}

	Matrix col_replacement(size_t a, size_t b, T scalar) {
		assert(a < cols && b < cols);
		for (size_t i = 0; i < rows; ++i) {
			_contents[i*cols+a] -= _contents[i*cols+b] * scalar;
		}
		return *this;
	}


	size_t rank() {
		if (_pivots == 0) {
			return RREF(*this)._pivots;
		}
		return _pivots;
	}
	
	size_t nullity() {
		return cols - rank();
	}

	T dot_cols(size_t a, size_t b) {
		assert(a < cols && b < cols);
		T sum = 0;
		for (size_t i = 0; i < rows; ++i) {
			sum += _contents[i*cols + a] * _contents[i*cols + b];
		}
		return sum;
	}

	T col_magnitude(size_t col) {
		assert(col < cols);
		return std::sqrt(dot_cols(col, col));
	}

	void normalize_col(size_t col) {
		assert(col < cols);
		T norm = col_magnitude(col);
		for (size_t i = 0; i < rows; ++i) {
			_contents[i*cols + col] = _contents[i*cols + col] / norm;
		}
	}

	// Normalizes the current matrix
	void normalize() {
		for (size_t i = 0; i < cols; ++i) {
			normalize_col(i);
		}
	}

	// Orthogonalizes the current matrix using stable Graham-Schmidt Process
	void orthogonalize() {
		// Ensure matrix is a basis
		assert(rank() == cols);

		for (size_t j = 0; j < cols; ++j) {
			// create u_j
			normalize_col(j);
			for (size_t k = j+1; k < cols; ++k) {
				col_replacement(k,j,dot_cols(k,j));
			}
		}
	}

	bool is_orthogonal() {
		for (size_t i = 0; i < cols-1; ++i) {
			for (size_t j = i+1; j < cols; ++j) {
				if (std::fabs(col_magnitude(i)-1) > 1e-12 ||  std::fabs(dot_cols(i,j)) > 1e-12) {
					return false;
				}
			}
		}
		return true;
	}

	bool is_upper_triangular() {
		if (rows != cols) return false;
		for (size_t i = 0; i < rows; ++i) {
			for (size_t j = 0; j < i; ++j) {
				if (_contents[i*cols + j]) return false;
			}
		}
		return true;
	}

	bool is_lower_triangular() {
		if (rows != cols) return false;
		for (size_t i = 0; i < rows; ++i) {
			for (size_t j = i+1; j < rows; ++j) {
				if (_contents[i*cols + j]) return false;
			}
		}
		return true;
	}

};

template<typename T>
Matrix<T> RREF(const Matrix<T>& mat) {
	// Using Guassian Elimination so not 100% numerically stable
	// Partial pivoting used though less likely to be unstable
	// If something breaks in 6 months then it might be here

	// Update two days later:
	// It is numerically unstable :D

	size_t rows = mat.rows;
	size_t cols = mat.cols;

	Matrix<T> res(rows, cols);
	res = mat;
	int p_row = 0;
	int p_col = 0;

	// Forward phase
	while (p_row < rows && p_col <cols) {			
		// find k-th pivot
		// argmax
		size_t i_max = 0;
		T highest = 0;
		for (size_t i = p_row; i <rows; ++i) {
			if (std::abs(res(i, p_col)) >= highest) {
				highest = std::abs(res(i, p_col));
				i_max = i;
			}
		}

		if (res(i_max, p_col) == 0) {
			++p_col;
		} else {
			res.row_interchange(p_row, i_max);

			for (size_t i = p_row + 1; i <rows; ++i) {
				T f = res(i, p_col) / res(p_row, p_col);
				res.row_replacement(i, p_row, f);
			}
			++p_row;
			++p_col;
		}
	}
	// zero any potential floating point errors
	for (size_t i = 0; i < res.size; ++i) {
		T x = res[i];
		if (std::fabs(x) <= 1e-12) {
			res[i] = 0;
		}
	}

	// Backward phase
	p_row =rows - 1;
	p_col =cols - 1;

	while (p_row >= 0 && p_col >= 0) {
		// find pivot col
		for (size_t i = 0; i <cols; ++i) {
			if (res(p_row, i)) {
				p_col = i;
				++res._pivots;
				res.pivot_cols.push_back(i);
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

	// zero any potential floating point errors
	res.round();
	
	return res;
}

template<typename T>
Matrix<T> transpose(const Matrix<T>& mat) {
	Matrix<T> res(mat.cols, mat.rows);
	for (size_t i = 0; i < mat.rows; ++i)
	{
		for (size_t j = 0; j < mat.cols; ++j)
		{
			res(j,i) = mat(i, j);
		}
		
	}
	return res;
}