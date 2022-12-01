#pragma once

#include <iostream>
#include <cassert>
#include <array>
#include <vector>
#include <cmath>
#include <limits>
#include "mat.hpp"

namespace fact {


    // Parameters:
    // Input Matrix
    // A: m x n matrix, rank(A) = dim(A)
    // Output Matrices:
    // Q: m x n orthogonal matrix
    // R: n x n upper triangular invertible matrix, diagonal strictly positive
    // Returns if successful
    template<typename T>
    bool QR(Matrix<T> A, Matrix<T> Q, Matrix<T> R) {
        if (Q.rows != A.rows || Q.cols != A.cols) return false;
        if (R.rows != A.cols || R.cols != A.cols) return false;

        Q = A;
        std::cout << Q;
        // set Q
        Q.orthogonalize();
        Matrix<double> Q_t = transpose(Q);

        // set R
        R = Q_t * A;

        R.round();
        Q.round();

        return true;
    }


};