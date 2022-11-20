#include <iostream>
#include <cassert>
#include <array>
#include "mat.hpp"
#include "vec.hpp"

template<typename T, size_t rows>
class Basis {
private:
    std::vector<Matrix<T, rows, 1>> _elems;
    size_t size = 0;

public:
    Basis() = default;

    Basis(std::vector<Vec<T, rows>> &elems) {
        size = elems.size();

    }

    Matrix join() {
        Matrix<T, rows, 

    }

};