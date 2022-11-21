#pragma once

#include <iostream>
#include <cassert>
#include <vector>
#include <array>
#include "vec.hpp"

template<typename T, size_t dim>
class Basis {
private:
    std::vector<Vec<T, dim>> _elems;
    size_t _size = 0;

public:
    Basis() = default;

    Basis(std::vector<Vec<T, dim>> &elems) {
        _size = elems.size();
        _elems = elems;
    }

     Vec<T, dim>& operator[](size_t ind) {
        assert(ind < _size);
        return _elems[ind];
    } 

    const Vec<T, dim>& operator[](size_t ind) const {
        assert(ind < _size);
        return _elems[ind];
    }

    size_t size() const {
        return _size;
    }

    friend std::ostream &operator<<(std::ostream &output, Basis<T, dim> b) {
        output << "{";
        for (size_t i = 0; i < b.size(); ++i) {
            output << b[i];
            if (i != b.size()-1)
                output << "\n";
        }
        output << "}\n";
        return output;
    }
};