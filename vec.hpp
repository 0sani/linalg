#pragma once

#include <iostream>
#include <cassert>
#include <array>
#include <cmath>


template<typename T, size_t size>
class Vec {
private:
    T _contents[size];
    T magnitude;

public:
    Vec() = default;

    Vec(const std::array<T, size>& arr) {
        T l = 0;
        for (size_t i = 0; i < size; ++i) {
            _contents[i] = arr[i];
            l += std::pow(arr[i],2);
        }
        magnitude = std::sqrt(l);
    }

    Vec operator=(const Vec<T, size>& v) {
		for (size_t i = 0; i < size; ++i) {
			_contents[i] = v[i];
		}
		return *this;
	}

    T& operator[](size_t ind) {
        assert(ind < size);
        return _contents[ind];
    } 

    const T& operator[](size_t ind) const {
        assert(ind < size);
        return _contents[ind];
    }

    Vec& operator*=(T scalar) {
        for (size_t i = 0; i < size; i++)
			_contents[i] *= scalar;
		return *this;
    }

    Vec& operator/=(T scalar) {
		for (size_t i = 0; i < size; i++)
			_contents[i] /= scalar;
		return *this;
	}

    Vec& operator+=(Vec<T, size> other) {
		for (size_t i = 0; i < size; i++)
			_contents[i] += other[i];
		return *this;
	}

    Vec& operator-=(Vec<T, size> other) {
		for (size_t i = 0; i < size; i++)
			_contents[i] += other[i];
		return *this;
	}

    Vec operator+(Vec<T, size> other) {
        Vec<T, size> res;
        for (size_t i = 0; i < size; ++i) {
            res[i] = _contents[i] + other[i];
        }
        return res;
    }

    Vec operator-(Vec<T, size> other) {
        Vec<T, size> res;
        for (size_t i = 0; i < size; ++i) {
            res[i] = _contents[i] - other[i];
        }
        return res;
    }

    // Dot product
    T operator*(Vec<T, size> other) {
        T total = 0;
        for (size_t i = 0; i < size; ++i) {
            total +=  _contents[i] * other[i];
        }
        return total;
    }

    // Scalar Multiplication
    Vec operator*(T scalar) {
        Vec<T, size> res;
        for (size_t i = 0; i < size; ++i) {
            res[i] = _contents[i] * scalar;
        }
        return res;
    }

    Vec operator/(T scalar) {
        Vec<T, size> res;
        for (size_t i = 0; i < size; ++i) {
            res[i] = _contents[i] / scalar;
        }
        return res;
    }

    friend std::ostream &operator<<(std::ostream &output, Vec<T, size> v) {
        output << "[";
        for (size_t i = 0; i < size; ++i) {
            output << v[i];
            if (i != size-1)
                output << " ";
            else
                output << "]";
        }
        return output;
    }

    // Projection onto the other vector 
    Vec proj_onto(Vec<T, size> other) {
        Vec<T, size> res = *this;
        return other * ((res*other) / other.norm());
    }

    T& norm() {
        return magnitude;
    }

};