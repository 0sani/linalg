#include <iostream>
#include <iostream>
#include <cassert>
#include <array>
#include <cmath>


template<typename T, size_t size>
class Vec {
private:
    T _contents[size];

public:
    Vec() = default;

    Vec(const std::array<T, size>& arr) {
        for (size_t i = 0; i < size; ++i)
            _contents[i] = arr[i];
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
            res[i] = _contents[i] + other[i]
        }
        return res;
    }

    Vec operator-(Vec<T, size> other) {
        Vec<T, size> res;
        for (size_t i = 0; i < size; ++i) {
            res[i] = _contents[i] - other[i]
        }
        return res;
    }

    // Dot product
    T& operator*(Vec<T, size> other) {
        T total = 0;
        for (size_t i = 0; i < size; ++i) {
            total +=  _contents[i] * other[i]
        }

    }

};