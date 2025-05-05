#include <mpi.h>
#include <omp.h>
#include <cstdint>
#include <iostream>
#include <cstring>

using namespace std;

const int n = 8; // Number of vertices in the graph


struct Vertex {
    uint8_t data[2 * n];        
    uint8_t rightmost_correct;

    Vertex() : rightmost_correct(0) {
        memset(data, 0, sizeof(data));
    }

    // Access values
    uint8_t& operator[](int idx) {
        return data[idx];
    }

    const uint8_t& operator[](int idx) const {
        return data[idx];
    }

    // Access inverse via ()
    uint8_t operator()(uint8_t value) const {
        return data[n + value - 1];  // inv[value] = position
    }

    // Access rightmost incorrect index
    uint8_t r() const {
        return rightmost_correct;
    }

    // Precompute inv[] and r()
    void compute_inverse() {
        uint8_t r_pos = 0;
        for (uint8_t i = 0; i < n; ++i) {
            data[n + data[i] - 1] = i;
            if (data[i] != i + 1) {
                r_pos = std::max(r_pos, i);
            }
        }
        rightmost_correct = r_pos;
    }

    // Equality checks
    bool operator==(const Vertex& other) const {
        return memcmp(data, other.data, n * sizeof(uint8_t)) == 0;
    }

    bool operator!=(const Vertex& other) const {
        return !(*this == other);
    }

    friend ostream& operator<<(ostream& os, const Vertex& v) {
        for (int i = 0; i < n; ++i) {
            os << (int)v[i];
        }
        return os;
    }
};


int main() {
  
}