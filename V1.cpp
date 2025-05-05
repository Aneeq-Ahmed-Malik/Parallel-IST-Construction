#include <mpi.h>
#include <omp.h>
#include <cstdint>
#include <iostream>
#include <cstring>
#include <set>
#include <queue>

using namespace std;

const int n = 8; // Number of vertices in the graph

struct Vertex
{
    uint8_t data[2 * n];
    uint8_t rightmost_correct;

    Vertex() : rightmost_correct(0){
        memset(data, 0, sizeof(data));
    }

    // Access values
    uint8_t &operator[](int idx)
    {
        return data[idx];
    }

    const uint8_t &operator[](int idx) const
    {
        return data[idx];
    }

    // Access inverse via ()
    uint8_t operator()(uint8_t value) const
    {
        return data[n + value - 1]; // inv[value] = position
    }

    // Access rightmost incorrect index
    uint8_t r() const
    {
        return rightmost_correct;
    }

    // Precompute inv[] and r()
    void compute_inverse()
    {
        uint8_t r_pos = 0;
        for (uint8_t i = 0; i < n; ++i)
        {
            data[n + data[i] - 1] = i;
            if (data[i] != i + 1)
            {
                r_pos = std::max(r_pos, i);
            }
        }
        rightmost_correct = r_pos;
    }

    // Equality checks
    bool operator==(const Vertex &other) const
    {
        return memcmp(data, other.data, n * sizeof(uint8_t)) == 0;
    }

    bool operator!=(const Vertex &other) const
    {
        return !(*this == other);
    }

    friend ostream &operator<<(ostream &os, const Vertex &v)
    {
        for (int i = 0; i < n; ++i)
        {
            os << (int)v[i];
        }
        return os;
    }
};

size_t factorial(int k);

Vertex Swap(Vertex &v, uint8_t x);
Vertex FindPosition(Vertex &v, Vertex &I_n, uint8_t t);
Vertex Parent1(Vertex &v, Vertex &I_n, uint8_t t);

void generateBubbleSortNetworkOptimized()
{
    int rank, size, namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Get_processor_name(processor_name, &namelen);

    int total_pairs = n * (n - 1);
    int pairs_per_rank = total_pairs / size;
    int remainder = total_pairs % size;
    int start_idx = rank * pairs_per_rank + min(rank, remainder);
    int num_pairs = pairs_per_rank + (rank < remainder ? 1 : 0);

    int rank_size = factorial(n - 2) * num_pairs;
    cout << "Rank " << rank << " on processor " << processor_name << " processing " << num_pairs << " pairs." << endl;

    for (int k = start_idx; k < start_idx + num_pairs; ++k)
    {
        uint8_t F = (k / (n - 1)) + 1;
        uint8_t S = (k % (n - 1)) + 1;
        if (S >= F)
            S++;

        Vertex initial;
        initial[0] = F;
        initial[1] = S;

        cout << "decoded pair (" << (int)F << ", " << (int)S << ") from index " << k << endl;
        cout << "Rank " << rank << " processing pair (" << (int)F << ", " << (int)S << ")" << endl;
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    generateBubbleSortNetworkOptimized();
    MPI_Finalize();
    return 0;
}

size_t factorial(int k)
{
    size_t res = 1;
    for (int i = 2; i <= k; ++i)
        res *= i;
    return res;
}

Vertex Swap(Vertex &v, uint8_t x)
{
    Vertex p = v;
    uint8_t i = v(x); // inverse
    swap(p[i], p[i + 1]);
    return p;
}

Vertex FindPosition(Vertex &v, Vertex &I_n, uint8_t t)
{
    Vertex p;
    if (t == 2 && Swap(v, t) == I_n)
        p = Swap(v, t - 1);
    else if (v[n - 2] == t || v[n - 2] == n - 1)
        p = Swap(v, v.r() + 1);
    else
        p = Swap(v, t);
    return p;
}

Vertex Parent1(Vertex &v, Vertex &I_n, uint8_t t)
{
    Vertex p;
    if (v[n - 1] == n)
    {
        if (t != n - 1)
            p = FindPosition(v, I_n, t);
        else
            p = Swap(v, v[n - 2]);
    }
    else
    {
        if (v[n - 1] == n - 1 && v[n - 2] == n && Swap(v, n) != I_n)
        {
            if (t == 1)
                p = Swap(v, n);
            else
                p = Swap(v, t - 1);
        }
        else
        {
            if (v[n - 1] == t)
                p = Swap(v, n);
            else
                p = Swap(v, t);
        }
    }
    return p;
}
