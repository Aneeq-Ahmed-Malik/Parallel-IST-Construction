#include <mpi.h>
#include <omp.h>
#include <cstdint>
#include <iostream>
#include <cstring>
#include <set>
#include <queue>
#include <algorithm>

using namespace std;

const int n = 8; // Number of vertices in the graph


struct Vertex {
    uint8_t data[2 * n];          // [ values | inv ]
    int lehmar;
    uint8_t rightmost_correct;

    Vertex() : lehmar(0), rightmost_correct(0) {
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

    void setLehmarCode(int l){
        this->lehmar = l;
    }

    int getLehmarCode() const {
        return this->lehmar;
    }

    friend ostream& operator<<(ostream& os, const Vertex& v) {
        for (int i = 0; i < n; ++i) {
            os << (int)v[i];
        }
        return os;
    }
};


struct ParentTable
{
    uint8_t data[n * (n - 1)]; // Stores up to (n - 1) parents of size n each

    ParentTable()
    {
        memset(data, 0, sizeof(data));
    }

    // Store a full permutation at index i
    void store(int i, const uint8_t *perm)
    {
        memcpy(&data[i * n], perm, n);
    }

    // Get a pointer to the i-th stored permutation
    const uint8_t *get(int i) const
    {
        return &data[i * n];
    }

    // Optional non-const accessor
    uint8_t *get(int i)
    {
        return &data[i * n];
    }

    friend ostream &operator<<(ostream &os, const ParentTable &pt)
    {
        for (int j = 0; j < n - 1; ++j)
        {
            for (int k = 0; k < n; ++k)
                cout << (int)pt.get(j)[k] << " ";
            cout << " -> t = " << j + 1 << endl;
        }
        return os;
    }
};


int permutation_rank(const uint8_t* perm_suffix, const uint8_t* value_to_index, int suffix_len) {
    size_t rank = 0;
    bool* used = new bool[suffix_len]();  // Initialize to false

    for (int i = 0; i < suffix_len; ++i) {
        uint8_t idx = value_to_index[perm_suffix[i]];
        int count = 0;
        for (int j = 0; j < idx; ++j) {
            if (!used[j]) count++;
        }
        rank = rank * (suffix_len - i) + count;
        used[idx] = true;
    }

    delete[] used;
    return rank;
}
    

int compute_rank(uint8_t* perm) {
    uint8_t* sorted_remaining = new uint8_t[n - 2];
    memcpy(sorted_remaining, &perm[2], (n - 2) * sizeof(uint8_t));
    sort(sorted_remaining, sorted_remaining + (n - 2));

    uint8_t* value_to_index = new uint8_t[n + 1](); // zero-initialized
    for (int i = 0; i < n - 2; ++i)
        value_to_index[sorted_remaining[i]] = i;

    int lehmar = permutation_rank(&perm[2], value_to_index, n - 2);

    delete[] sorted_remaining;
    delete[] value_to_index;
    return lehmar;
}

size_t factorial(int k);

Vertex Swap(Vertex &v, uint8_t x);
Vertex FindPosition(Vertex &v, Vertex &I_n, uint8_t t);
Vertex Parent1(Vertex &v, Vertex &I_n, uint8_t t);

void generateBubbleSortNetworkOptimized()
{

    double start_time = MPI_Wtime();

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
    int vertex_per_pair = factorial(n - 2);

    Vertex* local_vertices = new Vertex[rank_size]; // Allocate enough (overestimate, safe)

    Vertex I_n;
    for (int i = 0; i < n; ++i)
        I_n[i] = i + 1;

    omp_set_num_threads(2);

    #pragma omp parallel
    for (int k = start_idx; k < start_idx + num_pairs; ++k)
    {
        uint8_t F = (k / (n - 1)) + 1;
        uint8_t S = (k % (n - 1)) + 1;
        if (S >= F)
            S++;

        Vertex initial;
        initial[0] = F;
        initial[1] = S;


        uint8_t val = 1;
        for (int i = 2; i < n; ++i) {
            while (val == F || val == S) val++;
            initial[i] = val++;
        }

        // sorted_remaining using DMA
        uint8_t* sorted_remaining = new uint8_t[n - 2];
        memcpy(sorted_remaining, &initial[2], (n - 2) * sizeof(uint8_t));
        sort(sorted_remaining, sorted_remaining + (n - 2));

        // value_to_index using DMA
        uint8_t* value_to_index = new uint8_t[n + 1](); // zero-initialized
        for (int i = 0; i < n - 2; ++i) {
            value_to_index[sorted_remaining[i]] = i;
        }


        initial.compute_inverse();  // Precompute the inverse
        initial.setLehmarCode(permutation_rank(&initial[2], value_to_index, n - 2));

        int perm_size = n - 2;
        vector<uint8_t> perm(perm_size);
        memcpy(perm.data(), &initial[2], perm_size);

        vector<int> c(perm_size, 0);
        int base_index = (k - start_idx) * (rank_size / num_pairs);
        int count = 0;

        while (true) {
            Vertex v = initial;
            for (int i = 0; i < perm_size; ++i)
                v[i + 2] = perm[i];
            v.compute_inverse();
            v.setLehmarCode(permutation_rank(&v[2], value_to_index, perm_size));
            local_vertices[base_index + count++] = v;

            int i = 0;
            while (i < perm_size) {
                if (c[i] < i) {
                    if (i % 2 == 0)
                        swap(perm[0], perm[i]);
                    else
                        swap(perm[c[i]], perm[i]);
                    ++c[i];
                    i = 0;
                    break;
                } else {
                    c[i] = 0;
                    ++i;
                }
            }
            if (i == perm_size) break;
        }

        delete[] sorted_remaining;
        delete[] value_to_index;
    }

    // Synchronize all ranks
    MPI_Barrier(MPI_COMM_WORLD);


    ParentTable** parent_map = new ParentTable*[num_pairs];

    for (int i = 0; i < num_pairs; ++i)
        parent_map[i] = new ParentTable[vertex_per_pair];

    #pragma omp parallel for collapse(2)
    for(int i = 0; i < num_pairs; ++i) {
        for (int j = 0; j < vertex_per_pair; ++j) {
            int idx = i * vertex_per_pair + j;
            int lahmar_code = local_vertices[idx].getLehmarCode();
            for (int t = 1; t < n; ++t)
                parent_map[i][lahmar_code].store(t - 1, Parent1(local_vertices[idx], I_n, t).data);
        }
    }
   
    MPI_Barrier(MPI_COMM_WORLD);

    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;
    if (rank == 0)
        std::cout << "Elapsed time: " << elapsed_time << " seconds" << std::endl;

    delete[] local_vertices;

    if (rank == 0) {
        while (true) {
            cout << "Enter a permutation of size " << n << " (space-separated), or -1 to quit:\n";
            vector<uint8_t> query_perm(n);
            for (int i = 0; i < n; ++i) {
                int x;
                cin >> x;
                query_perm[i] = static_cast<uint8_t>(x);
            }
            
            // first compute the pair_idx
            uint8_t F = query_perm[0];
            uint8_t S = query_perm[1];
            if (S > F) S--;
            uint8_t pair_idx = (F - 1) * (n - 1) + (S - 1);
            int lahmar_code = compute_rank(query_perm.data());

            if (pair_idx < start_idx || pair_idx >= start_idx + num_pairs) {
                MPI_Bcast(&lahmar_code, 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Bcast(&pair_idx, 1, MPI_UINT8_T, 0, MPI_COMM_WORLD);

                // wait for response
                MPI_Recv(nullptr, 0, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                continue;
            }

            ParentTable*& pt = parent_map[pair_idx - start_idx];
    
            cout << "Parents of this permutation (from " << processor_name <<" rank : "<<rank << "):\n";
            cout << pt[lahmar_code] << endl;
            
        }
    }
    else {
        while (true) {
            int lahmar_code;
            uint8_t pair_idx;
            MPI_Bcast(&lahmar_code, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&pair_idx, 1, MPI_UINT8_T, 0, MPI_COMM_WORLD);
            if (pair_idx < start_idx || pair_idx >= start_idx + num_pairs)
                continue;

            ParentTable*& pt = parent_map[pair_idx - start_idx];
    
            cout << "Parents of this permutation (from " << processor_name <<" rank : "<<rank << "):\n";
            cout << pt[lahmar_code] << endl;

            MPI_Send(nullptr, 0, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
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