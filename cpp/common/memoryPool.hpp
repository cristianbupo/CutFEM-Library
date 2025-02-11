#ifndef CUTFEM_COMMON_MEMORYPOOL_HPP
#define CUTFEM_COMMON_MEMORYPOOL_HPP

#include <iostream>
#include <vector>
#include <span>
#include <list>

// class MemoryPool;
// class MemoryChunk {
//   public:
//     // Move constructor
//     // MemoryChunk(MemoryChunk &&other) noexcept : memory_pool(other.memory_pool), data(other.data),
//     index(other.index)
//     // {
//     //     other.memory_pool = nullptr;
//     // }

//     MemoryChunk(std::span<double> chunk_data) : data(chunk_data) {}

//     // MemoryChunk(MemoryPool *pool, std::span<double> chunk_data, std::size_t idx)
//     //     : memory_pool(pool), data(chunk_data), index(idx) {}

//     ~MemoryChunk() = default;

//     std::span<double> data;
// };

/// @brief Allocate memory that can be reused by different threads
class MemoryPool {

  public:
    MemoryPool() = default;

    MemoryPool(std::size_t n_chunk, std::size_t l_chunk) : chunk_size(l_chunk), number_chunks(n_chunk) {
        pool.resize(n_chunk * chunk_size);
    }

    std::size_t size() const { return pool.size(); }

    std::span<double> get_data(std::size_t i) { return std::span<double>(pool.data() + i * chunk_size, chunk_size); }

  private:
    std::size_t chunk_size{0};

    std::size_t number_chunks{0};

    std::vector<double> pool{};
};

// inline MemoryChunk::~MemoryChunk() {
//     //     if (memory_pool) {
//     // #pragma omp critical
//     //         memory_pool->available_memory_chunk.push_back(index);
//     //     }
// }

#endif