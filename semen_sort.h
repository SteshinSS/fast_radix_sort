#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <immintrin.h>
#include <mmintrin.h>
#include <cstring>

constexpr int ELEMENTS_TO_RUN_STD_SORT = 1000;
constexpr int ELEMENTS_TO_RUN_IN_CACHE = 10000;
constexpr uint32_t CACHE_LINE_SIZE = 64;

template <class RandomAccessIterator>
inline bool
is_sorted_or_find_extremes(RandomAccessIterator current, RandomAccessIterator end,
                           RandomAccessIterator& max, RandomAccessIterator& min)
{
    min = max = current;
    //This assumes we have more than 1 element based on prior checks.
    while (!(*(current + 1) < *current)) {
        //If everything is in sorted order, return
        if (++current == end - 1)
            return true;
    }

    //The maximum is the last sorted element
    max = current;
    //Start from the first unsorted element
    while (++current < end) {
        if (*max < *current)
            max = current;
        else if (*current < *min)
            min = current;
    }
    return false;
}

inline unsigned
rough_log_2_size(unsigned input)
{
    unsigned result = 0;
    while ((input >> result) && (result < (8*sizeof(int)))) ++result;
    return result;
}

template<typename ValueType>
inline unsigned GetKey(ValueType value, unsigned min, unsigned log_div) {
    value ^= 1U << (8 * sizeof(ValueType) - 1);
    return (value - min) >> log_div;
}

// read overloaded InCacheSort(RandomAccessIterator begin, RandomAccessIterator end) below first
template<class RandomAccessIterator, uint32_t log_total_buckets>
void InCacheSort(RandomAccessIterator begin, RandomAccessIterator end, unsigned u_min, unsigned u_max) {
    using ValueType = typename std::iterator_traits<RandomAccessIterator>::value_type;

    const uint32_t log_range = rough_log_2_size(u_max - u_min);
    const uint32_t log_div = log_range > log_total_buckets ? log_range - log_total_buckets : 0;

    std::vector<uint32_t> histogram(1U << log_total_buckets, 0);
    for (RandomAccessIterator item = begin; item != end; ++item) {
        histogram[GetKey(*item, u_min, log_div)]++;
    }
    std::vector<uint32_t> offsets(1U << log_total_buckets);
    offsets[0] = histogram[0];
    for (int i = 1; i < histogram.size(); ++i) {
        offsets[i] = offsets[i - 1] + histogram[i];
    }

    int bucket = 0;
    while (offsets[bucket] == 0) {
        ++bucket;
    }

    // filling buckets from right to left
    int position;
    int end_position = 0;
    do {
        ValueType current = *(begin + end_position);
        do {
            unsigned new_bucket = GetKey(current, u_min, log_div);
            position = --offsets[new_bucket];
            std::swap(*(begin + position), current);
        } while (position != end_position);
        do {
            end_position += histogram[bucket++];
        } while (bucket != histogram.size() && end_position == offsets[bucket]);
    } while (bucket != histogram.size());
    if (log_div == 0) {
        return;
    }

    for (unsigned i = 0; i + 1 < offsets.size(); ++i) {
        if (offsets[i + 1] - offsets[i] <= 1) {
            continue;
        }
        if (offsets[i + 1] - offsets[i] < ELEMENTS_TO_RUN_STD_SORT) {
            std::sort(begin + offsets[i], begin + offsets[i + 1]);
        } else {
            InCacheSort<RandomAccessIterator, log_total_buckets>(begin + offsets[i], begin + offsets[i + 1], u_min + (i << log_div), u_min + ((i + 1) << log_div) - 1);
        }
    }
    if (end - (begin + offsets.back()) > 1) {
        if (end - (begin + offsets.back()) < ELEMENTS_TO_RUN_STD_SORT) {
            std::sort(begin + offsets.back(), end);
        } else {
            InCacheSort<RandomAccessIterator, log_total_buckets>(begin + offsets.back(), end, u_min + ((offsets.size() - 1) << log_div), u_max);
        }
    }
}


template<class RandomAccessIterator, uint32_t log_total_buckets>
void InCacheSort(RandomAccessIterator begin, RandomAccessIterator end) {
    using ValueType = typename std::iterator_traits<RandomAccessIterator>::value_type;
    RandomAccessIterator min;
    RandomAccessIterator max;
    if (is_sorted_or_find_extremes(begin, end, max, min)) {
        return;
    }

    /*
     * Map extremes into unsigned is such way:
     * x -> x + (-int::min)
     * int::min -> 0
     * 0 -> int::max - 1
     * int::max -> unsigned::max
     */
    unsigned u_min = (*min) ^ (1U << (8 * sizeof(ValueType) - 1));
    unsigned u_max = (*max) ^ (1U << (8 * sizeof(ValueType) - 1));
    const uint32_t log_range = rough_log_2_size(u_max - u_min);
    const uint32_t log_div = log_range > log_total_buckets ? log_range - log_total_buckets : 0;

    std::vector<uint32_t> histogram(1U << log_total_buckets, 0);
    for (RandomAccessIterator item = begin; item != end; ++item) {
        histogram[GetKey(*item, u_min, log_div)]++;
    }
    std::vector<uint32_t> offsets(1U << log_total_buckets);
    offsets[0] = histogram[0];
    for (int i = 1; i < histogram.size(); ++i) {
        offsets[i] = offsets[i - 1] + histogram[i];
    }

    int bucket = 0;
    while (offsets[bucket] == 0) {
        ++bucket;
    }

    // filling buckets from right to left
    int position;
    int end_position = 0;
    do {
        ValueType current = *(begin + end_position);
        do {
            unsigned new_bucket = GetKey(current, u_min, log_div);
            position = --offsets[new_bucket];
            std::swap(*(begin + position), current);
        } while (position != end_position);
        do {
            end_position += histogram[bucket++];
        } while (bucket != histogram.size() && end_position == offsets[bucket]);
    } while (bucket != histogram.size());
    if (log_div == 0) {
        return;
    }
    for (unsigned i = 0; i + 1 < offsets.size(); ++i) {
        if (offsets[i + 1] - offsets[i] <= 1) {
            continue;
        }
        if (offsets[i + 1] - offsets[i] < ELEMENTS_TO_RUN_STD_SORT) {
            std::sort(begin + offsets[i], begin + offsets[i + 1]);
        } else {
            InCacheSort<RandomAccessIterator, log_total_buckets>(begin + offsets[i], begin + offsets[i + 1], u_min + (i << log_div), u_min + ((i + 1) << log_div) - 1);
        }
    }
    if (end - (begin + offsets.back()) > 1) {
        if (end - (begin + offsets.back()) > 1) {
            if (end - (begin + offsets.back()) < ELEMENTS_TO_RUN_STD_SORT) {
                std::sort(begin + offsets.back(), end);
            } else {
                InCacheSort<RandomAccessIterator, log_total_buckets>(begin + offsets.back(), end, u_min + ((offsets.size() - 1) << log_div), u_max);
            }
        }
    }
}

template<class RandomAccessIterator, uint32_t log_total_buckets>
void OutOfCacheSort(RandomAccessIterator begin, RandomAccessIterator end) {
    using ValueType = typename std::iterator_traits<RandomAccessIterator>::value_type;

    RandomAccessIterator min;
    RandomAccessIterator max;
    if (is_sorted_or_find_extremes(begin, end, max, min)) {
        return;
    }

    /*
     * Map extremes into unsigned is such way:
     * x -> x + (-int::min)
     * int::min -> 0
     * 0 -> int::max - 1
     * int::max -> unsigned::max
     */
    unsigned u_min = (*min) ^ (1U << (8 * sizeof(ValueType) - 1));
    unsigned u_max = (*max) ^ (1U << (8 * sizeof(ValueType) - 1));
    const uint32_t log_range = rough_log_2_size(u_max - u_min);
    const uint32_t log_div = log_range > log_total_buckets ? log_range - log_total_buckets : 0;

    const int total_buckets = 1U << log_total_buckets;
    std::vector<size_t> histogram(total_buckets, 0);
    for (RandomAccessIterator item = begin; item != end; ++item) {
        histogram[GetKey(*item, u_min, log_div)]++;
    }
    std::vector<size_t> buckets_ends(total_buckets);
    buckets_ends[0] = histogram[0];
    for (int i = 1; i < total_buckets; ++i) {
        buckets_ends[i] = buckets_ends[i - 1] + histogram[i];
    }

    /*
     * In this sort we will use software buffer for storing and permutating elements.
     * There is buffer (buffers[p]) for every bucket (p). Each buffer is size of cache line.
     * We will read cache line from initial array, permutate elements and write cache line back, when buffer is ready.
     */
    const size_t items_per_cache_line = CACHE_LINE_SIZE / sizeof(ValueType);
    //std::vector<std::vector<ValueType>> buffers(total_buckets, std::vector<ValueType>(items_per_cache_line));
    ValueType buffers[total_buckets][items_per_cache_line] __attribute__ ((aligned (64)));

    /*
     * We will insert correct elements from the buffer start, so the end will contain dirty elements,
     * which should be move into some other buffer. The position of first dirty element is in the dirty_positions[p].
     */
    std::vector<size_t> dirty_positions(total_buckets);

    /*
     * When buffer is full of correct elements, we will load it into the initial array. The position of where
     * we should write our cache line contains in the bucket_offsets[p]. Array's elements before offset are
     * either in other bucket or already in right position.
     */
    std::vector<size_t> bucket_offsets(total_buckets);

    /*
     * First and last cache lines of a bucket could be half-empty. In this case we will fill buffer starting from offset,
     * so start of buffer will be empty. Start of any (dirty or correct) elements contains in buffer_offsets[p].
     */
    std::vector<size_t> buffer_offsets(total_buckets);

    size_t start = 0;
    for (size_t p = 0; p < total_buckets; ++p) {
        /* Let's assume first bucket is big enough and first elements in the array look like this
         * 0 1 2 3 4 5 6 7 8 9
         * _ _ _ / \ _ _ _ _ _
         * Underscores mean placements in cache lines.
         *
         * In this case we want take only four first elements in starting buffers.
         * We will place them at the end of buffers's cache line, so buffers[p] looks like:
         * {..., 0, 1, 2, 3}
         * In case of 64 bytes per cache line, and 4-bytes ints, the buffers will starts at position 12 of 16.
         * buffer_offsets[p] contains this start position.
         */

        buffer_offsets[p] = (reinterpret_cast<size_t>(&(*(begin + start))) % CACHE_LINE_SIZE) / sizeof(ValueType);
        if (buckets_ends[p] - start < items_per_cache_line) {
            // We need special treatment in case of small buckets
            buffer_offsets[p] = std::max(buffer_offsets[p], items_per_cache_line - (buckets_ends[p] - start));
        }

        // Copying elements into the buffer
        for (int j = 0; j < (items_per_cache_line - buffer_offsets[p]); ++j) {
            buffers[p][j + buffer_offsets[p]] = *(begin + start + j);
        }
        dirty_positions[p] = buffer_offsets[p];
        bucket_offsets[p] = start;
        start += histogram[p];
    }

    ValueType current;
    size_t start_bucket = 0;
    size_t bucket = 1;
    size_t pos;

    size_t first_dirty_bucket = 0;
    // Find first non-empty bucket
    while (start_bucket < total_buckets && bucket_offsets[start_bucket] == buckets_ends[start_bucket]) {
        ++start_bucket;
    }
    first_dirty_bucket = start_bucket;
    current = buffers[start_bucket][dirty_positions[start_bucket]];
    bool is_done = false;

    /*
     * The main loop works in such way: we choose first dirty item to start. Next we permutate elements
     * until any of two events occur: either the permutation cycle is over or we filled the cache line. If we done a
     * permutation cycle, we choose next dirty element and repeat the cycle. In the second case, we place cache line
     * in the initial array and load next cache line into the buffer.
     */
    do {
        do {
            if (__builtin_expect((bucket == start_bucket), 0)) { // haven't check if the difference is significant
                // permutation cycle is over, so current is already placed into correct buffer
                while (first_dirty_bucket < total_buckets && bucket_offsets[first_dirty_bucket] == buckets_ends[first_dirty_bucket]) {
                    ++first_dirty_bucket;
                }
                if (first_dirty_bucket == total_buckets) {
                    // all buckets are full, all elements are in the right places
                    is_done = true;
                    break;
                }
                start_bucket = first_dirty_bucket;
                current = buffers[start_bucket][dirty_positions[start_bucket]];
            }
            bucket = GetKey(current, u_min, log_div);
            pos = dirty_positions[bucket]++;
            std::swap(buffers[bucket][pos], current);
        } while (pos != items_per_cache_line - 1);
        if (is_done) {
            break;
        }

        // write cache line into the initial array
        for (size_t i = buffer_offsets[bucket]; i < items_per_cache_line; ++i) {
            *(begin + bucket_offsets[bucket] + i - buffer_offsets[bucket]) = buffers[bucket][i];
        }

        /*
         * There is possibility to write buffer straight into memory, without caching. I tried _mm256_stream_si256,
         * but saw no improvements. There could be difference with AVX512, but unfortunately I haven't got these instructions.
         */

        /*
        if (__builtin_expect((buffer_offsets[bucket] == 0), 1)) {
            //_mm512_stream_si512(reinterpret_cast<__m512i *>(&(*(begin + bucket_offsets[bucket]))), *(reinterpret_cast<__m512i*>((buffers[bucket]))));
            _mm256_stream_si256(reinterpret_cast<__m256i*>(&(*(begin + bucket_offsets[bucket]))), *(reinterpret_cast<__m256i*>((buffers[bucket]))));
            _mm256_stream_si256(reinterpret_cast<__m256i*>(&(*(begin + bucket_offsets[bucket] + 32 / sizeof(ValueType)))), *(reinterpret_cast<__m256i*>((buffers[bucket] + 32 / sizeof(ValueType)))));
        }  else {
            for (size_t i = buffer_offsets[bucket]; i < items_per_cache_line; ++i) {
                *(begin + bucket_offsets[bucket] + i - buffer_offsets[bucket]) = buffers[bucket][i];
            }
        } */

        bucket_offsets[bucket] += items_per_cache_line - buffer_offsets[bucket];
        if (__builtin_expect((bucket_offsets[bucket] != buckets_ends[bucket]), 1)) {
            // load next cache line
            if (__builtin_expect((buckets_ends[bucket] - bucket_offsets[bucket] < items_per_cache_line), 0)) {
                buffer_offsets[bucket] = items_per_cache_line - (buckets_ends[bucket] - bucket_offsets[bucket]);
                for (int i = 0; i < items_per_cache_line - buffer_offsets[bucket]; ++i) {
                    buffers[bucket][i + buffer_offsets[bucket]] = *(begin + bucket_offsets[bucket] + i);
                }
                dirty_positions[bucket] = buffer_offsets[bucket];
            } else {
                buffer_offsets[bucket] = 0;
                for (int i = 0; i < items_per_cache_line; ++i) {
                    buffers[bucket][i] = *(begin + bucket_offsets[bucket] + i);
                }
                dirty_positions[bucket] = 0;
            }
        }
    } while (true);
    if (log_div == 0) {
        return;
    }

    if (bucket_offsets.front() > 1) {
        if (bucket_offsets.front() < ELEMENTS_TO_RUN_STD_SORT) {
            std::sort(begin, begin + bucket_offsets.front());
        } else if (bucket_offsets.front() < ELEMENTS_TO_RUN_IN_CACHE) {
            InCacheSort<std::vector<int>::iterator, log_total_buckets>(begin , begin + bucket_offsets.front(), u_min, u_min + (1 << log_div) - 1);
        }
        else {
            OutOfCacheSort<RandomAccessIterator, log_total_buckets>(begin, begin + bucket_offsets.front());
        }
    }
    for (size_t i = 0; i + 1 < total_buckets; ++i) {
        if (bucket_offsets[i + 1] - bucket_offsets[i] <= 1) {
            continue;
        }
        if (bucket_offsets[i + 1] - bucket_offsets[i] < ELEMENTS_TO_RUN_STD_SORT) {
            std::sort(begin + bucket_offsets[i], begin + bucket_offsets[i + 1]);
        } else if (bucket_offsets[i + 1] - bucket_offsets[i] < ELEMENTS_TO_RUN_IN_CACHE) {
            InCacheSort<std::vector<int>::iterator, log_total_buckets >(begin + bucket_offsets[i], begin + bucket_offsets[i + 1], u_min + ((i + 1) << log_div), u_min + ((i + 2) << log_div) - 1);
        }
        else {
            OutOfCacheSort<RandomAccessIterator, log_total_buckets>(begin + bucket_offsets[i], begin + bucket_offsets[i + 1]);
        }
    }
}