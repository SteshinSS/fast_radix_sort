#pragma once

#include <iostream>
#include <vector>
#include <algorithm>

constexpr int MINIMUM_ELEMENTS = 1000;

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

    int i;
    int i_end = 0;
    do {
        ValueType t = *(begin + i_end);
        do {
            unsigned new_bucket = GetKey(t, u_min, log_div);
            i = --offsets[new_bucket];
            std::swap(*(begin + i), t);
        } while (i != i_end);
        do {
            i_end += histogram[bucket++];
        } while (bucket != histogram.size() && i_end == offsets[bucket]);
    } while (bucket != histogram.size());

    if (log_div == 0) {
        return;
    }

    for (unsigned i = 0; i + 1 < offsets.size(); ++i) {
        if (offsets[i + 1] - offsets[i] <= 1) {
            continue;
        }
        if (offsets[i + 1] - offsets[i] < MINIMUM_ELEMENTS) {
            std::sort(begin + offsets[i], begin + offsets[i + 1]);
        } else {
            InCacheSort<RandomAccessIterator, log_total_buckets>(begin + offsets[i], begin + offsets[i + 1], u_min + (i << log_div), u_min - 1 + ((i + 1) << log_div));
        }
    }
    if (end - (begin + offsets.back()) > 1) {
        if (end - (begin + offsets.back()) < MINIMUM_ELEMENTS) {
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

    int i;
    int i_end = 0;
    do {
        ValueType t = *(begin + i_end);
        do {
            unsigned new_bucket = GetKey(t, u_min, log_div);
            i = --offsets[new_bucket];
            std::swap(*(begin + i), t);
        } while (i != i_end);
        do {
            i_end += histogram[bucket++];
        } while (bucket != histogram.size() && i_end == offsets[bucket]);
    } while (bucket != histogram.size());

    if (log_div == 0) {
        return;
    }
    for (unsigned i = 0; i + 1 < offsets.size(); ++i) {
        if (offsets[i + 1] - offsets[i] <= 1) {
            continue;
        }
        if (offsets[i + 1] - offsets[i] < MINIMUM_ELEMENTS) {
            std::sort(begin + offsets[i], begin + offsets[i + 1]);
        } else {
            InCacheSort<RandomAccessIterator, log_total_buckets>(begin + offsets[i], begin + offsets[i + 1], u_min + (i << log_div), u_min - 1 + ((i + 1) << log_div));
        }
    }
    if (end - (begin + offsets.back()) > 1) {
        if (end - (begin + offsets.back()) > 1) {
            if (end - (begin + offsets.back()) < MINIMUM_ELEMENTS) {
                std::sort(begin + offsets.back(), end);
            } else {
                InCacheSort<RandomAccessIterator, log_total_buckets>(begin + offsets.back(), end, u_min + ((offsets.size() - 1) << log_div), u_max);
            }
        }
    }
}

template<uint32_t log_total_buckets>
void OutOfCacheSort(const std::vector<int>::iterator& begin, const std::vector<int>::iterator& end) {
    const int total_buckets = 1U << log_total_buckets;
    std::vector<int>::iterator min;
    std::vector<int>::iterator max;
    if (is_sorted_or_find_extremes(begin, end, max, min)) {
        return;
    }
    unsigned u_min = (*min) ^ (1U << (8 * sizeof(int) - 1));
    unsigned u_max = (*max) ^ (1U << (8 * sizeof(int) - 1));
    const uint32_t log_range = rough_log_2_size(u_max - u_min);
    const uint32_t log_div = log_range > log_total_buckets ? log_range - log_total_buckets : 0;

    std::vector<uint32_t> histogram(total_buckets, 0);
    for (std::vector<int>::iterator item = begin; item != end; ++item) {
        histogram[GetKey(*item, u_min, log_div)]++;
    }
    std::vector<int> ends(total_buckets);
    ends[0] = histogram[0];
    for (int i = 1; i < total_buckets; ++i) {
        ends[i] = ends[i - 1] + histogram[i];
    }

    std::vector<int> last_item(total_buckets, 0);
    std::vector<int> offset(total_buckets, 0);
    std::vector<size_t> buffer_offset(total_buckets, 0);

    const size_t cacheline_size = 64 / sizeof(int);
    // const size_t cacheline_size = 4 * sizeof(int);
    std::vector<std::vector<int>> buffer(total_buckets, std::vector<int>(cacheline_size, 0));
    uint32_t start = 0;
    for (int i = 0; i < total_buckets; ++i) {
        buffer_offset[i] = (reinterpret_cast<size_t>(&(*(begin + start))) % 64) / sizeof(int);
        //buffer_offset[i] = (begin + start)
        if (ends[i] - start < cacheline_size) {
            buffer_offset[i] = std::max(buffer_offset[i], cacheline_size - (ends[i] - start));
        }

        for (int j = 0; j < (cacheline_size - buffer_offset[i]); ++j) {
            buffer[i][j + buffer_offset[i]] = *(begin + start + j);
        }
        last_item[i] = buffer[i][cacheline_size - 1];
        buffer[i][cacheline_size - 1] = buffer_offset[i];
        offset[i] = start;
        start += histogram[i];
    }
    int current = *begin;
    int bucket;
    int pos;
    
    do {
        std::cout << "hello" << std::endl;
        do {
            bucket = GetKey(current, u_min, log_div);
            std::cout << bucket << std::endl;
            pos = buffer[bucket][cacheline_size - 1]++;
            std::swap(buffer[bucket][pos], current);

        } while (pos != cacheline_size - 1);
        for (int i = buffer_offset[bucket]; i < cacheline_size; ++i) {
            *(begin + offset[bucket] + i - buffer_offset[bucket]) = buffer[bucket][i];
        }
        offset[bucket] += cacheline_size - buffer_offset[bucket];
        current = last_item[bucket];
        if (offset[bucket] == ends[bucket]) {
            bucket = 0;
            while (bucket < total_buckets && offset[bucket] == ends[bucket]) {
                ++bucket;
            }
            if (bucket == total_buckets) {
                break;
            } else {

                continue;
            }
        }

        if (ends[bucket] - offset[bucket] < cacheline_size) {
            buffer_offset[bucket] = cacheline_size - (ends[bucket] - offset[bucket]);
            for (int i = 0; i < cacheline_size - buffer_offset[bucket]; ++i) {
                buffer[bucket][i + buffer_offset[bucket]] = *(begin + offset[bucket] + i);
            }
            last_item[bucket] = buffer[bucket][cacheline_size - 1];
            buffer[bucket][cacheline_size - 1] = buffer_offset[bucket];
        } else {
            buffer_offset[bucket] = 0;
            for (int i = 0; i < cacheline_size; ++i) {
                buffer[bucket][i] = *(begin + offset[bucket] + i);
            }
            last_item[bucket] = buffer[bucket][cacheline_size - 1];
            buffer[bucket][cacheline_size - 1] = 0;
        }
    } while (true);
}