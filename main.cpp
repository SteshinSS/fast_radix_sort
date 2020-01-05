#include "semen_sort.h"
#include <boost/sort/spreadsort/spreadsort.hpp>
#include <iostream>
#include <vector>
#include <random>
#include <limits>
#include <algorithm>
#include <chrono>


void PrintArray(const std::vector<int>& vec) {
    for (int element : vec) {
        std::cout << element << " ";
    }
    std::cout << std::endl;
}


std::vector<int> GenerateArray(int n, int min, int max, std::mt19937& mersenne_engine) {
    std::uniform_int_distribution<int> dist {min, max};
    auto gen = [&dist, &mersenne_engine](){
        return dist(mersenne_engine);
    };

    std::vector<int> vec(n);
    std::generate(begin(vec), end(vec), gen);
    return vec;
}

void RunStressTest(int total_iterations) {
    std::random_device rnd_device;
    std::mt19937 mersenne_engine {rnd_device()};

    std::uniform_int_distribution<int> length_generator {0, 1000000};
    std::uniform_int_distribution<int> complexity_generator {std::numeric_limits<int>::min(), std::numeric_limits<int>::max()};

    for (int i = 0; i < total_iterations; ++i) {
        if (i % 10 == 0) {
            std::cout << "Iteration " << i << std::endl;
        }
        int length = length_generator(mersenne_engine);
        int a = complexity_generator(mersenne_engine);
        int b = complexity_generator(mersenne_engine);
        std::vector<int> vec = GenerateArray(length, std::min(a, b), std::max(a, b), mersenne_engine);
        std::vector<int> std_vec = vec;
        std::vector<int> semen_1_vec = vec;
        std::vector<int> semen_2_vec = vec;
        std::vector<int> semen_3_vec = vec;
        std::vector<int> semen_4_vec = vec;
        std::vector<int> semen_10_vec = vec;
        std::vector<int> semen_15_vec = vec;
        std::sort(std_vec.begin(), std_vec.end());

        InCacheSort<std::vector<int>::iterator, 1U>(semen_1_vec.begin(), semen_1_vec.end());
        InCacheSort<std::vector<int>::iterator, 2U>(semen_2_vec.begin(), semen_2_vec.end());
        InCacheSort<std::vector<int>::iterator, 3U>(semen_3_vec.begin(), semen_3_vec.end());
        InCacheSort<std::vector<int>::iterator, 4U>(semen_4_vec.begin(), semen_4_vec.end());
        InCacheSort<std::vector<int>::iterator, 10U>(semen_10_vec.begin(), semen_10_vec.end());
        InCacheSort<std::vector<int>::iterator, 15U>(semen_15_vec.begin(), semen_15_vec.end());

        if (std_vec != semen_1_vec || std_vec != semen_2_vec || std_vec != semen_3_vec ||
            std_vec != semen_4_vec || std_vec != semen_10_vec || std_vec != semen_15_vec) {
            std::cout << "Error: " << std::endl;
            PrintArray(vec);
        }
    }
}


void RunStressTest2(int total_iterations) {
    std::random_device rnd_device;
    std::mt19937 mersenne_engine {rnd_device()};

    std::uniform_int_distribution<int> length_generator {1, 1000000};
    std::uniform_int_distribution<int> complexity_generator {std::numeric_limits<int>::min(), std::numeric_limits<int>::max()};

    for (int i = 0; i < total_iterations; ++i) {
        if (i % 10 == 0) {
            std::cout << "Iteration " << i << std::endl;
        }
        int length = length_generator(mersenne_engine);
        int a = complexity_generator(mersenne_engine);
        int b = complexity_generator(mersenne_engine);
        std::vector<int> vec = GenerateArray(length, std::numeric_limits<int>::min(), std::numeric_limits<int>::max(), mersenne_engine);
        std::vector<int> std_vec = vec;
        std::vector<int> semen_1_vec = vec;
        std::vector<int> semen_2_vec = vec;
        std::vector<int> semen_3_vec = vec;
        std::vector<int> semen_4_vec = vec;
        std::vector<int> semen_10_vec = vec;
        std::vector<int> semen_15_vec = vec;
        std::sort(std_vec.begin(), std_vec.end());

        OutOfCacheSort<1U>(semen_1_vec.begin(), semen_1_vec.end());
        OutOfCacheSort<2U>(semen_2_vec.begin(), semen_2_vec.end());
        OutOfCacheSort<3U>(semen_3_vec.begin(), semen_3_vec.end());
        OutOfCacheSort<4U>(semen_4_vec.begin(), semen_4_vec.end());
        OutOfCacheSort<10U>(semen_10_vec.begin(), semen_10_vec.end());
        OutOfCacheSort<15U>(semen_15_vec.begin(), semen_15_vec.end());

        if (std_vec != semen_1_vec || std_vec != semen_2_vec || std_vec != semen_3_vec ||
            std_vec != semen_4_vec || std_vec != semen_10_vec || std_vec != semen_15_vec) {
            std::cout << "Error: " << std::endl;
            PrintArray(vec);
        }
    }
}

void RunBenchmarks() {
    std::random_device rnd_device;
    std::mt19937 mersenne_engine {rnd_device()};
    std::uniform_int_distribution<int> complexity_generator {std::numeric_limits<int>::min(), std::numeric_limits<int>::max()};

    int a = complexity_generator(mersenne_engine);
    int b = complexity_generator(mersenne_engine);
    std::vector<int> vec = GenerateArray(100, -2, 1, mersenne_engine);
    std::vector<int> copy = vec;
    std::vector<int> orig = vec;
    std::sort(copy.begin(), copy.end());
    //std::vector<int> vec = {2, 1, 1, -50, 60, 20, 500, -50};
    //PrintArray(vec);
    //std::sort(vec.begin(), vec.end());
    //std::reverse(vec.begin(), vec.end());
    auto start = std::chrono::high_resolution_clock::now();
    //InCacheSort<std::vector<int>::iterator, 10>(vec.begin(), vec.end());
    OutOfCacheSort<10>(vec.begin(), vec.end());
    //std::sort(vec.begin(), vec.end());
    //boost::sort::spreadsort::spreadsort(vec.begin(), vec.end());
    auto end = std::chrono::high_resolution_clock::now();
    PrintArray(vec);
    if (copy == vec) {
        std::cout << "OK" << std::endl;
    } else {
        std::cout << "Error: " << std::endl;
        PrintArray(orig);
    }
    std::cout << std::chrono::duration <double, std::micro>(end - start).count() << std::endl;
    //std::cout << std::min(a, b) << std::endl;
    //std::cout << std::max(a, b) << std::endl;
}

int main()
{
/*
    for (int i = 0; i < 1000; ++i) {
        RunBenchmarks();
    } */
    //RunBenchmarks();
    RunStressTest2(10000);
    return 0;
}