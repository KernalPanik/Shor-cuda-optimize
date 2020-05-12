#include <iostream>
#include <stdio.h>

#include <chrono>

#include "Shor.h"
#include "qft.h"
#include "quantumRegister.h"

int main(int argc, char** argv)
{
    //, 3901, 39211, 49949, 52939, 108673
    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();
    //int N = 2813;
   
    std::vector<int> Nval = { 575/*2813, 2813, 3901, 8927, 11387, 18419, 30883, 45173, 59431 */};
        
    double avgTime = 0;
    double avgTimeQFT = 0;

    for (auto N : Nval)
    {
        printf("Factoring %d\n", N);

        for (int i = 0; i < 1; i++)
        {
            __int64 timeElapsedForQFT = 0;
            start = chrono::steady_clock::now();
            int p = Shor(N, 100, timeElapsedForQFT);
            int q = N / p;
            end = chrono::steady_clock::now();
            auto time = chrono::duration_cast<chrono::milliseconds>(end - start).count();
            avgTime += time;
            avgTimeQFT += timeElapsedForQFT;
            printf("Found factor of %d to be %d in %ld ms QFT time: %ld \n ", N, p, time, timeElapsedForQFT);
        }
        avgTime /= 3;
        avgTimeQFT /= 3;

        printf("Average elapsed time for factoring: %lf ms\n", avgTime);

        printf("Average Time for QFT: %lf ms \n ", avgTimeQFT);

    }

    return 0;
}
