#include <iostream>
#include <stdio.h>

#include <chrono>

#include "Shor.h"
#include "qft.h"
#include "quantumRegister.h"
#include <random>
#include <fstream>

#define PI 3.14159265

template <class T>
std::vector <std::vector<T>> Multiply(std::vector <std::vector<T>>& a, std::vector <std::vector<T>>& b)
{
    const int n = a.size();     // a rows
    const int m = a[0].size();  // a cols
    const int p = b[0].size();  // b cols

    std::vector <std::vector<T>> c(n, std::vector<T>(p, 0));
    for (auto j = 0; j < p; ++j)
    {
        for (auto k = 0; k < m; ++k)
        {
            for (auto i = 0; i < n; ++i)
            {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return c;
}


std::vector<std::vector<amp>> generateNoiseValues(double micro, double sigma, double &noise);

int main(int argc, char** argv)
{
    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();
    //int N = 2813;
   
    std::vector<int> Nval = { 575};
        
    double avgTime = 0;
    double avgTimeQFT = 0;

    ofstream myfile;
    myfile.open("output.csv");
    for (auto N : Nval)
    {
        printf("Factoring %d\n", N);

        for (int i = 0; i < 1000; i++)
        {
            __int64 timeElapsedForQFT = 0;
            __int64 guessCount = 0;
            double noiseValue = 0;
            start = chrono::steady_clock::now();
            auto noises = generateNoiseValues(5, 2, noiseValue);
            int p = Shor(N, 100, timeElapsedForQFT, guessCount, noises);
            int q = N / p;
            end = chrono::steady_clock::now();
            auto time = chrono::duration_cast<chrono::milliseconds>(end - start).count();
            avgTime += time;
            avgTimeQFT += timeElapsedForQFT;
            std::string strToPrint = std::to_string(N) + ", " + std::to_string(time) + ", " + std::to_string(timeElapsedForQFT) + ", " + std::to_string(guessCount) + ", " + std::to_string(noiseValue) + "\n";
            myfile << strToPrint;
            printf("Found factor of %d to be %d in %ld ms in %d guesses. QFT time: %ld; noise value is %f \n ", N, p, time, guessCount, timeElapsedForQFT, noiseValue);
        }
        avgTime /= 3;
        avgTimeQFT /= 3;

        printf("Average elapsed time for factoring: %lf ms\n", avgTime);

        printf("Average Time for QFT: %lf ms \n ", avgTimeQFT);

    }
    myfile.close();

    return 0;
}

std::vector<std::vector<amp>> generateNoiseValues(double micro, double sigma, double &noise)
{
    std::normal_distribution<> d{ micro, sigma };
    std::random_device rd{};
    std::mt19937 gen{ rd() };

    auto signVal = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());

    noise = d(gen);
    bool b = signVal();
    if(b)
        noise *= -1;

    double angle1 = 45  + noise; // deg
    double angle2 = 180 + noise; // deg

    std::vector<std::vector<amp>> R;
    std::vector<amp> r1;
    std::vector<amp> r2;
    r1.push_back(cos(angle1 * PI / 180));
    r1.push_back(-sin(angle1 * PI / 180));
    r2.push_back(sin(angle1 * PI / 180));
    r2.push_back(cos(angle1 * PI / 180));
    
    R.push_back(r1);
    R.push_back(r2);

    amp i(0, 1);
    std::vector<std::vector<amp>> P;
    std::vector<amp> p1;
    std::vector<amp> p2;
    p1.push_back(1);
    p1.push_back(0);
    p2.push_back(0);
    p2.push_back(exp(i*(angle2* PI / 180)));

    P.push_back(p1);
    P.push_back(p2);

    std::vector<std::vector<amp>> res = Multiply(R, P);
    return res;
}

