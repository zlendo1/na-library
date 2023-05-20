#ifndef COMPRESSION_H
#define COMPRESSION_H

#include <stdexcept>
#include <cmath>
#include <vector>
#include <complex>

const double PI = 4 * std::atan(1);

std::complex<double> W(double N) {
    return std::complex<double>(std::cos(2 * PI / N), std::sin(2 * PI / N));
}

void FFT(std::vector<std::complex<double>> &original, std::vector<std::complex<double>> &transform, int N, int s = 0, int d = 0, int t = 1) {
    if (N == 1) {
        transform[d] = original[s];
    } else {
        std::complex<double> mi = 1;
        std::complex<double> w = 1. / W(N);

        for (int k = s; k < s + N / 2; ++k) {
            std::complex<double> u = original[k];

            original[k] = u + original[k + N / 2];

            original[k + N / 2] = mi * (u - original[k + N / 2]);

            mi *= w;
        }

        FFT(original, transform, N / 2, s, d, 2 * t);
        FFT(original, transform, N / 2, s + N / 2, d + t, 2 * t);
    }
}

void InverseFFT(std::vector<std::complex<double>> &transform, std::vector<std::complex<double>> &original, int N, int s = 0, int d = 0, int t = 1) {
    if (N == 1) {
        original[d] = transform[s];
    } else {
        InverseFFT(transform, original, N / 2, s, d, 2 * t);
        InverseFFT(transform, original, N / 2, s + t, d + N / 2, 2 * t);

        std::complex<double> mi = 1;
        std::complex<double> w = W(N);

        for(int k = d; k < d + N / 2; ++k) {
            std::complex<double> u = original[k];
            std::complex<double> v = mi * original[k + N / 2];

            original[k] = (u + v) / 2.;
            original[k + N / 2] = (u - v) / 2.;

            mi *= w;
        }
    }
}

bool FFTSizeCompatible(int N) {
    return std::log2(N) == std::floor(std::log2(N));
}

std::vector<double> LossyCompress(std::vector<double> data, int new_size) {
    int N = data.size(), M = new_size;

    if (!FFTSizeCompatible(N)) {
        throw std::range_error("Data size must be a power of two");
    }

    if (M <= 1 || new_size > N) {
        throw std::range_error("Bad new size");
    }

    std::vector<std::complex<double>> y(N);

    for (int i = 0; i < N / 2; ++i) {
        y[i] = data[2 * i];
    }

    for (int i = N / 2; i < N; ++i) {
        y[i] = data[2 * (N - i) - 1];
    }

    std::vector<std::complex<double>> yp(N);

    FFT(y, yp, N);

    std::vector<double> xp(M);

    std::complex<double> mi = 1;
    std::complex<double> w = std::pow(W(2 * N), -0.5);

    for (int k = 0; k < M - 1; ++k) {
        xp[k] = std::real(mi * yp[k]);

        mi *= w;
    }

    xp.back() = N;

    return xp;
}

std::vector<double> LossyDecompress(std::vector<double> compressed) {
    const int M = compressed.size(), N = compressed.back();

    compressed.pop_back();

    if (!FFTSizeCompatible(N) || N < M) {
        throw std::logic_error("Bad compressed sequence");
    }

    compressed.resize(N, 0);

    std::vector<std::complex<double>> yp(N);

    std::complex<double> mi = 1;
    std::complex<double> w = std::pow(W(2 * N), 0.5);

    yp[0] = compressed[0];
    
    for (int k = 1; k < M - 1; ++k) {
        mi *= w;

        yp[k] = 2. * mi * compressed[k];
    }

    std::vector<std::complex<double>> y(N);

    InverseFFT(yp, y, N);

    std::vector<double> x(N);

    for (int k = 0; k < N; ++k) {
        if (k % 2 == 0) {
            x[k] = std::real(y[k / 2]);
        } else {
            x[k] = std::real(y[N - (k + 1) / 2]);
        }
    }

    return x;
}

#endif