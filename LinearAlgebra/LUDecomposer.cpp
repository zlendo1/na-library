#pragma once

#include "Matrix.cpp"

class LUDecomposer {
    Matrix matrica;
    std::vector<int> permutacije;

public:
    LUDecomposer(Matrix m);

    void Solve(const Vector &b, Vector &x) const;
    Vector Solve(Vector b) const;
    void Solve(Matrix &b, Matrix &x) const;
    Matrix Solve(Matrix b) const;

    Matrix GetCompactLU() const;
    Matrix GetL() const;
    Matrix GetU() const;
    Vector GetPermuation() const;
};

LUDecomposer::LUDecomposer(Matrix m) : matrica(std::move(m)), permutacije(matrica.NRows()) {
    matrica.CheckSquareness();

    std::vector<std::vector<double>> &a = matrica.elementi;
    std::vector<int> &w = permutacije;
    const int n = matrica.NRows();
    const double epsilon = matrica.GetEpsilon();

    for (int j = 0; j < n; ++j) {
        long double s;
        
        for (int i = 0; i <= j; ++i) {
            s = a[i][j];

            for (int k = 0; k < i; ++k) {
                s -= (long double)a[i][k] * a[k][j];
            }

            a[i][j] = s;
        }

        int p = j;

        for (int i = j + 1; i < n; ++i) {
            s = a[i][j];

            for (int k = 0; k < j; ++k) {
                s -= (long double) a[i][k] * a[k][j];
            }

            a[i][j] = s;

            if (std::abs(s) > std::abs(a[p][j])) {
                p = i;
            }
        }

        if (std::abs(a[p][j]) < epsilon) {
            throw std::domain_error("Matrix is singular");
        }

        matrica.SwapRows(j, p);

        w[j] = p;
        
        long double mi = a[j][j];

        for (int i = j + 1; i < n; ++i) {
            a[i][j] /= mi;
        }
    }
}

void LUDecomposer::Solve(const Vector &b, Vector &x) const {
    if (b.NElems() != matrica.NRows() || x.NElems() != matrica.NCols()) {
        throw std::domain_error("Incompatible formats");
    }

    std::vector<double> new_b(b.elementi);
    std::vector<long double> y(b.NElems());
    const int n = matrica.NRows();

    for (int i = 0; i < n; ++i) {
        long double s = new_b[permutacije[i]];
        new_b[permutacije[i]] = new_b[i];

        for (int j = 0; j < i; ++j) {
            s -= matrica[i][j] * y[j];
        }

        y[i] = s;
    }

    for (int i = n - 1; i >= 0; --i) {
        long double s = y[i];

        for (int j = i + 1; j < n; ++j) {
            s -= (long double) matrica[i][j] * x[j];
        }

        x[i] = s / matrica[i][i];
    }
}

Vector LUDecomposer::Solve(Vector b) const {
    Solve(b, b);

    return b;
}

void LUDecomposer::Solve(Matrix &b, Matrix &x) const {
    Matrix::CheckMultiplicativity(matrica, x);
    Matrix::CheckComparability(x, b);

    Matrix y(b.NRows(), b.NCols()), new_b(b);
    const int n = b.NRows();
    const int m = b.NCols();

    for (int i = 0; i < n; ++i) {
        int p = permutacije[i];

        std::vector<double> s = new_b.elementi[p];
        new_b.elementi[p] = new_b.elementi[i];

        for (int k = 0; k < m; ++k) {
            for (int j = 0; j < i; ++j) {
                s[k] -= matrica[i][j] * y[j][k];
            }
        }

        y.elementi[i] = std::move(s);
    }

    for (int i = n - 1; i >= 0; --i) {
        std::vector<double> &s = y.elementi[i];

        for (int k = 0; k < m; ++k) {
            for (int j = i + 1; j < n; ++j) {
                s[k] -= matrica[i][j] * x[j][k];
            }

            x[i][k] = s[k] / matrica[i][i];
        }
    }
}

Matrix LUDecomposer::Solve(Matrix b) const {
    Solve(b, b);

    return b;
}

Matrix LUDecomposer::GetCompactLU() const {
    return matrica;
}

Matrix LUDecomposer::GetL() const {
    Matrix rez(matrica.NRows(), matrica.NRows());

    for (int i = 0; i < matrica.NRows(); ++i) {
        for (int j = i; j >= 0; --j) {
            if (j == i) {
                rez.elementi[i][j] = 1;

                continue;
            }

            rez.elementi[i][j] = matrica.elementi[i][j];
        }
    }

    return rez;
}

Matrix LUDecomposer::GetU() const {
    Matrix rez(matrica.NRows(), matrica.NRows());

    for (int i = 0; i < matrica.NRows(); ++i) {
        for (int j = i; j < matrica.NRows(); ++j) {
            rez.elementi[i][j] = matrica.elementi[i][j];
        }
    }

    return rez;
}

Vector LUDecomposer::GetPermuation() const {
    Vector rez(permutacije.size());

    for (int i = 0; i < permutacije.size(); ++i) {
        rez.elementi[i] = permutacije[i] + 1;
    }

    return rez;
}