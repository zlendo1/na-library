#pragma once

#include "Matrix.cpp"

class QRDecomposer {
    Matrix matrica;
    std::vector<double> RDijag;

public:
    QRDecomposer(Matrix m);

    void Solve(const Vector &b, Vector &x) const;
    Vector Solve(Vector b) const;
    void Solve(Matrix &b, Matrix &x) const;
    Matrix Solve(Matrix b) const;

    Vector MulQWith(Vector v) const;
    Matrix MulQWith(Matrix m) const;

    Vector MulQTWith(Vector v) const;
    Matrix MulQTWith(Matrix m) const;

    Matrix GetQ() const;
    Matrix GetR() const;
};

QRDecomposer::QRDecomposer(Matrix mat) : matrica(std::move(mat)), RDijag(matrica.NCols()) {
    if (matrica.NRows() < matrica.NCols()) {
        throw std::domain_error("Invalid matrix format");
    }

    std::vector<std::vector<double>> &a = matrica.elementi;
    const int n = matrica.NCols(), m = matrica.NRows();
    const double epsilon = matrica.GetEpsilon();

    for (int k = 0; k < n; ++k) {
        long double s = 0;

        for (int i = k; i < m; ++i) {
            s += (long double) a[i][k] * a[i][k];
        }

        s = std::sqrt(s);

        long double mi = std::sqrt(s * (s + std::abs(a[k][k])));

        if (mi < epsilon) {
            throw std::domain_error("Matrix is singular");
        }

        if (a[k][k] < 0) {
            s = -s;
        }

        a[k][k] = (a[k][k] + s) / mi;

        for (int i = k + 1; i < m; ++i) {
            a[i][k] /= mi;
        }

        RDijag[k] = -s;

        for (int j = k + 1; j < n; ++j) {
            s = 0;

            for (int i = k; i < m; ++i) {
                s += (long double) a[i][k] * a[i][j];
            }

            for (int i = k; i < m; ++i) {
                a[i][j] -= s * a[i][k];
            }
        }
    }
}

void QRDecomposer::Solve(const Vector &b, Vector &x) const {
    matrica.CheckSquareness();

    if (b.NElems() != matrica.NRows() || x.NElems() != matrica.NCols()) {
        throw std::domain_error("Incompatible formats");
    }

    Vector y = std::move(MulQTWith(b));
    const int n = b.NElems();

    for (int i = n - 1; i >= 0; --i) {
        long double s = y[i];

        for (int j = i + 1; j < n; ++j) {
            s -= (long double) matrica[i][j] * x[j];
        }

        x[i] = s / RDijag[i];
    }
}

Vector QRDecomposer::Solve(Vector b) const {
    matrica.CheckSquareness();

    Solve(b, b);

    return b;
}

void QRDecomposer::Solve(Matrix &b, Matrix &x) const {
    matrica.CheckSquareness();
    Matrix::CheckMultiplicativity(matrica, x);
    Matrix::CheckComparability(b, x);

    Matrix y = std::move(MulQTWith(b));
    const int n = b.NRows(), m = b.NCols();

    for (int i = n - 1; i >= 0; --i) {
        std::vector<double> &s = y.elementi[i];

        for (int k = 0; k < m; ++k) {
            for (int j = i + 1; j < n; ++j) {
                s[k] -= matrica[i][j] * x[j][k];
            }

            x[i][k] = s[k] / RDijag[i];
        }
    }
}

Matrix QRDecomposer::Solve(Matrix b) const {
    Solve(b, b);

    return b;
}

Vector QRDecomposer::MulQWith(Vector vec) const {
    if (matrica.NCols() != vec.NElems()) {
        throw std::domain_error("Incompatible formats");
    }

    if (!matrica.IsSquare()) { //Left as is with the explicit consent of professor Juric
        throw std::domain_error("Currently not defined for non-square decomposed matrices");
    }

    const std::vector<std::vector<double>> &a = matrica.elementi;
    std::vector<double> &v = vec.elementi;
    const int m = matrica.NRows(), n = matrica.NCols();

    for (int k = n - 1; k >= 0; --k) {
        long double s = 0;

        for (int i = k; i < m; ++i) {
            s += (long double) a[i][k] * v[i];
        }

        for (int i = k; i < m; ++i) {
            v[i] -= s * a[i][k];
        }
    }

    return vec;
}

Matrix QRDecomposer::MulQWith(Matrix mat) const {
    if (matrica.NCols() != mat.NRows()) {
        throw std::domain_error("Incompatible formats");
    }

    if (!matrica.IsSquare()) { //Left as is with the explicit consent of professor Juric
        throw std::domain_error("Currently not defined for non-square decomposed matrices");
    }

    const std::vector<std::vector<double>> &a = matrica.elementi;
    std::vector<std::vector<double>> &M = mat.elementi;
    const int m = matrica.NRows(), n = matrica.NCols(), matNCols = mat.NCols();

    for (int j = 0; j < matNCols; ++j) {
        for (int k = n - 1; k >= 0; --k) {
            long double s = 0;

            for (int i = k; i < m; ++i) {
                s += (long double) a[i][k] * M[i][j];
            }

            for (int i = k; i < m; ++i) {
                M[i][j] -= s * a[i][k];
            }
        }
    }

    return mat;
}

Vector QRDecomposer::MulQTWith(Vector vec) const {
    if (matrica.NRows() != vec.NElems()) {
        throw std::domain_error("Incompatible formats");
    }

    if (!matrica.IsSquare()) { //Left as is with the explicit consent of professor Juric
        throw std::domain_error("Currently not defined for non-square decomposed matrices");
    }

    const std::vector<std::vector<double>> &a = matrica.elementi;
    std::vector<double> &v = vec.elementi;
    const int m = matrica.NRows(), n = matrica.NCols();

    for (int k = 0; k < n; ++k) {
        long double s = 0;

        for (int i = k; i < m; ++i) {
            s += (long double) a[i][k] * v[i];
        }

        for (int i = k; i < m; ++i) {
            v[i] -= s * a[i][k];
        }
    }

    return vec;
}

Matrix QRDecomposer::MulQTWith(Matrix mat) const {
    if (matrica.NRows() != mat.NRows()) {
        throw std::domain_error("Incompatible formats");
    }

    if (!matrica.IsSquare()) { //Left as is with the explicit consent of professor Juric
        throw std::domain_error("Currently not defined for non-square decomposed matrices");
    }

    const std::vector<std::vector<double>> &a = matrica.elementi;
    std::vector<std::vector<double>> &M = mat.elementi;
    const int m = matrica.NRows(), n = matrica.NCols(), matNCols = mat.NCols();

    for (int j = 0; j < matNCols; ++j) {
        for (int k = 0; k < n; ++k) {
            long double s = 0;

            for (int i = k; i < m; ++i) {
                s += (long double) a[i][k] * M[i][j];
            }

            for (int i = k; i < m; ++i) {
                M[i][j] -= s * a[i][k];
            }
        }
    }

    return mat;
}

Matrix QRDecomposer::GetQ() const {
    Matrix Q(matrica.NRows(), matrica.NCols());

    const std::vector<std::vector<double>> &a = matrica.elementi;
    std::vector<std::vector<double>> &q = Q.elementi;
    const int m = Q.NRows(), n = Q.NCols();

    for (int j = 0; j < n; ++j) {
        q[j][j] =1;

        for (int k = n - 1; k >= 0; --k) {
            long double s = 0;

            for (int i = k; i < m; ++i) {
                s += (long double) a[i][k] * q[i][j];
            }

            for (int i = k; i < m; ++i) {
                q[i][j] -= s * a[i][k];
            }
        }
    }

    return Q;
}

Matrix QRDecomposer::GetR() const {
    Matrix R(matrica.NCols(), matrica.NCols());

    const std::vector<std::vector<double>> &a = matrica.elementi;
    std::vector<std::vector<double>> &r = R.elementi;
    const int n = R.NCols();

    for (int i = 0; i < n; ++i) {
        r[i][i] = RDijag[i];

        for (int j = i + 1; j < n; ++j) {
            r[i][j] = a[i][j];
        }
    }

    return R;
}