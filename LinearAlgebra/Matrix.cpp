#pragma once

#include "Vector.cpp"

class Matrix {
private:
    std::vector<std::vector<double>> elementi;

    static bool AreComparable(const Matrix &m1, const Matrix &m2);
    static bool AreMultiplyable(const Matrix &m1, const Matrix &m2);
    bool IsSquare() const;

    static void CheckComparability(const Matrix &m1, const Matrix &m2);
    static void CheckMultiplicativity(const Matrix &m1, const Matrix &m2);
    void CheckSquareness() const;

    int FindPivot(int from, int col) const;
    void SwapRows(int current_row, int new_row);

public:
    friend LUDecomposer;
    friend QRDecomposer;
    
    Matrix(int m, int n);
    Matrix(const Vector &v);
    Matrix(std::initializer_list<std::vector<double>> l);

    Matrix(const Matrix &m);
    Matrix(Matrix &&m);

    Matrix& operator=(const Matrix &m);
    Matrix& operator=(Matrix &&m);

    int NRows() const;
    int NCols() const;

    double* operator[](int i);
    const double* operator[](int i) const;
    double& operator()(int i, int j);
    double operator()(int i, int j) const;

    double Norm() const;
    friend double MatrixNorm(const Matrix &m);
    double GetEpsilon() const;

    void Print(int width, double eps) const;
    friend void PrintMatrix(const Matrix &m, int width, double eps);

    friend Matrix operator+(const Matrix &m1, const Matrix &m2);
    Matrix& operator+=(const Matrix &m);
    friend Matrix operator-(const Matrix &m1, const Matrix &m2);
    Matrix& operator-=(const Matrix &m);
    friend Matrix operator*(double s, const Matrix &m);
    friend Matrix operator*(const Matrix &m, double s);
    Matrix& operator*=(double s);
    friend Matrix operator *(const Matrix &m1, const Matrix &m2);
    Matrix& operator*=(const Matrix &m);
    friend Vector operator*(const Matrix &m, const Vector &v);

    friend Matrix Transpose(const Matrix &m);
    void Transpose();

    void Chop(double eps = -1);
    bool EqualTo(const Matrix &m, double eps = -1) const;

    friend Matrix LeftDiv(Matrix m1, Matrix m2);
    friend Vector LeftDiv(Matrix m, Vector v);
    friend Matrix operator/(const Matrix &m, double s);
    Matrix &operator/=(double s);
    friend Matrix operator/(Matrix m1, Matrix m2);
    Matrix &operator/=(Matrix m);

    double Det() const;
    friend double Det(Matrix m);

    void Invert();
    friend Matrix Inverse(Matrix m);

    void ReduceToRREF();
    friend Matrix RREF(Matrix m);

    int Rank() const;
    friend int Rank(Matrix m);
};

bool Matrix::AreComparable(const Matrix &m1, const Matrix &m2) {
    return m1.NRows() == m2.NRows() && m1.NCols() == m2.NCols();
}

bool Matrix::AreMultiplyable(const Matrix &m1, const Matrix &m2) {
    return m1.NCols() == m2.NRows();
}

bool Matrix::IsSquare() const {
    return NRows() == NCols();
}

void Matrix::CheckComparability(const Matrix &m1, const Matrix &m2) {
    if (!AreComparable(m1, m2)) {
        throw std::domain_error("Incompatible formats");
    }
}

void Matrix::CheckMultiplicativity(const Matrix &m1, const Matrix &m2) {
    if (!AreMultiplyable(m1, m2)) {
        throw std::domain_error("Incompatible formats");
    }
}

void Matrix::CheckSquareness() const {
    if (!IsSquare()) {
        throw std::domain_error("Matrix is not square");
    }
}

Matrix::Matrix(int m, int n) {
    if (m <= 0 || n <= 0) {
        throw std::range_error("Bad dimension");
    }

    elementi.resize(m, std::vector<double>(n, 0));
}

Matrix::Matrix(const Vector &v) : Matrix(v.NElems(), 1) {
    for (int i = 0; i < v.NElems(); ++i) {
        elementi[i][0] = v[i];
    }
}

Matrix::Matrix(std::initializer_list<std::vector<double>> l) : Matrix(l.size(), l.begin()->size()) {
    for (auto x = l.begin(); x != l.end(); x++) {
        if (x->size() != l.begin()->size()) {
            throw std::logic_error("Bad matrix");
        }
    }

    elementi = l;
}

int Matrix::NRows() const {
    return elementi.size();
}

int Matrix::NCols() const {
    return elementi.at(0).size();
}

double* Matrix::operator[](int i) {
    return &(elementi[i][0]);
}

const double* Matrix::operator[](int i) const {
    return &(elementi[i][0]);
}

double& Matrix::operator()(int i, int j) {
    try {
        return elementi.at(i - 1).at(j - 1);
    } catch (...) {
        throw std::range_error("Invalid index");
    }
}

double Matrix::operator()(int i, int j) const {
    try {
        return elementi.at(i - 1).at(j - 1);
    } catch (...) {
        throw std::range_error("Invalid index");
    }
}

double Matrix::Norm() const {
    long double sum = 0;

    for (int i = 0; i < NRows(); ++i) {
        for (int j = 0; j < NCols(); ++j) {
            sum += (long double) elementi[i][j] * elementi[i][j];
        }
    }

    return std::sqrt(sum);
}

double MatrixNorm(const Matrix &m) {
    return m.Norm();
}

double Matrix::GetEpsilon() const {
    return 10 * Norm() * std::numeric_limits<double>::epsilon();
}

void Matrix::Print(int width = 10, double eps = -1) const {
    double prag = eps < 0 ? GetEpsilon() : eps;

    for (int i = 0; i < NRows(); ++i) {
        for (int j = 0; j < NCols(); ++j) {
            const double &element = std::abs(elementi[i][j]) < prag ? 0 : elementi[i][j];
            std::cout << std::setw(element >= 0 ? width : width + 1) << element;
        }

        std::cout << std::endl;
    }
}

void PrintMatrix(const Matrix &m, int width = 10, double eps = -1) {
    m.Print(width, eps);
}

Matrix operator+(const Matrix &m1, const Matrix &m2) {
    Matrix m3(m1);

    return m3 += m2;
}

Matrix& Matrix::operator+=(const Matrix &m) {
    Matrix::CheckComparability(*this, m);

    for (int i = 0; i < NRows(); ++i) {
        for (int j = 0; j < NCols(); ++j) {
            elementi[i][j] += m.elementi[i][j];
        }
    }

    return *this;
}

Matrix operator-(const Matrix &m1, const Matrix &m2) {
    Matrix m3(m1);

    return m3 -= m2;
}

Matrix& Matrix::operator-=(const Matrix &m) {
    return *this += (-1) * m;
}

Matrix operator*(double s, const Matrix &m) {
    Matrix rez(m);

    return rez *= s;
}

Matrix operator*(const Matrix &m, double s) {
    return s * m;
}

Matrix& Matrix::operator*=(double s) {
    for (int i = 0; i < NRows(); ++i) {
        for (int j = 0; j < NCols(); ++j) {
            elementi[i][j] *= s;
        }
    }

    return *this;
}

Matrix operator*(const Matrix &m1, const Matrix &m2) {
    Matrix::CheckMultiplicativity(m1, m2);

    Matrix rez(m1.NRows(), m2.NCols());

    for (int i = 0; i < rez.NRows(); ++i) {
        for (int j = 0; j < rez.NCols(); ++j) {
            long double s = 0;
            
            for (int k = 0; k < m1.NCols(); ++k) {
                s += (long double) m1.elementi[i][k] * m2.elementi[k][j];
            }

            rez.elementi[i][j] = s;
        }
    }

    return rez;
}

Matrix& Matrix::operator*=(const Matrix &m) {
    return *this = *this * m;
}

Vector operator*(const Matrix &m, const Vector &v) {
    if (m.NCols() != v.NElems()) {
        throw std::domain_error("Incompatible formats");
    }

    Vector rez(m.NRows());

    for (int i = 0; i < m.NRows(); ++i) {
        long double s = 0;
        
        for (int j = 0; j < m.NCols(); ++j) {
            s += (long double) m.elementi[i][j] * v[j];
        }

        rez[i] = s;
    }

    return rez;
}

Matrix Transpose(const Matrix &m) {
    Matrix rez(m);

    rez.Transpose();

    return rez;
}

void Matrix::Transpose() {
    if (NCols() == NRows()) {
        for (int i = 0; i < NRows(); ++i) {
            for (int j = i; j < NCols(); ++j) {
                std::swap(elementi[i][j], elementi[j][i]);
            }
        }

        return;
    }

    Matrix rez(NCols(), NRows());

    for (int i = 0; i < NRows(); ++i) {
        for (int j = 0; j < NCols(); ++j) {
            rez.elementi[j][i] = elementi[i][j];
        }
    }

    *this = rez;
}

void Matrix::Chop(double eps) {
    double epsilon = eps < 0 ? GetEpsilon() : eps;

    for (auto vek : elementi) {
        for (double &x : vek) {
            if (std::abs(x) < epsilon) {
                x = 0;
            }
        }
    }
}

bool Matrix::EqualTo(const Matrix &m, double eps) const {
    if (this == &m) {
        return true;
    }

    if (NRows() != m.NRows() || NCols() != m.NCols()) {
        return false;
    }

    double epsilon = eps < 0 ? GetEpsilon() : eps;

    for (int i = 0; i < NRows(); ++i) {
        for (int j = 0; j < NCols(); ++j) {
            if (std::abs(elementi[i][j] - m.elementi[i][j]) > epsilon) {
                return false;
            }
        }
    }

    return true;
}

Matrix LeftDiv(Matrix m1, Matrix m2) {
    if (m1.NRows() != m1.NCols()) {
        throw std::domain_error("Divisor matrix is not square");
    }
    
    Matrix::CheckMultiplicativity(m1, m2);

    std::vector<std::vector<double>> &a = m1.elementi, &b = m2.elementi;
    const int n = m1.NRows(), m = m2.NCols();
    const double epsilon = m1.GetEpsilon();
    
    for (int k = 0; k < n; ++k) {
        int p = m1.FindPivot(k, k);

        if (std::abs(a[p][k]) < epsilon) {
            throw std::domain_error("Divisor matrix is singular");
        }

        m1.SwapRows(k, p);
        m2.SwapRows(k, p);

        for (int i = k + 1; i < n; ++i) {
            long double mi = (long double) a[i][k] / a[k][k];

            for (int j = k + 1; j < n; ++j) {
                a[i][j] -= mi * a[k][j];
            }

            for (int j = 0; j < m; ++j) {
                b[i][j] -= mi * b[k][j];
            }
        }
    }

    for (int k = 0; k < m; ++k) {
        for (int i = n - 1; i >= 0; --i) {
            long double s = b[i][k];

            for (int j = i + 1; j < n; ++j) {
                s -= (long double) a[i][j] * b[j][k];
            }

            b[i][k] = s / a[i][i];
        }
    }

    return m2;
}

Vector LeftDiv(Matrix m, Vector v) {
    if (m.NRows() != m.NCols()) {
        throw std::domain_error("Divisor matrix is not square");
    }
    
    Matrix::CheckMultiplicativity(m, v);

    std::vector<std::vector<double>> &a = m.elementi;
    const int n = m.NRows();
    const double epsilon = m.GetEpsilon();

    for (int k = 0; k < n; ++k) {
        int p = m.FindPivot(k, k);

        if (std::abs(a[p][k]) < epsilon) {
            throw std::domain_error("Divisor matrix is singular");
        }

        if (p != k) {
            std::swap(a[k], a[p]);
            std::swap(v[k], v[p]);
        }

        for (int i = k + 1; i < n; ++i) {
            long double mi = a[i][k] / a[k][k];

            for (int j = k + 1; j < n; ++j) {
                a[i][j] -= mi * a[k][j];
            }

            v[i] -= mi * v[k];
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        long double s = v[i];

        for (int j = i + 1; j < n; ++j) {
            s -= (long double) a[i][j] * v[j];
        }

        v[i] = s / a[i][i];
    }

    return v;
}

Matrix operator/(const Matrix &m, double s) {
    Matrix rez = m;

    return rez /= s;
}

Matrix &Matrix::operator/=(double s) {
    if (std::abs(s) < GetEpsilon()) {
        throw std::domain_error("Division by zero");
    }

    for (auto &vek : elementi) {
        for (double &x : vek) {
            x /= s;
        }
    }

    return *this;
}

Matrix operator/(Matrix m1, Matrix m2) {
    return m1 /= m2;
}

Matrix &Matrix::operator/=(Matrix m) {
    if (m.NRows() != m.NCols()) {
        throw std::domain_error("Divisor matrix is not square");
    }
    
    if (m.NRows() != NCols()) {
        throw std::domain_error("Incompatible formats");
    }

    std::vector<std::vector<double>> &a = m.elementi, &b = elementi;
    const int n_Row = NCols(), m_Col = NRows();
    const double epsilon = m.GetEpsilon();
    
    for (int k = 0; k < n_Row; ++k) {
        int p = k;

        for (int i = k + 1; i < n_Row; ++i) {
            if (std::abs(a[k][i]) > std::abs(a[k][p])) {
                p = i;
            }
        }

        if (std::abs(a[k][p]) < epsilon) {
            throw std::domain_error("Divisor matrix is singular");
        }

        if (p != k) {
            for (int i = 0; i < n_Row; ++i) {
                std::swap(a[i][k], a[i][p]);
            }

            for (int i = 0; i < m_Col; ++i) {
                std::swap(b[i][k], b[i][p]);
            }
        }

        for (int i = k + 1; i < n_Row; ++i) {
            long double mi = (long double) a[k][i] / a[k][k];

            for (int j = k + 1; j < n_Row; ++j) {
                a[j][i] -= mi * a[j][k];
            }

            for (int j = 0; j < m_Col; ++j) {
                b[j][i] -= mi * b[j][k];
            }
        }
    }

    for (int k = 0; k < m_Col; ++k) {
        for (int i = n_Row - 1; i >= 0; --i) {
            long double s = b[k][i];

            for (int j = i + 1; j < n_Row; ++j) {
                s -= (long double) a[j][i] * b[k][j];
            }

            b[k][i] = s / a[i][i];
        }
    }

    return *this;
}

double Det(Matrix m) {
    m.CheckSquareness();

    std::vector<std::vector<double>> &a = m.elementi;

    const int n = m.NRows();
    const double epsilon = m.GetEpsilon();
    
    long double d = 1;

    for (int k = 0; k < n; ++k) {
        int p = m.FindPivot(k, k);

        if (std::abs(a[p][k]) < epsilon) {
            return 0;
        }

        if (p != k) {
            std::swap(a[p], a[k]);
            d = -d;
        }

        d *= a[k][k];

        for (int i = k + 1; i < n; ++i) {
            long double mi = (long double) a[i][k] / a[k][k];

            for (int j = k + 1; j < n; ++j) {
                a[i][j] -= mi * a[k][j];
            }
        }
    }

    return d;
}

double Matrix::Det() const {
    return ::Det(*this);
}

void Matrix::Invert() {
    CheckSquareness();

    std::vector<std::vector<double>> &a = elementi;
    std::vector<int> w(NRows());

    const int n = NRows();
    const double epsilon = GetEpsilon();
    int p;

    for (int k = 0; k < n; ++k) {
        p = FindPivot(k, k);

        if (std::abs(a[p][k]) < epsilon) {
            throw std::domain_error("Matrix is singular");
        }

        SwapRows(k, p);

        w[k] = p;

        long double mi = a[k][k];
        a[k][k] = 1;

        for (int j = 0; j < n; ++j) {
            a[k][j] /= mi;
        }

        for (int i = 0; i < n; ++i) {
            if (i != k) {
                mi = a[i][k];
                a[i][k] = 0;

                for (int j = 0; j < n; ++j) {
                    a[i][j] -= mi * a[k][j];
                }
            }
        }
    }

    for (int j = n - 1; j >= 0; --j) {
        p = w[j];

        if (p != j) {
            for (int i = 0; i < n; ++i) {
                std::swap(a[i][p], a[i][j]);
            }
        }
    }
}

Matrix Inverse(Matrix m) {
    m.Invert();
    return m;
}

void Matrix::ReduceToRREF() {
    std::vector<std::vector<double>> &a = elementi;

    const int n = NRows(), m = NCols();
    int current_row = 0;
    double epsilon = GetEpsilon();

    for (int k = 0; k < m; ++k) {
        if (current_row == n) {
            break;
        }

        int p = FindPivot(current_row, k);

        if (std::abs(a[p][k]) < epsilon) {
            continue;
        }

        SwapRows(current_row, p);

        long double mi = a[current_row][k];

        for (int j = k; j < m; ++j) {
            a[current_row][j] /= mi;
        }

        for (int i = 0; i < n; ++i) {
            if (i == current_row) {
                continue;
            }
            
            if (std::abs(a[i][k]) < epsilon) {
                continue;
            }

            mi = a[i][k];
            
            for (int j = k; j < m; ++j) {
                a[i][j] -= mi * a[current_row][j];
            }
        }

        ++current_row;
    }
}

Matrix RREF(Matrix m) {
    m.ReduceToRREF();
    return m;
}

int Rank(Matrix m) {
    int rank = m.NRows();
    double epsilon = m.GetEpsilon();
    
    m.ReduceToRREF();

    for (int i = m.NRows() - 1; i >= 0; --i) {
        if (std::abs(m.elementi[i][m.NCols() - 1]) > epsilon) {
            break;
        }

        --rank;
    }

    return rank;
}

int Matrix::Rank() const {
    return ::Rank(*this);
}

int Matrix::FindPivot(int from_row, int col) const {
    int p = from_row;
    
    for (int i = from_row + 1; i < NRows(); ++i) {
        if (std::abs(elementi[i][col]) > std::abs(elementi[p][col])) {
            p = i;
        }
    }

    return p;
}

void Matrix::SwapRows(int current_row, int new_row) {
    if (current_row != new_row) {
        std::swap(elementi[current_row], elementi[new_row]);
    }
}

Matrix::Matrix(const Matrix &m) {
    *this = m;
}

Matrix::Matrix(Matrix &&m) {
    *this = std::move(m);
}

Matrix& Matrix::operator=(const Matrix &m) {
    elementi = m.elementi;

    return *this;
}

Matrix& Matrix::operator=(Matrix &&m) {
    elementi = std::move(m.elementi);

    return *this;
}