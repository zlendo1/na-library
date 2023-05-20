#pragma once

#include "LinearAlgebraHead.cpp"

class Vector {
private:
    std::vector<double> elementi;

    static void CheckCompatibility(const Vector &v1, const Vector &v2);

public:
    friend Matrix;
    friend LUDecomposer;
    friend QRDecomposer;
    
    explicit Vector(int n);
    Vector(std::initializer_list<double> l);

    Vector(const Vector &v);
    Vector(Vector &&v);

    Vector& operator=(const Vector &v);
    Vector& operator=(Vector &&v);

    int NElems() const;
    double& operator[](int i);
    double operator[](int i) const;
    double& operator()(int i);
    double operator()(int i) const;

    double Norm() const;
    friend double VectorNorm(const Vector &v);
    double GetEpsilon() const;

    void Chop(double eps = -1);
    bool EqualTo(const Vector &v, double eps = -1) const;

    void Print(char separator, double eps) const;
    friend void PrintVector(const Vector &v, char separator, double eps);

    friend Vector operator+(const Vector &v1, const Vector &v2);
    Vector &operator+=(const Vector &v);
    friend Vector operator-(const Vector &v1, const Vector &v2);
    Vector &operator-=(const Vector &v);
    friend Vector operator*(double s, const Vector &v);
    friend Vector operator*(const Vector &v, double s);
    Vector &operator*=(double s);
    friend double operator*(const Vector &v1, const Vector &v2);
    friend Vector operator/(const Vector &v, double s);
    Vector &operator/=(double s);
};

void Vector::CheckCompatibility(const Vector &v1, const Vector &v2) {
    if (v1.NElems() != v2.NElems()) {
        throw std::domain_error("Incompatible formats");
    }
}

Vector::Vector(int n) {
    if (n <= 0) {
        throw std::range_error("Bad dimension");
    }

    elementi.resize(n);
}

Vector::Vector(std::initializer_list<double> l) {
    if (!l.size()) {
        throw std::range_error("Bad dimension");
    }

    elementi = l;
}

int Vector::NElems() const {
    return elementi.size();
}

double& Vector::operator[](int i) {
    return elementi[i];
}

double Vector::operator[](int i) const {
    return elementi[i];
}

double& Vector::operator()(int i) {
    try {
        return elementi.at(i - 1);
    } catch (...) {
        throw std::range_error("Invalid index");
    }
}

double Vector::operator()(int i) const {
    try {
        return elementi.at(i - 1);
    } catch (...) {
        throw std::range_error("Invalid index");
    }
}

double Vector::Norm() const {
    long double suma = 0L;

    for (double x : elementi) {
        suma += (long double) x * x;
    }

    return std::sqrt(suma);
}

double VectorNorm(const Vector &v) {
    return v.Norm();
}

double Vector::GetEpsilon() const {
    return 10 * Norm() * std::numeric_limits<double>::epsilon();
}

void Vector::Print(char separator = '\n', double eps = -1) const {
    double prag = eps >= 0 ? eps : GetEpsilon();

    for (int i = 0; i < elementi.size(); ++i) {
        std::cout << (std::abs(elementi[i]) < prag ? 0 : elementi[i]);

        if (i != elementi.size() - 1 || separator == '\n') {
            std::cout << separator;
        }
    }
}

void PrintVector(const Vector &v, char separator = '\n', double eps = -1) {
    v.Print(separator, eps);
}

Vector operator+(const Vector &v1, const Vector &v2) {
    Vector::CheckCompatibility(v1, v2);

    Vector v3(v1.NElems());

    for (int i = 0; i < v3.NElems(); ++i) {
        v3.elementi[i] = v1.elementi[i] + v2.elementi[i];
    }

    return v3;
}

Vector& Vector::operator+=(const Vector &v) {
    Vector::CheckCompatibility(*this, v);

    for (int i = 0; i < NElems(); ++i) {
        elementi[i] += v.elementi[i];
    }

    return *this;
}

Vector operator-(const Vector &v1, const Vector &v2) {
    return v1 + ((-1) * v2);
}

Vector& Vector::operator-=(const Vector &v) {
    return *this += (-1) * v;
}

Vector operator*(double s, const Vector &v) {
    Vector rez = v;

    for (int i = 0; i < rez.NElems(); ++i) {
        rez.elementi[i] *= s;
    }

    return rez;
}

Vector operator*(const Vector &v, double s) {
    return s * v;
}

Vector& Vector::operator*=(double s) {
    for (int i = 0; i < NElems(); ++i) {
        elementi[i] *= s;
    }

    return *this;
}

double operator*(const Vector &v1, const Vector &v2) {
    Vector::CheckCompatibility(v1, v2);

    long double sum = 0L;

    for (int i = 0; i < v1.NElems(); ++i) {
        sum += (long double) v1.elementi[i] * v2.elementi[i];
    }

    return sum;
}

Vector operator/(const Vector &v, double s) {
    if (std::abs(s) < v.GetEpsilon()) {
        throw std::domain_error("Division by zero");
    }
    
    Vector rez = v;

    for (int i = 0; i < rez.NElems(); ++i) {
        rez.elementi[i] /= s;
    }

    return rez;
}

Vector& Vector::operator/=(double s) {
    if (std::abs(s) < GetEpsilon()) {
        throw std::domain_error("Division by zero");
    }

    for (int i = 0; i < NElems(); ++i) {
        elementi[i] /= s;
    }

    return *this;
}

void Vector::Chop(double eps) {
    double epsilon_value = eps < 0 ? GetEpsilon() : eps;
    
    for (double &x : elementi) {
        if (std::abs(x) < epsilon_value) {
            x = 0;
        }
    }
}

bool Vector::EqualTo(const Vector &v, double eps) const {
    if (this == &v) {
        return true;
    }
    
    if (elementi.size() != v.elementi.size()) {
        return false;
    }

    double epsilon_value = eps < 0 ? GetEpsilon() : eps;

    for (int i = 0; i < elementi.size(); ++i) {
        if (std::abs(elementi[i] - v.elementi[i]) > epsilon_value) {
            return false;
        }
    }

    return true;
}

Vector::Vector(const Vector &v) {
    *this = v;
}

Vector::Vector(Vector &&v) {
    *this = std::move(v);
}

Vector& Vector::operator=(const Vector &v) {
    elementi = v.elementi;

    return *this;
}

Vector& Vector::operator=(Vector &&v) {
    elementi = std::move(v.elementi);

    return *this;
}