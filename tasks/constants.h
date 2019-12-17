#ifndef CMA_LABORATORY_WORK_2__CONSTANTS_H_
#define CMA_LABORATORY_WORK_2__CONSTANTS_H_
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <complex>
#include <random>
#include <map>

typedef std::vector<std::vector<double>> Matrix;
typedef std::vector<double> Vector;
typedef std::complex<double> Complex;

const double DIF_EPS = 1E-7;
const double DIF_EPS_QR = 1E-14;
const double EPS = 1E-9;
const long long ROUND_CONST_LL = 10000000;
const double ROUND_CONST_LD = 10000000.0;
const double PI = 3.1415926535;

struct Compare {
  bool operator()(const Complex &lhs, const Complex &rhs) {
    return abs(lhs) < abs(rhs);
  }
};

template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &data) {
  for (const auto &el : data) {
    out << el << std::endl;
  }
  return out;
}

template<typename T>
std::istream &operator>>(std::istream &in, std::vector<T> &data) {
  for (auto &el : data) {
    in >> el;
  }
  return in;
}

template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<std::vector<T>> &data) {
  for (const auto &row : data) {
    for (const auto &el : row) {
      out << el << " ";
    }
    out << std::endl;
  }
  return out;
}

template<typename T>
std::istream &operator>>(std::istream &in, std::vector<std::vector<T>> &data) {
  for (auto &el : data) {
    std::cin >> el;
  }
  return in;
}

std::ostream &operator<<(std::ostream &out, const std::complex<double> &data);

void FFT(std::vector<std::complex<double>> &data, bool invert);
void Multiply(const std::vector<double> &first_pol,
              std::vector<double> &second_pol, size_t cur_size);

std::complex<double> RoundEigenValue(const std::complex<double> &value);
int GetRandomNum(size_t mod);
int Sign(double value);
void PrintPolynom(const Vector &polynom);
Complex MyPow(Complex val, int n);

#endif //CMA_LABORATORY_WORK_2__CONSTANTS_H_
