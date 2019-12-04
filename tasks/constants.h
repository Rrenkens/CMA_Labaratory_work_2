#ifndef CMA_LABORATORY_WORK_2__CONSTANTS_H_
#define CMA_LABORATORY_WORK_2__CONSTANTS_H_
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <complex>

typedef std::complex<double> Base;

const double EPS = 1E-9;
const double PI = 3.1415926535;

template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &data) {
  for (const auto &el : data) {
    out << el << " ";
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
  for (const auto &el : data) {
    out << el << "\n";
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

void FFT(std::vector<Base> &data, bool invert);
void Multiply(const std::vector<double> &first_pol,
              std::vector<double> &second_pol, size_t cur_size);

#endif //CMA_LABORATORY_WORK_2__CONSTANTS_H_
