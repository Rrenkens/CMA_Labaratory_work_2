#ifndef CMA_LABORATORY_WORK_2_TASKS_MATRIX_H_
#define CMA_LABORATORY_WORK_2_TASKS_MATRIX_H_
#include "constants.h"


typedef std::vector<std::vector<double>> Matrix;
typedef std::vector<double> Vector;

struct Compare {
  bool operator ()(const std::complex<double>& lhs, const std::complex<double>& rhs) {
    return abs(lhs) < abs(rhs);
  }
};


Matrix CreateMatrix(size_t row, size_t col);
Matrix CreateEMatrix(size_t row);
Matrix operator*(const Matrix &lhs, const Matrix &rhs);
Matrix operator+(const Matrix &lhs, const Matrix &rhs);
Matrix operator-(const Matrix &lhs, const Matrix &rhs);
Matrix operator*=(Matrix &lhs, const Matrix &rhs);
Matrix TransposeMatrix(const Matrix &a);
double NormOfVector(const Vector &data);
Vector VectorRotationWithEuclideanNorm(const Vector &data);
double NormOfMatrix(const Matrix &a);
double BinPow(double a, int n);
std::vector<std::vector<std::complex<double>>> MatrixToComplex(const Matrix& matrix);

template <typename T>
std::vector<T> VectorRotationWithCubNorm(std::vector<T>& data) {
  if(data.empty()) {
    throw std::invalid_argument("Normalization of empty vector");
  }
  T norm = data.front();
  for(const auto & el : data) {
    if(abs(el) > abs(norm)) {
      norm = el;
    }
  }
  if (abs(norm) < EPS) {
    throw std::invalid_argument("Normalization of a vector with norm == 0");
  }
  for(auto &el : data) {
    el /= norm;
  }
  return data;
}

#endif //CMA_LABORATORY_WORK_2_TASKS_MATRIX_H_
