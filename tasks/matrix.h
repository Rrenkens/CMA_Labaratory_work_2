#ifndef CMA_LABORATORY_WORK_2_TASKS_MATRIX_H_
#define CMA_LABORATORY_WORK_2_TASKS_MATRIX_H_
#include "constants.h"

Matrix CreateMatrix(size_t row, size_t col);
Matrix CreateEMatrix(size_t row);
Matrix operator*(const Matrix &lhs, const Matrix &rhs);
Matrix operator+(const Matrix &lhs, const Matrix &rhs);
Matrix operator-(const Matrix &lhs, const Matrix &rhs);
Matrix operator*=(Matrix &lhs, const Matrix &rhs);
Matrix TransposeMatrix(const Matrix &a);
double EuclideanNormOfVector(const Vector &data);
Vector VectorRotationWithEuclideanNorm(const Vector &data);
std::vector<std::vector<std::complex<double>>> MatrixToComplex(const Matrix &matrix);
std::vector<std::complex<double>> VectorToComplex(const Vector& data);

template<typename T>
std::vector<T> operator*(const std::vector<std::vector<T>> &lhs, const std::vector<T> &rhs) {
  if (lhs.front().size() != rhs.size()) {
    throw std::invalid_argument("These matrices and vector are not compatible");
  }
  std::vector<T> ret(rhs.size(), 0);
  for (size_t row = 0; row < lhs.size(); row++) {
    for (size_t col = 0; col < rhs.size(); col++) {
      ret[row] += lhs[row][col] * rhs[col];
    }
  }
  return ret;
}

template<typename T>
std::vector<T> operator*(const T val, const std::vector<T> &rhs) {
  std::vector<T> ret(rhs.size());
  for (size_t row = 0; row < rhs.size(); row++) {
    ret[row] = val * rhs[row];
  }
  return ret;
}

template<typename T>
std::vector<T> operator+(const std::vector<T> &lhs, const std::vector<T> &rhs) {
  if (lhs.size() != rhs.size()) {
    throw std::invalid_argument("You cannot add two vectors of different orders");
  }
  Vector ret(rhs.size(), 0);
  for (size_t row = 0; row < lhs.size(); row++) {
    ret[row] = lhs[row] + rhs[row];
  }
  return ret;
}

template<typename T>
std::vector<T> operator-(const std::vector<T> &lhs, const std::vector<T> &rhs) {
  if (lhs.size() != rhs.size()) {
    throw std::invalid_argument("You cannot add two vectors of different orders");
  }
  std::vector<T> ret(rhs.size(), 0);
  for (size_t row = 0; row < lhs.size(); row++) {
    ret[row] = lhs[row] - rhs[row];
  }
  return ret;
}

template<typename T>
T CubNormOfVector(const std::vector<T> &data) {
  if (data.empty()) {
    throw std::invalid_argument("Normalization of empty vector");
  }
  T norm = fabs(data.front());
  for (const auto &el : data) {
    norm = std::max(fabs(el), fabs(norm));
  }
  if (fabs(norm) < EPS) {
    throw std::invalid_argument("Normalization of a vector with norm == 0");
  }
  return norm;
}

template<typename T>
std::vector<T> VectorRotationWithCubNorm(const std::vector<T> &data) {
  std::vector<T> temp = data;
  T norm = CubNormOfVector(data);
  for (auto &el : temp) {
    el /= norm;
  }
  return temp;
}

#endif //CMA_LABORATORY_WORK_2_TASKS_MATRIX_H_
