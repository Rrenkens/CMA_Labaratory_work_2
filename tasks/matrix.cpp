#include "matrix.h"

Matrix CreateMatrix(size_t row, size_t col) {
  return Matrix(row, Vector(col));
}

Matrix CreateEMatrix(size_t row) {
  Matrix temp(row, Vector(row, 0));
  for (size_t cur_row = 0; cur_row < row; cur_row++) {
    temp[cur_row][cur_row] = 1;
  }
  return temp;
}

Matrix operator*(const Matrix &lhs, const Matrix &rhs) {
  if (lhs.front().size() != rhs.size()) {
    throw std::invalid_argument("These matrices are not compatible");
  }
  Matrix ret(lhs.size(), Vector(rhs.front().size()));
  for (size_t i = 0; i < lhs.size(); i++) {
    for (size_t j = 0; j < rhs.front().size(); j++) {
      for (size_t k = 0; k < rhs.size(); k++) {
        ret[i][j] += lhs[i][k] * rhs[k][j];
      }
    }
  }
  return ret;
}

Matrix operator*=(Matrix &lhs, const Matrix &rhs) {
  lhs = std::move(lhs * rhs);
  return lhs;
}

Matrix operator+(const Matrix &lhs, const Matrix &rhs) {
  Matrix ret = lhs;
  for (size_t row = 0; row < lhs.size(); row++) {
    for (size_t col = 0; col < rhs.front().size(); col++) {
      ret[row][col] += rhs[row][col];
    }
  }
  return ret;
}

Matrix operator-(const Matrix &lhs, const Matrix &rhs) {
  Matrix ret = lhs;
  for (size_t row = 0; row < lhs.size(); row++) {
    for (size_t col = 0; col < rhs.front().size(); col++) {
      ret[row][col] -= rhs[row][col];
    }
  }
  return ret;
}

Matrix TransposeMatrix(const Matrix &a) {
  Matrix ret(a[0].size(), Vector(a.size()));
  for (size_t row = 0; row < a.size(); row++) {
    for (size_t col = 0; col < a[row].size(); col++) {
      ret[col][row] = a[row][col];
    }
  }
  return ret;
}

double NormOfVector(const Vector &data) {
  double norm = 0;
  for (const auto &el : data) {
    norm += el * el;
  }
  return sqrt(norm);
}

Vector VectorRotationWithEuclideanNorm(const Vector &data) {
  double norm = NormOfVector(data);
  Vector ans = data;
  if (norm < EPS) {
    throw std::invalid_argument("Normalization of a vector with norm == 0");
  }
  for (auto &el : ans) {
    el /= norm;
  }
  return ans;
}

double NormOfMatrix(const Matrix &a) {
  double norm = 0;
  for (const auto &i:a) {
    for (const auto &j:i) {
      norm += j * j;
    }
  }
  return sqrt(norm);
}
double BinPow(double a, int n) {
  if (n == 0) {
    return 1;
  }
  if (n % 2 == 1) {
    return BinPow(a, n - 1) * a;
  } else {
    double b = BinPow(a, n / 2);
    return b * b;
  }
}

std::vector<std::vector<std::complex<double>>>
MatrixToComplex(const Matrix &matrix) {
  std::vector < std::vector < std::complex < double >>
      > ret(matrix.size(), std::vector < std::complex < double >> (matrix.size()));
  for (size_t row = 0; row < matrix.size(); row++) {
    for (size_t col = 0; col < matrix.size(); col++) {
      ret[row][col] = {matrix[row][col], 0};
    }
  }
  return ret;
}

