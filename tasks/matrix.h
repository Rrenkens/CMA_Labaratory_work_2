#ifndef CMA_LABORATORY_WORK_2_TASKS_MATRIX_H_
#define CMA_LABORATORY_WORK_2_TASKS_MATRIX_H_
#include "constants.h"

typedef std::vector<std::vector<double>> Matrix;
typedef std::vector<double> Vector;

Matrix CreateMatrix(size_t row, size_t col);
Matrix operator*(const Matrix &lhs, const Matrix &rhs);
Matrix operator+(const Matrix &lhs, const Matrix &rhs);
Matrix operator-(const Matrix &lhs, const Matrix &rhs);
Matrix operator*=(Matrix &lhs, const Matrix &rhs);
Matrix TransposeMatrix(const Matrix &a);
double NormOfVector(const Vector &data);
Vector NormalizeVector(const Vector &data);
double NormOfMatrix(const Matrix &a);

#endif //CMA_LABORATORY_WORK_2_TASKS_MATRIX_H_
