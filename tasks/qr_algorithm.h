#ifndef CMA_LABORATORY_WORK_2__QR_ALGORITHM_H_
#define CMA_LABORATORY_WORK_2__QR_ALGORITHM_H_
#include "constants.h"
#include "matrix.h"

void MulBlockOfMatrix(Matrix &matrix, size_t start_row, size_t end_row,
                      size_t start_col, size_t end_col,
                      const Matrix &mul_matrix, bool left, bool right) {
  Matrix temp_matrix(end_row - start_row, Vector(end_col - start_col));
  for (size_t cur_row = start_row, i = 0; cur_row < end_row; cur_row++, i++) {
    for (size_t cur_col = start_col, j = 0; cur_col < end_col; cur_col++, j++) {
      temp_matrix[i][j] = matrix[cur_row][cur_col];
    }
  }

  if (left) {
    temp_matrix = mul_matrix * temp_matrix;
  }
  if (right) {
    temp_matrix *= mul_matrix;
  }

  for (size_t cur_row = start_row, i = 0; cur_row < end_row; cur_row++, i++) {
    for (size_t cur_col = start_col, j = 0; cur_col < end_col; cur_col++, j++) {
      matrix[cur_row][cur_col] = temp_matrix[i][j];
    }
  }
}

//Reduction matrix to the Hessenberg form by the reflection method
//The total complexity of the Hessenberg reduction = O(n^3)
Matrix ReductionToHessenberg(Matrix &matrix) {
  Matrix Q(matrix.size(), Vector(matrix.size(), 0));
  for (size_t pos = 0; pos < Q.size(); pos++) {
    Q[pos][pos] = 1;
  }
  for (size_t col = 0; col < static_cast<int>(matrix.size()) - 2; col++) {
    Vector w(matrix.size() - 1 - col);

    for (size_t row = col + 1, pos = 0; row < matrix.size(); row++, pos++) {
      w[pos] = matrix[row][col];
    }
    w[0] = w[0] <= 0 ? w[0] - NormOfVector(w) : w[0] + NormOfVector(w);
    w = NormalizeVector(w);

    Matrix H(matrix.size() - col - 1, Vector(matrix.size() - col - 1, 0));

    //H = E - 2ww^{T}
    for (size_t row = 0; row < w.size(); row++) {
      for (size_t cur_col = row; cur_col < w.size(); cur_col++) {
        if (row != cur_col) {
          H[cur_col][row] = H[row][cur_col] = -2 * w[row] * w[cur_col];
        } else {
          H[row][cur_col] = -2 * w[row] * w[cur_col] + 1;
        }
      }
    }

    MulBlockOfMatrix(matrix, 0, col + 1,
                     col + 1, matrix.size(), H, false, true);
    MulBlockOfMatrix(matrix, col + 1, matrix.size(),
                     0, col + 1, H, true, false);
    MulBlockOfMatrix(matrix, col + 1, matrix.size(),
                     col + 1, matrix.size(), H, true, true);
    MulBlockOfMatrix(Q, 0, matrix.size(), col + 1,
                     matrix.size(), H, false, true);
  }

  return Q;
}
#endif //CMA_LABORATORY_WORK_2__QR_ALGORITHM_H_
