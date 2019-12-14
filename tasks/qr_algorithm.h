#ifndef CMA_LABORATORY_WORK_2__QR_ALGORITHM_H_
#define CMA_LABORATORY_WORK_2__QR_ALGORITHM_H_
#include "constants.h"
#include "matrix.h"
#include <exception>
#include <map>

//Function for multiplying matrix blocks
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
    w = VectorRotationWithEuclideanNorm(w);

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

//Find QR decomposition of matrix
Matrix QrDecompos(Matrix &matrix) {
  Matrix Q = CreateEMatrix(matrix.size());
  for (int col = 0; col < (int) matrix.size() - 1; col++) {
    if (matrix[col + 1][col] * matrix[col + 1][col] + matrix[col][col] * matrix[col][col] < DIF_EPS) {
      continue;
    }
    double cos = matrix[col][col]
        / (sqrt(matrix[col + 1][col] * matrix[col + 1][col] + matrix[col][col] * matrix[col][col]));
    double sin =
        -matrix[col + 1][col]
            / (sqrt(matrix[col + 1][col] * matrix[col + 1][col] + matrix[col][col] * matrix[col][col]));
    Matrix T_i = CreateMatrix(2, 2);
    T_i[0][0] = T_i[1][1] = cos;
    T_i[0][1] = sin;
    T_i[1][0] = -sin;
    MulBlockOfMatrix(Q, 0, matrix.size(), col, col + 2, T_i, false, true);
    std::swap(T_i[0][1], T_i[1][0]);
    MulBlockOfMatrix(matrix, col, col + 2, 0, matrix.size(), T_i, true, false);
  }
  return Q;
}

//Function for finding the eigenvalue of block 2x2
std::vector<std::complex<double>> Eigenvalues(Matrix &matrix) {
  std::vector<std::complex<double>> eigenvalues;
  for (size_t col = 0; col < matrix.size(); col++) {
    if (col + 1 != matrix.size() && fabs(matrix[col + 1][col]) >= EPS) {
      double D = (matrix[col][col] + matrix[col + 1][col + 1]) * (matrix[col][col] + matrix[col + 1][col + 1])
          - 4 * (matrix[col][col] * matrix[col + 1][col + 1] - matrix[col][col + 1] * matrix[col + 1][col]);
      if (D < 0) {
        D = -D;
        eigenvalues.emplace_back((matrix[col][col] + matrix[col + 1][col + 1]) / 2, -sqrt(D) / 2);
        eigenvalues.emplace_back((matrix[col][col] + matrix[col + 1][col + 1]) / 2, sqrt(D) / 2);
      } else {
        eigenvalues.emplace_back((matrix[col][col] + matrix[col + 1][col + 1] - sqrt(D)) / 2, 0);
        eigenvalues.emplace_back((matrix[col][col] + matrix[col + 1][col + 1] + sqrt(D)) / 2, 0);
      }
      col++;
    } else {
      eigenvalues.emplace_back(matrix[col][col], 0);
    }
  }
  return eigenvalues;
}

//Check to stop QR algorithm
bool CheckEigenValues(std::vector<std::complex<double>> &cur, std::vector<std::complex<double>> &next) {
  if (next.size() != cur.size()) {
    return false;
  }
  for (size_t i = 0; i < cur.size(); i++) {
    if (abs(cur[i] - next[i]) >= DIF_EPS) {
      return false;
    }
  }
  return true;
}

//Eigen vector = Q * Hessenberg eigen vector
std::vector<std::complex<double>> EigenVector(const std::vector<std::vector<std::complex<double>>> &matrix,
                                              const std::vector<std::complex<double>> &vector) {
  std::vector<std::complex<double>> temp(vector.size(), 0);
  for (size_t row = 0; row < matrix.size(); row++) {
    for (size_t col = 0; col < matrix.size(); col++) {
      temp[row] += matrix[row][col] * vector[col];
    }
  }
  return temp;
}

//Find eigen vector for eigenvalue
std::map<std::complex<double>, std::vector<std::vector<std::complex<double>>>, Compare>
EigenValuesWithEigenVector(const std::vector<std::vector<std::complex<double>>> &matrix,
                           const std::vector<std::complex<double>> &eigen_values,
                           const std::vector<std::vector<std::complex<double>>> &q) {
  std::map<std::complex<double>, std::vector<std::vector<std::complex<double>>>, Compare> ans{};
  for (const auto &el : eigen_values) {
    std::complex<double> temp_el = RoundEigenValue(el);
    if (ans.find(temp_el) == ans.end()) {
      std::vector<std::vector<std::complex<double>>> gauss_matrix = matrix;

      // A - \lambda * E
      for (size_t row = 0; row < gauss_matrix.size(); row++) {
        gauss_matrix[row][row] -= el;
      }
      std::vector<int> independent_val;

      //Gauss
      for (size_t row = 0, col = 0; row < gauss_matrix.size() && col < gauss_matrix.size(); col++) {
        size_t pos = row;
        //Choose main element
        for (size_t cur_row = row + 1; cur_row < std::min(col + 2, gauss_matrix.size()); cur_row++) {
          if (abs(gauss_matrix[cur_row][col]) > abs(gauss_matrix[pos][col])) {
            pos = cur_row;
          }
        }
        if (row != pos) {
          std::swap(gauss_matrix[row], gauss_matrix[pos]);
        }
        if (abs(gauss_matrix[row][col]) < EPS) {
          independent_val.push_back(col);
          continue;
        }

        std::complex<double> diff = gauss_matrix[row][col];
        for (size_t cur_col = col; cur_col < gauss_matrix.size(); cur_col++) {
          gauss_matrix[row][cur_col] /= diff;
        }
        if (row + 1 != gauss_matrix.size()) {
          diff = gauss_matrix[row + 1][col];
          for (size_t cur_col = col; cur_col < gauss_matrix.size(); cur_col++) {
            gauss_matrix[row + 1][cur_col] -= diff * gauss_matrix[row][cur_col];
          }
        }
        row++;
      }

      for (int row = (int) gauss_matrix.size() - 1; row > 0; row--) {
        if (abs(gauss_matrix[row][row]) >= EPS) {
          for (int cur_row = row - 1; cur_row >= 0; cur_row--) {
            std::complex<double> dif = gauss_matrix[cur_row][row];
            gauss_matrix[cur_row][row] = 0;
            for (const auto &pos : independent_val) {
              gauss_matrix[cur_row][pos] -= dif * gauss_matrix[row][pos];
            }
          }
        }
      }

      //Create eigen vectors
      for (const auto &pos : independent_val) {
        std::vector<std::complex<double>> eigen_vector;
        for (size_t row = 0; row < matrix.size(); row++) {
          if (row == pos) {
            eigen_vector.emplace_back(1, 0);
          } else if ((abs(gauss_matrix[row][pos]) >= EPS &&
              (abs(gauss_matrix[row][row]) >= EPS))) {
            eigen_vector.push_back(-gauss_matrix[row][pos] / gauss_matrix[row][row]);
          } else {
            eigen_vector.emplace_back(0, 0);
          }
        }
        //Find eigen vector of Required matrix
        eigen_vector = EigenVector(q, eigen_vector);

        //Vector rationing
        ans[temp_el].push_back(VectorRotationWithCubNorm(eigen_vector));
      }
    }
  }
  return ans;
}

void QRAlgorithm(Matrix &matrix) {
  //We need them to reconstruct eigenvectors
  std::vector<std::vector<std::complex<double>>> q = MatrixToComplex(ReductionToHessenberg(matrix));
  std::vector<std::vector<std::complex<double>>> hessenberg = MatrixToComplex(matrix);

  //new_matrix matrix on the last QR iteration,
  //new_eigen_values eigen values of this matrix
  Matrix new_matrix, cur_matrix;
  new_matrix = cur_matrix = matrix;

  //QR algorithm
  std::vector<std::complex<double>> cur_eigen_values, new_eigen_values;
  cur_eigen_values = new_eigen_values = Eigenvalues(matrix);
  do {
    cur_matrix = matrix = new_matrix;
    cur_eigen_values = new_eigen_values;
    Matrix Q = QrDecompos(matrix);
    new_matrix = matrix * Q;
    new_eigen_values = Eigenvalues(new_matrix);
  } while (!CheckEigenValues(new_eigen_values, cur_eigen_values));
  std::map<std::complex<double>, std::vector<std::vector<std::complex<double>>>,
           Compare> ans = EigenValuesWithEigenVector(hessenberg, new_eigen_values, q);
  for (const auto &el : ans) {
    std::cout << "Eigen value: " << el.first << std::endl << "Eigen vectors for this value:" << std::endl;
    for (const auto &ell : el.second) {
      std::cout << ell << std::endl;
    }
  }
}

#endif //CMA_LABORATORY_WORK_2__QR_ALGORITHM_H_
