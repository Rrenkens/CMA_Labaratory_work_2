#ifndef CMA_LABORATORY_WORK_2__QR_ALGORITHM_H_
#define CMA_LABORATORY_WORK_2__QR_ALGORITHM_H_
#include "../matrix.h"

//Function for multiplying matrix blocks
void MulBlockOfMatrix(Matrix &matrix, size_t start_row, size_t end_row,
                      size_t start_col, size_t end_col,
                      const Matrix &mul_matrix, bool left, bool right) {
  Matrix temp_matrix(end_row - start_row, Vector(end_col - start_col));
  for (size_t cur_row = start_row, i = 0; cur_row < end_row; cur_row++, i++) {
    for (size_t cur_col = start_col, j = 0;
         cur_col < end_col; cur_col++, j++) {
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
    for (size_t cur_col = start_col, j = 0;
         cur_col < end_col; cur_col++, j++) {
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
    w[0] = w[0] <= 0 ? w[0] - EuclideanNormOfVector(w) :
           w[0] + EuclideanNormOfVector(w);
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

    //A^{i + 1} =   A_{i}        U_{i} * H
    //            H * L_{i}    H * B_{i} * H
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
Matrix QrDecompos(Matrix &matrix, Matrix &q, bool is_symmetric) {
  Matrix Q = CreateEMatrix(matrix.size());
  for (int col = 0; col < (int) matrix.size() - 1; col++) {

    if (matrix[col + 1][col] * matrix[col + 1][col] +
        matrix[col][col] * matrix[col][col] < DIF_EPS_QR) {
      continue;

    }
    double cos = matrix[col][col] / (sqrt(matrix[col + 1][col] *
        matrix[col + 1][col] + matrix[col][col] * matrix[col][col]));
    double sin = -matrix[col + 1][col] / (sqrt(matrix[col + 1][col] *
        matrix[col + 1][col] + matrix[col][col] * matrix[col][col]));

    Matrix T_i = CreateMatrix(2, 2);
    T_i[0][0] = T_i[1][1] = cos;
    T_i[0][1] = sin;
    T_i[1][0] = -sin;

    MulBlockOfMatrix(Q, 0, matrix.size(), col, col + 2, T_i, false, true);

    if (is_symmetric) {
      MulBlockOfMatrix(q, 0, matrix.size(), col, col + 2, T_i, false, true);
    }

    std::swap(T_i[0][1], T_i[1][0]);
    MulBlockOfMatrix(matrix, col, col + 2, 0, matrix.size(), T_i, true, false);
  }
  return Q;
}

//Function for finding the eigenvalue of block 2x2
std::vector<Complex> Eigenvalues(Matrix &matrix) {
  std::vector<Complex> eigenvalues;
  for (size_t col = 0; col < matrix.size(); col++) {
    if (col + 1 != matrix.size() && fabs(matrix[col + 1][col]) >= DIF_EPS) {
      double D = (matrix[col][col] + matrix[col + 1][col + 1]) * (matrix[col][col] + matrix[col + 1][col + 1])
          - 4 * (matrix[col][col] * matrix[col + 1][col + 1] - matrix[col][col + 1] * matrix[col + 1][col]);
      if (D > EPS) {
        return {};
      } else if (fabs(D) < EPS) {
        eigenvalues.emplace_back((matrix[col][col] + matrix[col + 1][col + 1]) / 2, 0);
        eigenvalues.emplace_back((matrix[col][col] + matrix[col + 1][col + 1]) / 2, 0);
      } else {
        D = -D;
        eigenvalues.emplace_back((matrix[col][col] + matrix[col + 1][col + 1]) / 2, -sqrt(D) / 2);
        eigenvalues.emplace_back((matrix[col][col] + matrix[col + 1][col + 1]) / 2, sqrt(D) / 2);
      }
      col++;
    } else {
      eigenvalues.emplace_back(matrix[col][col], 0);
    }
  }
  return eigenvalues;
}

//Check to stop QR algorithm
bool CheckEigenValues(std::vector<Complex> &cur, std::vector<Complex> &next) {
  if (next.size() != cur.size() || next.empty() || cur.empty()) {
    return false;
  }
  for (size_t i = 0; i < cur.size(); i++) {
    if (abs(cur[i] - next[i]) >= DIF_EPS_QR) {
      return false;
    }
  }
  return true;
}

//Eigen vector = Q * Hessenberg eigen vector
std::vector<Complex> EigenVector(
    const std::vector<std::vector<Complex>> &matrix,
    const std::vector<Complex> &vector) {
  std::vector<Complex> temp(vector.size(), 0);
  for (size_t row = 0; row < matrix.size(); row++) {
    for (size_t col = 0; col < matrix.size(); col++) {
      temp[row] += matrix[row][col] * vector[col];
    }
  }
  return temp;
}

//Find eigen vector for eigenvalue using Gauss
std::map<Complex, std::vector<std::vector<Complex>>, Compare>
EigenValuesWithEigenVector(const std::vector<std::vector<Complex>> &matrix,
                           const std::vector<Complex> &eigen_values,
                           const std::vector<std::vector<Complex>> &q) {
  std::map<Complex,
           std::vector<std::vector<Complex>>, Compare> ans;
  for (const auto &el : eigen_values) {
    Complex temp_el = RoundEigenValue(el);
    if (ans.find(temp_el) == ans.end()) {
      std::vector<std::vector<Complex>> gauss_matrix = matrix;

      // A - \lambda * E
      for (size_t row = 0; row < gauss_matrix.size(); row++) {
        gauss_matrix[row][row] -= el;
      }
      std::vector<int> independent_val;

      //Forward running of Gauss algorithm
      for (size_t row = 0, col = 0; row < gauss_matrix.size() &&
          col < gauss_matrix.size(); col++) {
        size_t pos = row;

        //Choose main element
        for (size_t cur_row = row + 1;
             cur_row < std::min(col + 2, gauss_matrix.size()); cur_row++) {
          if (abs(gauss_matrix[cur_row][col]) > abs(gauss_matrix[pos][col])) {
            pos = cur_row;
          }
        }

        if (row != pos) {
          std::swap(gauss_matrix[row], gauss_matrix[pos]);
        }

        //Declare an independent variable
        if (abs(gauss_matrix[row][col]) < DIF_EPS) {
          independent_val.push_back(col);
          continue;
        }

        Complex diff = gauss_matrix[row][col];
        for (size_t cur_col = col; cur_col < gauss_matrix.size(); cur_col++) {
          gauss_matrix[row][cur_col] /= diff;
        }
        if (row + 1 != gauss_matrix.size()) {
          diff = gauss_matrix[row + 1][col];
          for (size_t cur_col = col;
               cur_col < gauss_matrix.size(); cur_col++) {
            gauss_matrix[row + 1][cur_col] -=
                diff * gauss_matrix[row][cur_col];
          }
        }
        row++;
      }

      //Reverse running of Gauss algorithm
      for (int row = (int) gauss_matrix.size() - 1; row > 0; row--) {
        if (abs(gauss_matrix[row][row]) >= DIF_EPS) {
          for (int cur_row = row - 1; cur_row >= 0; cur_row--) {
            Complex dif = gauss_matrix[cur_row][row];
            gauss_matrix[cur_row][row] = 0;
            for (const auto &pos : independent_val) {
              gauss_matrix[cur_row][pos] -= dif * gauss_matrix[row][pos];
            }
          }
        }
      }

      //Create eigen vectors
      for (const auto &pos : independent_val) {
        std::vector<Complex> eigen_vector;
        for (size_t row = 0; row < matrix.size(); row++) {
          if (row == pos) {
            eigen_vector.emplace_back(1, 0);
          } else if ((abs(gauss_matrix[row][pos]) >= DIF_EPS &&
              (abs(gauss_matrix[row][row]) >= DIF_EPS))) {
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

//Check for the symmetry of the matrix
bool CheckForSymmetric(const Matrix &matrix) {
  for (size_t row = 0; row < matrix.size(); row++) {
    for (size_t col = row + 1; col < matrix.size(); col++) {
      if (fabs(matrix[row][col] - matrix[col][row]) >= EPS) {
        return false;
      }
    }
  }
  return true;
}

//QR algorithm
void QRAlgorithm(Matrix &matrix) {

  size_t count_of_iterations = 0;

  //We need them to reconstruct eigenvectors
  const std::vector<std::vector<Complex>> copy_matrix =
      MatrixToComplex(matrix);
  bool is_symmetric = CheckForSymmetric(matrix);
  Matrix q = ReductionToHessenberg(matrix);
  std::vector<std::vector<Complex>> hessenberg
      = MatrixToComplex(matrix);

  //new_matrix matrix on the last QR iteration,
  //new_eigen_values eigen values of this matrix
  Matrix new_matrix, cur_matrix;
  new_matrix = cur_matrix = matrix;

  //QR algorithm
  std::vector<Complex> cur_eigen_values, new_eigen_values;
  cur_eigen_values = new_eigen_values = Eigenvalues(matrix);
  do {
    count_of_iterations++;
    cur_matrix = matrix = new_matrix;
    cur_eigen_values = new_eigen_values;
    Matrix Q = QrDecompos(matrix, q, is_symmetric);
    new_matrix = matrix * Q;
    new_eigen_values = Eigenvalues(new_matrix);
  } while (!CheckEigenValues(new_eigen_values, cur_eigen_values));

  std::map<Complex,
           std::vector<std::vector<Complex>>, Compare> ans;

  //Search for eigenvectors.
  //If this matrix was symmetric, then we search
  //eigenvectors using the qr algorithm
  if (!is_symmetric) {
    ans = EigenValuesWithEigenVector(hessenberg,
                                     new_eigen_values, MatrixToComplex(q));
  } else {
    q = TransposeMatrix(q);
    for (size_t i = 0; i < new_eigen_values.size(); i++) {
      ans[RoundEigenValue(new_eigen_values[i])].push_back(
          VectorRotationWithCubNorm(VectorToComplex(q[i])));
    }
  }

  //Print only eigenvalues
//  for (const auto &el : ans) {
//    std::cout << el.first << std::endl;
//  }

//  std::cout << "Count of iterations = "
//            << count_of_iterations << std::endl;

  //Print all eigenvalue with eigenvectors (without attached vectors)
//  for (const auto &eigenvalue : ans) {
//    std::cout << "Eigen value: " << eigenvalue.first << std::endl
//              << "Eigen vectors for this value:" << std::endl;
//    for (const auto &eigenvectors: eigenvalue.second) {
//      std::cout << eigenvectors << std::endl;
////      std::cout << copy_matrix * eigenvectors << std::endl;
////      std::cout << eigenvalue.first * eigenvectors << std::endl;
//    }
//  }
}
#endif //CMA_LABORATORY_WORK_2__QR_ALGORITHM_H_
