#ifndef CMA_LABORATORY_WORK_2_TASKS_DANILEVSKIY_DANILEVSKIY_H_
#define CMA_LABORATORY_WORK_2_TASKS_DANILEVSKIY_DANILEVSKIY_H_

#include "../matrix.h"
#include "newton.h"

void SwapCol(Matrix &matrix, size_t first_col, size_t second_col, size_t pos) {
  for (size_t row = 0; row <= pos; row++) {
    std::swap(matrix[row][first_col], matrix[row][second_col]);
  }
}

void MulMatrix(Matrix &matrix, const Vector &right, const Vector &left, size_t pos) {
  for (size_t row = 0; row < right.size(); row++) {
    Vector temp_vector = matrix[row];
    for (size_t col = 0; col < right.size(); col++) {
      temp_vector[col] = col == pos ? matrix[row][col] * right[col] :
                         matrix[row][col] + matrix[row][pos] * right[col];
    }
    matrix[row] = temp_vector;
  }

  for (size_t col = 0; col < left.size(); col++) {
    double temp = 0;
    for (size_t row = 0; row < left.size(); row++) {
      temp += matrix[row][col] * left[row];
    }
    matrix[pos][col] = temp;
  }
}

void MulSimilarity(Matrix &matrix, Vector &data, size_t pos) {
  for (size_t col = 0; col < matrix.size(); col++) {
    matrix[pos][col] = data[col];
  }

  for (size_t row = pos + 1; row + 1 < matrix.size(); row++) {
    Vector temp_vector = matrix[row];
    for (size_t col = 0; col < matrix.size(); col++) {
      temp_vector[col] = col == pos ? matrix[row][col] * data[col] :
                         matrix[row][pos] * data[col] + matrix[row][col];
    }
    matrix[row] = temp_vector;
  }
}

//Finding the roots of equations using Aberth method
std::vector<Complex> FindAllRoot(Vector p) {
  if (p.size() <= 1) {
    return {};
  }
  int n = static_cast<int>(p.size()) - 1;
  Vector d = Derivative(p);
  std::vector<Complex> ans(n);
  std::vector<Complex> w(n);
  double r = fabs(n * p[0] / 2 * p[1]) + fabs(p[n - 1] / (2 * n * p[n]));
  double theta = 2 * PI / n, c = theta / (n + 1);
  for (size_t i = 0; i < ans.size(); ++i) {
    ans[i] = Complex(r * cos(i * theta + c), r * sin(i * theta + c));
  }
  double norm = 1;
  while (norm > EPS) {
    norm = 0;
    for (size_t i = 0; i < w.size(); ++i) {
      Complex coef = Compute(p, ans[i]) / Compute(d, ans[i]);
      Complex sum = 0;
      for (size_t j = 0; j < ans.size(); ++j) {
        if (j != i) {
          sum += static_cast<Complex>(1.0) / (ans[i] - ans[j]);
        }
      }
      w[i] = coef / (static_cast<Complex>(1) - coef * sum);
    }
    for (size_t i = 0; i < w.size(); ++i) {
      ans[i] -= w[i];
      norm += abs(w[i] / ans[i]);
    }
  }
  return ans;
}

//Do so that the coefficient at x^{n} is 1
void ChangeSign(Vector &polynom) {
  if (polynom.back() > 0) {
    return;
  }
  for (auto &el : polynom) {
    el = -el;
  }
}

//Find eigenvectors of Frobenius matrix
//y = (x^{n - 1}, x^{n - 2}, ... , x, 1)
std::map<Complex, std::vector<std::vector<Complex>>, Compare>
FindEigenVectorsWithFrobenius(const std::vector<std::vector<Complex>> &matrix,
                              const std::vector<Complex> &eigenvalues) {
  std::map<std::complex<double>, std::vector<std::vector<std::complex<double>>>, Compare> ans;
  for (const auto &el : eigenvalues) {
    Complex temp_el = RoundEigenValue(el);
    if (ans.find(temp_el) == ans.end()) {
      std::vector<Complex> frabenius_eigen_vector(matrix.size());
      for (int i = 0; i < matrix.size(); i++) {
        frabenius_eigen_vector[i] = MyPow(el, static_cast<int>(matrix.size()) - 1 - i);
      }
      std::vector<Complex> eigen_vector = matrix * frabenius_eigen_vector;
      ans[temp_el].push_back(VectorRotationWithCubNorm(eigen_vector));
    }
  }
  return ans;
}

//Find eigenvectors using Gauss algorithm
std::map<Complex, std::vector<std::vector<Complex>>, Compare>
FindEigenVectorsWithGauss(const std::vector<std::vector<Complex>> &matrix,
                          const std::vector<Complex> &eigenvalues) {
  std::map<std::complex<double>, std::vector<std::vector<std::complex<double>>>, Compare> ans;
  for (const auto &el : eigenvalues) {
    Complex temp_el = RoundEigenValueD(el);
    if (ans.find(temp_el) == ans.end()) {
      std::vector<std::vector<std::complex<double>>> gauss_matrix = matrix;

      // A - \lambda * E
      for (size_t row = 0; row < gauss_matrix.size(); row++) {
        gauss_matrix[row][row] -= el;
      }
      std::vector<int> independent_val;

      //Forward running of Gauss algorithm
      for (size_t row = 0, col = 0; row < gauss_matrix.size() && col < gauss_matrix.size(); col++) {
        size_t pos = row;
        //Choose main element
        for (size_t cur_row = row + 1; cur_row < gauss_matrix.size(); cur_row++) {
          if (abs(gauss_matrix[cur_row][col]) > abs(gauss_matrix[pos][col])) {
            pos = cur_row;
          }
        }

        if (row != pos) {
          std::swap(gauss_matrix[row], gauss_matrix[pos]);
        }

        //Declare an independent variable
        if (abs(gauss_matrix[row][col]) < DIF_D) {
          independent_val.push_back(col);
          continue;
        }

        Complex diff = gauss_matrix[row][col];
        for (size_t cur_col = col; cur_col < gauss_matrix.size(); cur_col++) {
          gauss_matrix[row][cur_col] /= diff;
        }
        for (size_t cur_row = row + 1; cur_row < gauss_matrix.size(); cur_row++) {
          diff = gauss_matrix[cur_row][col];
          for (size_t cur_col = col; cur_col < gauss_matrix.size(); cur_col++) {
            gauss_matrix[cur_row][cur_col] -= diff * gauss_matrix[row][cur_col];
          }
        }
        row++;
      }

      //Reverse running of Gauss algorithm
      for (int row = (int) gauss_matrix.size() - 1; row > 0; row--) {
        if (abs(gauss_matrix[row][row]) >= DIF_D) {
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
          } else if ((abs(gauss_matrix[row][pos]) >= DIF_D &&
              (abs(gauss_matrix[row][row]) >= DIF_D))) {
            eigen_vector.push_back(-gauss_matrix[row][pos] / gauss_matrix[row][row]);
          } else {
            eigen_vector.emplace_back(0, 0);
          }
        }

        //Vector rationing
        ans[temp_el].push_back(VectorRotationWithCubNorm(eigen_vector));
      }
    }
  }
  return ans;
}

void Danilevskiy(Matrix &matrix) {
  Matrix copy_of_matrix = matrix;
  std::vector<std::vector<Complex>> complex_matrix = MatrixToComplex(matrix);
  Matrix similarity_transformation_matrix = CreateEMatrix(matrix.size());

  std::vector<int> pos_of_frobenius_matrix_ = {static_cast<int>(matrix.size())};

  //We bring to normal form of Frobenius
  for (int row = static_cast<int>(copy_of_matrix.size()) - 1; row >= 0; row--) {
    if (row - 1 < 0 || std::fabs(copy_of_matrix[row][row - 1]) < EPS) {
      bool is_need_to_divide_into_frobenius = true;
      for (int col = row - 2; col >= 0; col--) {
        if (std::fabs(copy_of_matrix[row][col]) >= EPS) {
          std::swap(copy_of_matrix[row - 1], copy_of_matrix[col]);
          SwapCol(copy_of_matrix, row - 1, col, row);
          SwapCol(similarity_transformation_matrix, row - 1, col, matrix.size() - 1);
          is_need_to_divide_into_frobenius = false;
          break;
        }
      }
      if (is_need_to_divide_into_frobenius) {
        pos_of_frobenius_matrix_.push_back(row);
        continue;
      }
    }

    Vector similarity_transformation(pos_of_frobenius_matrix_.back()),
        inverse_similarity_transformation(pos_of_frobenius_matrix_.back());

    for (size_t col = 0; col < pos_of_frobenius_matrix_.back(); col++) {
      similarity_transformation[col] = col == row - 1 ? 1 / copy_of_matrix[row][col] :
                                       -copy_of_matrix[row][col] / copy_of_matrix[row][row - 1];
      inverse_similarity_transformation[col] = copy_of_matrix[row][col];
    }
    MulMatrix(copy_of_matrix, similarity_transformation,
              inverse_similarity_transformation, row - 1);
    MulMatrix(similarity_transformation_matrix, similarity_transformation, {}, row - 1);
  }

  std::reverse(pos_of_frobenius_matrix_.begin(), pos_of_frobenius_matrix_.end());

  //We find the characteristic polynomial of the matrix
  std::vector<double> cur_polynom = {1};
  for (size_t i = 1; i < pos_of_frobenius_matrix_.size(); i++) {
    std::vector<double> temp_polynom(
        copy_of_matrix[pos_of_frobenius_matrix_[i - 1]].begin() + pos_of_frobenius_matrix_[i - 1],
        copy_of_matrix[pos_of_frobenius_matrix_[i - 1]].begin() +
            pos_of_frobenius_matrix_[i]);
    std::reverse(temp_polynom.begin(), temp_polynom.end());
    temp_polynom.push_back(-1);
    Multiply(temp_polynom, cur_polynom, pos_of_frobenius_matrix_[i] + 1);
  }

  ChangeSign(cur_polynom);
  PrintPolynom(cur_polynom);

    std::vector<Complex> eigenvalues = VectorToComplex(FindRoot(cur_polynom));
  //std::vector<Complex> eigenvalues = FindAllRoot(cur_polynom);

  //Print only eigenvalues
  for (const auto &el : eigenvalues) {
    std::cout << el << std::endl;
  }

  std::map<Complex, std::vector<std::vector<Complex>>, Compare> ans;

  //Choose the method of searching for eigenvectors depending on whether
  //the matrix was divided into Frobenius cells.
  if (pos_of_frobenius_matrix_.size() == 2) {
    ans = FindEigenVectorsWithFrobenius(MatrixToComplex(similarity_transformation_matrix), eigenvalues);
  } else {
    ans = FindEigenVectorsWithGauss(MatrixToComplex(matrix), eigenvalues);
  }

  for (const auto &el : ans) {
    std::cout << "Eigen value: " << el.first << std::endl << "Eigen vectors for this value:" << std::endl;
    for (const auto &ell : el.second) {
      std::cout << ell << std::endl;
//      std::cout << complex_matrix * ell << std::endl;
//      std::cout << el.first * ell << std::endl;
    }
  }
}

#endif //CMA_LABORATORY_WORK_2_TASKS_DANILEVSKIY_DANILEVSKIY_H_
