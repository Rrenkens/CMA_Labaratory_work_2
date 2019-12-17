#ifndef CMA_LABORATORY_WORK_2_TASKS_DANILEVSKIY_DANILEVSKIY_H_
#define CMA_LABORATORY_WORK_2_TASKS_DANILEVSKIY_DANILEVSKIY_H_

#include "../matrix.h"

/*
⠄⠄⠄⣾⣿⠿⠿⠶⠿⢿⣿⣿⣿⣿⣦⣤⣄⢀⡅⢠⣾⣛⡉⠄⠄⠄⠸⢀⣿⠄
⠄⠄⢀⡋⣡⣴⣶⣶⡀⠄⠄⠙⢿⣿⣿⣿⣿⣿⣴⣿⣿⣿⢃⣤⣄⣀⣥⣿⣿⠄
⠄⠄⢸⣇⠻⣿⣿⣿⣧⣀⢀⣠⡌⢻⣿⣿⣿⣿⣿⣿⣿⣿⣿⠿⠿⠿⣿⣿⣿⠄
⠄⢀⢸⣿⣷⣤⣤⣤⣬⣙⣛⢿⣿⣿⣿⣿⣿⣿⡿⣿⣿⡍⠄⠄⢀⣤⣄⠉⠋⣰
⠄⣼⣖⣿⣿⣿⣿⣿⣿⣿⣿⣿⢿⣿⣿⣿⣿⣿⢇⣿⣿⡷⠶⠶⢿⣿⣿⠇⢀⣤
⠘⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣽⣿⣿⣿⡇⣿⣿⣿⣿⣿⣿⣷⣶⣥⣴⣿⡗
⢀⠈⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡟⠄
⢸⣿⣦⣌⣛⣻⣿⣿⣧⠙⠛⠛⡭⠅⠒⠦⠭⣭⡻⣿⣿⣿⣿⣿⣿⣿⣿⡿⠃⠄
⠘⣿⣿⣿⣿⣿⣿⣿⣿⡆⠄⠄⠄⠄⠄⠄⠄⠄⠹⠈⢋⣽⣿⣿⣿⣿⣵⣾⠃⠄
⠄⠘⣿⣿⣿⣿⣿⣿⣿⣿⠄⣴⣿⣶⣄⠄⣴⣶⠄⢀⣾⣿⣿⣿⣿⣿⣿⠃⠄⠄
⠄⠄⠈⠻⣿⣿⣿⣿⣿⣿⡄⢻⣿⣿⣿⠄⣿⣿⡀⣾⣿⣿⣿⣿⣛⠛⠁⠄⠄⠄
⠄⠄⠄⠄⠈⠛⢿⣿⣿⣿⠁⠞⢿⣿⣿⡄⢿⣿⡇⣸⣿⣿⠿⠛⠁⠄⠄⠄⠄⠄
⠄⠄⠄⠄⠄⠄⠄⠉⠻⣿⣿⣾⣦⡙⠻⣷⣾⣿⠃⠿⠋⠁⠄⠄⠄⠄⠄⢀⣠⣴
⣿⣿⣿⣶⣶⣮⣥⣒⠲⢮⣝⡿⣿⣿⡆⣿⡿⠃⠄⠄⠄⠄⠄⠄⠄⣠⣴⣿⣿⣿
*/

void SwapCol(Matrix &matrix, size_t first_col, size_t second_col, size_t pos) {
  for (size_t row = 0; row <= pos; row++) {
    std::swap(matrix[row][first_col], matrix[row][second_col]);
  }
}

void MulMatrix(Matrix &matrix, Vector &right, Vector &left, size_t pos) {
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

Complex Compute(const Vector &polynom, Complex val) {
  Complex ans = 0;
  for (int i = 0; i < polynom.size(); i++) {
    ans += polynom[i] * MyPow(val, static_cast<int>(polynom.size()) - i - 1);
  }
  return ans;
}

std::vector<Complex> DurandKerner(const Vector &polynom) {
  if (polynom.empty()) {
    return {};
  }
  std::vector<Complex> ans(polynom.size() - 1);
  ans[0] = 1;
  for (size_t i = 1; i < ans.size(); ++i) {
    ans[i] = ans[i - 1] * Complex(0.4, 0.9);
  }
  double norm = 1;
  while (norm > EPS) {
    norm = 0;
    std::vector<Complex> new_ans(ans.size());
    for (size_t i = 0; i < new_ans.size(); ++i) {
      Complex diff = Compute(polynom, ans[i]);
      for (size_t j = 0; j < new_ans.size(); ++j) {
        if (j != i) {
          diff /= (ans[i] - ans[j]);
        }
      }
      new_ans[i] = ans[i] - diff;
      norm += abs(diff);
    }
    ans = new_ans;
  }
  return ans;
}

void ChangeSign(Vector &polynom) {
  std::reverse(polynom.begin(), polynom.end());
  double sign = 1;
  if (polynom.front() < 0) {
    sign = -1;
  }
  for (auto &el : polynom) {
    el = sign * el;
  }
}

std::map<Complex, std::vector<std::vector<Complex>>, Compare>
FindEigenVectorsWithFrabenius(const std::vector<std::vector<Complex>> &matrix,
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
      ans[el].push_back(VectorRotationWithCubNorm(eigen_vector));
    }
  }
  return ans;
}

std::map<Complex, std::vector<std::vector<Complex>>, Compare>
FindEigenVectorsWithGauss(const std::vector<std::vector<Complex>> &matrix,
                          const std::vector<Complex> &eigenvalues) {
  std::map<std::complex<double>, std::vector<std::vector<std::complex<double>>>, Compare> ans;
  for (const auto &el : eigenvalues) {
    Complex temp_el = RoundEigenValue(el);
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
        for (size_t cur_row = row + 1; cur_row < gauss_matrix.size(); cur_row++) {
          if (abs(gauss_matrix[cur_row][col]) > abs(gauss_matrix[pos][col])) {
            pos = cur_row;
          }
        }
        if (row != pos) {
          std::swap(gauss_matrix[row], gauss_matrix[pos]);
        }
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
          for (size_t cur_col = col; cur_col < gauss_matrix.size(); cur_col++) {
            gauss_matrix[row + 1][cur_col] -= diff * gauss_matrix[row][cur_col];
          }
        }
        row++;
      }

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

        //Vector rationing
        ans[temp_el].push_back(VectorRotationWithCubNorm(eigen_vector));
      }
    }
  }
}

void Danilevskiy(Matrix &matrix) {
  Matrix copy_of_matrix = matrix;
  Matrix similarity_transformation_matrix = CreateEMatrix(matrix.size());

  std::vector<int> pos_of_frobenius_matrix_ = {static_cast<int>(matrix.size())};

  for (int row = static_cast<int>(copy_of_matrix.size()) - 1; row >= 0; row--) {
    if (row - 1 < 0 || std::fabs(copy_of_matrix[row][row - 1]) < EPS) {
      bool is_need_to_divide_into_frobenius = true;
      for (int col = row - 2; col >= 0; col--) {
        if (std::fabs(copy_of_matrix[row][col]) >= EPS) {
          std::swap(copy_of_matrix[row - 1], copy_of_matrix[col]);
          SwapCol(copy_of_matrix, row - 1, col, row);
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
    MulSimilarity(similarity_transformation_matrix, similarity_transformation, row - 1);
  }

  std::reverse(pos_of_frobenius_matrix_.begin(), pos_of_frobenius_matrix_.end());

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
  std::vector<Complex> eigenvalues = DurandKerner(cur_polynom);
  std::map<Complex, std::vector<std::vector<Complex>>, Compare> ans;

  if (pos_of_frobenius_matrix_.size() == 2) {
    ans = FindEigenVectorsWithFrabenius(MatrixToComplex(similarity_transformation_matrix), eigenvalues);
  } else {
    ans = FindEigenVectorsWithGauss(MatrixToComplex(matrix), eigenvalues);
  }
  for (const auto &el : ans) {
    std::cout << "Eigen value: " << el.first << std::endl << "Eigen vectors for this value:" << std::endl;
    for (const auto &ell : el.second) {
      std::cout << ell << std::endl;
    }
  }
}

#endif //CMA_LABORATORY_WORK_2_TASKS_DANILEVSKIY_DANILEVSKIY_H_
