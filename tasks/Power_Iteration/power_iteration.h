#ifndef CMA_LABORATORY_WORK_2_TASKS_POWER_ITERATION_POWER_ITERATION_H_
#define CMA_LABORATORY_WORK_2_TASKS_POWER_ITERATION_POWER_ITERATION_H_

#include "../matrix.h"

struct Data {
  std::vector<Vector> v_i;
  std::vector<Vector> u_i;
  Vector norm_of_v_i;
};

Vector CreateRandomInitialApproach(size_t size) {
  Vector initial_approach(size);
  do {
    for (auto &el : initial_approach) {
      el = GetRandomNum(10);
    }
  } while (fabs(EuclideanNormOfVector(initial_approach)) < EPS);
  return initial_approach;
}

size_t PosOfMaxElement(const Vector &data, double norm) {
  size_t pos = 0;
  for (size_t i = 0; i < data.size(); i++) {
    if (fabs(fabs(data[i]) - norm) < EPS) {
      pos = i;
      break;
    }
  }
  return pos;
}

double CalculationOfSimpleEigenvalue(const Data &data, size_t cur_iteration) {
  return data.norm_of_v_i[cur_iteration % 4]
      * Sign(data.u_i[(cur_iteration + 3) % 4]
             [PosOfMaxElement(data.v_i[cur_iteration % 4],
                              data.norm_of_v_i[cur_iteration % 4])]);
}

double CalculationOfEigenvalueWithDifferentSign(const Data &data,
                                                size_t cur_iteration) {
  return sqrt(data.norm_of_v_i[cur_iteration % 4] *
      data.norm_of_v_i[(cur_iteration + 3) % 4]);
}

Complex CalculationOfComplexEigenvalue(const Data &data, size_t cur_iteration, bool &flag) {
  size_t pos = PosOfMaxElement(data.v_i[(cur_iteration + 2) % 4],
                               data.norm_of_v_i[(cur_iteration + 2) % 4]);
  double r_num = data.v_i[cur_iteration % 4][pos] *
      data.v_i[(cur_iteration + 2) % 4][pos] *
      data.norm_of_v_i[(cur_iteration + 3) % 4] -
      data.norm_of_v_i[(cur_iteration + 2) % 4] *
          data.v_i[(cur_iteration + 3) % 4][pos] *
          data.v_i[(cur_iteration + 3) % 4][pos];

  double r_den = data.v_i[(cur_iteration + 3) % 4][pos] *
      data.u_i[(cur_iteration + 1) % 4][pos] -
      data.u_i[(cur_iteration + 2) % 4][pos] *
          data.u_i[(cur_iteration + 2) % 4][pos] *
          data.norm_of_v_i[(cur_iteration + 2) % 4];

  double r = sqrt(fabs(r_num / r_den));

  double cos = (data.norm_of_v_i[(cur_iteration + 3) % 4] *
      data.v_i[cur_iteration % 4][pos] + r * r *
      data.u_i[(cur_iteration + 2) % 4][pos]) /
      (2 * r * data.v_i[(cur_iteration + 3) % 4][pos]);
  if (fabs(cos) > 1) {
    flag = false;
    return {0, 0};
  }
  return Complex{r * cos, r * sqrt(1 - cos * cos)};
}

bool Check(const Vector &u) {
  for (const auto &el : u) {
    if (fabs(el) >= EPS) {
      return false;
    }
  }
  return true;
}

bool CheckForComplex(std::vector<Complex> &u) {
  for (const auto &el : u) {
    if (abs(el) >= EPS) {
      return false;
    }
  }
  return true;
}

void PowerIteration(Matrix &matrix) {
  std::vector<std::vector<Complex>> complex_matrix = MatrixToComplex(matrix);
  Data data;
  size_t cur_iteration = 0;
  data.u_i.resize(4), data.v_i.resize(4), data.norm_of_v_i.resize(4);

  double lambda_1, lambda_2;
  Complex lambda_31;
  bool complex_flag = true;
  Vector diff_1, u_21, u_22, diff_2;
  std::vector<Complex> v_complex, u_complex, u_31, u_32, diff_3;

  Vector initial_vector;

  while (true) {
    if (cur_iteration == 0) {
      initial_vector = CreateRandomInitialApproach(matrix.size());
      data.v_i[cur_iteration % 4] = initial_vector;
    } else {
      data.v_i[cur_iteration % 4] = matrix * data.u_i[(cur_iteration + 3) % 4];
    }

    data.norm_of_v_i[cur_iteration % 4] = CubNormOfVector(data.v_i[cur_iteration % 4]);

    if (cur_iteration >= 5) {
      lambda_1 = CalculationOfSimpleEigenvalue(data, cur_iteration);
      lambda_2 = CalculationOfEigenvalueWithDifferentSign(data, cur_iteration);
      lambda_31 = CalculationOfComplexEigenvalue(data, cur_iteration, complex_flag);

      diff_1 = matrix * data.u_i[(cur_iteration + 3) % 4] - lambda_1 * data.u_i[(cur_iteration + 3) % 4];

      u_21 = data.v_i[cur_iteration % 4] + lambda_2 * data.u_i[(cur_iteration + 3) % 4];
      u_22 = data.v_i[cur_iteration % 4] - lambda_2 * data.u_i[(cur_iteration + 3) % 4];
      diff_2 = matrix * u_21 - lambda_2 * u_21;

      if (complex_flag) {
        v_complex = VectorToComplex(data.v_i[(cur_iteration + 3) % 4]),
            u_complex = VectorToComplex(data.u_i[(cur_iteration + 2) % 4]);
        Complex lambda_32 = {lambda_31.real(), -lambda_31.imag()};
        u_31 = v_complex - lambda_32 * u_complex;
        u_32 = v_complex - lambda_31 * u_complex;
        diff_3 = complex_matrix * u_31 - lambda_31 * u_31;
      }

      if (Check(diff_1)) {
        std::cout << "Case with simple eigenvalue" << std::endl;
        std::cout << "Eigen value = " << lambda_1 << std::endl <<
                  "Eigen vector for this eigenvalue:" << std::endl <<
                  VectorRotationWithCubNorm(data.u_i[(cur_iteration + 3) % 4]) << std::endl;
        break;
      } else if (Check(diff_2)) {
        std::cout << "Case with eigenvalue with different sign" << std::endl;
        std::cout << "Eigen value = " << lambda_2 << ", " << -lambda_2 << std::endl <<
                  "Eigen vector for this eigenvalue:" << std::endl << VectorRotationWithCubNorm(u_21)
                  << std::endl << VectorRotationWithCubNorm(u_22) << std::endl;
        break;
      } else if (complex_flag && CheckForComplex(diff_3)) {
        std::cout << "Case with complex eigenvalue" << std::endl;
        std::cout << "Eigen value = " << lambda_31 << ", " << lambda_31 << std::endl <<
                  "Eigen vector for this eigenvalue:" << std::endl << VectorRotationWithCubNorm(u_31)
                  << std::endl << VectorRotationWithCubNorm(u_32) << std::endl;
        break;
      }
    }
    complex_flag = true;
    data.u_i[cur_iteration % 4] = VectorRotationWithCubNorm(data.v_i[cur_iteration % 4]);

    if (cur_iteration >= 1e5 || Check(initial_vector - data.u_i[cur_iteration % 4])) {
      cur_iteration = 0;
      continue;
    }
    cur_iteration++;
  }
}

#endif //CMA_LABORATORY_WORK_2_TASKS_POWER_ITERATION_POWER_ITERATION_H_
