#ifndef CMA_LABORATORY_WORK_2_TASKS_POWER_ITERATION_POWER_ITERATION_H_
#define CMA_LABORATORY_WORK_2_TASKS_POWER_ITERATION_POWER_ITERATION_H_

#include "../matrix.h"

struct Data {
  std::vector<Vector> v_i; //The vectors obtained by the power method.
  std::vector<Vector> u_i; //Normalized vectors v_i
  Vector norm_of_v_i; //Norms of vectors v_i
};

//Create random initial approach, each element from range [-10, 10]
Vector CreateRandomInitialApproach(size_t size) {
  Vector initial_approach(size);
  do {
    for (auto &el : initial_approach) {
      el = GetRandomNum(10);
    }
  } while (fabs(EuclideanNormOfVector(initial_approach)) < EPS);
  return initial_approach;
}

//Find position of max element in vector(abs of this element == norm of vector)
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

//Сalculation eigenvalue for case, if simple real
//eigenvalue or eigenvalue real and multiple.
double CalculationOfSimpleEigenvalue(const Data &data, size_t cur_iteration) {
  return data.norm_of_v_i[cur_iteration % 4]
      * Sign(data.u_i[(cur_iteration + 3) % 4]
             [PosOfMaxElement(data.v_i[cur_iteration % 4],
                              data.norm_of_v_i[cur_iteration % 4])]);
}

//Calculation eigenvalue for case, if two largest abs eigenvalues
//are real and opposite in sign.
double CalculationOfEigenvalueWithDifferentSign(const Data &data,
                                                size_t cur_iteration) {
  return sqrt(data.norm_of_v_i[cur_iteration % 4] *
      data.norm_of_v_i[(cur_iteration + 3) % 4]);
}

//Calculation eigenvalue for complex case. If abs(r_den) < EPS or
//abs(cos) > 1 we note that it certainly does not converge to a complex case.
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

  if (fabs(cos) > 1 || fabs(r_den) < EPS) {
    flag = false;
    return {0, 0};
  }

  return Complex{r * cos, r * sqrt(1 - cos * cos)};
}

//Check that all vector coordinates < EPS
bool Check(const Vector &u) {
  for (const auto &el : u) {
    if (fabs(el) >= EPS) {
      return false;
    }
  }
  return true;
}

//Check that all vector coordinates < EPS for complex case
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
  Complex lambda_31, lambda_32;
  bool complex_flag = true;
  Vector diff_1, u_21, u_22, diff_21, diff_22;
  std::vector<Complex> v_complex, u_complex, u_31, u_32, diff_31, diff_32;

  Vector initial_vector;

  while (true) {
    //We skip the first 4 iterations, because 4 vectors are needed to calculate the complex case
    if (cur_iteration == 0) {
      initial_vector = CreateRandomInitialApproach(matrix.size());
      data.v_i[cur_iteration % 4] = initial_vector;
    } else {
      data.v_i[cur_iteration % 4] = matrix * data.u_i[(cur_iteration + 3) % 4];
    }

    data.norm_of_v_i[cur_iteration % 4] = CubNormOfVector(data.v_i[cur_iteration % 4]);

    if (cur_iteration >= 4) {
      lambda_1 = CalculationOfSimpleEigenvalue(data, cur_iteration);
      lambda_2 = CalculationOfEigenvalueWithDifferentSign(data, cur_iteration);
      lambda_31 = CalculationOfComplexEigenvalue(data, cur_iteration, complex_flag);

      //We calculate the eigenvectors for the resulting eigenvalues.
      //Then we find the vector diff = AU - lambda * U.
      //If all diff coordinates < EPS, this means case has converged.
      diff_1 = matrix * data.u_i[(cur_iteration + 3) % 4] - lambda_1 * data.u_i[(cur_iteration + 3) % 4];

      u_21 = data.v_i[cur_iteration % 4] + lambda_2 * data.u_i[(cur_iteration + 3) % 4];
      u_22 = data.v_i[cur_iteration % 4] - lambda_2 * data.u_i[(cur_iteration + 3) % 4];
      diff_21 = matrix * u_21 - lambda_2 * u_21;
      diff_22 = matrix * u_22 + lambda_2 * u_22;

      if (complex_flag) {
        v_complex = VectorToComplex(data.v_i[(cur_iteration + 3) % 4]),
            u_complex = VectorToComplex(data.u_i[(cur_iteration + 2) % 4]);
        lambda_32 = {lambda_31.real(), -lambda_31.imag()};
        u_31 = v_complex - lambda_32 * u_complex;
        u_32 = v_complex - lambda_31 * u_complex;
        diff_31 = complex_matrix * u_31 - lambda_31 * u_31;
        diff_32 = complex_matrix * u_32 - lambda_32 * u_32;
      }

      //Check if any case has converged.
      if (Check(diff_1) &&
          CubNormOfVector(data.u_i[(cur_iteration + 3) % 4]) >= DIF_EPS) {

        std::cout << matrix * data.u_i[(cur_iteration + 3) % 4] << std::endl;
        std::cout << lambda_1 * data.u_i[(cur_iteration + 3) % 4] << std::endl;

        std::cout << "Count of iteration = " << cur_iteration << std::endl;
        std::cout << "Case with simple eigenvalue" << std::endl;
        std::cout << "Eigen value = " << lambda_1 << std::endl <<
                  "Eigen vector for this eigenvalue:" << std::endl <<
                  VectorRotationWithCubNorm(data.u_i[(cur_iteration + 3) % 4]) << std::endl;

        break;
      } else if (Check(diff_21) && Check(diff_22) &&
          CubNormOfVector(u_21) >= DIF_EPS &&
          CubNormOfVector(u_22) >= DIF_EPS) {

        std::cout << "Count of iteration = " << cur_iteration << std::endl;
        std::cout << "Case with eigenvalue with different sign" << std::endl;
        std::cout << "Eigen value = " << lambda_2 << std::endl <<
                  "Eigen vector for this eigenvalue:" <<
                  std::endl << VectorRotationWithCubNorm(u_21) << std::endl;
        std::cout << "Eigen value = " << -lambda_2 << std::endl <<
                  "Eigen vector for this eigenvalue:" << std::endl <<
                  VectorRotationWithCubNorm(u_22) << std::endl;

        break;
      } else if (complex_flag &&
          CheckForComplex(diff_31) &&
          CheckForComplex(diff_32) &&
          abs(CubNormOfVector(u_31)) >= DIF_EPS &&
          abs(CubNormOfVector(u_32)) >= DIF_EPS) {

        std::cout << "Count of iteration = " << cur_iteration << std::endl;
        std::cout << "Case with complex eigenvalue" << std::endl;
        std::cout << "Eigen value = " << lambda_31 << std::endl <<
                  "Eigen vector for this eigenvalue:" << std::endl
                  << VectorRotationWithCubNorm(u_31) << std::endl;
        std::cout << "Eigen value = " << lambda_32 << std::endl <<
                  "Eigen vector for this eigenvalue:" << std::endl
                  << VectorRotationWithCubNorm(u_32) << std::endl;
        break;
      }
    }

    complex_flag = true;
    data.u_i[cur_iteration % 4] = VectorRotationWithCubNorm(data.v_i[cur_iteration % 4]);

    //If the current vector coincides with the initial vector
    //or the number of iterations is too large, we specify a
    //different initial approximation.
    if (cur_iteration >= 1e5 || Check(initial_vector - data.u_i[cur_iteration % 4])) {
      cur_iteration = 0;
      continue;
    }

    cur_iteration++;
  }
}
#endif //CMA_LABORATORY_WORK_2_TASKS_POWER_ITERATION_POWER_ITERATION_H_
