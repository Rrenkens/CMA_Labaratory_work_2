#ifndef CMA_LABORATORY_WORK_2_TASKS_DANILEVSKIY_NEWTON_H_
#define CMA_LABORATORY_WORK_2_TASKS_DANILEVSKIY_NEWTON_H_

#include "../matrix.h"

//Calculating the polynomial value at a given point
template<typename T>
T Compute(const Vector &polynom, T val) {
  T ans = 0;
  for (int i = 0; i < polynom.size(); i++) {
    ans += static_cast<T>(polynom[i]) * MyPow(val, i);
  }
  return ans;
}

//Finds the derivative of a linear function
std::vector<double> Derivative(const std::vector<double> &pol) {
  std::vector<double> derivative(pol.size() - 1);
  for (size_t i = 1; i < pol.size(); i++) {
    derivative[i - 1] = pol[i] * i;
  }
  return derivative;
}

//Check Fourier condition
bool FourierCondition(const std::vector<double> &pol,
                      const std::vector<double> &second_derivative,
                      double val) {
  return Sign(Compute(pol, val)) * Sign(Compute(second_derivative, val)) > 0;
}

//Bisection method
double BisectionMethod(const std::vector<double> &pol,
                       const std::vector<double> &second_derivative,
                       double l, double r) {
  int sign_l = Sign(Compute(pol, l));
  while (fabs(r - l) > EPS) {
    double m = (l + r) / 2;
    int sign_m = Sign(Compute(pol, m));
    if (sign_m != sign_l) {
      r = m;
    } else {
      l = m;
      sign_l = sign_m;
    }
    if (FourierCondition(pol, second_derivative, l)) {
      return l;
    } else if (FourierCondition(pol, second_derivative, r)) {
      return r;
    }
  }
  return l;
}

//Newthon method
double NewthonMethod(const std::vector<double> &pol,
                     std::vector<double> &derivative_pol,
                     double initial_approx) {
  while (true) {
    double cur_approx = initial_approx -
        Compute(pol, initial_approx) / Compute(derivative_pol, initial_approx);
    if (fabs(initial_approx - cur_approx) < EPS) {
      return cur_approx;
    }
    initial_approx = cur_approx;
  }
}

//Finding the roots of equations using Newton and bisection method
std::vector<double> FindRoot(const std::vector<double> &pol) {
  if (pol.size() == 2) {
    return {-pol.front() / pol.back()};
  } else {
    std::vector<double> derivative_pol = Derivative(pol);
    std::vector<double> second_derivative_pol = Derivative(derivative_pol);

    std::vector<double> root_of_derivative = FindRoot(derivative_pol);
    root_of_derivative.push_back(INF),
        root_of_derivative.insert(root_of_derivative.begin(), -INF);
    std::vector<double> root_of_polynom;

    for (size_t i = 0; i + 1 < root_of_derivative.size(); i++) {
      if (Sign(Compute(pol, root_of_derivative[i])) ==
          Sign(Compute(pol, root_of_derivative[i + 1]))) {
        continue;
      }
      double initial_approach = BisectionMethod(pol, second_derivative_pol,
                                                root_of_derivative[i],
                                                root_of_derivative[i + 1]);
      root_of_polynom.push_back(NewthonMethod(pol,
                                              derivative_pol,
                                              initial_approach));
    }

    return root_of_polynom;
  }
}

#endif //CMA_LABORATORY_WORK_2_TASKS_DANILEVSKIY_NEWTON_H_
