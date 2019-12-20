#include <iostream>
#include "tasks/Qr_Algorithm/qr_algorithm.h"
#include "tasks/Power_Iteration/power_iteration.h"
#include "tasks/Danilevskiy/danilevskiy.h"
#include "tasks/Danilevskiy/newton.h"
#include <iomanip>

int main() {
  std::cout << std::fixed << std::setprecision(9);
  size_t size_;
  std::cout << "Input size of matrix" << std::endl;
  std::cin >> size_;

  Matrix matrix = CreateMatrix(size_, size_);

  std::cin >> matrix;

  //PowerIteration(matrix);
  //Danilevskiy(matrix);
  //QRAlgorithm(matrix);
  return 0;
}
