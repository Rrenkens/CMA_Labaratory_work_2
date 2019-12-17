#include <iostream>
#include "tasks/Qr_Algorithm/qr_algorithm.h"
#include "tasks/Power_Iteration/power_iteration.h"
#include "tasks/Danilevskiy/danilevskiy.h"
#include <iomanip>

int main() {
  std::cout << std::fixed << std::setprecision(9);
  Matrix matrix = CreateMatrix(20, 20);

  std::cin >> matrix;

  //PowerIteration(matrix);
  //Danilevskiy(matrix);
  //QRAlgorithm(matrix);
  return 0;
}
