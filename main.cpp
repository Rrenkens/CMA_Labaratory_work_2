#include <iostream>
#include "tasks/Qr_Algorithm/qr_algorithm.h"
#include "tasks/Power_Iteration/power_iteration.h"
#include <iomanip>

int main() {
  std::cout << std::fixed << std::setprecision(20);
  Matrix matrix = CreateMatrix(10, 10);
  std::cin >> matrix;
  //QRAlgorithm(matrix);
  PowerIteration(matrix);
  return 0;
}
