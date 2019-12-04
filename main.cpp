#include <iostream>
#include "tasks/matrix.h"
#include "tasks/qr_algorithm.h"
#include <iomanip>

int main() {
  std::cout << std::fixed << std::setprecision(8);
  Matrix matrix = CreateMatrix(4, 4);
  std::cin >> matrix;
  ReductionToHessenberg(matrix);
  return 0;
}
