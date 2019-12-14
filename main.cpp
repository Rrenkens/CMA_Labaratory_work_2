#include <iostream>
#include "tasks/matrix.h"
#include "tasks/qr_algorithm.h"
#include <iomanip>


int main() {
  std::cout << std::fixed << std::setprecision(9);
  Matrix matrix = CreateMatrix(3, 3);
  std::cin >> matrix;
  QRAlgorithm(matrix);
  return 0;
}
