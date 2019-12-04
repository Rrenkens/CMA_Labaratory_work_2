#include "constants.h"

void FFT(std::vector<Base> &data, bool invert) {
  int n = data.size();
  if (n == 1) {
    return;
  }
  std::vector<Base> data_0(n / 2), data_1(n / 2);
  for (int i = 0, j = 0; i < n; i += 2, ++j) {
    data_0[j] = data[i];
    data_1[j] = data[i + 1];
  }
  FFT(data_0, invert), FFT(data_1, invert);

  double ang = 2 * PI / n * (invert ? -1 : 1);
  Base w(1), wn(cos(ang), sin(ang));
  for (int i = 0; i < n / 2; ++i) {
    data[i] = data_0[i] + w * data_1[i];
    data[i + n / 2] = data_0[i] - w * data_1[i];
    if (invert) {
      data[i] /= 2, data[i + n / 2] /= 2;
    }
    w *= wn;
  }
}

void Multiply(const std::vector<double> &first_pol,
              std::vector<double> &second_pol, size_t cur_size) {
  std::vector<Base> first_fft(first_pol.begin(), first_pol.end()),
      second_fft(second_pol.begin(), second_pol.end());
  size_t n = 1;
  while (n < std::max(first_pol.size(), second_pol.size())) {
    n <<= 1;
  }
  n <<= 1;
  first_fft.resize(n), second_fft.resize(n);

  FFT(first_fft, false), FFT(second_fft, false);
  for (size_t i = 0; i < n; i++) {
    first_fft[i] *= second_fft[i];
  }
  FFT(first_fft, true);

  second_pol.resize(cur_size);
  for (size_t i = 0; i < cur_size; i++) {
    second_pol[i] = first_fft[i].real();
  }
}