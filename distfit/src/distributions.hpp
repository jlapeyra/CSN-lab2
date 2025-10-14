#include <cmath>
#include <vector>
#include "optimizer.hpp"

inline float pdf_geometric(float k, const std::vector<float>& theta) {
  float p = theta[0];
  if (p <= 0.0f || p >= 1.0f || k < 1.0f) return 0.0f;
  return std::pow(1.0f - p, k - 1) * p;
}

inline float pdf_poisson(float k, const std::vector<float>& theta) {
  float lambda = theta[0];
  if (lambda <= 0.0f || k < 0.0f) return 0.0f;
  return std::exp(k * std::log(lambda) - lambda - std::lgamma(k + 1));
}

// Helper for Altmann-Zeta
inline float hurwitz_zeta(float s, float q, int max_terms = 100) {
  float sum = 0.0f;
  for (int n = 0; n < max_terms; ++n) sum += std::pow(n + q, -s);
  return sum;
}

inline float pdf_altmann_zeta(float k, float s, float q, float norm) {
  if (k < 0.0f) return 0.0f;
  return std::pow(k + q, -s) / norm;
}

inline float pdf_truncated_zeta(float k, float alpha, int k_max, float norm) {
  if (k < 1.0f || k > k_max) return 0.0f;
  return std::pow(k, -alpha) / norm;
}

template <typename PDF>
inline auto makePDF(const std::string& name, VectorFieldConfiguration cfg, PDF&& pdf, const std::vector<float>& data) {
  auto fn = [data, pdf](const std::vector<float>& theta) {
    float logL = 0.0f;
    for (float x : data) {
      float p = pdf(x, theta);
      if (p <= 1e-20f) p = 1e-20f;
      logL += std::log(p);
    }
    return -logL;
  };
  return FunctionDefinition<decltype(fn)>{cfg, fn, name, 0.0f};
}

inline auto makePDF_AltmannZeta(const std::string& name, VectorFieldConfiguration cfg, const std::vector<float>& data) {
  auto fn = [data](const std::vector<float>& theta) {
    float s = theta[0], q = theta[1];
    if (s <= 1.0f || q < 0.0f) return INFINITY;
    float norm = hurwitz_zeta(s, 1.0f + q, 100);
    float logL = 0.0f;
    for (float x : data) {
      float p = pdf_altmann_zeta(x, s, q, norm);
      if (p <= 1e-20f) p = 1e-20f;
      logL += std::log(p);
    }
    return -logL;
  };
  return FunctionDefinition<decltype(fn)>{cfg, fn, name, 0.0f};
}

inline auto makePDF_TruncatedZeta(const std::string& name, VectorFieldConfiguration cfg, const std::vector<float>& data) {
  auto fn = [data](const std::vector<float>& theta) {
    float alpha = theta[0];
    int   k_max = static_cast<int>(theta[1]);
    if (alpha <= 1.0f) return INFINITY;
    float norm = 0.0f;
    for (int n = 1; n <= k_max; ++n) norm += std::pow(n, -alpha);
    float logL = 0.0f;
    for (float x : data) {
      float p = pdf_truncated_zeta(x, alpha, k_max, norm);
      if (p <= 1e-20f) p = 1e-20f;
      logL += std::log(p);
    }
    return -logL;
  };
  return FunctionDefinition<decltype(fn)>{cfg, fn, name, 0.0f};
}

// Convenience wrappers
inline auto pdfGeometric(const std::vector<float>& data) {
  return makePDF("Geometric", {{{0.0001f, 0.9999f}}}, pdf_geometric, data);
}

inline auto pdfPoisson(const std::vector<float>& data) {
  return makePDF("Poisson", {{{0.0001f, 1000.0f}}}, pdf_poisson, data);
}

inline auto pdfAltmannZeta(const std::vector<float>& data) {
  return makePDF_AltmannZeta("AltmannZeta", {{{1.01f, 10.0f}, {0.0f, 100.0f}}}, data);
}

inline auto pdfTruncatedZeta(const std::vector<float>& data) {
  return makePDF_TruncatedZeta("TruncatedZeta", {{{1.01f, 10.0f}, {1.0f, 1000.0f}}}, data);
}
