#pragma once
#include <functional>
#include "html_generator.hpp"
#include "distributions.hpp"


// --- Data computation ---

struct CDFComparisonData {
  std::vector<std::string>         labels;
  std::vector<std::vector<double>> datasets; // [empirical_tail, model_tail]
  std::string                      model_name;
};

// Compute empirical and theoretical CDF tails (no logs, just raw values)
inline CDFComparisonData compute_cdf_fit_data(
  const std::vector<float>&                              data,
  const std::vector<float>&                              theta,
  const std::string&                                     model_name,
  std::function<float(float, const std::vector<float>&)> pdf) {


  std::vector<float> unique_k;
  auto               empirical_cdf   = compute_empirical_cdf(data, unique_k);
  auto               theoretical_cdf = compute_theoretical_cdf(pdf, theta, 1, (int)unique_k.back());

  std::vector<std::string> labels;
  std::vector<double>      empirical_tail, model_tail;

  for (size_t i = 0; i < unique_k.size(); ++i) {
    float k    = unique_k[i];
    float tail = 1.0f - empirical_cdf[i];
    if (tail > 0.0f) {
      labels.push_back(std::to_string((int)k));
      empirical_tail.push_back(tail);
    }
  }

  for (size_t k = 0; k < theoretical_cdf.size(); ++k) {
    float tail = 1.0f - theoretical_cdf[k];
    if (tail > 0.0f)
      model_tail.push_back(tail);
  }

  return CDFComparisonData{
    .labels     = std::move(labels),
    .datasets   = {empirical_tail, model_tail},
    .model_name = model_name};
}

// --- Plotting ---

inline void plot_cdf_fit(
  html_generator&          html,
  const CDFComparisonData& cdf_data,
  const std::string&       id) {
  html_generator::chart_begin chart{
    .id             = id,
    .labels         = cdf_data.labels,
    .datasets       = cdf_data.datasets,
    .dataset_labels = {"Empirical", cdf_data.model_name}};

  html << html_generator::row_begin{}
       << "<h3>CDF Fit Comparison (" << cdf_data.model_name << ") [log-log]</h3>"
       << html_generator::row_end{}
       << chart
       << html_generator::chart_end{};
}
