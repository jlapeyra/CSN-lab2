#pragma once
#include <iostream>
#include <vector>
#include <fstream>

struct ScopedRedirect {
  std::streambuf* old;
  std::ofstream   file;
  ScopedRedirect(const std::string& path) {
    file.open(path);
    old = std::cout.rdbuf(file.rdbuf());
  }
  ~ScopedRedirect() {
    std::cout.rdbuf(old);
  }
};

struct html_generator {
  enum class CellType { None,
                        Data,
                        Header };

  struct html_begin {};
  struct html_end {};
  struct table_begin {};
  struct table_end {};
  struct row_begin {};
  struct row_end {};
  struct header_row_begin {};
  struct image {
    std::string path;
  };

  // Internal state
  CellType current_cell = CellType::None;
  bool     in_row       = false;

  html_generator& operator<<(const table_begin&) {
    std::cout << "<table border='1'>\n";
    return *this;
  }

  html_generator& operator<<(const table_end&) {
    std::cout << "</table>\n";
    return *this;
  }

  html_generator& operator<<(const row_begin&) {
    std::cout << "  <tr>";
    in_row       = true;
    current_cell = CellType::Data;
    return *this;
  }

  html_generator& operator<<(const header_row_begin&) {
    std::cout << "  <tr>";
    in_row       = true;
    current_cell = CellType::Header;
    return *this;
  }

  html_generator& operator<<(const image& image) {
    std::cout << "<img src=\"" << image.path << "\">\n";
    return *this;
  }

  html_generator& operator<<(const row_end&) {
    std::cout << "</tr>\n";
    in_row       = false;
    current_cell = CellType::None;
    return *this;
  }

  html_generator& operator<<(const html_begin&) {
    std::cout << "<html><head>"
              << "<link rel=\"stylesheet\" href=\"nexus.css\">"
              << "<script src=\"https://cdn.jsdelivr.net/npm/chart.js\"></script>"
              << "</head><body>" << std::endl;
    return *this;
  }

  html_generator& operator<<(const html_end&) {
    std::cout << "</body></html>" << std::endl;
    return *this;
  }

  // Generic cell output
  template <typename T>
  html_generator& operator<<(const T& value) {
    if (!in_row) return *this; // Ignore if not in a row

    if (current_cell == CellType::Header)
      std::cout << "<th>" << value << "</th>";
    else if (current_cell == CellType::Data)
      std::cout << "<td>" << value << "</td>";

    return *this;
  }

  // Chart support
  struct chart_begin {
    std::string                      id;
    std::vector<std::string>         labels;
    std::vector<std::vector<double>> datasets;
    std::vector<std::string>         dataset_labels;
  };
  struct chart_end {};

  // Chart generator (multi-line plot)
  html_generator& operator<<(const chart_begin& chart) {
    std::cout << "<canvas id='" << chart.id << "' height='300'></canvas>\n";
    std::cout << "<script>\n"
                 "const ctx = document.getElementById('"
              << chart.id << "').getContext('2d');\n"
              << "new Chart(ctx, {\n"
                 "    type: 'line',\n"
                 "    data: {\n"
                 "        labels: [";

    // X-axis labels
    for (size_t i = 0; i < chart.labels.size(); ++i) {
      std::cout << "\"" << chart.labels[i] << "\"";
      if (i + 1 < chart.labels.size()) std::cout << ",";
    }

    std::cout << "],\n        datasets: [\n";

    // Predefined color palette (orange-based and readable on dark background)
    static const std::vector<std::string> colors = {
      "rgba(255,127,0,1)",   // orange
      "rgba(255,200,0,1)",   // amber
      "rgba(255,255,255,1)", // white
      "rgba(128,200,255,1)", // blue
      "rgba(255,100,100,1)", // red
      "rgba(128,255,128,1)"  // green
    };

    // Each dataset = one optimizer
    for (size_t i = 0; i < chart.datasets.size(); ++i) {
      std::string color = colors[i % colors.size()];
      std::cout << "{ label: '" << chart.dataset_labels[i] << "', data: [";
      for (size_t j = 0; j < chart.datasets[i].size(); ++j) {
        std::cout << chart.datasets[i][j];
        if (j + 1 < chart.datasets[i].size()) std::cout << ",";
      }
      std::cout << "], fill: false, borderColor: '" << color
                << "', borderWidth: 2, tension: 0.25, pointRadius: 3, "
                   "pointHoverRadius: 5, pointBackgroundColor: '"
                << color
                << "', pointBorderColor: '#000000' }";
      if (i + 1 < chart.datasets.size()) std::cout << ",";
      std::cout << "\n";
    }

    // Chart options
    std::cout
      << "]},\n"
         "    options: {\n"
         "        responsive: true,\n"
         "        plugins: {\n"
         "            legend: {\n"
         "                labels: {\n"
         "                    color: '#ff7f00',\n"
         "                    font: { family: 'JetBrains Mono, monospace', size: 12 }\n"
         "                }\n"
         "            },\n"
         "            title: {\n"
         "                display: true,\n"
         "                text: 'Optimizer Cost Convergence',\n"
         "                color: '#ff7f00',\n"
         "                font: { family: 'JetBrains Mono, monospace', size: 16, weight: 'bold' }\n"
         "            }\n"
         "        },\n"
         "        scales: {\n"
         "            x: {\n"
         "                ticks: { color: '#e0e0e0' },\n"
         "                grid: { color: 'rgba(255,127,0,0.15)' }\n"
         "            },\n"
         "            y: {\n"
         "                ticks: { color: '#e0e0e0' },\n"
         "                grid: { color: 'rgba(255,127,0,0.1)' },\n"
         "                beginAtZero: true\n"
         "            }\n"
         "        }\n"
         "    }\n"
         "});\n"
         "</script>\n";
    return *this;
  }

  html_generator& operator<<(const chart_end&) {
    return *this;
  }
};

template <typename T>
std::vector<std::vector<T>> transpose(const std::vector<std::vector<T>>& mat) {
  if (mat.empty()) return {};

  size_t                      rows = mat.size();
  size_t                      cols = mat[0].size();
  std::vector<std::vector<T>> result(cols, std::vector<T>(rows));

  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      result[j][i] = mat[i][j];
    }
  }

  return result;
}
