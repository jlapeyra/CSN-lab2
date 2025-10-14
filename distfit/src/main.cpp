#include "distributions.hpp"
#include <stdlib.h>
#include <iomanip>
#include <fstream>
#include <future>
#include <filesystem>

namespace fs = std::filesystem;

int ITERATIONS = 5000;

std::vector<float> read_distribution(const std::string& str) {
  std::ifstream      ifs(str);
  std::vector<float> degreees;
  float              degree;

  while (ifs >> degree) degreees.push_back(degree);
  return degreees;
}


void print_distribution(const std::vector<float>& distribution) {
  for (float f : distribution) std::cout << f << " ";
  std::cout << std::endl;
}

// --- Print table header ---
void print_header() {
  std::cout << std::setw(15) << "Distribution"
            << std::setw(15) << "Param1"
            << std::setw(15) << "Param2"
            << std::setw(15) << "NegLogL" << "\n";
  std::cout << std::string(60, '-') << "\n";
}


template <typename PDF>
OptimizerResult run_optimizer(const std::string& name, const std::vector<float>& data, PDF&& costFn) {
  RandomOptimizer<decltype(costFn)> optimizer;
  OptimizerArguments                args;
  args.maxIterations = ITERATIONS;
  return optimizer(args, costFn);
}

template <typename PDF>
void fit_distribution_single(const std::string& name, const std::vector<float>& data, PDF&& costFn) {
  RandomOptimizer<decltype(costFn)> optimizer;
  OptimizerArguments                args;

  args.maxIterations = ITERATIONS;

  auto result = optimizer(args, costFn);

  std::cout << std::setw(15) << name;
  if (!result.args.empty()) {
    std::cout << std::setw(15) << result.args[0];
    if (result.args.size() > 1) std::cout << std::setw(15) << result.args[1];
    else
      std::cout << std::setw(15) << "-";
  } else {
    std::cout << std::setw(15) << "-" << std::setw(15) << "-";
  }
  std::cout << std::setw(15) << result.cost << "\n";
}

int run_single(const std::string& path) {

  auto degrees = read_distribution(path);

  std::vector<std::future<void>> jobs;

  jobs.push_back(std::async(std::launch::async, [degrees]() { fit_distribution_single("Geometric", degrees, pdfGeometric(degrees)); }));
  jobs.push_back(std::async(std::launch::async, [degrees]() { fit_distribution_single("Poisson", degrees, pdfPoisson(degrees)); }));
  jobs.push_back(std::async(std::launch::async, [degrees]() { fit_distribution_single("AltmannZeta", degrees, pdfAltmannZeta(degrees)); }));
  jobs.push_back(std::async(std::launch::async, [degrees]() { fit_distribution_single("TruncZeta", degrees, pdfTruncatedZeta(degrees)); }));
  return 0;
}

int run_combined(const std::string& dir) {

  fs::path dirPath(dir);
  if (!fs::exists(dirPath) || !fs::is_directory(dirPath)) {
    std::cerr << "Invalid directory: " << dirPath << "\n";
    return 1;
  }

  std::vector<fs::path> files;
  for (auto& p : fs::directory_iterator(dirPath))
    if (fs::is_regular_file(p)) files.push_back(p.path());

  std::vector<std::vector<OptimizerResult>> results;

  results.resize(files.size());

  int fileId = 0;

  std::vector<std::string> models = {"Geometric", "Poisson", "AltmannZeta", "TruncZeta"};

  std::cerr << "-> Going to fit models: (Geometric, Poission AltmannZeta TruncZeta) for files: " << files.size() << std::endl;
  for (int i = 0; i < files.size(); i++)
    std::cerr << files[i] << std::endl;

  std::cerr << "..." << std::endl;

  std::atomic<int> remaining = files.size();

#pragma omp parallel for
  for (int i = 0; i < files.size(); i++) {
    auto& file = files[i];

    auto degrees = read_distribution(file.string());
    if (degrees.empty()) continue;

    results[i].push_back(run_optimizer("Geometric", degrees, pdfGeometric(degrees)));
    results[i].push_back(run_optimizer("Poisson", degrees, pdfPoisson(degrees)));
    results[i].push_back(run_optimizer("AltmannZeta", degrees, pdfAltmannZeta(degrees)));
    results[i].push_back(run_optimizer("TruncZeta", degrees, pdfTruncatedZeta(degrees)));

    std::cerr << files[i] << " Done" << std::endl;
    std::cerr << (remaining--) << "/" << files.size() << std::endl;
  }

  std::cerr << "Results: " << std::endl;
  // --- Print NegLogL table ---
  std::cout << std::setw(25) << "File";
  for (auto& m : models) std::cout << std::setw(15) << m;
  std::cout << "\n"
            << std::string(25 + models.size() * 15, '-') << "\n";

  for (int i = 0; i < files.size(); ++i) {
    std::cout << std::setw(25) << files[i].filename().string();
    for (int j = 0; j < models.size(); ++j)
      std::cout << std::setw(15) << results[i][j].cost;
    std::cout << "\n";
  }

  // --- Print optimized parameters table ---
  std::cout << "\nOptimized Parameters:\n";
  std::cout << std::setw(25) << "File";
  for (auto& m : models) std::cout << std::setw(30) << m;
  std::cout << "\n"
            << std::string(25 + models.size() * 30, '-') << "\n";

  for (int i = 0; i < files.size(); ++i) {
    std::cout << std::setw(25) << files[i].filename().string();
    for (int j = 0; j < models.size(); ++j) {
      auto& params = results[i][j].args;
      if (!params.empty()) {
        std::cout << std::setw(15) << params[0];
        if (params.size() > 1) std::cout << std::setw(15) << params[1];
        else
          std::cout << std::setw(15) << "-";
      } else {
        std::cout << std::setw(15) << "-" << std::setw(15) << "-";
      }
    }
    std::cout << "\n";
  }

  // --- Model selection using AIC ---
  std::cout << "\nModel Selection (based on AIC):\n";
  std::cout << std::setw(25) << "File"
            << std::setw(20) << "Best Model"
            << std::setw(30) << "Parameters"
            << std::setw(15) << "AIC" << "\n";
  std::cout << std::string(25 + 20 + 30 + 15, '-') << "\n";

  for (int i = 0; i < files.size(); ++i) {
    double bestAIC      = std::numeric_limits<double>::max();
    int    bestModelIdx = -1;

    for (int j = 0; j < models.size(); ++j) {
      int    k   = static_cast<int>(results[i][j].args.size()); // number of parameters
      double aic = 2.0 * k + 2.0 * results[i][j].cost;          // AIC = 2k - 2 lnL, cost = -lnL
      if (aic < bestAIC) {
        bestAIC      = aic;
        bestModelIdx = j;
      }
    }

    // Print best model and parameters
    std::cout << std::setw(25) << files[i].filename().string();
    if (bestModelIdx >= 0) {
      std::cout << std::setw(20) << models[bestModelIdx];

      auto&              params = results[i][bestModelIdx].args;
      std::ostringstream paramStr;
      for (size_t p = 0; p < params.size(); ++p) {
        if (p > 0) paramStr << ", ";
        paramStr << params[p];
      }
      std::cout << std::setw(30) << paramStr.str();
      std::cout << std::setw(15) << bestAIC;
    } else {
      std::cout << std::setw(20) << "-" << std::setw(30) << "-" << std::setw(15) << "-";
    }
    std::cout << "\n";
  }

  return 0;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <-d> [file/directory]\n";
    return 1;
  }
  bool        isDir = false;
  std::string pathArg;

  if (argc == 3 && std::string(argv[1]) == "-d") {
    isDir   = true;
    pathArg = argv[2];
  } else if (argc == 2) {
    pathArg = argv[1];
  } else {
    std::cerr << "Invalid arguments.\n";
    return 1;
  }

  if (isDir)
    return run_combined(pathArg);
  else
    return run_single(pathArg);
}
