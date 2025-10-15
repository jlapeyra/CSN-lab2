#include "distributions.hpp"
#include <stdlib.h>
#include <iomanip>
#include <fstream>
#include <future>
#include <filesystem>
#include "util.hpp"

namespace fs = std::filesystem;

int  ITERATIONS      = 5000;
bool FAST            = false;
bool COMPARE         = false;
bool USE_ALTERNATIVE = false;

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
  OptimizerArguments args;
  args.maxIterations = ITERATIONS;

  if (FAST) {
    // Only run RandomStallOptimizer for speed
    RandomStallOptimizer<decltype(costFn)> stallOpt;
    OptimizerResult                        resStall = stallOpt(args, costFn);

    if (COMPARE) {
      std::cout << "--------------------------------------------\n";
      std::cout << "Model: " << name << "\n";
      std::cout << "RandomStallOptimizer: cost = " << resStall.cost
                << ", iterations = " << resStall.iterations << "\n";
    }
    return resStall;
  } else {
    // Run both optimizers
    RandomOptimizer<decltype(costFn)> randomOpt;
    OptimizerResult                   resRandom = randomOpt(args, costFn);

    RandomStallOptimizer<decltype(costFn)> stallOpt;
    OptimizerResult                        resStall = stallOpt(args, costFn);

    if (COMPARE) {
      float deltaCost = resRandom.cost - resStall.cost;
      int   deltaIter = resRandom.iterations - resStall.iterations;
      float metric    = static_cast<float>(deltaIter) / (deltaCost + 1e-6f);

      std::cout << "--------------------------------------------\n";
      std::cout << "Model: " << name << "\n";
      std::cout << "RandomOptimizer:      cost = " << resRandom.cost
                << ", iterations = " << resRandom.iterations << "\n";
      std::cout << "RandomStallOptimizer: cost = " << resStall.cost
                << ", iterations = " << resStall.iterations << "\n";
      std::cout << "Delta cost = " << deltaCost
                << ", Delta iterations = " << deltaIter << "\n";
      std::cout << "Metric (iters saved per cost unit): " << metric << "\n";
    }

    // Return the better result
    return (resStall.cost < resRandom.cost) ? resStall : resRandom;
  }
}

template <typename PDF>
void fit_distribution_single(const std::string& name, const std::vector<float>& data, const PDF& costFn) {
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

  optimizeReportHTML(optimizer, costFn, name + ".html");
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

struct ModelFitter {
  fs::path                                  dirPath;
  std::vector<fs::path>                     files;
  std::vector<std::string>                  models = {"Geometric", "Poisson"}; //, "AltmannZeta", "TruncZeta"};
  std::vector<std::vector<OptimizerResult>> results;
  std::vector<std::vector<double>>          resultsAIC;

  ModelFitter(const std::string& dir) :
    dirPath(dir) {
    if (!fs::exists(dirPath) || !fs::is_directory(dirPath)) {
      throw std::runtime_error("Invalid directory: " + dirPath.string());
    }

    for (auto& p : fs::directory_iterator(dirPath))
      if (fs::is_regular_file(p)) files.push_back(p.path());

    results.resize(files.size());

    if (USE_ALTERNATIVE) {
      models.push_back("AltmannZeta");
      models.push_back("TruncZeta");
    }
  }
  void runOptimizers() {
    std::atomic<int> remaining = static_cast<int>(files.size());

    std::cerr << "Starting optimizer runs for " << files.size() << " files...\n";

#pragma omp parallel for
    for (int i = 0; i < files.size(); ++i) {
      auto& file    = files[i];
      auto  degrees = read_distribution(file.string());
      if (degrees.empty()) {
        std::cerr << file << " skipped (empty distribution)\n";
        continue;
      }

      results[i].push_back(run_optimizer("Geometric", degrees, pdfGeometric(degrees)));
      results[i].push_back(run_optimizer("Poisson", degrees, pdfPoisson(degrees)));
      if (USE_ALTERNATIVE) {
        results[i].push_back(run_optimizer("AltmannZeta", degrees, pdfAltmannZeta(degrees)));
        results[i].push_back(run_optimizer("TruncZeta", degrees, pdfTruncatedZeta(degrees)));
      }

      std::cerr << file.filename() << " done (" << --remaining << "/" << files.size() << " remaining)\n";
    }

    std::cerr << "All optimizers completed.\n";
  }

  void computeAIC() {
    resultsAIC = mapf(results, [](const OptimizerResult& result) {
      int k = static_cast<int>(result.args.size());
      return 2.0 * k + 2.0 * result.cost; // AIC = 2k + 2*NegLogL
    });
  }

  void printNegLogL() const {
    const int fileWidth = 25;
    const int colWidth  = 15;

    std::cout << std::setw(fileWidth) << "File";
    for (auto& m : models) std::cout << std::setw(colWidth) << m;
    std::cout << "\n"
              << std::string(fileWidth + models.size() * colWidth, '-') << "\n";

    for (size_t i = 0; i < files.size(); ++i) {
      std::cout << std::setw(fileWidth) << files[i].filename().string();
      for (size_t j = 0; j < models.size(); ++j)
        std::cout << std::setw(colWidth) << results[i][j].cost;
      std::cout << "\n";
    }

    std::cout << "\n";
  }

  void printParameters() const {
    const int fileWidth  = 25;
    const int paramWidth = 30;

    std::cout << "\nOptimized Parameters:\n";
    std::cout << std::setw(fileWidth) << "File";
    for (auto& m : models) std::cout << std::setw(paramWidth) << m;
    std::cout << "\n"
              << std::string(fileWidth + models.size() * paramWidth, '-') << "\n";

    for (size_t i = 0; i < files.size(); ++i) {
      std::cout << std::setw(fileWidth) << files[i].filename().string();
      for (size_t j = 0; j < models.size(); ++j) {
        const auto&        params = results[i][j].args;
        std::ostringstream paramStr;

        if (!params.empty()) {
          for (size_t p = 0; p < params.size(); ++p) {
            if (p > 0) paramStr << ", ";
            paramStr << params[p];
          }
          std::cout << std::setw(paramWidth) << paramStr.str();
        } else {
          std::cout << std::setw(paramWidth) << "-";
        }
      }
      std::cout << "\n";
    }

    std::cout << "\n";
  }


  void printAICTable() const {
    std::cout << "\nAIC Table (ΔAIC relative to best model per file):\n";
    std::cout << std::setw(25) << "File";
    for (auto& m : models) std::cout << std::setw(15) << m;
    std::cout << "\n"
              << std::string(25 + models.size() * 15, '-') << "\n";

    for (int i = 0; i < files.size(); ++i) {
      double minAIC = *std::min_element(resultsAIC[i].begin(), resultsAIC[i].end());
      std::cout << std::setw(25) << files[i].filename().string();

      for (double aic : resultsAIC[i]) {
        double delta = aic - minAIC;
        if (std::abs(delta) < 1e-8)
          std::cout << std::setw(15) << ("*" + std::to_string(delta));
        else
          std::cout << std::setw(15) << delta;
      }
      std::cout << "\n";
    }

    std::cout << "\n";
  }

  void printAkaikeWeights() const {
    std::cout << "\nAkaike Weights (probability each model is best):\n";
    std::cout << std::setw(30) << "File";
    for (auto& m : models) std::cout << std::setw(15) << m;
    std::cout << "\n"
              << std::string(30 + models.size() * 15, '-') << "\n";

    for (int i = 0; i < files.size(); ++i) {
      const auto& aics   = resultsAIC[i];
      double      minAIC = *std::min_element(aics.begin(), aics.end());

      double              sumW = 0.0;
      std::vector<double> weights;
      weights.reserve(aics.size());

      for (double aic : aics) {
        double w = std::exp(-0.5 * (aic - minAIC));
        weights.push_back(w);
        sumW += w;
      }

      std::cout << std::setw(30) << files[i].filename().string();
      for (double w : weights) std::cout << std::setw(15) << (w / sumW);
      std::cout << "\n";
    }

    std::cout << "\n";
  }

  void printBestModelSummary() const {
    std::cout << "\nModel Selection (Best Model per File):\n";
    std::cout << std::setw(25) << "File"
              << std::setw(20) << "Best Model"
              << std::setw(30) << "Parameters"
              << std::setw(15) << "AIC" << "\n";
    std::cout << std::string(25 + 20 + 30 + 15, '-') << "\n";

    for (int i = 0; i < files.size(); ++i) {
      double bestAIC      = std::numeric_limits<double>::max();
      int    bestModelIdx = -1;
      for (int j = 0; j < models.size(); ++j) {
        if (resultsAIC[i][j] < bestAIC) {
          bestAIC      = resultsAIC[i][j];
          bestModelIdx = j;
        }
      }

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

    std::cout << "\n";
  }

  void writeFittedCDFs() const {
    fs::path outputDir = "fitted_models";
    fs::create_directories(outputDir);

    for (int i = 0; i < files.size(); ++i) {
      const auto& degrees = read_distribution(files[i].string());
      if (degrees.empty()) continue;

      // Identify best model
      double bestAIC      = std::numeric_limits<double>::max();
      int    bestModelIdx = -1;
      for (int j = 0; j < models.size(); ++j) {
        if (resultsAIC[i][j] < bestAIC) {
          bestAIC      = resultsAIC[i][j];
          bestModelIdx = j;
        }
      }
      if (bestModelIdx < 0) continue;

      const auto&        theta     = results[i][bestModelIdx].args;
      const std::string& modelName = models[bestModelIdx];

      std::vector<float> unique_k;

      auto empirical_cdf = compute_empirical_cdf(degrees, unique_k);

      std::vector<float> theoretical_cdf;
      theoretical_cdf.reserve(unique_k.size());

      for (float k : unique_k) {
        if (modelName == "Geometric")
          theoretical_cdf.push_back(cdf_geometric(k, theta));
        else if (modelName == "Poisson")
          theoretical_cdf.push_back(cdf_poisson(k, theta));
        else if (modelName == "AltmannZeta" && USE_ALTERNATIVE)
          theoretical_cdf.push_back(cdf_altmann_zeta(k, theta[0], theta[1]));
        else if (modelName == "TruncZeta" && USE_ALTERNATIVE)
          theoretical_cdf.push_back(cdf_truncated_zeta(k, theta[0], static_cast<int>(theta[1])));

        else
          theoretical_cdf.push_back(0.0f);
      }

      fs::path      outFile = outputDir / (files[i].filename().stem().string() + "_" + modelName + "_fit.txt");
      std::ofstream out(outFile);
      if (!out.is_open()) {
        std::cerr << "Failed to open output: " << outFile << "\n";
        continue;
      }

      out << "K EmpiricalCDF TheoreticalCDF\n";
      for (size_t n = 0; n < unique_k.size() && n < theoretical_cdf.size(); ++n)
        out << unique_k[n] << " " << empirical_cdf[n] << " " << theoretical_cdf[n] << "\n";

      std::cout << " → " << outFile.filename().string() << " written.\n";
    }
  }

  void
  generatePlots() {
    fs::path outputDir  = "fitted_models";
    fs::path plotScript = "plot_cdf.R";
    if (!fs::exists(plotScript)) {
      std::cerr << "R script not found: " << plotScript << "\n";
      return;
    }

    for (const auto& file : fs::directory_iterator(outputDir)) {
      if (file.path().extension() != ".txt") continue;

      std::string inputFile  = file.path().string();
      std::string outputFile = file.path().string() + ".png";

      std::string stem      = file.path().stem().string();
      std::string firstWord = stem.substr(0, stem.find('_'));

      // Build a descriptive title
      std::string title = "CDF Comparison for " + firstWord;

      // Construct system call
      std::string cmd = "Rscript \"" + plotScript.string() + "\" \"" +
        inputFile + "\" \"" + outputFile + "\" \"" + title + "\"";

      std::cout << " → Generating plot: " << outputFile << "\n";
      int ret = std::system(cmd.c_str());
      if (ret != 0) {
        std::cerr << "Failed to generate plot for " << inputFile << "\n";
      }
    }
  }
};

// Usage:
int run_combined(const std::string& dir) {
  try {
    ModelFitter fitter(dir);
    std::cout << "=== Running optimizers on all datasets ===\n";
    fitter.runOptimizers();

    std::cout << "=== Computing AIC values ===\n";
    fitter.computeAIC();

    std::cout << "=== Printing negative log-likelihoods ===\n";
    fitter.printNegLogL();

    std::cout << "=== Printing fitted parameters ===\n";
    fitter.printParameters();

    std::cout << "=== Printing AIC table ===\n";
    fitter.printAICTable();

    std::cout << "=== Printing Akaike weights ===\n";
    fitter.printAkaikeWeights();

    std::cout << "=== Printing best model summary ===\n";
    fitter.printBestModelSummary();

    std::flush(std::cout);

    std::cout << "=== Writing fitted CDFs to output files ===\n";
    fitter.writeFittedCDFs();

    std::flush(std::cout);

    std::cout << "=== Generating CDFs plots===\n";
    fitter.generatePlots();

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }
  return 0;
}


int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " [-fast] [-compare] [-alternative] <-d> [file/directory] | [file]\n";
    return 1;
  }

  bool        isDir = false;
  std::string pathArg;

  // Global flags
  FAST            = false;
  COMPARE         = false;
  USE_ALTERNATIVE = false;

  // Parse arguments
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "-d") {
      if (i + 1 < argc) {
        isDir   = true;
        pathArg = argv[++i];
      } else {
        std::cerr << "Error: -d requires a directory argument.\n";
        return 1;
      }
    } else if (arg == "-fast") {
      FAST = true;
    } else if (arg == "-compare") {
      COMPARE = true;
    } else if (arg == "-alternative") {
      USE_ALTERNATIVE = true;
    } else if (pathArg.empty()) {
      pathArg = arg;
    } else {
      std::cerr << "Unknown argument: " << arg << "\n";
      return 1;
    }
  }

  if (pathArg.empty()) {
    std::cerr << "No input file or directory specified.\n";
    return 1;
  }

  if (isDir)
    return run_combined(pathArg);
  else
    return run_single(pathArg);
}
