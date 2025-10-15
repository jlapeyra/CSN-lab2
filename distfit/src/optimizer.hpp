#pragma once
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <vector>
#include "html_generator.hpp"

// ---------------------------------------------------------------
// 🧩 Utilities
// ---------------------------------------------------------------

inline float randf() {
  static thread_local std::mt19937                          gen(std::random_device{}());
  static thread_local std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  return dist(gen);
}

inline float vectorDistance(const std::vector<float>& a, const std::vector<float>& b) {
  float d = 0.f;
  for (size_t i = 0; i < a.size(); ++i)
    d += (a[i] - b[i]) * (a[i] - b[i]);
  return std::sqrt(d);
}

// ---------------------------------------------------------------
// 🌌 Vector Field Configuration
// ---------------------------------------------------------------

struct VectorFieldConfiguration {
  struct Field {
    float min;
    float max;
  };

  std::vector<Field> fields;

  std::vector<float> randomVector() const {
    std::vector<float> v(fields.size());
    for (size_t i = 0; i < fields.size(); ++i)
      v[i] = randf() * (fields[i].max - fields[i].min) + fields[i].min;
    return v;
  }

  size_t dimensions() const { return fields.size(); }
};

// ---------------------------------------------------------------
// 🧠 Function Definition (templated, callable-based)
// ---------------------------------------------------------------

template <typename Callable> struct FunctionDefinition {
  VectorFieldConfiguration config;
  Callable                 func;
  std::string              name;
  float                    expectedOptimal = 0.0f;

  float operator()(const std::vector<float>& args) const { return func(args); }
};

// ---------------------------------------------------------------
// ⚙️ Optimizer Base Types
// ---------------------------------------------------------------

struct OptimizerArguments {
  int                          maxIterations = 1000;
  std::map<std::string, float> hyper;
};

struct OptimizerResult {
  std::vector<float> args;
  float              cost       = std::numeric_limits<float>::infinity();
  int                iterations = 0;
};


// ---------------------------------------------------------------
// 🧩 Example Function Builders
// ---------------------------------------------------------------

inline auto makeOptimizationFunctionVoronoi(size_t dim, int points) {
  VectorFieldConfiguration cfg;
  cfg.fields.resize(dim, {0.0f, 1.0f});

  std::vector<std::vector<float>> centers(points);
  for (auto& c : centers)
    c = cfg.randomVector();

  auto fn = [centers](const std::vector<float>& args) -> float {
    float minDist = std::numeric_limits<float>::infinity();
    for (auto& c : centers)
      minDist = std::min(minDist, vectorDistance(c, args));
    return minDist;
  };

  return FunctionDefinition<decltype(fn)>{cfg, fn, "Voronoi", 0.0f};
}

// ---------------------------------------------------------------
// 🧬 Random Optimizer
// ---------------------------------------------------------------

template <typename Function> struct RandomOptimizer {
  std::string name = "Random";

  OptimizerResult operator()(const OptimizerArguments& args, const Function& f) const {
    const auto&     cfg = f.config;
    OptimizerResult r;
    for (int i = 0; i < args.maxIterations; ++i) {
      auto  candidate = cfg.randomVector();
      float score     = f(candidate);
      if (score < r.cost) {
        r.cost = score;
        r.args = candidate;
      }
    }
    r.iterations = args.maxIterations;
    return r;
  }
};

// ---------------------------------------------------------------
// 🧬 Random Stall Optimizer with Early Stopping
// ---------------------------------------------------------------

template <typename Function>
struct RandomStallOptimizer {
  std::string name = "RandomStall";

  struct Result : public OptimizerResult {
    bool earlyStopped = false;
  };

  Result operator()(const OptimizerArguments& args, const Function& f) const {
    const auto& cfg = f.config;

    Result r;
    r.cost = std::numeric_limits<float>::infinity();

    float epsilon    = 1e-3;
    int   stallLimit = std::min(100, args.maxIterations / 10);

    int stallCounter = 0;

    for (int i = 0; i < args.maxIterations; ++i) {
      auto  candidate = cfg.randomVector();
      float score     = f(candidate);

      if (score < (r.cost - epsilon))
        stallCounter = 0;
      else
        stallCounter++;

      if (score < r.cost) {
        r.cost = score;
        r.args = candidate;
      }

      if (stallCounter >= stallLimit) {
        r.earlyStopped = true;
        r.iterations   = i + 1;
        return r;
      }
    }

    r.iterations = args.maxIterations;
    return r;
  }
};

// ---------------------------------------------------------------
// 🔥 Simulated Annealing
// ---------------------------------------------------------------

template <typename Function> struct SimulatedAnnealingOptimizer {
  std::string name = "SimulatedAnnealing";

  OptimizerResult operator()(const OptimizerArguments& args,
                             const Function&           f) const {
    const auto& cfg          = f.config;
    auto        current      = cfg.randomVector();
    float       currentScore = f(current);
    auto        best         = current;
    float       bestScore    = currentScore;

    float temperature =
      args.hyper.count("temperature") ? args.hyper.at("temperature") : 1.0f;
    const float coolingRate =
      args.hyper.count("cooling") ? args.hyper.at("cooling") : 0.995f;
    const float minTemperature = 0.0001f;

    int iterations = 0;
    while (temperature > minTemperature && iterations < args.maxIterations) {
      auto   neighbor = current;
      size_t idx      = rand() % neighbor.size();
      neighbor[idx]   = cfg.fields[idx].min +
        randf() * (cfg.fields[idx].max - cfg.fields[idx].min);
      float neighborScore = f(neighbor);

      if (neighborScore < currentScore ||
          randf() < std::exp((currentScore - neighborScore) / temperature)) {
        current      = neighbor;
        currentScore = neighborScore;
      }

      if (neighborScore < bestScore) {
        best      = neighbor;
        bestScore = neighborScore;
      }

      temperature *= coolingRate;
      ++iterations;
    }

    return {best, bestScore, iterations};
  }
};

// ---------------------------------------------------------------
// 🧗 Hill Climbing
// ---------------------------------------------------------------

template <typename Function> struct HillClimbingOptimizer {
  std::string name = "HillClimbing";

  OptimizerResult operator()(const OptimizerArguments& args,
                             const Function&           f) const {
    const auto& cfg          = f.config;
    auto        current      = cfg.randomVector();
    float       currentScore = f(current);

    auto  best      = current;
    float bestScore = currentScore;

    const float stepSize =
      args.hyper.count("step") ? args.hyper.at("step") : 0.1f;
    const int neighbors = 10;

    for (int i = 0; i < args.maxIterations; ++i) {
      std::vector<float> bestNeighbor      = current;
      float              bestNeighborScore = currentScore;

      for (int j = 0; j < neighbors; ++j) {
        auto neighbor = current;
        for (size_t k = 0; k < neighbor.size(); ++k) {
          float range = cfg.fields[k].max - cfg.fields[k].min;
          float delta = (randf() * 2.0f - 1.0f) * stepSize * range;
          neighbor[k] = std::clamp(neighbor[k] + delta, cfg.fields[k].min, cfg.fields[k].max);
        }
        float neighborScore = f(neighbor);
        if (neighborScore < bestNeighborScore) {
          bestNeighbor      = neighbor;
          bestNeighborScore = neighborScore;
        }
      }

      current      = bestNeighbor;
      currentScore = bestNeighborScore;

      if (currentScore < bestScore) {
        best      = current;
        bestScore = currentScore;
      }
    }

    return {best, bestScore, args.maxIterations};
  }
};

// ---------------------------------------------------------------
// 🚶 Random Walk
// ---------------------------------------------------------------

template <typename Function> struct RandomWalkOptimizer {
  std::string name = "RandomWalk";

  OptimizerResult operator()(const OptimizerArguments& args,
                             const Function&           f) const {
    const auto& cfg          = f.config;
    auto        current      = cfg.randomVector();
    float       currentScore = f(current);

    auto  best      = current;
    float bestScore = currentScore;

    const float step = args.hyper.count("step") ? args.hyper.at("step") : 0.05f;

    for (int i = 0; i < args.maxIterations; ++i) {
      auto next = current;
      for (size_t j = 0; j < next.size(); ++j) {
        float range = cfg.fields[j].max - cfg.fields[j].min;
        float delta = (randf() * 2.0f - 1.0f) * step * range;
        next[j] =
          std::clamp(next[j] + delta, cfg.fields[j].min, cfg.fields[j].max);
      }
      float score = f(next);
      if (score < currentScore) {
        current      = next;
        currentScore = score;
      }

      if (score < bestScore) {
        best      = next;
        bestScore = score;
      }
    }

    return {best, bestScore, args.maxIterations};
  }
};

// ---------------------------------------------------------------
// 🌻 Inverse Sunflower (quasi-random spiral sampling)
// ---------------------------------------------------------------

template <typename Function> struct InverseSunflowerOptimizer {
  std::string name = "InverseSunflower";

  OptimizerResult operator()(const OptimizerArguments& args,
                             const Function&           f) const {
    const auto& cfg         = f.config;
    const float goldenAngle = 3.14159265359f * (3.0f - std::sqrt(5.0f));

    OptimizerResult r;

    for (int i = 0; i < args.maxIterations; ++i) {
      float radius = std::sqrt((float)i / args.maxIterations);
      float theta  = i * goldenAngle;

      auto candidate = cfg.randomVector();
      if (cfg.dimensions() >= 2) {
        candidate[0] = 0.5f + radius * std::cos(theta);
        candidate[1] = 0.5f + radius * std::sin(theta);
      }

      for (size_t j = 0; j < candidate.size(); ++j)
        candidate[j] =
          std::clamp(candidate[j], cfg.fields[j].min, cfg.fields[j].max);

      float score = f(candidate);
      if (score < r.cost) {
        r.cost = score;
        r.args = candidate;
      }
    }

    r.iterations = args.maxIterations;
    return r;
  }
};

// ---------------------------------------------------------------
// 🧬 Differential Evolution
// ---------------------------------------------------------------

template <typename Function> struct DifferentialEvolutionOptimizer {
  std::string name = "DifferentialEvolution";

  OptimizerResult operator()(const OptimizerArguments& args,
                             const Function&           f) const {
    const auto& cfg = f.config;
    const int   popSize =
      args.hyper.count("population") ? (int)args.hyper.at("population") : 20;
    const float F  = args.hyper.count("F") ? args.hyper.at("F") : 0.8f;
    const float CR = args.hyper.count("CR") ? args.hyper.at("CR") : 0.9f;

    std::vector<std::vector<float>> population(popSize);
    std::vector<float>              scores(popSize);
    for (int i = 0; i < popSize; ++i) {
      population[i] = cfg.randomVector();
      scores[i]     = f(population[i]);
    }

    for (int iter = 0; iter < args.maxIterations; ++iter) {
      for (int i = 0; i < popSize; ++i) {
        int a, b, c;
        do {
          a = rand() % popSize;
        } while (a == i);
        do {
          b = rand() % popSize;
        } while (b == i || b == a);
        do {
          c = rand() % popSize;
        } while (c == i || c == a || c == b);

        std::vector<float> trial = population[i];
        int                R     = rand() % cfg.dimensions();
        for (size_t j = 0; j < cfg.dimensions(); ++j) {
          if (randf() < CR || j == (size_t)R) {
            trial[j] =
              population[a][j] + F * (population[b][j] - population[c][j]);
            trial[j] =
              std::clamp(trial[j], cfg.fields[j].min, cfg.fields[j].max);
          }
        }

        float trialScore = f(trial);
        if (trialScore < scores[i]) {
          population[i] = trial;
          scores[i]     = trialScore;
        }
      }
    }

    int bestIndex =
      std::min_element(scores.begin(), scores.end()) - scores.begin();
    return {population[bestIndex], scores[bestIndex], args.maxIterations};
  }
};

// ---------------------------------------------------------------
// 🌻 Inverse Sunflower Finder (Penalized + Multi-Candidate)
// ---------------------------------------------------------------

template <typename Function>
struct InverseSunflowerFinderOptimizer {
  std::string name = "InverseSunflowerFinder";

  OptimizerResult operator()(const OptimizerArguments& conf, const Function& f) const {
    const auto&  config = f.config;
    const size_t dims   = config.dimensions();

    const size_t K            = conf.hyper.count("K") ? (size_t)conf.hyper.at("K") : 2;
    const size_t N            = conf.hyper.count("N") ? (size_t)conf.hyper.at("N") : 20;
    float        alpha        = conf.hyper.count("alpha") ? conf.hyper.at("alpha") : 0.5f;
    const float  lambda       = conf.hyper.count("lambda") ? conf.hyper.at("lambda") : 0.99f;
    const float  delta        = conf.hyper.count("delta") ? conf.hyper.at("delta") : 0.05f;
    const float  penaltyConst = conf.hyper.count("penalty") ? conf.hyper.at("penalty") : 0.04f;
    const float  minAlpha     = conf.hyper.count("minAlpha") ? conf.hyper.at("minAlpha") : 0.0001f;

    std::mt19937                    gen(std::random_device{}());
    std::normal_distribution<float> dist(0.0f, 1.0f);

    // Initial root point(s): center of search space
    std::vector<std::vector<float>> rootPoints(1, std::vector<float>(dims));
    for (size_t i = 0; i < dims; i++)
      rootPoints[0][i] = (config.fields[i].min + config.fields[i].max) * 0.5f;

    float              bestScore = f(rootPoints[0]);
    std::vector<float> best      = rootPoints[0];

    int iterations = 0;
    while (alpha > minAlpha && iterations < conf.maxIterations) {
      const size_t                    total = N * rootPoints.size();
      std::vector<std::vector<float>> points(total, std::vector<float>(dims));
      std::vector<float>              costs(total, 0.0f);

      // --- Sample points around current roots ---
      for (size_t k = 0; k < rootPoints.size(); ++k) {
        for (size_t i = 0; i < N; ++i) {
          std::vector<float> direction(dims);
          float              norm = 0.0f;
          for (size_t j = 0; j < dims; ++j) {
            direction[j] = dist(gen);
            norm += direction[j] * direction[j];
          }
          norm = std::sqrt(norm);
          if (norm > 0)
            for (float& d : direction) d /= norm;

          size_t pidx = i + N * k;
          for (size_t j = 0; j < dims; ++j) {
            points[pidx][j] = rootPoints[k][j] + direction[j] * alpha;
            points[pidx][j] = std::clamp(points[pidx][j],
                                         config.fields[j].min,
                                         config.fields[j].max);
          }
          costs[pidx] = f(points[pidx]);
        }
      }

      // --- Apply penalized proximity cost ---
      std::vector<float> penalized = costs;
      for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = i + 1; j < points.size(); ++j) {
          float d = vectorDistance(points[i], points[j]);
          if (d < delta) {
            float penalty = penaltyConst * (delta - d);
            penalized[i] += penalty;
            penalized[j] += penalty;
          }
        }
      }

      // --- Sort and pick K best ---
      std::vector<size_t> indices(points.size());
      std::iota(indices.begin(), indices.end(), 0);
      std::sort(indices.begin(), indices.end(),
                [&](size_t a, size_t b) { return penalized[a] < penalized[b]; });

      rootPoints.resize(K);
      for (size_t i = 0; i < K; ++i)
        rootPoints[i] = points[indices[i]];

      best      = rootPoints[0];
      bestScore = costs[indices[0]];

      // --- Step decay ---
      alpha *= lambda;
      ++iterations;
    }

    return {best, bestScore, iterations};
  }
};

// --- Optimizer factory ---
template <typename Function>
auto getOptimizer(const std::string& optimizerID) {
  if (optimizerID == "Random")
    return RandomOptimizer<Function>{};
  else if (optimizerID == "HillClimb")
    return HillClimbingOptimizer<Function>{};
  else if (optimizerID == "SA")
    return SimulatedAnnealingOptimizer<Function>{};
  else if (optimizerID == "RW")
    return RandomWalkOptimizer<Function>{};
  else if (optimizerID == "InverseSunflower")
    return InverseSunflowerOptimizer<Function>{};
  else if (optimizerID == "DifferentialEvolution")
    return DifferentialEvolutionOptimizer<Function>{};
  else
    throw std::runtime_error("Unknown optimizer: " + optimizerID);
}

void printOptimizerList() {
  std::vector<std::string> optimizers = {
    "Random",
    "HillClimb",
    "SA",
    "RW",
    "InverseSunflower",
    "DifferentialEvolution"};

  std::cout << "Available optimizers:\n";
  std::cout << "--------------------\n";
  for (const auto& opt : optimizers) {
    std::cout << " - " << opt << "\n";
  }
  std::cout << std::endl;
}

// ---------------------------------------------------------------
// 📊 Runner
// ---------------------------------------------------------------

template <typename Optimizer, typename Function>
auto optimizeReportData(const Optimizer& opt, const Function& fn) {
  std::cout << "Running optimizer [" << opt.name << "] on function [" << fn.name << "]\n";

  int                             k    = 10;
  float                           best = std::numeric_limits<float>::infinity();
  std::vector<float>              bestArgs;
  std::vector<int>                iterationSteps;
  std::vector<float>              costHistory;
  std::vector<std::vector<float>> argHistory; // each argument dimension evolution

  for (int i = 0; i < 4; ++i) {
    OptimizerArguments args;
    args.maxIterations = k;
    auto result        = opt(args, fn);

    std::cout << "Iterations: " << k << " → cost: " << result.cost << "\n";
    iterationSteps.push_back(k);
    costHistory.push_back(result.cost);

    if (!result.args.empty()) {
      if (argHistory.empty()) argHistory.resize(result.args.size());
      for (size_t j = 0; j < result.args.size(); ++j)
        argHistory[j].push_back(result.args[j]);
    }

    if (result.cost < best) {
      best     = result.cost;
      bestArgs = result.args;
    }
    k *= 10;
  }

  std::cout << "Best cost found: " << best << "\n\n";

  return std::tuple{
    iterationSteps, costHistory, argHistory, best, bestArgs};
}

// ----------------------------------------------------------------------
// 📊 HTML Optimizer Evolution Report
// ----------------------------------------------------------------------
template <typename Optimizer, typename Function>
void optimizeReportHTML(const Optimizer& opt, const Function& fn, const std::string& htmlFile = "optimizer_report.html") {
  auto [iterations, costs, argHist, bestCost, bestArgs] = optimizeReportData(opt, fn);

  ScopedRedirect out(htmlFile);

  html_generator html;
  html << html_generator::html_begin{};

  std::cout << "<h1 style='color:#ff7f00;font-family:JetBrains Mono, monospace;'>"
            << "Optimizer Report: " << opt.name << "</h1>\n";

  html << html_generator::table_begin{};
  html << html_generator::header_row_begin{} << "Iterations" << "Cost" << html_generator::row_end{};
  for (size_t i = 0; i < iterations.size(); ++i)
    html << html_generator::row_begin{} << iterations[i] << costs[i] << html_generator::row_end{};
  html << html_generator::table_end{};

  std::vector<std::string> labels;
  for (auto k : iterations) labels.push_back(std::to_string(k));
  html << html_generator::chart_begin{
    .id             = "cost_chart",
    .labels         = labels,
    .datasets       = {std::vector<double>(costs.begin(), costs.end())},
    .dataset_labels = {opt.name}};
  html << html_generator::chart_end{};

  if (!argHist.empty()) {
    std::vector<std::vector<double>> argData(argHist.size());
    for (size_t i = 0; i < argHist.size(); ++i)
      argData[i] = std::vector<double>(argHist[i].begin(), argHist[i].end());

    std::vector<std::string> argLabels;
    for (size_t i = 0; i < argData.size(); ++i)
      argLabels.push_back("arg" + std::to_string(i));

    html << html_generator::chart_begin{
      .id             = "arg_evolution",
      .labels         = labels,
      .datasets       = argData,
      .dataset_labels = argLabels};
    html << html_generator::chart_end{};
  }

  std::cout << "<h3>Best cost: " << bestCost << "</h3>\n";
  html << html_generator::html_end{};

  std::cout.rdbuf(nullptr); // flush
  std::cerr << "✅ HTML report written to " << htmlFile << "\n";
}
