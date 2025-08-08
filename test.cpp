#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <random>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <memory>
#include <limits>
#include <unordered_map>

class Parameters {
public:
    std::unordered_map<std::string, std::string> param_map;

    Parameters(const std::string& file) {
        std::ifstream in(file);
        std::string key, value;
        while (in >> key >> value) param_map[key] = value;
    }

    std::string get_string(const std::string& key) const {
        return param_map.at(key);
    }

    double get_double(const std::string& key) const {
        return std::stod(param_map.at(key));
    }

    int get_int(const std::string& key) const {
        return std::stoi(param_map.at(key));
    }
};

Parameters* g_params;

class SimulationConfig {
public:
    static std::string filename() { return g_params->get_string("filename"); }
    static double temperature;
};
double SimulationConfig::temperature = 0.0;

class SimulationConstants {
public:
    static constexpr double kB = 8.617333262145e-5;
    static double jump_distance() { return g_params->get_double("jump_distance"); }
    static double r_recombine() { return g_params->get_double("r_recombine"); }
    static double t_max;
    static int n_runs() { return g_params->get_int("n_runs"); }
};
double SimulationConstants::t_max = 0.0;

// RANDOMISER
std::mt19937 rng(std::random_device{}());
std::uniform_real_distribution<> uniform(0.0, 1.0);

struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;

    double distance(const Vec3& other) const {
        return std::sqrt((x - other.x) * (x - other.x) +
                         (y - other.y) * (y - other.y) +
                         (z - other.z) * (z - other.z));
    }
    
    bool random_jump() {
        double theta = std::acos(1.0 - 2.0 * uniform(rng));
        double phi = 2.0 * M_PI * uniform(rng);
        x += SimulationConstants::jump_distance() * std::sin(theta) * std::cos(phi);
        y += SimulationConstants::jump_distance() * std::sin(theta) * std::sin(phi);
        z += SimulationConstants::jump_distance() * std::cos(theta);

        if (x > 214.0) x = 0.0;
        if (x < 0.0)   x = 214.0;
        if (y > 214.0) y = 0.0;
        if (y < 0.0)   y = 214.0;
        if (z < -107.0) z = -106.0;

        return z > 213.0;
    }
};

class Defect {
public:
    Vec3 pos;
    double rate = 0.0;
    int jump_count = 0;
    bool alive = true;

    Defect(const Vec3& p) : pos(p) {}
    virtual ~Defect() = default;

    virtual void compute_rate(double T) = 0;
    virtual bool is_interstitial() const = 0;
    virtual char symbol() const = 0;

    bool do_jump_and_update() {
        bool escaped = pos.random_jump();
        if (escaped) {
            alive = false;
            return true;
        }
        ++jump_count;
        return false;
    }
};

class Interstitial : public Defect {
public:
    Interstitial(const Vec3& p) : Defect(p) {}
    void compute_rate(double T) override {
        double wi = g_params->get_double("wi");
        double Em_i = g_params->get_double("Em_i");
        rate = wi * std::exp(-Em_i / (SimulationConstants::kB * T));
    }
    bool is_interstitial() const override { return true; }
    char symbol() const override { return 'I'; }
};

class Vacancy : public Defect {
public:
    Vacancy(const Vec3& p) : Defect(p) {}
    void compute_rate(double T) override {
        double wv = g_params->get_double("wv");
        double Em_v = g_params->get_double("Em_v");
        rate = wv * std::exp(-Em_v / (SimulationConstants::kB * T));
    }
    bool is_interstitial() const override { return false; }
    char symbol() const override { return 'V'; }
};

class Impurity : public Defect {
public:
    Impurity(const Vec3& p) : Defect(p) {}
    void compute_rate(double) override {
        rate = 0.0;
    }
    bool is_interstitial() const override { return false; }
    char symbol() const override { return 'X'; }
    static double Em_recomb() { return g_params->get_double("Em_recomb"); }
};

std::vector<std::unique_ptr<Defect>> read_defects(const std::string& filename, double T) {
    std::ifstream infile(filename);
    if (!infile) throw std::runtime_error("Failed to open file: " + filename);

    std::vector<std::unique_ptr<Defect>> defects;
    std::string line;
    std::getline(infile, line); 
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        char ch;
        Vec3 pos;
        int occ1 = 0, occ2 = 0;
        double dist = 0.0;

        iss >> ch >> pos.x >> pos.y >> pos.z >> ch;
        iss >> ch >> occ1 >> occ2 >> ch;
        iss >> dist;

        std::unique_ptr<Defect> dptr;
        if (occ1 == 2 && occ2 == 0) {
            dptr = std::make_unique<Interstitial>(pos);
        } else if (occ1 == 0 && occ2 == 0) {
            dptr = std::make_unique<Vacancy>(pos);
        } else {
            dptr = std::make_unique<Impurity>(pos);
        }

        dptr->compute_rate(T);
        defects.push_back(std::move(dptr));
    }

    return defects;
}

// NEW FUNCTION: Attempt recombination after each jump event
void recombine_defects(std::vector<std::unique_ptr<Defect>>& defects) {
    double r_recomb = SimulationConstants::r_recombine();
    // Simple O(N^2) pairwise check - might be slow for large systems
    for (size_t i = 0; i < defects.size(); ++i) {
        if (!defects[i]->alive || !defects[i]->is_interstitial()) continue;

        for (size_t j = 0; j < defects.size(); ++j) {
            if (i == j) continue;
            if (!defects[j]->alive || defects[j]->is_interstitial()) continue;

            // Check distance
            if (defects[i]->pos.distance(defects[j]->pos) <= r_recomb) {
                // Recombine: kill both defects
                defects[i]->alive = false;
                defects[j]->alive = false;
                break;  // Move to next interstitial
            }
        }
    }
}

void simulate(const std::string& filename, double T,
              double& surviving_frac, double& jump_ratio,
              int& total_jumps_i, int& total_jumps_v) {
    std::vector<std::unique_ptr<Defect>> defects = read_defects(filename, T);
    double t = 0.0;
    while (t < SimulationConstants::t_max) {
        double total_rate = 0.0;
        for (auto& d : defects) {
            if (d->alive) total_rate += d->rate;
        }
        if (total_rate == 0.0) break;
        double dt = -std::log(uniform(rng)) / total_rate;
        t += dt;

        double r = uniform(rng) * total_rate;
        double cumulative = 0.0;
        for (auto& d : defects) {
            if (!d->alive) continue;
            cumulative += d->rate;
            if (cumulative >= r) {
                d->do_jump_and_update();

                // After a defect jumps, attempt recombination
                recombine_defects(defects);

                break;
            }
        }
    }

    int total_i = 0, alive_i = 0;
    int total_v = 0, alive_v = 0;
    total_jumps_i = 0;
    total_jumps_v = 0;

    for (const auto& d : defects) {
        if (d->is_interstitial()) {
            ++total_i;
            if (d->alive) ++alive_i;
            total_jumps_i += d->jump_count;
        } else if (d->symbol() == 'V') {
            ++total_v;
            if (d->alive) ++alive_v;
            total_jumps_v += d->jump_count;
        }
    }

    surviving_frac = static_cast<double>(alive_i + alive_v) / (total_i + total_v);
    jump_ratio = total_jumps_v > 0 ? static_cast<double>(total_jumps_i) / total_jumps_v : 0.0;
}

void run_stats(const std::string& filename, double T) {
    double surviving_frac = 0.0;
    double jump_ratio = 0.0;
    int total_jumps_i = 0;
    int total_jumps_v = 0;

    simulate(filename, T, surviving_frac, jump_ratio, total_jumps_i, total_jumps_v);

    std::cout << "Surviving fraction: " << surviving_frac << "\n";
    std::cout << "Jump ratio: " << jump_ratio << "\n";
    std::cout << "Interstitial jumps: " << total_jumps_i << "\n";
    std::cout << "Vacancy jumps: " << total_jumps_v << "\n";
}

int main(int argc, char* argv[]) {
    try {
        Parameters params("parameters.txt");
        g_params = &params;

        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg.substr(0, 2) == "-T") {
                SimulationConfig::temperature = std::stod(arg.substr(2));
            } else if (arg.substr(0, 2) == "-t") {
                SimulationConstants::t_max = std::stod(arg.substr(2));
            }
        }

        run_stats(SimulationConfig::filename(), SimulationConfig::temperature);
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
    return 0;
}

