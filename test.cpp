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
//recombination radius for Sn is 2.520, for Ge 3.090

class Parameters {
public:
    std::unordered_map<std::string, std::string> param_map;

    Parameters(const std::string& file) {
        std::ifstream in(file);
        std::string key, value;
        while (in >> key >> value) param_map[key] = value;
    }
	//mthod to retrieve a parameter value as a std::string.
    std::string get_string(const std::string& key) const {
        return param_map.at(key);
    }
	// First gets the string from param_map, then uses std::stod (string-to-double) to convert it
    double get_double(const std::string& key) const {
        return std::stod(param_map.at(key));
    }

    int get_int(const std::string& key) const {
        return std::stoi(param_map.at(key));
    }
};
 //Declares a global pointer named g_params to aParameters object
Parameters* g_params;

class SimulationConfig {
public:
    static std::string filename() { return g_params->get_string("filename"); }
    static double temperature;
};
double SimulationConfig::temperature = 0.0; // initialised outside the class because static members need to be defined separately

class SimulationConstants {
public:
    static constexpr double kB = 8.617333262145e-5;
    static double jump_distance() { return g_params->get_double("jump_distance"); }
    static double r_recombine() { return g_params->get_double("r_recombine"); }
    static double t_max;
    static int n_runs() { return g_params->get_int("n_runs"); }
};
double SimulationConstants::t_max = 0.0;

std::mt19937 rng(std::random_device{}());
std::uniform_real_distribution<> uniform(0.0, 1.0);

struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;
	//CALCULATES EUCLEDEAN DISTANCE
    double distance(const Vec3& other) const {
        return std::sqrt((x - other.x) * (x - other.x) +
                         (y - other.y) * (y - other.y) +
                         (z - other.z) * (z - other.z));
    }
	//GENERATES A RANDOM DIURECTION IN 3D
    bool random_jump() {
        double theta = std::acos(1.0 - 2.0 * uniform(rng));
        double phi = 2.0 * M_PI * uniform(rng);
        x += SimulationConstants::jump_distance() * std::sin(theta) * std::cos(phi);
        y += SimulationConstants::jump_distance() * std::sin(theta) * std::sin(phi);
        z += SimulationConstants::jump_distance() * std::cos(theta);
        
	//CELL BOUNDARY CONDITIONS
        if (x > 214.0) x = 0.0;
        if (x < 0.0)   x = 214.0;
        if (y > 214.0) y = 0.0;
        if (y < 0.0)   y = 214.0;
        if (z < -107.0) z = -106.0;

        return z > 213.0; //IF ESCAPED 
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
    virtual bool is_interstitial() const { return false; }
    virtual bool is_vacancy() const { return false; }
    virtual bool is_impurity() const { return false; }

    virtual char symbol() const = 0; //this is for the defect name when writing the file
	// akes the defect perform a random jump. Return true if escaped
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
	//initializes a defect with a given position p.
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
    
    bool is_vacancy() const override { return true; }
    char symbol() const override { return 'V'; }
};

class Impurity : public Defect {
public:
    Impurity(const Vec3& p) : Defect(p) {}
    void compute_rate(double) override {
        rate = 0.0;
    }
    
    bool is_impurity() const override { return true; }
    char symbol() const override { return 'X'; }
    static double r_vx() { return g_params->get_double("r_vx"); }
};

std::vector<std::unique_ptr<Defect>> read_defects(const std::string& filename, double T) {
    std::ifstream infile(filename);
    if (!infile) throw std::runtime_error("Failed to open file: " + filename);

    std::vector<std::unique_ptr<Defect>> defects; //store all defects read
    std::string line; //read and discard the first line
    std::getline(infile, line);
    while (std::getline(infile, line)) { //read all lines, skip if empty
        if (line.empty()) continue;
        std::istringstream iss(line);
        char ch; //ch used to ignore certain characters 
        Vec3 pos;
        int occ1 = 0, occ2 = 0;
        double dist = 0.0;

        iss >> ch >> pos.x >> pos.y >> pos.z >> ch;
        iss >> ch >> occ1 >> occ2 >> ch;
        iss >> dist;
	//identify the defect based on the occupancy in the data file
        std::unique_ptr<Defect> dptr;
        if (occ1 == 2 && occ2 == 0) {
            dptr = std::make_unique<Interstitial>(pos);
        } else if (occ1 == 0 && occ2 == 0) {
            dptr = std::make_unique<Vacancy>(pos);
        } else {
            dptr = std::make_unique<Impurity>(pos);
        }

        dptr->compute_rate(T);
        defects.push_back(std::move(dptr)); //Moves the unique pointer into the vector.
    }

    return defects; //Returns the vector of unique pointers to all the created defects.
}
//void == appplies changes directry, no return value
void recombine_defects(std::vector<std::unique_ptr<Defect>>& defects) { 
    double r_recomb = SimulationConstants::r_recombine();
    double r_vx = Impurity::r_vx();

    for (size_t i = 0; i < defects.size(); ++i) {
        if (!defects[i]->alive) continue;

        for (size_t j = i + 1; j < defects.size(); ++j) {
            if (!defects[j]->alive) continue;

            double dist = defects[i]->pos.distance(defects[j]->pos);
            if (dist > r_recomb) continue;

            // Immediate recombination V + I 
            if ((defects[i]->is_interstitial() && defects[j]->is_vacancy()) ||
                (defects[j]->is_interstitial() && defects[i]->is_vacancy())) {

                defects[i]->alive = false;
                defects[j]->alive = false;
                break;
            }
            
            //X + V immediate recombination
            if ((defects[i]->is_vacancy() && defects[j]->is_impurity()) ||
                (defects[j]->is_vacancy() && defects[i]->is_impurity())) {

                defects[i]->alive = false;
                defects[j]->alive = false;

                break;  
            }
        }
    }

    // make sure X and V are not alivve and the distance is less than r_recomb, calculate X-V midpoint 
    for (size_t i = 0; i < defects.size(); ++i) {
        if (defects[i]->alive) continue;
        if (!defects[i]->is_vacancy()) continue;

        for (size_t j = 0; j < defects.size(); ++j) {
            if (i == j) continue;
            if (defects[j]->alive) continue;
            if (!defects[j]->is_impurity()) continue;

            double dist = defects[i]->pos.distance(defects[j]->pos);
            if (dist < r_recomb) continue;

            Vec3 midpoint{
                0.5 * (defects[i]->pos.x + defects[j]->pos.x),
                0.5 * (defects[i]->pos.y + defects[j]->pos.y),
                0.5 * (defects[i]->pos.z + defects[j]->pos.z)
            };

            // Check if I is within r_vx of this midpoint
            for (size_t k = 0; k < defects.size(); ++k) {
                if (!defects[k]->alive) continue;
                if (!defects[k]->is_interstitial()) continue;

                double dist_ik = defects[k]->pos.distance(midpoint);
                if (dist_ik <= r_vx) {
                    // Kill interstitial and vacancy, revive X
                    defects[k]->alive = false;
                    defects[i]->alive = false;
                    defects[j]->alive = true; 
                    break;
                }
            }
        }
    }
}


void simulate(const std::string& filename, double T,
              double& surviving_frac, double& jump_ratio,
              int& total_jumps_i, int& total_jumps_v) {
    std::vector<std::unique_ptr<Defect>> defects = read_defects(filename, T);
    std::ofstream traj("trajectory.xyz");

    double t = 0.0;
    while (t < SimulationConstants::t_max) {
        int alive_count = 0;
        for (const auto& d : defects) {
            if (d->alive) ++alive_count;
        }
        traj << alive_count << "\n";
	traj << "Lattice=\"214.438 0.0 0.0 0.0 214.438 0.0 0.0 0.0 346.657\" "
	     << "Origin=\"-0.104531 -0.104531 -107.324\" "
	     << "Properties=species:S:1:pos:R:3 Time=" << t << "\n";
        for (const auto& d : defects) {
            if (!d->alive) continue;
            traj << d->symbol() << " " << d->pos.x << " " << d->pos.y << " " << d->pos.z << "\n";
        }

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
//STATISTICS
void run_stats(const std::string& filename, double T) {
    int n_runs = SimulationConstants::n_runs();
    double sum_surviving_frac = 0.0;
    double sum_jump_ratio = 0.0;
    int sum_total_jumps_i = 0;
    int sum_total_jumps_v = 0;

    for (int run = 0; run < n_runs; ++run) {
    
	    double surviving_frac = 0.0;
	    double jump_ratio = 0.0;
	    int total_jumps_i = 0;
	    int total_jumps_v = 0;

	    simulate(filename, T, surviving_frac, jump_ratio, total_jumps_i, total_jumps_v);
	    sum_surviving_frac += surviving_frac;
            sum_jump_ratio += jump_ratio;
            sum_total_jumps_i += total_jumps_i;
            sum_total_jumps_v += total_jumps_v;
    }
    std::cout << "Total runs: " << (n_runs) << "\n";
    std::cout << "Average surviving fraction: " << (sum_surviving_frac / n_runs) << "\n";
    std::cout << "Average jump ratio: " << (sum_jump_ratio / n_runs) << "\n";
    std::cout << "Average interstitial jumps: " << sum_total_jumps_i / n_runs << "\n";
    std::cout << "Average vacancy jumps: " << sum_total_jumps_v / n_runs << "\n";
}
	    
//make it possibele to define temperature and time without making changes to the script 
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

