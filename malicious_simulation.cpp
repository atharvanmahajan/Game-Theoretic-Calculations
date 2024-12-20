#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <numeric>
#include <algorithm>
#include <map>
#include <cstdlib>
#include <limits>

// Function to generate power law distribution
std::vector<int> generatePowerLawDistribution(int N, double r_0, double alpha, std::mt19937& gen, std::uniform_real_distribution<>& dis) {
    std::vector<int> distribution;
    distribution.reserve(N);
    int N_0 = static_cast<int>(std::floor(N * r_0));
    int remaining = N - N_0;
    distribution.push_back(N_0);
    
    while (remaining > 0) {
        double u = dis(gen);
        int value = std::ceil(std::exp(-std::log(u) / alpha));
        if (value > remaining || value > std::numeric_limits<int>::max() / 2) {
            value = remaining;
        }
        distribution.push_back(value);
        remaining -= value;
    }
    
    return distribution;
}

// Function to perform the simulation
void simulate(int n, double r_0, int N_iter, int n_iter, double alpha) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);

    // Generate the group sizes using power law distribution
    std::vector<int> N = generatePowerLawDistribution(n, r_0, alpha, gen, dis);
    int idx_largest = std::max_element(N.begin() + 1, N.end()) - N.begin();
    N[idx_largest] -= 1;

    std::vector<double> exp_m_arr;
    std::vector<double> prob_q_arr;

    // Create the array of group sizes weighting
    std::vector<int> array;
    for (size_t j = 0; j < N.size(); ++j) {
        array.insert(array.end(), N[j], j);
    }

    // Simulation iterations
    for (int i = 0; i < N_iter; ++i) {
        std::map<int, std::pair<int, int>> m_dict;

        for (int j = 0; j < n_iter; ++j) {
            std::shuffle(array.begin(), array.end(), gen);
            std::vector<int> init(N.size(), 0);
            init[idx_largest] = 1;
            int idx = 0;
            int m = 1;
            double cutoff = 0.5;
            std::string status = "lost";

            while (true) {
                m += 1;
                if (idx >= array.size()) break;
                int selected_number = array[idx];
                idx += 1;
                init[selected_number] += 1;
                cutoff += 0.5;
                if (m >= 3 && init[idx_largest] >= cutoff) {
                    status = "won";
                    break;
                }
                else if (m >= 3 && std::any_of(init.begin(), init.end(), [cutoff](int i){ return i > cutoff; })) {
                    status = "lost";
                    break;
                }
            }

            m_dict[m].first++;
            if (status == "won") {
                m_dict[m].second++;
            }
        }

        double exp_m = 0;
        int total_m = 0;
        for (auto& p : m_dict) {
            exp_m += p.first * p.second.first;
            total_m += p.second.first;
        }
        exp_m /= total_m;

        double prob_q = static_cast<double>(std::accumulate(m_dict.begin(), m_dict.end(), 0,
            [](int value, const std::map<int, std::pair<int, int>>::value_type& p) { return value + p.second.second; })) / n_iter;

        exp_m_arr.push_back(exp_m);
        prob_q_arr.push_back(prob_q);
    }

    // Output results
    std::cout << "Distribution: ";
    for (int num : N) {
        std::cout << num << " ";
    }
    std::cout << std::endl;

    // Output results
    std::cout << "Average termination length: " << std::accumulate(exp_m_arr.begin(), exp_m_arr.end(), 0.0) / exp_m_arr.size() << std::endl;
    std::cout << "Average collusion probability: " << std::accumulate(prob_q_arr.begin(), prob_q_arr.end(), 0.0) / prob_q_arr.size() << std::endl;
}

int main(int argc, char *argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " n r_0 N_iter n_iter alpha" << std::endl;
        return 1;
    }

    int n = std::atoi(argv[1]);
    double r_0 = std::atof(argv[2]);
    int N_iter = std::atoi(argv[3]);
    int n_iter = std::atoi(argv[4]);
    double alpha = std::atof(argv[5]);

    simulate(n, r_0, N_iter, n_iter, alpha);

    return 0;
}