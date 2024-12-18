#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <locale>
#include <cmath>
#include <utility>
#include "default.hpp"

using namespace std;

vector<double> generateSample(size_t size, double mean, double stddev) {
    vector<double> sample(size);
    double z;
    for (size_t i = 0; i < size; ++i) { 
        z = norm_ppf(rand() / double(RAND_MAX));
        sample[i] = mean + z * stddev;
    }

    return sample;
}

pair<double, double> wilcoxonTwoSampleTest(const vector<double>& sample1, const vector<double>& sample2) {
    size_t n = sample1.size();
    size_t m = sample2.size();

    if (n != m) {
        cerr << "Размеры выборок должны быть одинаковыми!" << endl;
        return { -1, -1 };
    }

    vector<double> differences;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            differences.push_back(sample1[i] - sample2[j]);
        }
    }

    vector<double> absDifferences;
    for (double diff : differences) {
        absDifferences.push_back(abs(diff));
    }

    vector<double> sortedAbsDifferences = absDifferences;
    sort(sortedAbsDifferences.begin(), sortedAbsDifferences.end());

    vector<double> ranks(absDifferences.size());
    for (size_t i = 0; i < absDifferences.size(); ++i) {
        ranks[i] = find(sortedAbsDifferences.begin(), sortedAbsDifferences.end(), absDifferences[i]) - sortedAbsDifferences.begin() + 1;
    }

    double W = 0;
    for (size_t i = 0; i < differences.size(); ++i) {
        if (differences[i] > 0) {
            W += ranks[i];
        }
    }
    double alpha = 0.975;
    double expected = (n * (m + n + 1)) / 2.0;
    double variance = (n * m * (m + n + 1)) / 12.0;
    double stddev = sqrt(variance);


    double t = t_ppf(m + n - 2, alpha);
    double z = norm_ppf(alpha);

    double w_right = 0.5 * (t + z);
    double w_left = -w_right;
    //cout << w_left << " " << w_right;

    double W1 = (W - expected + 0.5) / stddev;

    double W_star_1 = (W1 / 2) * (1 + sqrt(abs((m + n - 2) / (m + n - 1 - W1 * W1))));
    if (W_star_1 > w_right && W_star_1 < w_left) {
        //cout << "\nH0+";
    }
    else {
        //cout << "\nH0-";
    }
    return { W1, W_star_1 };
}


