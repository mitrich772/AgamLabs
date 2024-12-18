#pragma once

#include"default.hpp";
#include<string>
#include<cmath>
#include<fstream>
#include<iostream>
#include<ctime>
#include <algorithm>
#include <cstring>
#include <vector>
#include <random>
#include <math.h>
#include <iomanip>
#include <stdexcept>


std::vector<double> generateSampleWeibull(int n, double scale = 12045.03, double shape = 1.668) {
    if (n <= 0 || scale <= 0 || shape <= 0) {
        throw std::invalid_argument("n, scale, and shape must be positive.");
    }

    std::vector<double> result;
    result.reserve(n);

   
    std::random_device rd;                              
    std::mt19937 gen(rd());                             
    std::weibull_distribution<> dist(shape, scale);     

    // Генерация выборки
    for (int i = 0; i < n; ++i) {
        result.push_back(dist(gen));
    }

    return result;
}
/**
 * @brief Функция для проверки соответствия выборки нормальному распределению
 * @param sample выборка по вейболовскому распределению
 * @param alpha Уровень значимости
 */
std::string smirnovCriterion(std::vector<double> sample, double alpha = 0.05) {

    const int size = sample.size();
    // Сортируем выборку в порядке возрастания
    std::vector<double> sortedSample = sample;
    std::sort(sortedSample.begin(), sortedSample.end());

    // Вычисляем накопленные частоты W(x_i)
    std::vector<double> W;

    W.reserve(size);
    for (unsigned int i = 0; i < size; ++i) W.push_back((i + 0.5) / size);

    // Вычисляем статистику критерия Смирнова
    double omegaSquared = 1.0 / (12 * size);
    for (unsigned int i = 0; i < size; ++i) {
        double F_xi = norm_cdf(sortedSample[i]);
        omegaSquared += std::pow(F_xi - W[i], 2);
    }

    // Критическое значение
    double omegaAdjusted = omegaSquared * (1.0 + 0.5 / size);

    // Критические значения для разных уровней значимости
    double criticalValue;
    if (alpha == 0.15) {
        criticalValue = 0.091;
    } else if (alpha == 0.10) {
        criticalValue = 0.104;
    } else if (alpha == 0.05) {
        criticalValue = 0.126;
    } else if (alpha == 0.01) {
        criticalValue = 0.178;
    } else {
        std::cerr << "Unsupported significance level." << std::endl;
    }
    //std::cout << "omegaAdjusted = " << omegaAdjusted << std::endl;
    // Проверяем гипотезу
    if (omegaAdjusted <= criticalValue) {
        return("Null hypothesis is accepted at significance level alpha = " + to_string(alpha));
    } else {
         return("Null hypothesis is rejected at significance level alpha = " + to_string(alpha));
    }
}

/**
 * @brief Функция для проверки соответствия выборки нормальному распределению
 * @param sample выборка по вейболовскому распределению
 * @param alpha Уровень значимости
 */
std::string andersonDarlingCriterion(std::vector<double> sample, double alpha = 0.05) {
    const int size = sample.size();

    // Сортируем выборку в порядке возрастания
    std::vector<double> sortedSample = sample;
    std::sort(sortedSample.begin(), sortedSample.end());

    // Вычисляем статистику критерия Андерсона-Дарлинга
    double A2 = 0.0;
    for (int i = 1; i <= size; ++i) {
        const double F_xi = norm_cdf(sortedSample[i - 1]);
        A2 += (2 * i - 1) * std::log(F_xi) + (2 * size - 2 * i + 1) * std::log(1 - norm_cdf(sortedSample[size - i]));
    }
    A2 = -size - (A2 / size);

    // Корректируем статистику
    const double A2_adjusted = (A2 - 0.7 / size) * (1.0 + 3.6 / size - 8.0 / (size * size));

    // Критические значения для разных уровней значимости
    double criticalValue;
    if (alpha == 0.15) {
        criticalValue = 0.576;
    } else if (alpha == 0.10) {
        criticalValue = 0.656;
    } else if (alpha == 0.05) {
        criticalValue = 0.787;
    } else if (alpha == 0.01) {
        criticalValue = 1.092;
    } else {
        std::cerr << "Unsupported significance level." << std::endl;
    }

    //std::cout << "A2_adjusted = " << A2_adjusted << std::endl;
    // Проверяем гипотезу
    if (A2_adjusted <= criticalValue) {
        return("Null hypothesis is accepted at significance level alpha = " + to_string(alpha));
    } else {
         return("Null hypothesis is rejected at significance level alpha = " + to_string(alpha));
    }
}