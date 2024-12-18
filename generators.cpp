#pragma once
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





/**
 * @brief Функция для генерации выборки из распределения Вейбулла
 * @param size Размер выборки
 * @param scale Параметр масштаба
 * @param shape Параметр формы
 * @param filename Имя файла для сохранения
 */
void generateWeibull(const std::string& filename, int size ,double scale = 12045.03, double shape = 1.668,double beta = 0.95) {
    std::ofstream outFile(filename);
    std::vector<double> sample;
    sample.reserve(size); // Генератор случайных чисел

    std::random_device rd;  // Получаем случайное число из устройства
    std::mt19937 gen(rd()); // Инициализируем генератор
    std::weibull_distribution<> d(shape, scale); // Вейбулла с параметрами shape и scale
    // Генерация данных
    for (int i = 0; i < size; ++i) {
        sample.push_back(d(gen)); // Генерируем данные и добавляем в вектор
    }
    
    if (!outFile.is_open()) {
        std::cerr << "Error opening file for writing sample: " << filename << std::endl;
        if (errno == ENOENT) {
            std::cerr << "File or directory does not exist." << std::endl;
        } else if (errno == EACCES) {
            std::cerr << "Permission denied." << std::endl;
        } else {
            std::cerr << "Unknown error: " << strerror(errno) << std::endl;
        }
        return;
    }

    // Приведение вывода под формат функции MLE_Normal()
    // Запись размера выборки
    outFile << "Samples_size\n";
    outFile << size << "\n";

    // Запись параметра beta
    outFile << "beta\n";
    outFile << beta << "\n";

    // Запись шага минимизации
    // outFile << "step_of_minimization\n";
    // outFile << 0.5 << "\n";

    // Запись точности
    outFile << "eps_output\n";
    outFile << "1.e-15\n";

    // Лимит итераций
    // outFile << "lim_of_iteration\n";
    // outFile << 500 << "\n";

    // Генерация выборки и запись данных
    outFile << "Data\n";

    std::sort(sample.begin(), sample.end());

    for (double value : sample) {

        outFile << value << " ";
    }
    outFile << "\n";

    // Запись информации о цензурировании
    outFile << "Censorizes\n";
    for (int i = 0; i < size; ++i) {
        int censor = (i % 2 == 0) ? 1 : 0; // Пример: цензурирование каждого второго значения
        outFile << censor;
        if (i < size - 1) {
            outFile << " ";
        }
    }
    outFile << "\n";

    // Запись значения kp
    outFile << "kp\n";
    outFile << 22 << "\n";

    // Запись значений P
    outFile << "P\n";
    outFile << "0.025 0.075 0.125 0.175 0.225 0.275 0.325 0.375 0.425 0.475 0.525\n";
    outFile << "0.575 0.625 0.675 0.725 0.775 0.825 0.875 0.925 0.975 0.99 0.995\n";

    outFile.close();
    std::cout << "Sample successfully saved to file: " << filename << std::endl;
}
/**
 * @brief Функция для генерации выборки из распределения Вейбулла
 * @param filename Имя файла для сохранения
 */
void generateMLE_Weibull(const std::string& filename,double scale = 12045.03, double shape = 1.668) {
    int size; // Размер выборки
    //double scale; // Параметр масштаба
    //double shape; // Параметр формы


    std::cout << "Enter sample size "<< filename << ": " ;
    std::cin >> size;
    // std::cout << "Enter scale parameter: ";
    // std::cin >> scale;
    // std::cout << "Enter form parameter: ";
    // std::cin >> shape;

    generateWeibull("input/" + filename + ".inp", size, scale, shape);
}
/**
 * @brief Функция для чтения данных из файла
 * @param filename Имя файла с выборкой
 * @param sample Выборка
 * @param censor Цензурирование
 * @param beta Параметр beta
 * @param step Шаг минимизации
 * @param eps Точность
 * @param lim Лимит итераций
 * @param kp Количество значений P
 * @param p Массив значений P
 */
void readSampleFromFile(const std::string& filename, std::vector<double>& sample, std::vector<int>& censor, double& beta, double& step, double& eps, double& lim, int& kp, double*& p) {
    std::string s1; // Мусорная переменная для считывания текста
    unsigned int size; // Размер выборки
    double z;

    // Чтение файла с выборкой
    std::ifstream input("Inp/" + filename + ".inp");
    if (!input.is_open()) {
        std::cerr << "Error opening file for reading sample: " << "Inp/" + filename + ".inp" << std::endl;
        if (errno == ENOENT) {
            std::cerr << "File or directory does not exist." << std::endl;
        } else if (errno == EACCES) {
            std::cerr << "Permission denied." << std::endl;
        } else {
            std::cerr << "Unknown error: " << strerror(errno) << std::endl;
        }
        return;
    }
    std::getline(input, s1);      //Text
    input >> size;                       //Sample size
    std::getline(input, s1);
    input >> beta;                      //Beta
    std::getline(input, s1);
    input >> step;                      //Step of minimization
    std::getline(input, s1);
    input >> eps;                       //Eps output
    std::getline(input, s1);
    input >> lim;                       //Lim of iteration
    std::getline(input, s1);

    sample.reserve(size);
    for (unsigned int i = 0; i < size; ++i) {
        input >> z;
        sample.push_back(z);
    }
    std::getline(input, s1); //text

    censor.reserve(size);
    for (unsigned int i = 0; i < size; ++i) {
        input >> z;
        censor.push_back(static_cast<int>(z));
    }
    std::getline(input, s1);

    input >> kp;
    p = new double[kp];

    std::getline(input, s1);
    for (unsigned int i = 0; i < kp; ++i) input >> p[i];
    input.close();
}



/**
 * @brief Функция для генерации выборки из нормального распределения и записи в файл
 * @param filename Имя файла для записи
 * @param sampleSize Размер выборки
 */
void generateNormal(const std::string& filename, int sampleSize, double beta = 0.95 ,double mean = 5.0, double stddev = 1.0) {
    std::vector<double> sample;
    sample.reserve(sampleSize);

    std::random_device rd;  // Используется для инициализации генератора случайных чисел (настоящее случайное значение)
    std::mt19937 gen(rd()); // Инициализируем генератор
    
    std::normal_distribution<> d(mean, stddev);
    for (int i = 0; i < sampleSize; ++i) {
        sample.push_back(d(gen)); // Генерируем данные и добавляем в вектор
    }
    // Сортировка данных по возрастанию
    std::sort(sample.begin(), sample.end());

    //std::uniform_real_distribution<> dis(0.0, 1.0);    // Определяем равномерное распределение от 0 до 1
    // for (unsigned int i = 0; i < sampleSize; ++i) {
    //     double u1 = dis(gen);
    //     double u2 = dis(gen);
    //     double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(2 * 3.1415 * u2);
    //     sample.push_back(mean + z * stddev);
    // }

    // Запись выборки в файл
    std::ofstream output(filename);
    if (!output.is_open()) {
        std::cerr << "Error opening file for writing: " <<  filename  << std::endl;
        return;
    }

    output << "Samples_size\n";
    output << sampleSize << "\n";
    output << "beta\n";
    output << beta <<"\n";
    // output << "step_of_minimization\n";
    // output << "0.5\n";
    output << "eps_output\n";
    output << "1.e-15\n";
    // output << "lim_of_iteration\n";
    // output << "500\n";
    output << "Data\n";
    for (const auto& value : sample) {
        output << value << " ";
    }
    output << "\n";
    output << "Censorizes\n";
    for (unsigned int i = 0; i < sampleSize; ++i) {
        output << (i % 2) << " ";  // Чередуем 1 и 0 для примера цензурирования
    }
    output << "\n";
    output << "kp\n";
    output << "15\n";
    output << "P\n";
    output << "0.005 0.01 0.025 0.05 0.1 0.2 0.3 0.5 0.7 0.8 0.9 0.95 0.975 0.99 0.995\n";

    output.close();
    std::cout << "Sample successfully saved to file: " << filename << std::endl;
}
void write_to_file_MLS(const std::string& filename, int Samples_size, const std::vector<double>& Data, const std::vector<int>& Censorizes) {
    std::string filepath = "input/" + filename + ".inp";

    // Удаляем файл, если он существует
    if (std::remove(filepath.c_str()) != 0) {
        std::cerr << "Ошибка при удалении файла: " << filepath << std::endl;
    }


    std::ofstream out1(filepath); // Открываем файл для записи
    if (!out1) {
        std::cerr << "Ошибка при открытии файла!" << std::endl;
        return;
    }

    // Запись параметров
    out1 << "Samples_size\n" << Samples_size << "\n";

    // Запись данных
    out1 << "Data\n";
    for (const auto& value : Data) {
        out1 << std::fixed << std::setprecision(9) << value << " ";
    }
    out1 << "\n";

    // Запись цензоров
    out1 << "Censorizes\n";
    for (const auto& censor : Censorizes) {
        out1 << censor << " ";
    }
    out1 << "\n";

    // Запись дополнительных параметров (kp и P)
    int kp = 15; // Пример значения kp
    out1 << "kp\n" << kp << "\n";
    out1 << "P\n";
    std::vector<double> P = {0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.975, 0.99, 0.995};
    for (const auto& p_value : P) {
        out1 << p_value << " ";
    }
    out1 << "\n";
    std::cout << "Sample successfully saved to file: " << filepath << std::endl;

    out1.close(); // Закрываем файл
}
/**
 * @brief Функция для получения количества экспереметов, генерации выборки из нормального распределения и записи в файл 
 * @param filename Имя файла для записи
 */
void generateMLE_Normal(const std::string &filename) {
    int size; // Размер выборки
    //double mean; // Мат ожидание
    //double stddev; // Стандартное отклонение


    std::cout << "Enter sample size " << filename << ": "; ;
    std::cin >> size;
    // std::cout << "Enter mean: ";
    // std::cin >> mean;
    // std::cout << "Enter standard deviation: ";
    // std::cin >> stddev;

    generateNormal("input/" + filename + ".inp", size);
}
void generateMLS_Normal(const std::string& filename, double mean = 5.0, double stddev = 1.0) {
    // Генератор случайных чисел
    int sampleSize; // Размер выборки
    //double mean; // Мат ожидание
    //double stddev; // Стандартное отклонение


    std::cout << "Enter sample size " << filename << ": "; ;
    std::cin >> sampleSize;
    std::vector<double> Data;
    std::vector<int> Censorizes;

    std::random_device rd;  // Получаем случайное число из устройства
    std::mt19937 gen(rd()); // Инициализируем генератор
    std::normal_distribution<> d(mean, stddev); // Нормальное распределение с средним 5.0 и стандартным отклонением 1.0
    // Здесь был белка

    
    // Генерация данных
    for (int i = 0; i < sampleSize; ++i) {
        Data.push_back(d(gen)); // Генерируем данные и добавляем в вектор
    }

    // Сортировка данных по возрастанию
    std::sort(Data.begin(), Data.end());

    // Генерация цензоров (например, случайно 0 или 1)
    for (int i = 0; i < sampleSize; ++i) {
        Censorizes.push_back(rand() % 2); // Генерируем случайные цензоры
    }
    write_to_file_MLS(filename, sampleSize, Data, Censorizes);
}
void generateMLS_Weibull(const std::string& filename, double scale = 12045.03, double shape = 1.668){
    // Генератор случайных чисел
    int sampleSize; // Размер выборки
    //double mean; // Мат ожидание
    //double stddev; // Стандартное отклонение


    std::cout << "Enter sample size " << filename << ": "; ;
    std::cin >> sampleSize;
    std::vector<double> Data;
    std::vector<int> Censorizes;

    std::random_device rd;  // Получаем случайное число из устройства
    std::mt19937 gen(rd()); // Инициализируем генератор
    std::weibull_distribution<> d(shape, scale); // Нормальное распределение с средним 5.0 и стандартным отклонением 1.0
    // Здесь был белка

    
    // Генерация данных
    for (int i = 0; i < sampleSize; ++i) {
        Data.push_back(d(gen)); // Генерируем данные и добавляем в вектор
    }

    // Сортировка данных по возрастанию
    std::sort(Data.begin(), Data.end());

    // Генерация цензоров (например, случайно 0 или 1)
    for (int i = 0; i < sampleSize; ++i) {
        Censorizes.push_back(rand() % 2); // Генерируем случайные цензоры
    }
    write_to_file_MLS(filename, sampleSize, Data, Censorizes);
}
