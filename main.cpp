
#include"MLE_NORMAL.cpp"
#include"MLE_WEIBULL.cpp"
#include"MLS_WEIBULL_AND_NORMAL.cpp"
#include"generators.cpp"
#include"AncdersonSmirCriter.cpp"
#include"Fisher_Student.cpp"
#include"wilcoxon.cpp"
#include<iostream>
#include <initializer_list>
#include <numeric>




void writeToFile(const std::string& filename, std::initializer_list<std::string> args) {
    std::ofstream file("output/" + filename + ".out");

    if (!file.is_open()) {
        throw std::runtime_error("Не удалось открыть файл: " + filename);
    }

    for (const auto& str : args) {
        file << str << "\n"; 
    }
    std::cout << "Strings successfully saved to file: " << filename << std::endl;
    file.close(); // Закрытие файла
}
void wilcoxon_do(int sample_size = 20,double stddev = 1,double mean1 = 50,double mean2 = 52.5) {
    using namespace std;
    cout << "6 лаба:" << endl;
    std::ofstream file("output/Wilcoxon.out");
    vector<double> sample1 = generateSample(sample_size, mean1, stddev);
    vector<double> sample2 = generateSample(sample_size, mean2, stddev);
    file << "size:" << sample_size << endl;
    file << "samples 1: ";
    for (double x : sample1) file << x << " ";
    file << endl;
    file << "samples 2: ";
    for (double x : sample2) file << x << " ";
    file << endl;

    pair<double, double> statistics = wilcoxonTwoSampleTest(sample1, sample2);
    double W1 = statistics.first;
    double W_star_1 = statistics.second;
    
    file << "Нормализованная статистика W1: " << W1 << endl;
    file << "Нормализованная статистика W*[1]: " << W_star_1 << endl;

    boost::math::normal dist_normal(0, 1); // Стандартное нормальное распределение
    double p_value = 2 * (1 - boost::math::cdf(dist_normal, abs(W1))); // Двусторонний тест



    file << "p-значение: " << p_value << endl;
    if (p_value < 0.05) {
        file << "Нулевая гипотеза отвергается: выборки статистически различны." << endl;
    }
    else {
        file << "Нулевая гипотеза не отвергается: выборки статистически не различны." << endl;
    }
    file.close();

}
int main(){
    setlocale(LC_ALL, "ru_RU");
    generateMLE_Normal("MLE_Normal");
    generateMLE_Weibull("MLE_Weibull",2,1.5);
    generateMLS_Normal("MLS_Normal");
    generateMLS_Weibull("MLS_Weibull",2,1.5);

    mle_normal();
    mle_weibull();
    mls(0); // normal
    mls(1); // weibull

    //---5
    vector<double> weibullSamples = generateSampleWeibull(20,2,1.5); 
    writeToFile("Anderson&SmirnovCriterions", {
       "smirnovCriterion: " + std::string(smirnovCriterion(weibullSamples)),
       "andersonDarlingCriterion: " + std::string(andersonDarlingCriterion(weibullSamples))
        });  

    //---4
    int    n1 = 30;
    double a1 = 4.975;
    double s1 = 0.1641;

    std::vector<double> x1 = generateSampleNormal( n1, a1, s1 );

    int    n2 = 20;
    double a2 = 3.5;
    double s2 = 0.777;

    double alpha = 0.05;

    std::vector<double> x2 = generateSampleNormal( n2, a2, s2 );

    std::vector<std::vector<double>> samples = {
        {x1},
        {x2}
    };

    writeToFile("FisherStudent", { // да да...
            "Fisher: " + std::string(Fisher(x1, x2) ? "Accepted" : "Rejected"),
            "Student: " + std::string(Student(x1, x2) ? "Accepted" : "Rejected"),
            "Grabs: " + std::string(Grabs(x2) ? "Accepted" : "Rejected"),
            "checkHomogeneity: " + std::string(checkHomogeneity(x1, x2, alpha) ? "Accepted" : "Rejected"),
            "bartlettTest: " + std::string(bartlettTest(samples, alpha) ? "Accepted" : "Rejected")
        });

    //---6
    wilcoxon_do(30);
}