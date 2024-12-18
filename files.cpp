#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>
// #include "default.hpp"
#include <random>
using namespace std;

int main()
{
    // Параметры распределения Вейбулла
    double shape = 2.1231445; // параметр формы k
    double scale = 4.1236753; // параметр масштаба λ

    // Создаем генератор случайных чисел с использованием стандартного генератора
    random_device rd;          // Источник случайности
    mt19937       gen( rd() ); // Стандартный генератор (Mersenne Twister)

    // Определяем распределение Вейбулла
    weibull_distribution<double> weibull_dist( shape, scale );

    // Количество элементов в выборке
    size_t         sample_size = 20;
    vector<double> sample( sample_size );

    // Генерируем выборку
    cout << "Сгенерированные данные (распределение Вейбулла):\n";
    for ( size_t i = 0; i < sample_size; ++i ) {
        sample[i] = weibull_dist( gen );
    }
    sort( sample.begin(), sample.end() );
    cout << endl;
    for ( auto i = 0lu; i < sample_size; i++ ) {
        cout << fixed << setprecision( 4 ) << sample[i] << " ";
    }
    return 0;
}
