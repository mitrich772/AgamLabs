#include <algorithm>

#include "default.hpp"

static constexpr double alpha = 0.05;

void printvector( const vector<double>& y )
{
    int n = y.size();
    for ( int i = 0; i < n; i++ ) {
        cout << i << ". " << y[i] << endl;
    }
}

std::vector<double> generateSampleNormal( int n, double a, double s )
{
    std::vector<double> znorm( n );
    std::vector<double> result( n );

    for ( int i = 0; i < n; i++ ) {
        double z = std::rand() / double( RAND_MAX );
        znorm.push_back( norm_ppf( z ) );
    }

    for ( int i = 0; i < n; i++ ) {
        result.push_back( a + znorm[i] * s );
    }

    return result;
}

bool Fisher( const std::vector<double>& s1, const std::vector<double>& s2 )
{
    double avg1 = std::accumulate( s1.begin(), s1.end(), 0 ) / static_cast<double>( s1.size() );
    double avg2 = std::accumulate( s2.begin(), s2.end(), 0 ) / static_cast<double>( s2.size() );

    double disp1 = 0.0;
    for ( auto i = 0lu; i < s1.size(); i++ ) {
        disp1 += pow( ( s1[i] - avg1 ), 2 );
    }
    disp1 /= s1.size() - 1;

    double disp2 = 0.0;
    for ( auto i = 0lu; i < s2.size(); i++ ) {
        disp2 += pow( ( s2[i] - avg2 ), 2 );
    }
    disp2 /= s2.size() - 1;

    double F = disp1 / disp2;

    double f1 = s1.size() - 1;
    double f2 = s2.size() - 1;

    if ( F < 1 ) {
        F  = disp2 / disp1;
        f1 = s2.size() - 1;
        f2 = s2.size() - 1;

        std::swap( disp1, disp2 );
    }

    double CritF = f_ppf( 1 - alpha, f1, f2 );

    if ( F <= CritF ) {
        return true;
    }
    else {
        return false;
    }
}

bool Student( const std::vector<double>& s1, const std::vector<double>& s2 )
{
    double avg1 = std::accumulate( s1.begin(), s1.end(), 0 ) / static_cast<double>( s1.size() );
    double avg2 = std::accumulate( s2.begin(), s2.end(), 0 ) / static_cast<double>( s2.size() );

    double disp1 = 0.0;
    for ( auto i = 0lu; i < s1.size(); i++ ) {
        disp1 += pow( ( s1[i] - avg1 ), 2 );
    }
    disp1 /= s1.size() - 1;

    double disp2 = 0.0;
    for ( auto i = 0lu; i < s2.size(); i++ ) {
        disp2 += pow( ( s2[i] - avg2 ), 2 );
    }
    disp2 /= s2.size() - 1;

    double f1 = s1.size() - 1;
    double f2 = s2.size() - 1;

    double t, f, c, CritT;
    c = ( disp1 / ( f1 + 1 ) ) / ( disp1 / ( f1 + 1 ) + disp2 / ( f2 + 1 ) );

    f = 1. / ( c * c / f1 + ( 1 - c ) * ( 1 - c ) / f2 );

    t = ( avg1 - avg2 ) / sqrt( ( disp1 / ( f1 + 1 ) ) + ( disp2 / ( f2 + 1 ) ) );

    CritT = t_ppf( 1 - alpha / 2, f );

    if ( fabs( t ) <= CritT ) {
        return true;
    }
    else {
        return false;
    }
}

bool Grabs( std::vector<double>& y )
{
    double avg = std::accumulate( y.begin(), y.end(), 0 ) / static_cast<double>( y.size() );

    double cko = 0.0;
    for ( auto i = 0lu; i < y.size(); i++ ) {
        cko += pow( ( y[i] - avg ), 2 );
    }
    cko /= y.size() - 1;
    cko = sqrt( cko );

    cout << "avg = " << avg << " cko = " << cko << endl;

    y[y.size() - 10] = ( 8.23 );

    std::vector<double> un;

    for ( auto i = 0lu; i < y.size(); i++ ) {
        un.push_back( fabs( avg - y[i] ) / cko );
    }

    double u = 0.0;
    for ( auto i = 0lu; i < y.size(); i++ ) {
        double temp = fabs( y[i] - avg ) / cko;
        u           = max( u, temp );
    }

    std::cout << "u = " << u << '\n';

    int n = y.size();

    students_t_distribution<double> dist( n - 2 );
    double                          t  = quantile( dist, 1 - alpha / ( 2 * n ) );
    double                          ua = ( n - 1 ) * sqrt( t / ( n * ( n - 2 + t ) ) );

    for ( int i = 0; i < n; i++ ) {
        if ( un[i] <= ua ) {
        }
        else {
            std::cout << "есть выброс - " << y[i] << '\n';
            return false;
        }
    }

    return true;
}

// Функция для вычисления выборочной дисперсии
double calculateVariance(const std::vector<double>& data) {
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    double variance = 0.0;
    for (double value : data) {
        variance += (value - mean) * (value - mean);
    }
    return variance / (data.size() - 1);
}

// Функция для проверки однородности дисперсий
bool checkHomogeneity(const std::vector<double>& sample1, const std::vector<double>& sample2, double alpha) {
    // Вычисляем выборочные дисперсии
    double variance1 = calculateVariance(sample1);
    double variance2 = calculateVariance(sample2);

    // Вычисляем дисперсионное отношение F
    double F = (variance1 > variance2) ? (variance1 / variance2) : (variance2 / variance1);

    // Вычисляем степени свободы
    int df1 = sample1.size() - 1; 
    int df2 = sample2.size() - 1; 

    // Получаем критическое значение F
    boost::math::fisher_f_distribution<double> F_dist(df1, df2);
    double F_critical = boost::math::quantile(F_dist, 1 - alpha); 

    // Сравниваем F с критическим значением
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Calculated F: " << F << std::endl;
    std::cout << "Critical F: " << F_critical << std::endl;

    if (F <= F_critical) {
        return true;
    } else {
        return false;
    }
}

// Функция для проверки однородности выборок по критерию Бартлетта
bool 
bartlettTest(const std::vector<std::vector<double>>& samples, double alpha) {
    int k = samples.size();
    int N = 0;
    std::vector<double> variances(k);
    std::vector<int> sampleSizes(k);

    for (int i = 0; i < k; ++i) {
        sampleSizes[i] = samples[i].size();
        N += sampleSizes[i];
        variances[i] = calculateVariance(samples[i]);
    }

    double S2_p = 0.0;
    for (int i = 0; i < k; ++i) {
        S2_p += (sampleSizes[i] - 1) * variances[i];
    }
    S2_p /= (N - k);

    double T = 0.0;
    for (int i = 0; i < k; ++i) {
        T += (sampleSizes[i] - 1) * std::log(variances[i]);
    }
    T = (N - k) * std::log(S2_p) - T;

    int df = k - 1;
    boost::math::chi_squared_distribution<double> chi_dist(df);
    double chi_critical = boost::math::quantile(chi_dist, 1 - alpha);

    return (T > chi_critical); 
}




// int main()
// {
//     int    n1 = 30;
//     double a1 = 4.975;
//     double s1 = 0.1641;

//     std::vector<double> x1 = generateSampleNormal( n1, a1, s1 );

//     int    n2 = 20;
//     double a2 = 3.5;
//     double s2 = 0.777;

//     double alpha = 0.05;

//     std::vector<double> x2 = generateSampleNormal( n2, a2, s2 );

//     std::vector<std::vector<double>> samples = {
//         {x1},
//         {x2}
//     };

//     bartlettTest(samples, alpha);

//     std::cout << ( Fisher( x1, x2 ) ? "гипотеза принята" : "гипотеза отклонена" ) << '\n';
//     std::cout << ( Student( x1, x2 ) ? "гипотеза принята" : "гипотеза отклонена" ) << '\n';
//     std::cout << ( Grabs( x2 ) ? "гипотеза принята" : "гипотеза отклонена" ) << '\n';
//     std::cout << (checkHomogeneity(x1, x2, alpha) ? "гипотеза принята" : "гипотеза отклонена" ) << '\n';
//     std::cout << (bartlettTest(samples, alpha) ? "гипотеза принята" : "гипотеза отклонена" ) << '\n';


//     return 0;
// }
