#include "default.hpp"

void mle_weibull()
{
    string   str;
    int      tempcens, pcount, paramscount, xcount;
    double   beta, eps, tempx, avg, cko, delta, tup, tlow, s1, s2, c, temp;
    double * p, *fcum, *ycum;
    double** v;

    ifstream input( "input/MLE_Weibull.inp" );
    ofstream output( "output/MLE_Weibull.out" );

    input >> str;
    input >> nesm.n;
    input >> str;
    input >> beta;
    input >> str;
    input >> eps;
    input >> str;
    for ( int i = 0; i < nesm.n; i++ ) {
        input >> tempx;
        nesm.x.push_back( tempx ); //
    }
    input >> str;
    for ( int i = 0; i < nesm.n; i++ ) {
        input >> tempcens;
        nesm.r.push_back( tempcens );
    }
    input >> str;
    input >> pcount;
    input >> str;
    p = new double[pcount];
    for ( int i = 0; i < pcount; i++ ) {
        input >> p[i];
    }

    input.close();

    paramscount = 2;
    vector<double> params;
    vector<double> xcens;
    xcount = 0;
    avg    = 0;
    cko    = 0;
    for ( int i = 0; i < nesm.n; i++ ) {
        xcount += ( 1 - nesm.r[i] );
        if ( xcount ) {
            xcens.push_back( log( nesm.x[i] ) ); //
        }
        avg += ( 1 - nesm.r[i] ) * ( log( nesm.x[i] ) ); //
    }

    avg /= xcount;

    for ( int i = 0; i < nesm.n; i++ ) {
        cko += ( 1 - nesm.r[i] ) * ( log( nesm.x[i] ) - avg ) * ( log( nesm.x[i] ) - avg ); //
    }

    cko /= xcount;

    params.push_back( 1. / cko );

    int q = 0;
    int iterations;

    if ( xcount != nesm.n ) {
        neldermead( params, eps, WeibullMinFunction );
        q = WeibullMinFunction( params );
    }

    s1 = 0;
    s2 = 0;
    c  = 0;
    for ( int i = 0; i < nesm.n; i++ ) {
        s1 += ( 1 - nesm.r[i] ) * pow( ( nesm.x[i] ), params[0] ); //
        s2 += nesm.r[i] * pow( ( nesm.x[i] ), params[0] );         //
    }
    temp = ( s1 + s2 ) / xcount;
    cout << temp;
    c = exp( log( temp ) / params[0] );
    params.push_back( c );
    // Ковариационная матрица

    v = new double*[paramscount];
    for ( int i = 0; i < paramscount; i++ ) {
        v[i] = new double[paramscount];
    }
    for ( int i = 0; i < paramscount; i++ ) {
        for ( int j = 0; j < paramscount; j++ ) {
            v[i][j] = 0;
        }
    }

    CovMatrixMleW( nesm.n, nesm.x, nesm.r, params[1], params[0], v );

    double* x = new double[nesm.n];
    int*    r = new int[nesm.n];
    for ( int i = 0; i < nesm.n; i++ ) {
        x[i] = nesm.x[i]; //exp
    }
    for ( int i = 0; i < nesm.n; i++ ) {
        r[i] = nesm.r[i];
    }
    fcum = new double[xcount];
    ycum = new double[xcount];

    cum( nesm.n, x, r, xcount, fcum, ycum );

    //Доверительные границы

    double* t     = new double[pcount];
    double* xplow = new double[pcount];
    double* xpup  = new double[pcount];
    double* zp    = new double[pcount];
    double* xp    = new double[pcount];

    // v[0][0] = 0.0579915;
    // v[0][1] = -0.00846277;
    // v[1][0] = -0.00846277;
    // v[1][1] = 0.0395316;

    for ( int i = 0; i < pcount; i++ ) {
        zp[i] = log( log( 1. / ( 1. - p[i] ) ) );
        delta = zp[i] * sqrt( nesm.n );

        nctWeibull( beta, nesm.n - 1, v[0][0], v[1][1], v[0][1], delta, tlow, tup );

        xp[i]    = log( params[1] ) + 1. / params[0] * zp[i];
        xplow[i] = log( params[1] ) + 1. / params[0] * tlow / sqrt( nesm.n );
        xpup[i]  = log( params[1] ) + 1. / params[0] * tup / sqrt( nesm.n );
    }

    output << "Кол-во наблюдений" << endl;
    output << nesm.n << endl;
    output << "Наблюдения" << endl;
    for ( int i = 0; i < nesm.n; i++ ) {
        output << nesm.x[i] << "; ";
    }
    output << "\n";
    output << "Цензурирование" << endl;
    for ( int i = 0; i < nesm.n; i++ ) {
        output << nesm.r[i] << "; ";
    }
    output << "\n";
    output << "cp*=" << setprecision( 12 ) << fixed << avg << endl;
    output << "cko*=" << fixed << cko << endl;
    output << "Q = " << q << endl;
    output << "Iterations = " << iterations << endl;
    output << "b=" << params[0] << endl;
    output << "c=" << params[1] << endl;
    output << "P" << endl;
    output << setprecision( 3 );
    for ( int i = 0; i < pcount; i++ ) {
        output << p[i] << "; ";
    }
    output << setprecision( 12 );
    output << "\n";

    output << "v11=" << v[0][0] << endl;
    output << "v12=" << v[0][1] << endl;
    output << "v21=" << v[1][0] << endl;
    output << "v22=" << v[1][1] << endl;

    output << "XpLOW" << endl;
    for ( int i = 0; i < pcount; i++ ) {
        output << xplow[i] << "; ";
    }
    output << "\n";
    output << "XP" << endl;
    for ( int i = 0; i < pcount; i++ ) {
        output << xp[i] << "; ";
    }
    output << "\n";
    output << "XpUP" << "\n";
    for ( int i = 0; i < pcount; i++ ) {
        output << xpup[i] << "; ";
    }
    output << "\n";
    output << "ZP+5.0" << endl;
    for ( int i = 0; i < pcount; i++ ) {
        output << zp[i] + 5. << "; ";
    }
    output << "\n";

    output << "X cens" << endl;
    for ( int i = 0; i < xcount; i++ ) {
        output << log( ycum[i] ) << "; ";
    }
    output << "\n";

    output << "ZX" << endl;
    for ( int i = 0; i < xcount; i++ ) {
        output << 5.0 + log( log( 1. / ( 1. - fcum[i] ) ) ) << "; ";
    }
    output << "\n";


    output.close();
    delete[] xplow, xpup, zp, xp;
    nesm.r.clear(), nesm.x.clear(), xcens.clear(), params.clear();
}

// int main()
// {
//     mle_weibull();
// }
