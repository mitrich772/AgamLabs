#include "default.hpp"

void mle_normal()
{
    string   str;
    int      tempcens, pcount, paramscount, xcount;
    double   beta, eps, tempx, avg, cko, delta, tup, tlow;
    double * p, *fcum, *ycum;
    double** v;

    ifstream input( "input/MLE_Normal.inp" );
    ofstream output( "output/MLE_Normal.out" );

    input >> str;
    input >> nesm.n;
    input >> str;
    input >> beta;
    input >> str;
    input >> eps;
    input >> str;
    for ( int i = 0; i < nesm.n; i++ ) {
        input >> tempx;
        nesm.x.push_back( tempx );
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

    avg    = 0;
    cko    = 0;
    xcount = 0;
    for ( int i = 0; i < nesm.n; i++ ) {
        xcount += ( 1 - nesm.r[i] );
        if ( xcount ) {
            xcens.push_back( nesm.x[i] );
        }
        avg += ( 1 - nesm.r[i] ) * nesm.x[i];
        cko += ( 1 - nesm.r[i] ) * nesm.x[i] * nesm.x[i];
    }

    avg /= xcount;
    cko = sqrt( ( cko - avg * avg * xcount ) / ( xcount - 1 ) );

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

    params.push_back( avg );
    params.push_back( cko );

    int q = 0;
    int iterations;

    if ( xcount != nesm.n ) {
        iterations = neldermead( params, eps, NormalMinFunction );
        q          = NormalMinFunction( params );
        CovMatrixMleN( nesm.n, nesm.x, nesm.r, params[0], params[1], v );
    }

    //Доверительные границы

    double* t     = new double[pcount];
    double* xplow = new double[pcount];
    double* xpup  = new double[pcount];
    double* zp    = new double[pcount];
    double* xp    = new double[pcount];

    for ( int i = 0; i < pcount; i++ ) {
        // ZP - квантиль уровня P нормального распределениия
        zp[i] = norm_ppf( p[i] );
        // XP - квантиль нормального распределения, получается сдвигом X = s * Z + a
        xp[i] = params[0] + zp[i] * params[1];
        // Параметр нецентральности (2.82)
        delta = zp[i] * sqrt( nesm.n );

        if ( xcount == nesm.n ) {
            t[i] = nct_ppf( beta, nesm.n - 1, delta );
            //noncentralt_inv(beta,nesm.n-1,delta);
        }
        else {
            lmtaprn( beta, nesm.n - 1, v[0][0], v[1][1], v[0][1], delta, tlow, tup );
            xpup[i]  = params[0] + tup * params[1] / sqrt( nesm.n );
            xplow[i] = params[0] + tlow * params[1] / sqrt( nesm.n );
        }
    }

    if ( xcount == nesm.n ) {
        for ( int i = 0; i < xcount; i++ ) {
            xpup[i] = params[0] + t[i] * params[1] / sqrt( nesm.n );
        }
        for ( int i = 0; i < xcount; i++ ) {
            xplow[i] = params[0] - t[xcount - i - 1] * params[1] / sqrt( nesm.n );
        }
    }


    fcum = new double[nesm.n];
    ycum = new double[nesm.n];
    double* xq;
    int*    rq;
    xq = new double[nesm.n];
    rq = new int[nesm.n];
    for ( int i = 0; i < nesm.n; i++ ) {
        xq[i] = nesm.x[i];
        rq[i] = nesm.r[i];
    }
    cum( nesm.n, xq, rq, xcount, fcum, ycum );

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
    output << "cp=" << params[0] << endl;
    output << "cko=" << params[1] << endl;
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
        output << ycum[i] << "; ";
    }
    output << "\n";

    output << "ZX" << endl;
    for ( int i = 0; i < xcount; i++ ) {
        output << 5 + norm_ppf( fcum[i] ) << "; ";
    }
    output << "\n";


    output.close();
    delete[] xplow, xpup, zp, xp;
    nesm.r.clear(), nesm.x.clear(), xcens.clear(), params.clear();
}

// int main()
// {
//     mle_normal();
// }
