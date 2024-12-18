
#include "default.hpp"
void MleastSquare_weight( int n, int k, double** x, double** y, double** v, double**& db, double**& b, double*& yr )
{
    int    i, j;
    double s, **xt;

    xt = new double*[n];
    for ( i = 0; i < n; i++ )
        xt[i] = new double[n];
    v  = InverseMatrix( v, n );                                   // (n x n)
    xt = MultiplyMatrix( k, n, n, n, TransMatrix( n, k, x ), v ); // (k x n)

    db = InverseMatrix( MultiplyMatrix( k, n, n, k, xt, x ), k );               //covariance matrix factors (k x k)
    b  = MultiplyMatrix( k, n, n, 1, MultiplyMatrix( k, k, k, n, db, xt ), y ); // coef b[k][0]
    for ( i = 0; i < n; i++ ) {
        s = 0;
        for ( j = 0; j < k; j++ )
            s += b[j][0] * x[i][j];
        yr[i] = s;
    }
    delete[] xt;
}
void mls(int a)
{

    int     i, j, n, k, ku, kp, km, kmm, nc;
    string  ff;
    string  s1;
    double *xz, **x, **y, *xplow, *xpup, **v, *yr, **b, **db, *zp, *p, *xp, *fcum, *ycum, tlow, tup, delta;
    double  cp, cko, vrs, er, beta;
    int*    r;

    if ( a == 1 ) {
        ff = "MLS_Weibull";
    }
    if ( a == 0 ) {
        ff = "MLS_Normal";
    }

    ifstream inp( "input/" + ff + ".inp" );
    inp >> s1;
    inp >> n;
    inp >> s1;
    xz = new double[n];
    for ( i = 0; i < n; i++ )
        inp >> xz[i];
    inp >> s1;
    r = new int[n];
    for ( i = 0; i < n; i++ )
        inp >> r[i];
    inp >> s1;
    inp >> kp;
    p  = new double[kp];
    zp = new double[kp];
    xp = new double[kp];
    inp >> s1;
    for ( i = 0; i < kp; i++ )
        inp >> p[i];
    inp.close();

    for ( i = 0; i < n; i++ ) {
        if ( a == 1 ) xz[i] = log( xz[i] );
    }

    km = 0;
    for ( i = 0; i < n; i++ )
        if ( r[i] == 0 ) km += 1;
    fcum = new double[km];
    ycum = new double[km];
    cum( n, xz, r, km, fcum, ycum );

    k  = 2;
    yr = new double[km];
    x  = new double*[km];
    y  = new double*[km];
    v  = new double*[km];
    db = new double*[k];
    b  = new double*[k];


    for ( i = 0; i < km; i++ ) {
        x[i] = new double[km];
        y[i] = new double[km];
        v[i] = new double[km];
    }
    for ( i = 0; i < k; i++ ) {
        db[i] = new double[k];
        b[i]  = new double[k];
    }
    for ( i = 0; i < k; i++ ) {
        for ( j = 0; j < k; j++ ) {
            b[i][j]  = 0;
            db[i][j] = 0;
        }
    }

    j = 0;
    for ( i = 0; i < n; i++ ) {
        if ( r[i] == 0 ) {
            y[j][0] = xz[i];
            j++;
        }
    }

    standart( km, ycum, cp, cko );


    for ( i = 0; i < km; i++ ) {
        for ( j = i; j < km; j++ ) {
            if ( a == 0 ) ordern( n, fcum[i], fcum[j], er, vrs );
            if ( a == 1 ) orderw( n, fcum[i], fcum[j], er, vrs );
            v[j][i] = vrs;
            v[i][j] = vrs;
        }
        x[i][0] = 1;
        x[i][1] = er;
    }


    MleastSquare_weight( km, k, x, y, v, db, b, yr );

    nc      = 4;
    int m[] = { km, kp, kp, kp };

    beta = 0.95;

    xplow = new double[kp];
    xpup  = new double[kp];

    for ( i = 0; i < kp; i++ ) {
        if ( a == 0 ) {
            zp[i] = norm_ppf( p[i] );
            delta = zp[i] * sqrt( n );
            lmtaprn( beta, n - 1, db[0][0], db[1][1], db[0][1], delta, tlow, tup );
        }
        if ( a == 1 ) {
            zp[i] = log( log( 1. / ( 1. - p[i] ) ) );
            delta = zp[i] * sqrt( n );
            nctWeibull( beta, n - 1, db[0][0], db[1][1], db[0][1], delta, tlow, tup );
        }
        xp[i]    = b[0][0] + b[1][0] * zp[i];
        xplow[i] = b[0][0] + b[1][0] * tlow / sqrt( n );
        xpup[i]  = b[0][0] + b[1][0] * tup / sqrt( n );
    }



    //########################Out##############################################
    //std::cout << "AAAAAAAAAAAAAA|" << ff << endl; 
    ofstream out( "output/" + ff + ".out" );
    out << "Количество: " << n << "\n";
    out << "Наблюдения" << "\n";
    for ( i = 0; i < n; i++ )
        out << xz[i] << "; ";
    out << "\n";
    out << "R" << "\n";
    for ( i = 0; i < n; i++ )
        out << r[i] << "; ";
    out << "\n";

    out << "cp* =" << cp << "\n";
    out << "cko* =" << cko << "\n";
    out << "A(B) =" << b[0][0] << "\n";
    out << "S(C) =" << b[1][0] << "\n";

    out << "P" << "\n";
    for ( i = 0; i < kp; i++ )
        out << p[i] << "; ";
    out << "\n";

    out << "t11=" << db[0][0] << "\n";
    out << "t12=" << db[0][1] << "\n";
    out << "t22=" << db[1][1] << "\n";

    out << "X cens" << "\n";
    for ( int i = 0; i < km; i++ ) {
        out << ycum[i] << "; ";
    }
    out << "\n";
    out << "ZX" << "\n";

    if ( a == 0 ) {
        for ( i = 0; i < km; i++ ) {
            out << ( 5.0 + norm_ppf( fcum[i] ) ) << "; ";
        }
    }

    if ( a == 1 ) {
        for ( i = 0; i < km; i++ ) {
            out << ( 5.0 + log( log( 1. / ( 1. - fcum[i] ) ) ) ) << "; ";
        }
    }
    out << endl;

    out << "XpLOW" << "\n";
    for ( int i = 0; i < kp; i++ )
        out << xplow[i] << "; ";
    out << "\n";
    out << "XP" << "\n";
    for ( int i = 0; i < kp; i++ )
        out << xp[i] << "; ";
    out << "\n";
    out << "XpUP" << "\n";
    for ( int i = 0; i < kp; i++ )
        out << xpup[i] << "; ";
    out << "\n";
    out << "ZP+5.0" << endl;
    for ( int i = 0; i < kp; i++ ) {
        out << ( zp[i] + 5.0 ) << "; ";
    }
    out << endl;


    out.close();

    delete[] xz, x, y, yr, v, m, zp, p, xplow, xpup, xp;
}

// int main(){
//     mls_weibull();
// }