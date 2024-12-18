#ifndef DEFAULT_H
#define DEFAULT_H



#include <string.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip> 
#include <vector> 
#include <map> 
#include <numeric>
#include <functional>
#include <ctime>

#include <boost/math/distributions/weibull.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/non_central_t.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/distributions/non_central_f.hpp>
#include <boost/math/distributions/binomial.hpp>
#include<boost/math/distributions.hpp>

using namespace std;
using namespace boost::math;

double norm_cdf(double x);
double norm_ppf(double p);
double norm_pdf(double x);
double t_cdf(double x,double f);
double t_ppf(double p,double f);
double t_pdf(double x,double f);
double chi_cdf(double x,double f);
double chi_ppf(double p,double f);
double chi_pdf(double x,double f);
double f_cdf(double x,double f1,double f2);
double f_ppf(double p,double f1,double f2);
double f_pdf(double x,double f1,double f2);
double nct_cdf(double x,double f,double delta);
double nct_ppf(double p,double f,double delta);
double nct_pdf(double x,double f,double delta);
double nchi_cdf(double x,double f,double delta);
double nchi_ppf(double p,double f,double delta);
double nchi_pdf(double x,double f,double delta);
double ncf_cdf(double x,double f1,double f2,double delta);
double ncf_ppf(double p,double f1,double f2,double delta);
double ncf_pdf(double x,double f1,double f2,double delta);

void CovMatrixMleW(int n,vector<double> x,vector<int>r,double c, double b, double **&v);
void CovMatrixMleN(int n,vector<double> x,vector<int>r,double a, double s, double **&v);
double NormalMinFunction(vector<double>xsimpl);
double WeibullMinFunction(vector<double>xsimpl);
void MleastSquare_weight(int n,int k,double **x,double **y,double **v,double **&db,double **&b,double *&yr);
void cum(int n,double* x,int* r,int &k,double* fcum,double* ycum);
void nctWeibull(double beta,int f,double t1,double t2,double t12,double d,double &tlow,double &tup);

int neldermead(vector<double>&x0, double eps,double(*func)(vector<double>));

double** TransMatrix(int m, int n, double** a);
double** MultiplyMatrix(int rowsa, int colsa, int rowsb, int colsb, double** a, double** b);
double **InverseMatrix(double **a,int n);

///////////////////////////////////////////////

struct ne_simp {
 int n;
 string ts;
 vector <double>p;
 vector <double>x;
 vector<int>r;
 vector<int>nsample;
};
ne_simp nesm;

double weibull_ppf(double p, double b, double c)
{
    if(p<=0 || p>=1) return 0;
    weibull_distribution<>d(b, c);  
    return(quantile(d, p));
}
//#############Normal Distribution############################

double norm_cdf(double x) {
   normal_distribution<>d(0,1);  
   return(cdf(d,x)); 
}
double norm_ppf(double p) {
   if(p<=0 || p>=1) return 0;
   normal_distribution<>d(0,1);  
   return(quantile(d,p)); 
}
double norm_pdf(double x) {
   normal_distribution<>d(0,1);  
   return(pdf(d,x)); 
}

//#############Student Distribution############################

double t_cdf(double x,double f) {
   students_t_distribution<>d(f);  
   return(cdf(d,x)); 
}
double t_ppf(double p,double f) {
   if(p<=0 || p>=1) return 0;
   students_t_distribution<>d(f); 
   return(quantile(d,p)); 
}
double t_pdf(double x,double f) {
   students_t_distribution<>d(f);  
   return(pdf(d,x)); 
}
//#############Chi-Squared Distribution############################
double chi_cdf(double x,double f) {
   chi_squared_distribution<>d(f);  
   return(cdf(d,x)); 
}
double chi_ppf(double p,double f) {
   if(p<=0 || p>=1) return 0;
   chi_squared_distribution<>d(f); 
   return(quantile(d,p)); 
}
double chi_pdf(double x,double f) {
   chi_squared_distribution<>d(f);  
   return(pdf(d,x)); 
}
//#############F-Distribution############################
double f_cdf(double x,double f1,double f2) {
   fisher_f_distribution<>d(f1,f2);  
   return(cdf(d,x)); 
}
double f_ppf(double p,double f1,double f2) {
   if(p<=0 || p>=1) return 0;
   fisher_f_distribution<>d(f1,f2); 
   return(quantile(d,p)); 
}
double f_pdf(double x,double f1,double f2) {
   fisher_f_distribution<>d(f1,f2);  
   return(pdf(d,x)); 
}
//#############Non Central t-Distribution############################

double nct_cdf(double x,double f,double delta) {
  non_central_t_distribution<>d(f,delta);
  return(cdf(d,x));
}
double nct_ppf(double p,double f,double delta) {
  if(p<=0 || p>=1) return 0;
  non_central_t_distribution<>d(f,delta);
  return(quantile(d,p));
}
double nct_pdf(double x,double f,double delta) {
  non_central_t_distribution<>d(f,delta);
  return(pdf(d,x));
}
//#############Non Central chi-squared-Distribution############################

double nchi_cdf(double x,double f,double delta) {
  non_central_chi_squared_distribution<>d(f,delta);
  return(cdf(d,x));
}
double nchi_ppf(double p,double f,double delta) {
  if(p<=0 || p>=1) return 0;
  non_central_chi_squared_distribution<>d(f,delta);
  return(quantile(d,p));
}
double nchi_pdf(double x,double f,double delta) {
  non_central_chi_squared_distribution<>d(f,delta);
  return(pdf(d,x));
}
//#############Non Central F-Distribution############################

double ncf_cdf(double x,double f1,double f2,double delta) {
  non_central_f_distribution<>d(f1,f2,delta);
  return(cdf(d,x));
}
double ncf_ppf(double p,double f1,double f2,double delta) {
  if(p<=0 || p>=1) return 0;
  non_central_f_distribution<>d(f1,f2,delta);
  return(quantile(d,p));
}
double ncf_pdf(double x,double f1,double f2,double delta) {
  non_central_f_distribution<>d(f1,f2,delta);
  return(pdf(d,x));
}

//########################################################################
void standart(int n,double *x,double &mean,double &s){
	 	int i;
	 	mean=0;s=0;
		for(i=0;i<n;i++){
			mean+=x[i];
			s+=x[i]*x[i];
		}
		mean/=n;
		s=(s-pow(mean,2)*n)/(n-1);
		s=pow(s,0.5);
	}


//###################Kaplan-Meier####################################
void cum(int n,double* x,int* r,int &k,double* fcum,double* ycum) {

 int i,j;
 double s;
 
  k=0;
  for (i=1;i<=n;i++) {
      s = 1.0;
      for (j=1;j<=i;j++)  if (r[j-1]==0) s=s*(n*1.0-j*1.0)/(n*1.0- j*1.0+1.0);
         if (r[i-1]==0) {
          //fcum[k]= 1.-s;
          //fcum[k]=((1.-s)*n-0.375)/(n+0.25);
          //fcum[k]=((1.-s)*n-0.5)/n;
          fcum[k]=(1.-s)*n*1.0/(n+1.0);
          if (s==1 || s==0) fcum[k] = ((1-s)*n-0.375)/(n+0.25);
          ycum[k]=x[i-1];k++;
      }
  }
}

/////////////////////////////////////////////////////////////////

void ordern(int n, double pr, double ps, double& er, double& vrs) {
    double p, pr1, ps1, xr, xr1, xr2, xr3, xr4, xr5, xr6, dr, qr, qs;
    double xs1, xs2, xs3, xs4, xs5, xs6, ds, xs;
    double z1, z2, z3, z4, z5, z6, z7;

    p = 1;
    //xr =normal_inv(pr,0,1);
    xr=norm_ppf(pr);
    //xs =normal_inv(ps,0,1);
    xs=norm_ppf(ps);
    qr = 1. - pr;
    qs = 1. - ps;
    pr1 = pr * p; ps1 = ps * p;
    //dr = normal_pdf(xr,0,1);
    //ds = normal_pdf(xs,0,1);

    dr = norm_pdf(xr);
    ds = norm_pdf(xs);

    xr1 = p / dr; xr2 = xr * xr1 * xr1;
    xr3 = (2. * xr * xr + 1.) * pow((p / dr), 3);
    xr4 = (6. * xr * xr * xr + 7. * xr) * pow((p / dr), 4);
    xr5 = (24. * pow(xr, 4) + 46. * xr * xr + 7.) * pow((p / dr), 5);
    xr6 = (120. * pow(xr, 5) + 326. * xr * xr * xr + 127. * xr) * pow((p / dr), 6);
    xs1 = p / ds; xs2 = xs * xs1 * xs1;
    xs3 = (2. * xs * xs + 1.) * pow((p / ds), 3);
    xs4 = (6. * xs * xs * xs + 7. * xs) * pow((p / ds), 4);
    xs5 = (24. * pow(xs, 4) + 46. * xs * xs + 7.) * pow((p / ds), 5);
    xs6 = (120. * pow(xs, 5) + 326. * xs * xs * xs + 127. * xs) * pow(p / ds, 6);

    er = xr + pr * qr * xr2 / (2. * (n + 2.)) + pr * qr * ((qr - pr) * xr3 / 3. + pr * qr * xr4 / 8.) / pow((n + 2.), 2) + pr * qr * (-(qr - pr) * xr3 / 3. + (pow((qr - pr), 2) - pr * qr) * xr4 / 4. + qr * pr * (qr - pr) * xr5 / 6. + pow((qr * pr), 2) * xr6 / 48.) / pow((n + 2.), 3);

    z1 = (qr - pr) * xr2 * xs1 + (qs - ps) * xr1 * xs2 + pr * qr * xr3 * xs1 / 2. + ps * qs * xr1 * xs3 / 2. + pr * qs * xr2 * xs2 / 2.;
    z1 = z1 * pr * qs / pow((n + 2.), 2);
    z2 = -(qr - pr) * xr2 * xs1 - (qs - ps) * xr1 * xs2 + (pow((qr - pr), 2) - pr * qr) * xr3 * xs1;
    z3 = (pow((qs - ps), 2) - ps * qs) * xr1 * xs3 + (1.5 * (qr - pr) * (qs - ps) + 0.5 * ps * qr - 2. * pr * qs) * xr2 * xs2;
    z4 = (5. / 6.) * pr * qr * (qr - pr) * xr4 * xs1 + (5. / 6.) * ps * qs * (qs - ps) * xr1 * xs4 + (pr * qs * (qr - pr) + .5 * pr * qr * (qs - ps)) * xr3 * xs2;
    z5 = (pr * qs * (qs - ps) + 0.5 * ps * qs * (qr - pr)) * xr2 * xs3 + (1. / 8.) * pow((pr * qr), 2) * xr5 * xs1 + (1. / 8.) * pow((ps * qs), 2) * xr1 * xs5;
    z6 = 0.25 * pr * pr * qr * qs * xr4 * xs2 + 0.25 * pr * ps * qs * qs * xr2 * xs4 + (2. * (pr * pr * qs * qs) + 3. * pr * qr * ps * qs) * xr3 * xs3 / 12.;
    z7 = z2 + z3 + z4 + z5 + z6;
    vrs = z1 + pr * qs * z7 / pow((n + 2.), 3) + pr * qs * xr1 * xs1 / (n + 2.);
}

////////////////////////////////////////////////////////////////

void orderw(int n, double pr, double ps, double &er, double &vrs) {
    double  qr, qs,xr, xr1, xr2, xr3, xr4, xr5, xr55, xr6;
    double xs1, xs2, xs3, xs4, xs5;
    double z1, z2, z3, z4, z5, z6, z7, a1, b1, c1, d1;

    xr = log(log(1. / (1. - pr)));
    qr = 1. - pr;
    qs = 1. - ps;
    xr1 = 1. / (log(1. / (1. - pr)) * (1. - pr));
    xr2 = xr1 * (1. / (1. - pr) - xr1);
    xr3 = xr2 * xr2 / xr1 + xr1 * (1. / pow((1. - pr), 2) - xr2);
    xr4 = (3. * xr1 * xr2 * xr3 - 2. * pow(xr2, 3)) / pow(xr1, 2) + xr1 * (2. / pow((1. - pr), 3) - xr3);
    xr55 = (-12. * xr1 * xr2 * xr2 * xr3 + 3. * xr1 * xr1 * xr3 * xr3 + 4. * xr1 * xr1 * xr2 * xr4 + 6. * pow(xr2, 4));
    xr5 = xr55 / pow(xr1, 3) + xr1 * (6. / pow(1. - pr, 4) - xr4);
    a1 = -12. * pow(xr2, 3) * xr3 - 12. * xr1 * (2. * xr2 * xr3 * xr3 + xr2 * xr2 * xr4);
    b1 = 6. * xr1 * xr2 * xr3 * xr3 + 6. * xr1 * xr1 * xr3 * xr4;
    c1 = 8. * xr1 * xr2 * xr2 * xr4 + 4. * xr1 * xr1 * (xr3 * xr4 + xr2 * xr5);
    d1 = 24. * pow(xr2, 3) * xr3;
    xr6 = (pow(xr1, 3) * (a1 + b1 + c1 + d1) - 3. * xr1 * xr1 * xr2 * xr55) / pow(xr1, 6) + xr2 * (6. / pow((1. - pr), 4) - xr4) + xr1 * (24. / pow((1. - pr), 5) - xr5);
    xs1 = 1. / (log(1. / (1. - ps)) * (1. - ps));
    xs2 = xs1 * (1. / (1. - ps) - xs1);
    xs3 = xs2 * xs2 / xs1 + xs1 * (1. / pow((1. - ps), 2) - xs2);
    xs4 = (3. * xs1 * xs2 * xs3 - 2. * pow(xs2, 3)) / (xs1 * xs1) + xs1 * (2. / pow((1. - ps), 3) - xs3);
    xs5 = (-12. * xs1 * xs2 * xs2 * xs3 + 3. * xs1 * xs1 * xs3 * xs3 + 4. * xs1 * xs1 * xs2 * xs4 + 6. * pow(xs2, 4)) / pow(xs1, 3) + xs1 * (6. / pow((1. - ps), 4) - xs4);

    z1 = pr * qr * xr2 / (2. * (n + 2.));
    z2 = pr * qr * ((qr - pr) * xr3 / 3. + pr * qr * xr4 / 8.) / pow((n + 2.), 2);
    z3 = -(qr - pr) * xr3 / 3. + (pow((qr - pr), 2) - pr * qr) * xr4 / 4.;
    z4 = qr * pr * (qr - pr) * xr5 / 6. + qr * qr * pr * pr * xr6 / 48.;
    er = xr + z1 + z2 + pr * qr * (z3 + z4) / pow((n + 2.), 3);

    z1 = (qr - pr) * xr2 * xs1 + (qs - ps) * xr1 * xs2 + pr * qr * xr3 * xs1 / 2. + ps * qs * xr1 * xs3 / 2. + pr * qs * xr2 * xs2 / 2.;
    z1 = z1 * pr * qs / pow((n + 2.), 2);
    z2 = -(qr - pr) * xr2 * xs1 - (qs - ps) * xr1 * xs2 + (pow((qr - pr), 2) - pr * qr) * xr3 * xs1;
    z3 = (pow((qs - ps), 2) - ps * qs) * xr1 * xs3 + (1.5 * (qr - pr) * (qs - ps) + 0.5 * ps * qr - 2. * pr * qs) * xr2 * xs2;
    z4 = (5. / 6.) * pr * qr * (qr - pr) * xr4 * xs1 + (5. / 6.) * ps * qs * (qs - ps) * xr1 * xs4 + (pr * qs * (qr - pr) + .5 * pr * qr * (qs - ps)) * xr3 * xs2;
    z5 = (pr * qs * (qs - ps) + .5 * ps * qs * (qr - pr)) * xr2 * xs3 + (1. / 8.) * pr * pr * qr * qr * xr5 * xs1 + (1. / 8.) * ps * ps * qs * qs * xr1 * xs5;
    z6 = 0.25 * pr * pr * qr * qs * xr4 * xs2 + 0.25 * pr * ps * qs * qs * xr2 * xs4 + (2. * pr * pr * qs * qs + 3. * pr * qr * ps * qs) * xr3 * xs3 / 12.;
    z7 = z2 + z3 + z4 + z5 + z6;
    vrs = z1 + pr * qs * z7 / pow((n + 2.), 3) + pr * qs * xr1 * xs1 / (n + 2.);
}


///////////////////////////////////////////////////////

//Приближенные доверительные граница для квантиля
void lmtaprn(double beta,int f,double t1,double t2,double t12,double d,double &tlow,double &tup){

double zb,f1x,f2x,f4x,e11,e2x,e3x,e44;

zb=norm_ppf(beta);
f1x=t2/f;
f2x=2*t12/sqrt(f+1.);
f4x=1.0-f1x/2;
e3x=f4x*f4x-zb*zb*f1x;
e11=f4x*d+zb*zb*f2x/2;
e2x=d*d-zb*zb*t1;
e44=sqrt(abs(e11*e11-e2x*e3x));
tlow=(e11-e44)/e3x;
tup=(e11+e44)/e3x;
}
//************************************************************************
double invnontap(double beta,int f,double d) {
	double zb, z, f4x;
 zb=norm_ppf(beta);
 f4x=1.0-1.0/(4.0*f);
 z=(f4x*d+zb*sqrt(f4x*f4x-zb*zb/(2.*f)+d*d/(2.*f)))/(f4x*f4x-zb*zb/(2.*f));
 return z;
}

/////////////////////////////////////////////

//Доверительные граница для квантиля Вейбулла
void nctWeibull(double beta,int f,double t1,double t2,double t12,double d,double &tlow,double &tup){

/*
 input:
 beta-доверительная вероятность
 f-число степеней свободы
 t1,t22-диагональные элементы ковариационной матрица оценок параметров
 t12,t21-недиагональные элементы ковариационной матрица оценок параметров
 output:
 tlow,tup-нижняя и верхняя квантили (процентные точки) нецентрального распределения Стьюдента
*/

  double zb,f1x,f2x,f4x,e11,e2x,e3x,e44;

   zb=norm_ppf(beta);
   f1x=t2/f;
   f2x=2*t12/sqrt(f+1);
   f4x=1.0-f1x/2;
   e3x=f4x*f4x-zb*zb*f1x;
   e11=f4x*d+zb*zb*f2x/2;
   e2x=d*d-zb*zb*t1;
   e44=sqrt(abs(e11*e11-e2x*e3x));
   tlow=(e11-e44)/e3x;
   tup=(e11+e44)/e3x;

}

//#############################################################

int neldermead(vector<double>&x0, double eps,double(*func)(vector<double>)) {

    double rho,chi, psi, sigma, nonzdelt, zdelt;
    double maxfun,fval,fxr,fxe,fxc,fxcc;
    int i,j,k,N,maxiter,iterations,doshrink;

    rho = 1;chi = 2;psi = 0.5;sigma = 0.5;nonzdelt = 0.05;zdelt = 0.00025;
    N = x0.size();maxiter = N * 200; maxfun = 1.e20;

    vector<vector<double>> sim(N + 1,vector<double>(N,0.0));
    vector<double> fsim(N + 1);
    sim[0]= x0;
    vector<double>y = x0;

    for (k = 0; k < N; k++) {
        if (y[k] != 0) {
            y[k] = (1 + nonzdelt) * y[k];
        } else {
            y[k] = zdelt;
        }
        sim[k + 1] = y;
    }
  
    for (i = 0; i < N + 1; i++)   fsim[i] = func(sim[i]);

    iterations = 0;
    fval = maxfun;

//###########################################################

    while (true) {
        if (fval <= eps || iterations >= maxiter) break;
        vector<double> xbar(N,0.0);
        for (i = 0; i < N; i++) {
            for ( j = 0; j < N; j++) xbar[j] += sim[i][j];
        }
        for (i = 0; i < N; i++) xbar[i] /= N;
        
        vector<double> xr(N);
        for (i = 0; i < N; i++) xr[i] = (1 + rho) * xbar[i] - rho * sim[N][i];
        fxr = func(xr);

        doshrink = 0;
        if (fxr < fsim[0]) {
            vector<double> xe(N);
            for (i = 0; i < N; i++)  xe[i] = (1 + rho * chi) * xbar[i] - rho * chi * sim[N][i];
            fxe = func(xe);
            if (fxe < fxr) {
                sim[N] = xe;fsim[N] = fxe;
            } else {
                sim[N] = xr;fsim[N] = fxr;
            }
        } else {
            if (fxr < fsim[N - 1]) {
                sim[N] = xr;fsim[N] = fxr;
            } else {
                if (fxr < fsim[N]) {
                    vector<double> xc(N);
                    for (i = 0; i < N; i++) xc[i] = (1 + psi * rho) * xbar[i] - psi * rho * sim[N][i];
                    fxc = func(xc);
                    if (fxc <= fxr) {
                        sim[N] = xc; fsim[N] = fxc;
                    } else {
                        doshrink = 1;
                    }
                } else {
                    vector<double> xcc(N);
                    for (i = 0; i < N; i++)  xcc[i] = (1 - psi) * xbar[i] + psi * sim[N][i];
                    fxcc = func(xcc);
                    if (fxcc < fsim[N]) {
                        sim[N] = xcc;fsim[N] = fxcc;
                    } else {
                        doshrink = 1;
                    }
                }
                if (doshrink) {
                    for (j = 1; j < N + 1; j++) {
                        for (i = 0; i < N; i++)     sim[j][i] = sim[0][i] + sigma * (sim[j][i] - sim[0][i]);
                        fsim[j] = func(sim[j]);
                    }
                }
            }
        }

        iterations++;
        fval = func(sim[0]);
        vector<int> ind(N + 1);
        for (i = 0; i < N + 1; i++) ind[i] = i;
        sort(ind.begin(), ind.end(), [&](int i, int j) { return fsim[i] < fsim[j]; });
        vector<vector<double>> sim_sorted(N + 1, vector<double>(N));
        vector<double> fsim_sorted(N + 1);
        for (i = 0; i < N + 1; i++) {
            sim_sorted[i] = sim[ind[i]];
            fsim_sorted[i] = fsim[ind[i]];
        }
        sim = sim_sorted; fsim = fsim_sorted;

    }  //end while

//########################################################################
    x0=sim[0];
    return iterations;

}

void CovMatrixMleW(int n,vector<double> x,vector<int>r,double c, double b, double **&v) {

    int i, k;
    double s1, s2, z, cpw, ckow;
    cpw = log(c); ckow = 1 / b; s1 = 0; s2 = 0; k = 0;
    for (i = 0; i < n; i++) {
        z = (log(x[i]) - cpw) / ckow;
        s1 += (1 -r[i]) * z;
        s2 += z * z * exp(z);
        k += (1 - r[i]);
    }
    v[0][0] =double(k)/double(n); v[0][1] = (k + s1) / n; v[1][0] = (k + s1) / n; v[1][1] = (k + s2) / n;
    v=InverseMatrix(v,2);
 
}


//############################################

void CovMatrixMleN(int n,vector<double> x,vector<int>r,double a, double s, double **&v) {

    double z, p, d, s1, s2, s3, psi;
    int j, k;
    s1 = 0; s2 = 0; s3 = 0; k = 0;

    for (j = 0; j < n; j++) {
        z = (x[j] - a) / s;
        //p = normal_cdf(z,0,1);
        //d =normal_pdf(-0.5*z*z,0,1); 

        p=norm_cdf(z);
        d=norm_pdf(-0.5*z*z); 

        psi = d / (1 - p);
        s1 += r[j] * psi * (psi - z);
        s2 += r[j] * psi * z * (z * (psi - z) - 1);
        s3 += r[j] * psi * (z * (psi - z) - 1);
        k += (1 - r[j]);
    }

    v[0][0] = (k + s1) / n; v[0][1] = s3 / n;
    v[1][0] = s3/n; v[1][1] = (2 * k + s2) / n;
    v=InverseMatrix(v,2);
   
}

//**************************Операции с матрицами****************************
double** TransMatrix(int m, int n, double** a) {
	int i, j;
	double** b;
	b = new double* [n];
	for (i = 0; i < n; i++) b[i] = new double[m];
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) b[i][j] = a[j][i];
	}
	return b;
}

//Умножение матриц
double** MultiplyMatrix(int rowsa, int colsa, int rowsb, int colsb, double** a, double** b) {
	int i, j, k;
	double t;
	double** c;
	c = new double* [rowsa];
	for (i = 0; i < rowsa; i++) c[i] = new double[colsb];

	if (colsa != rowsb) return 0;
	for (k = 0; k < colsb; k++) {
		for (i = 0; i < rowsa; i++) {
			t = 0;
			for (j = 0; j < rowsb; j++) t += a[i][j] * b[j][k];
			c[i][k] = t;
		}
	}
	return c;
}

//***********************************************************************
double **InverseMatrix(double **a,int n) {
        double temp;
        int i,j,k;

        double **e;
	e = new double* [n];
	for (i = 0; i < n; i++) e[i] = new double[n];
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                e[i][j] = 0;
                 if (i == j)
                    e[i][j] = 1;
            }
 
        for (k = 0; k < n; k++) {
            temp = a[k][k];
            for (j = 0; j < n; j++) {
                a[k][j] /= temp;
                e[k][j] /= temp;
            }
 
            for (i = k + 1; i < n; i++) {
                temp = a[i][k];
                for (j = 0; j < n; j++) {
                    a[i][j] -= a[k][j] * temp;
                    e[i][j] -= e[k][j] * temp;
                }
            }
        }
       for (k = n - 1; k > 0; k--) {
            for (i = k - 1; i >= 0; i--)  {
                temp = a[i][k];
                 for (j = 0; j < n; j++)  {
                    a[i][j] -= a[k][j] * temp;
                    e[i][j] -= e[k][j] * temp;
                }
            }
        }
    return e; 
    }

double  NormalMinFunction(vector<double>xsimpl) {
    double s1,s2,s3,s4,z,psi,p,d,c1,c2;
    int i,kx;
    s1 = 0; s2 = 0; s3 = 0; s4 = 0; kx = 0;
    if (xsimpl[0]<=0) return 10000;
    if (xsimpl[1]<=0) return 10000;

    for (i=0;i<nesm.n;i++) {
            z = (nesm.x[i]-xsimpl[0])/xsimpl[1];
            d = norm_pdf(z);
            p = norm_cdf(z);
            psi = d / (1. - p);
            s1 +=(1.-nesm.r[i])*(nesm.x[i]-xsimpl[0]);
            s2 += (1.-nesm.r[i])*pow(nesm.x[i]-xsimpl[0],2);
            s3 += nesm.r[i]*psi;
            s4 += nesm.r[i]*psi*z;
            kx += 1-nesm.r[i];
    }
    c1=s1+xsimpl[1]*s3;
    c2=s2+pow(xsimpl[1],2)*(s4-kx);
    z=c1*c1+c2*c2;
    return z;
}

double WeibullMinFunction(vector<double>xsimpl)
{
    double s1,s2,s3,s4, s5, temp, temp2, temp3, c1, c2;
    int kx;
    double z;
    s1 = 0; s2 = 0; s3 = 0; s4 = 0; s5 = 0; kx = 0;

    if (xsimpl[0]<=0) return 10000;

    for (int i = 0; i < nesm.n; i++) 
    {
        s1 += (1 - nesm.r[i]) * log(nesm.x[i]);
        s2 += (1 - nesm.r[i]) * pow(nesm.x[i], xsimpl[0]);
        s3 += nesm.r[i] * pow(nesm.x[i], xsimpl[0]);
        s4 += (1 - nesm.r[i]) * pow(nesm.x[i], xsimpl[0]) * log(nesm.x[i]);
        s5 += nesm.r[i] * pow(nesm.x[i], xsimpl[0]) * log(nesm.x[i]);
        kx += 1 - nesm.r[i];
    }

    s1 = s1 + (double(kx)/xsimpl[0]);
    temp = s2 + s3;
    c1 = s1 * temp;
    temp2 = s4+s5;
    c2 = temp2 * kx;
    temp3 = (c1 - c2);
    z = pow(temp3, 2);
    return z;
}
#endif