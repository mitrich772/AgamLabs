#include <string.h>
#include <fstream>
#include <cmath>
#include <iomanip> 
#include <vector> 
#include <map> 
#include <numeric>
#include <functional>
#include <ctime>
#include<aclui.h>
#include<iostream>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/non_central_t.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/distributions/non_central_f.hpp>
#include <boost/math/distributions/binomial.hpp>

using namespace std;
using namespace boost::math;

double  NormalMinFunction(vector<double>xsimpl);
double norm_cdf(double x);
double norm_ppf(double p);
double norm_pdf(double x);
double nct_cdf(double x,double f,double delta);
double nct_ppf(double p,double f,double delta);
double nct_pdf(double x,double f,double delta);
int neldermead(vector<double>&x0, double eps,double(*func)(vector<double>)); 
void CovMatrixMleN(int n,vector<double> x,vector<int>r,double a, double s, double **&v);
double** TransMatrix(int m, int n, double** a);
double** MultiplyMatrix(int rowsa, int colsa, int rowsb, int colsb, double** a, double** b);
double **InverseMatrix(double **a,int n);

struct ne_simp {
 int n;
 vector <double>p;
 vector <double>x;
 vector<int>r;
};
ne_simp nesm;

//#############Normal Distribution############################
double norm_cdf(double x) { //cumulative distribution function
   normal_distribution<>d(0,1);  
   return(cdf(d,x)); 
}
double norm_ppf(double p) { //percent point function
   if(p<=0 || p>=1) return 0;
   normal_distribution<>d(0,1);  
   return(quantile(d,p)); 
}
double norm_pdf(double x) { //probability density function
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
//#############Non Central chi-squared-Distribution######################

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
//#######################Kaplan-Meier#################################
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

//#################################################################

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


//#####################################################################

void CovMatrixMleN(int n,vector<double> x,vector<int>r,double a, double s, double **&v) {

    double z, p, d, s1, s2, s3, psi;
    int j, k;
    s1 = 0; s2 = 0; s3 = 0; k = 0;

    for (j = 0; j < n; j++) {
        z = (x[j] - a) / s;
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
//##########################################################

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

//##############################################################

int main() {
    int i,j,k,kx,icount,kp,nx,lim;
    string s1;
    double **v,*fcum,*ycum,*xplow,*xpup,*zp,*p,tlow,tup,*xp;
    double cp,cko,q,eps,beta,delta,step,*t,tp,z;

//#################################################################
    
    ifstream inp("MLE_Normal.inp");
    ofstream out("Out/MLE_Normal.out");

    inp >> s1; //text
    inp >>nesm.n; //sample size //количество тестов эксперемнетов

    inp >> s1;
    inp >> beta; //beta amount

    inp >> s1;
    inp >> eps;
    
    inp >> s1; //text
    for (i = 0; i <nesm.n; i++) {
        inp>>z;nesm.x.push_back(z);
    }

    inp >> s1; //text
    for (i=0;i<nesm.n; i++) {
        inp>>j;nesm.r.push_back(j);
    }

    inp >> s1;
    inp>>kp; // кол-во вер-тей квантилей

    p=new double[kp];
    zp=new double[kp];
    xp=new double[kp];

    inp >> s1;
    for (i = 0; i < kp; i++){
        inp >> p[i]; //считываем каждый квантиль
    };

    inp.close();

//###########################################################

     nx=2;   //количество переменных
     vector<double>xsimpl; //вектор переменных подлежащих оценке
     v = new double* [nx]; //ковариационная матрица
    

    for (i = 0; i < nx; i++) v[i] = new double[nx];
    for (i = 0; i < nx; i++) {
        for (j = 0; j < nx; j++) v[i][j] = 0;
    }

    cp=0;cko=0;k=0;
    for(i=0;i<nesm.n;i++) {
        k+=(1-nesm.r[i]); // количество наблюдений
        cp+=(1-nesm.r[i])*nesm.x[i];
        cko+=(1-nesm.r[i])*nesm.x[i]*nesm.x[i];
    }

    cp/=k; //выборочное среднее по наблюдениям
    cko=sqrt((cko-cp*cp*k)/(k-1)); //выборочное ско по наблюдениям

    q=0;icount=0; 
    xsimpl.push_back(cp);xsimpl.push_back(cko); //первые приближения
    
    v[0][0]=1.;v[1][1]=0.5;v[0][1]=0.;v[1][0]=0.;

     q=0;
      
    if(k!=nesm.n) {
     // Simplex
     icount=neldermead(xsimpl,eps,NormalMinFunction);
     q=NormalMinFunction(xsimpl);
     CovMatrixMleN(nesm.n,nesm.x,nesm.r, xsimpl[0],xsimpl[1],v);
    }

    //a=xsimpl[0]; s=xsimpl[1]

  // Конец расчета оценок параметров a,s

    // Для графика
    kx=kp;
    int nc;
    nc=4;
    int m[]={k,kx,kx,kx};

   //Доверительные границы

    t=new double[kx];
    xplow=new double[kx];
    xpup=new double[kx];
     
    for (i = 0; i < kx; i++) {
      //квантиль нормированного нормального распределения
      zp[i]=norm_ppf(p[i]);  
      //квантиль нормального распределения
      xp[i]=xsimpl[0]+zp[i]*xsimpl[1]; // xp=a+zp*s
      delta = zp[i] * sqrt(nesm.n);  //параметр нецентральности
      t[i]= nct_ppf(beta,nesm.n-1,delta);
      // верхняя доверительная граница 
      xpup[i]=xsimpl[0]+t[i]*xsimpl[1]/sqrt(nesm.n); 
      // нижняя доверительная граница
      xplow[i]=xsimpl[0]-t[kx-i-1]*xsimpl[1]/sqrt(nesm.n);
   }      


//###################График#####################################

    vector<double>xx;vector<double>y; 

    fcum = new double[nesm.n];ycum = new double[nesm.n];
    double *xq;
    int *rq;
    xq=new double[nesm.n];rq=new int[nesm.n];
    for(i=0;i<nesm.n;i++) {xq[i]=nesm.x[i];rq[i]=nesm.r[i];}
    cum(nesm.n,xq,rq,k,fcum,ycum); // размер выборки, время наблюдений, цензурирование,кол выж,массив для оценок функции выж,массив для соотв знач

    for (i = 0; i < k; i++) {
       xx.push_back(ycum[i]);
       y.push_back( 5.0 + norm_ppf(fcum[i]));
    }

    for (i = 0; i < kx; i++) {
       xx.push_back(xplow[i]);
       y.push_back(zp[i]+5.0);
    }

    for (i = 0; i < kx; i++) {
       xx.push_back(xp[i]);
       y.push_back(zp[i]+5.0);
    }
    
    for (i = 0; i < kx; i++) {
       xx.push_back(xpup[i]);
       y.push_back(zp[i]+5.0);
    }
 
    ofstream out1("Out/MLE_Normal.csv");
    out1<<nc<<endl;
    for(i=0;i<nc;i++) out1<<m[i]<<";";
    out1<<endl;

    int km=0;
    for(i=0;i<nc;i++) {
      for(j=0;j<m[i];j++) out1<<xx[j+km]<<";";
      km+=m[i];
      out1<<endl;
    }
     
     km=0;
     for(i=0;i<nc;i++) {
      for(j=0;j<m[i];j++) out1<<y[j+km]<<";";
      km+=m[i];
      out1<<endl;
    } 
     out1.close();

 //#################Вывод############################# 

 
 out << "n=" <<nesm.n << "\n";
 out << "X" << "\n";
 for (i = 0; i <nesm.n; i++)  out <<nesm.x[i] << " , ";
 out << "\n";
 out << "R" << "\n";
 for (i = 0; i<nesm.n; i++)  out <<nesm.r[i] << " , ";
 out << "\n";
 out << "cp*=" <<setprecision(12) << fixed <<cp << "\n";
 out << "cko*="<<setprecision(12) << fixed <<cko << "\n";
 out << "Q="<<setprecision(12) << q << "\n";
 out << "icount="<<icount<<endl;
 out << "cp="<<setprecision(12) << fixed << xsimpl[0]<<"\n";
 out << "cko="<<setprecision(12) << fixed <<xsimpl[1]<<"\n";
 out << "P" << "\n";
 for (i = 0; i < kx; i++) out << p[i] << " ; ";
 out << "\n";
 out << "Xplow" << "\n";
 for (i = 0; i < kx; i++) out << xplow[i] << " ; ";
 out << "\n";
 out << "Xp" << "\n";
 for (i = 0; i < kx; i++) out << xp[i] << " ; ";
 out << "\n";
 out << "Xpup" << "\n";
 for (i = 0; i < kx; i++) out << xpup[i] << " ; ";
 out << "\n";
 out << "v11=" << v[0][0] << "\n";
 out << "v12=" << v[0][1] << "\n";
 out << "v21=" << v[1][0] << "\n";
 out << "v22=" << v[1][1] << "\n";
 out.close();

 delete[] fcum, ycum, v,xplow,xpup,zp,p,xp;
 nesm.r.clear();nesm.x.clear();xsimpl.clear();
 return(0);
}




//######################################################

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
            kx+=1-nesm.r[i];
    }
    c1=s1+xsimpl[1]*s3;
    c2=s2+pow(xsimpl[1],2)*(s4-kx);
    z=c1*c1+c2*c2;
    return z;
}

