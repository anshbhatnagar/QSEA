#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <complex>

using namespace std;

double m12 = 0.014;
double m22 = 0.014;
double lambda1 = 0.1775;
double lambda2 = 0.185;
double lambda12 = 0.1;
double V4[2][2] = {{lambda1, lambda12},{lambda12, lambda2}};
double alpha = -0.09;
const double pi = 3.14159265358979323846;

double phi1Range[2] = {-0.2,1.5};
double phi2Range[2] = {-0.2,0.2};
double kStart=2;

int phi1Steps = 100;
int phi2Steps = 20;
int kSteps = 1000;

double T = 0;

int N = 2;

int kStop=kSteps;

double dphi1 = abs(phi1Range[1]-phi1Range[0])/phi1Steps;
double dphi2 = abs(phi2Range[1]-phi2Range[0])/phi2Steps;
double dk = -kStart/kSteps;

double A[6] = {0, 2./9, 1./3, 3./4, 1, 5./6};
double B[6][5] = {{0, 0, 0, 0, 0},
                {2./9, 0, 0, 0, 0},
                {1./12, 1./4, 0, 0, 0},
                {69./128, -243./128, 135./64, 0, 0},
                {-17./12, 27./4, -27./5, 16./15, 0},
                {65./432, -5./16, 13./16, 4./27, 5./144}};
double C[6] = {1./9, 0, 9./20, 16./45, 1./12, 0};
double CH[6] = {47./450, 0, 12./25, 32./225, 1./30, 6./25};
double CT[6] = {1./150, 0, -3./100, 16./75, 1./20, -6./25};

double releps = pow(10, -6);
double abseps = pow(10, -7);

double sF = 0.9;

double minKstep = pow(10, -8);

double zeroCutoff = pow(10,-13);

bool zeroTemp = false;
bool err = false;

complex<double> operator* (const int& lhs, const complex<double>& rhs) {
    return lhs*real(rhs)+(lhs*imag(rhs))*1i;
}

complex<double> operator/ (const int& lhs, const complex<double>& rhs) {
    double x = real(rhs);
    double y = imag(rhs);
    double divx = x/(x*x+y*y);
    double divy = y/(x*x+y*y);
    return lhs*(divx-divy*1i);
}

int heaviside(double theta){
    if(theta<0){
        return 0;
    }else{
        return 1;
    }
}

double omegaN(int n){
    return 2*pi*n*T;
}

double matsubaraSum(double kt2){
    double sum=0;
    int n = 1;

    if(kt2>0){
        sum+=pow(kt2, 3./2);
    }

    while(kt2>pow(omegaN(n),2)){
        sum+=2*pow(kt2-pow(omegaN(n),2),3./2);
        n+=1;
    }

    return sum;
}

double initV(double phi1, double phi2){
    double freeV1 = m12/2*phi1*phi1+alpha*pow(phi1,3)/6+lambda1*pow(phi1,4)/24;
    double freeV2 = m22/2*phi2*phi2+lambda2*pow(phi2,4)/24;
    double int12 = lambda12/4*pow(phi1*phi2,2);
    return freeV1+freeV2+int12;
}

double V2(double phi1, double phi2, int phiIndex){
    if(phiIndex==0){
        return m12+alpha*phi1+lambda1*phi1*phi1/2+lambda12/2*phi2*phi2;
    }else{
        return m22+lambda2*phi2*phi2/2+lambda12/2*phi1*phi1;
    } 
}

bool checkNan(complex<double> Gk1, complex<double> Gk2){
    return isnan(real(Gk1)) or isnan(real(Gk2));
}

double uk(double k, double phi1, double phi2, double upp1, double upp2){
    double deriv = 0;
    double upp[2]={upp1,upp2};
    double kt2[2]={0,0};
    double vA[2]={0,0};
    double mB[2][2]={{0,0},{0,0}};
    complex<double> Gk[2]={{0,0},{0,0}};

    for(int a=0;a<N;a++){
        kt2[a]=k*k-V2(phi1, phi2, a);
        vA[a]=kt2[a]+upp[a];
    }

    for(int a=0;a<N;a++){
        for(int b=0;b<N;b++){
            if(zeroTemp){
                mB[a][b]=pow(kt2[b],2)/(64*pi*pi)*V4[a][b];
            }else{
                mB[a][b]=T/(12*pi*pi)*V4[a][b]*matsubaraSum(kt2[b]);
            }
        }
    }

    complex<double> A1 = vA[0];
    complex<double> A2 = vA[1];
    complex<double> B11 = mB[0][0];
    complex<double> B12 = mB[0][1];
    complex<double> B21 = mB[1][0];
    complex<double> B22 = mB[1][1];

    Gk[0]=-0.25*(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22)/(B11*(-(B12*B21) + B11*B22)) + 
   sqrt(-((A1*A2*B12 + pow(B12,2) - B12*B21 - pow(A1,2)*B22 + 2*B11*B22)/(B11*(B12*B21 - B11*B22))) + 
      pow(-(A2*B11*B12) - A1*B12*B21 + 2.*A1*B11*B22,2)/(4.*pow(B11,2)*pow(-(B12*B21) + B11*B22,2)) + 
      (-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2.*B11*B22)/(3.*(-(B11*B12*B21) + pow(B11,2)*B22)) + 
      (pow(2,0.3333333333333333)*(pow(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2.*B11*B22,2) - 
           3.*(A2*B12 - 2.*A1*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2.*A1*B11*B22) + 12.*B22*(-(B11*B12*B21) + pow(B11,2)*B22)))/
       (3.*B11*(-(B12*B21) + B11*B22)*pow(2*pow(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2.*B11*B22,3) - 
           9.*(A2*B12 - 2*A1*B22)*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2.*B11*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2.*A1*B11*B22) + 
           27.*B22*pow(-(A2*B11*B12) - A1*B12*B21 + 2.*A1*B11*B22,2) + 27.*pow(A2*B12 - 2.*A1*B22,2)*(-(B11*B12*B21) + pow(B11,2)*B22) - 
           72.*B22*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2.*B11*B22)*(-(B11*B12*B21) + pow(B11,2)*B22) + 
           sqrt(-4.*pow(pow(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2.*B11*B22,2) - 3.*(A2*B12 - 2*A1*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2.*A1*B11*B22) + 
                12.*B22*(-(B11*B12*B21) + pow(B11,2)*B22),3) + pow(2*pow(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2.*B11*B22,3) - 
               9.*(A2*B12 - 2*A1*B22)*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2.*B11*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2.*A1*B11*B22) + 
               27.*B22*pow(-(A2*B11*B12) - A1*B12*B21 + 2.*A1*B11*B22,2) + 27.*pow(A2*B12 - 2.*A1*B22,2)*(-(B11*B12*B21) + pow(B11,2)*B22) - 
               72.*B22*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2.*B11*B22)*(-(B11*B12*B21) + pow(B11,2)*B22),2)),0.3333333333333333)) + 
      pow(2.*pow(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22,3) - 
         9*(A2*B12 - 2*A1*B22)*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22) + 
         27*B22*pow(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22,2) + 27*pow(A2*B12 - 2*A1*B22,2)*(-(B11*B12*B21) + pow(B11,2)*B22) - 
         72*B22*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)*(-(B11*B12*B21) + pow(B11,2)*B22) + 
         sqrt(-4*pow(pow(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22,2) - 3*(A2*B12 - 2*A1*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22) + 
              12*B22*(-(B11*B12*B21) + pow(B11,2)*B22),3) + pow(2*pow(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22,3) - 
             9*(A2*B12 - 2*A1*B22)*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22) + 
             27*B22*pow(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22,2) + 27*pow(A2*B12 - 2*A1*B22,2)*(-(B11*B12*B21) + pow(B11,2)*B22) - 
             72*B22*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)*(-(B11*B12*B21) + pow(B11,2)*B22),2)),0.3333333333333333)/
       (3.*pow(2,0.3333333333333333)*B11*(-(B12*B21) + B11*B22)))/2. - 
   sqrt(-((A1*A2*B12 + pow(B12,2) - B12*B21 - pow(A1,2)*B22 + 2*B11*B22)/(B11*(B12*B21 - B11*B22))) + 
      pow(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22,2)/(2.*pow(B11,2)*pow(-(B12*B21) + B11*B22,2)) - 
      (-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)/(3.*(-(B11*B12*B21) + pow(B11,2)*B22)) - 
      (pow(2,0.3333333333333333)*(pow(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22,2) - 
           3*(A2*B12 - 2*A1*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22) + 12*B22*(-(B11*B12*B21) + pow(B11,2)*B22)))/
       (3.*B11*(-(B12*B21) + B11*B22)*pow(2*pow(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22,3) - 
           9*(A2*B12 - 2*A1*B22)*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22) + 
           27*B22*pow(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22,2) + 27*pow(A2*B12 - 2*A1*B22,2)*(-(B11*B12*B21) + pow(B11,2)*B22) - 
           72*B22*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)*(-(B11*B12*B21) + pow(B11,2)*B22) + 
           sqrt(-4*pow(pow(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22,2) - 3*(A2*B12 - 2*A1*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22) + 
                12*B22*(-(B11*B12*B21) + pow(B11,2)*B22),3) + pow(2*pow(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22,3) - 
               9*(A2*B12 - 2*A1*B22)*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22) + 
               27*B22*pow(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22,2) + 27*pow(A2*B12 - 2*A1*B22,2)*(-(B11*B12*B21) + pow(B11,2)*B22) - 
               72*B22*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)*(-(B11*B12*B21) + pow(B11,2)*B22),2)),0.3333333333333333)) - 
      pow(2*pow(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22,3) - 
         9*(A2*B12 - 2*A1*B22)*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22) + 
         27*B22*pow(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22,2) + 27*pow(A2*B12 - 2*A1*B22,2)*(-(B11*B12*B21) + pow(B11,2)*B22) - 
         72*B22*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)*(-(B11*B12*B21) + pow(B11,2)*B22) + 
         sqrt(-4*pow(pow(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22,2) - 3*(A2*B12 - 2*A1*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22) + 
              12*B22*(-(B11*B12*B21) + pow(B11,2)*B22),3) + pow(2*pow(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22,3) - 
             9*(A2*B12 - 2*A1*B22)*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22) + 
             27*B22*pow(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22,2) + 27*pow(A2*B12 - 2*A1*B22,2)*(-(B11*B12*B21) + pow(B11,2)*B22) - 
             72*B22*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)*(-(B11*B12*B21) + pow(B11,2)*B22),2)),0.3333333333333333)/
       (3.*pow(2,0.3333333333333333)*B11*(-(B12*B21) + B11*B22)) + ((8*(-(A2*B12) + 2*A1*B22))/(B11*(-(B12*B21) + B11*B22)) + 
         (4*(A1*A2*B12 + pow(B12,2) - B12*B21 - pow(A1,2)*B22 + 2*B11*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22))/
          (pow(B11,2)*(B12*B21 - B11*B22)*(-(B12*B21) + B11*B22)) - pow(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22,3)/(pow(B11,3)*pow(-(B12*B21) + B11*B22,3)))/
       (4.*sqrt(-((A1*A2*B12 + pow(B12,2) - B12*B21 - pow(A1,2)*B22 + 2*B11*B22)/(B11*(B12*B21 - B11*B22))) + 
           pow(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22,2)/(4.*pow(B11,2)*pow(-(B12*B21) + B11*B22,2)) + 
           (-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)/(3.*(-(B11*B12*B21) + pow(B11,2)*B22)) + 
           (pow(2,0.3333333333333333)*(pow(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22,2) - 
                3*(A2*B12 - 2*A1*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22) + 12*B22*(-(B11*B12*B21) + pow(B11,2)*B22)))/
            (3.*B11*(-(B12*B21) + B11*B22)*pow(2*pow(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22,3) - 
                9*(A2*B12 - 2*A1*B22)*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22) + 
                27*B22*pow(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22,2) + 27*pow(A2*B12 - 2*A1*B22,2)*(-(B11*B12*B21) + pow(B11,2)*B22) - 
                72*B22*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)*(-(B11*B12*B21) + pow(B11,2)*B22) + 
                sqrt(-4*pow(pow(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22,2) - 
                     3*(A2*B12 - 2*A1*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22) + 12*B22*(-(B11*B12*B21) + pow(B11,2)*B22),3) + 
                  pow(2*pow(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22,3) - 
                    9*(A2*B12 - 2*A1*B22)*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22) + 
                    27*B22*pow(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22,2) + 27*pow(A2*B12 - 2*A1*B22,2)*(-(B11*B12*B21) + pow(B11,2)*B22) - 
                    72*B22*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)*(-(B11*B12*B21) + pow(B11,2)*B22),2)),0.3333333333333333)) + 
           pow(2*pow(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22,3) - 
              9*(A2*B12 - 2*A1*B22)*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22) + 
              27*B22*pow(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22,2) + 27*pow(A2*B12 - 2*A1*B22,2)*(-(B11*B12*B21) + pow(B11,2)*B22) - 
              72*B22*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)*(-(B11*B12*B21) + pow(B11,2)*B22) + 
              sqrt(-4*pow(pow(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22,2) - 
                   3*(A2*B12 - 2*A1*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22) + 12*B22*(-(B11*B12*B21) + pow(B11,2)*B22),3) + 
                pow(2*pow(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22,3) - 
                  9*(A2*B12 - 2*A1*B22)*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)*(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22) + 
                  27*B22*pow(-(A2*B11*B12) - A1*B12*B21 + 2*A1*B11*B22,2) + 27*pow(A2*B12 - 2*A1*B22,2)*(-(B11*B12*B21) + pow(B11,2)*B22) - 
                  72*B22*(-(A1*A2*B12) - pow(B12,2) + B12*B21 + pow(A1,2)*B22 - 2*B11*B22)*(-(B11*B12*B21) + pow(B11,2)*B22),2)),0.3333333333333333)/
            (3.*pow(2,0.3333333333333333)*B11*(-(B12*B21) + B11*B22)))))/2.;

    Gk[1]=(-A2 - B21*Gk[0] + sqrt(4*B22 + pow(A2 + B21*Gk[0],2)))/(2.*B22);
    
    if(checkNan(Gk[0],Gk[1])){
        if(real(B11)>zeroCutoff and real(B22)>zeroCutoff){
            if(real(B21)<zeroCutoff and real(B12)>zeroCutoff){
                Gk[0]=-0.25*(B12*sqrt(pow(A2,2) + 4*B22))/(B11*B22) - (-(A2*B12) + 2*A1*B22)/(4.*B11*B22) + 
   sqrt((2*pow(A1,2))/pow(B11,2) + 8/B11 + (pow(A2,2)*pow(B12,2))/(pow(B11,2)*pow(B22,2)) - 
      (2*A1*A2*B12)/(pow(B11,2)*B22) + (2*pow(B12,2))/(pow(B11,2)*B22) + 
      (8*A1*B12)/(pow(B11,2)*sqrt(pow(A2,2) + 4*B22)) - 
      (pow(A2,3)*pow(B12,2))/(pow(B11,2)*pow(B22,2)*sqrt(pow(A2,2) + 4*B22)) + 
      (2*A1*pow(A2,2)*B12)/(pow(B11,2)*B22*sqrt(pow(A2,2) + 4*B22)) - 
      (4*A2*pow(B12,2))/(pow(B11,2)*B22*sqrt(pow(A2,2) + 4*B22)))/(2.*sqrt(2));
            }else if(real(B21)>zeroCutoff and real(B12)<zeroCutoff){
                Gk[0]=(-A1 + sqrt(pow(A1,2) + 4*B11))/(2.*B11);
            }else if(real(B21)<zeroCutoff and real(B12)<zeroCutoff){
                Gk[0]=(-A1 + sqrt(pow(A1,2) + 4*B11))/(2.*B11);
            }
            Gk[1]=(-A2 - B21*Gk[0] + sqrt(4*B22 + pow(A2 + B21*Gk[0],2)))/(2.*B22);
        }else if(real(B11)<zeroCutoff and real(B22)>zeroCutoff){
            if(real(B21)>zeroCutoff and real(B12)>zeroCutoff){
                Gk[0]=-((B12 - B21 - A2*B12*(-0.3333333333333333*(A2*B12 + A1*B22)/(B12*B22) - 
          (pow(2,0.3333333333333333)*(3*B12*(A1*A2 - B12 + B21)*B22 - pow(A2*B12 + A1*B22,2)))/
           (3.*B12*B22*pow(-2*pow(A2,3)*pow(B12,3) + 3*A1*pow(A2,2)*pow(B12,2)*B22 - 9*A2*pow(B12,3)*B22 + 
               9*A2*pow(B12,2)*B21*B22 + 3*pow(A1,2)*A2*B12*pow(B22,2) + 18*A1*pow(B12,2)*pow(B22,2) + 
               9*A1*B12*B21*pow(B22,2) - 2*pow(A1,3)*pow(B22,3) + 
               sqrt(pow(-2*pow(A2,3)*pow(B12,3) + 3*A1*pow(A2,2)*pow(B12,2)*B22 - 9*A2*pow(B12,3)*B22 + 
                   9*A2*pow(B12,2)*B21*B22 + 3*pow(A1,2)*A2*B12*pow(B22,2) + 18*A1*pow(B12,2)*pow(B22,2) + 
                   9*A1*B12*B21*pow(B22,2) - 2*pow(A1,3)*pow(B22,3),2) + 
                 4*pow(3*B12*(A1*A2 - B12 + B21)*B22 - pow(A2*B12 + A1*B22,2),3)),0.3333333333333333)) + 
          pow(-2*pow(A2,3)*pow(B12,3) + 3*A1*pow(A2,2)*pow(B12,2)*B22 - 9*A2*pow(B12,3)*B22 + 
             9*A2*pow(B12,2)*B21*B22 + 3*pow(A1,2)*A2*B12*pow(B22,2) + 18*A1*pow(B12,2)*pow(B22,2) + 
             9*A1*B12*B21*pow(B22,2) - 2*pow(A1,3)*pow(B22,3) + 
             sqrt(pow(-2*pow(A2,3)*pow(B12,3) + 3*A1*pow(A2,2)*pow(B12,2)*B22 - 9*A2*pow(B12,3)*B22 + 
                 9*A2*pow(B12,2)*B21*B22 + 3*pow(A1,2)*A2*B12*pow(B22,2) + 18*A1*pow(B12,2)*pow(B22,2) + 
                 9*A1*B12*B21*pow(B22,2) - 2*pow(A1,3)*pow(B22,3),2) + 
               4*pow(3*B12*(A1*A2 - B12 + B21)*B22 - pow(A2*B12 + A1*B22,2),3)),0.3333333333333333)/
           (3.*pow(2,0.3333333333333333)*B12*B22)) - 
       B12*B22*pow(-0.3333333333333333*(A2*B12 + A1*B22)/(B12*B22) - 
          (pow(2,0.3333333333333333)*(3*B12*(A1*A2 - B12 + B21)*B22 - pow(A2*B12 + A1*B22,2)))/
           (3.*B12*B22*pow(-2*pow(A2,3)*pow(B12,3) + 3*A1*pow(A2,2)*pow(B12,2)*B22 - 9*A2*pow(B12,3)*B22 + 
               9*A2*pow(B12,2)*B21*B22 + 3*pow(A1,2)*A2*B12*pow(B22,2) + 18*A1*pow(B12,2)*pow(B22,2) + 
               9*A1*B12*B21*pow(B22,2) - 2*pow(A1,3)*pow(B22,3) + 
               sqrt(pow(-2*pow(A2,3)*pow(B12,3) + 3*A1*pow(A2,2)*pow(B12,2)*B22 - 9*A2*pow(B12,3)*B22 + 
                   9*A2*pow(B12,2)*B21*B22 + 3*pow(A1,2)*A2*B12*pow(B22,2) + 18*A1*pow(B12,2)*pow(B22,2) + 
                   9*A1*B12*B21*pow(B22,2) - 2*pow(A1,3)*pow(B22,3),2) + 
                 4*pow(3*B12*(A1*A2 - B12 + B21)*B22 - pow(A2*B12 + A1*B22,2),3)),0.3333333333333333)) + 
          pow(-2*pow(A2,3)*pow(B12,3) + 3*A1*pow(A2,2)*pow(B12,2)*B22 - 9*A2*pow(B12,3)*B22 + 
             9*A2*pow(B12,2)*B21*B22 + 3*pow(A1,2)*A2*B12*pow(B22,2) + 18*A1*pow(B12,2)*pow(B22,2) + 
             9*A1*B12*B21*pow(B22,2) - 2*pow(A1,3)*pow(B22,3) + 
             sqrt(pow(-2*pow(A2,3)*pow(B12,3) + 3*A1*pow(A2,2)*pow(B12,2)*B22 - 9*A2*pow(B12,3)*B22 + 
                 9*A2*pow(B12,2)*B21*B22 + 3*pow(A1,2)*A2*B12*pow(B22,2) + 18*A1*pow(B12,2)*pow(B22,2) + 
                 9*A1*B12*B21*pow(B22,2) - 2*pow(A1,3)*pow(B22,3),2) + 
               4*pow(3*B12*(A1*A2 - B12 + B21)*B22 - pow(A2*B12 + A1*B22,2),3)),0.3333333333333333)/
           (3.*pow(2,0.3333333333333333)*B12*B22),2))/(A1*B21));
            }else if(real(B21)<zeroCutoff and real(B12)>zeroCutoff){
                Gk[0]=(-(A2*B12) + 2*A1*B22 - B12*sqrt(pow(A2,2) + 4*B22))/(2.*(-(A1*A2*B12) - pow(B12,2) + pow(A1,2)*B22));
            }else if(real(B21)>zeroCutoff and real(B12)<zeroCutoff){
                Gk[0]=1./A1;
            }else if(real(B21)<zeroCutoff and real(B12)<zeroCutoff){
                Gk[0]=1./A1;
            }
            Gk[1]=(-A2 - B21*Gk[0] + sqrt(4*B22 + pow(A2 + B21*Gk[0],2)))/(2.*B22);
        }else if(real(B11)>zeroCutoff and real(B22)<zeroCutoff){
            if(real(B21)>zeroCutoff and real(B12)>zeroCutoff){
                Gk[0]=-0.3333333333333333*(A2*B11 + A1*B21)/(B11*B21) - (pow(2,0.3333333333333333)*
      (3*B11*(A1*A2 + B12 - B21)*B21 - pow(A2*B11 + A1*B21,2)))/
    (3.*B11*B21*pow(-2*pow(A2,3)*pow(B11,3) + 3*A1*pow(A2,2)*pow(B11,2)*B21 + 9*A2*pow(B11,2)*B12*B21 + 
        3*pow(A1,2)*A2*B11*pow(B21,2) + 18*A2*pow(B11,2)*pow(B21,2) + 9*A1*B11*B12*pow(B21,2) - 
        2*pow(A1,3)*pow(B21,3) - 9*A1*B11*pow(B21,3) + 
        sqrt(pow(-2*pow(A2,3)*pow(B11,3) + 3*A1*pow(A2,2)*pow(B11,2)*B21 + 9*A2*pow(B11,2)*B12*B21 + 
            3*pow(A1,2)*A2*B11*pow(B21,2) + 18*A2*pow(B11,2)*pow(B21,2) + 9*A1*B11*B12*pow(B21,2) - 
            2*pow(A1,3)*pow(B21,3) - 9*A1*B11*pow(B21,3),2) + 
          4*pow(3*B11*(A1*A2 + B12 - B21)*B21 - pow(A2*B11 + A1*B21,2),3)),0.3333333333333333)) + 
   pow(-2*pow(A2,3)*pow(B11,3) + 3*A1*pow(A2,2)*pow(B11,2)*B21 + 9*A2*pow(B11,2)*B12*B21 + 
      3*pow(A1,2)*A2*B11*pow(B21,2) + 18*A2*pow(B11,2)*pow(B21,2) + 9*A1*B11*B12*pow(B21,2) - 
      2*pow(A1,3)*pow(B21,3) - 9*A1*B11*pow(B21,3) + 
      sqrt(pow(-2*pow(A2,3)*pow(B11,3) + 3*A1*pow(A2,2)*pow(B11,2)*B21 + 9*A2*pow(B11,2)*B12*B21 + 
          3*pow(A1,2)*A2*B11*pow(B21,2) + 18*A2*pow(B11,2)*pow(B21,2) + 9*A1*B11*B12*pow(B21,2) - 
          2*pow(A1,3)*pow(B21,3) - 9*A1*B11*pow(B21,3),2) + 
        4*pow(3*B11*(A1*A2 + B12 - B21)*B21 - pow(A2*B11 + A1*B21,2),3)),0.3333333333333333)/
    (3.*pow(2,0.3333333333333333)*B11*B21);
            }else if(real(B21)<zeroCutoff and real(B12)>zeroCutoff){
                Gk[0]=(-(A1*A2) - B12 + sqrt(4*pow(A2,2)*B11 + pow(A1*A2 + B12,2)))/(2.*A2*B11);
            }else if(real(B21)>zeroCutoff and real(B12)<zeroCutoff){
                Gk[0]=(-A1 + sqrt(pow(A1,2) + 4*B11))/(2.*B11);
            }else if(real(B21)<zeroCutoff and real(B12)<zeroCutoff){
                Gk[0]=(-A1 + sqrt(pow(A1,2) + 4*B11))/(2.*B11);
            }
            Gk[1]=1./(A2+B21*Gk[0]);
        }else if(real(B11)<zeroCutoff and real(B22)<zeroCutoff){
            if(real(B21)>zeroCutoff and real(B12)>zeroCutoff){
                Gk[0]=(-(A1*A2) - B12 + B21 + sqrt(pow(A1*A2 + B12 - B21,2) + 4*A1*A2*B21))/(2.*A1*B21);
                Gk[1]=(-(A1*A2) + B12 - B21 + sqrt(pow(A1*A2 + B12 - B21,2) + 4*A1*A2*B21))/(2.*A2*B12);
            }else if(real(B21)<zeroCutoff and real(B12)>zeroCutoff){
                Gk[0]=A2/(A1*A2+B12);
                Gk[1]=1./A2;
            }else if(real(B21)>zeroCutoff and real(B12)<zeroCutoff){
                Gk[0]=1./A1;
                Gk[1]=A1/(A1*A2+B21);
            }else if(real(B21)<zeroCutoff and real(B12)<zeroCutoff){
                Gk[0]=1./A1;
                Gk[1]=1./A2;
            }
        }
    }

    for(int a=0;a<N;a++){
        if(zeroTemp){
            deriv += k/(32*pi*pi)*pow(kt2[a],2)*heaviside(kt2[a])*real(Gk[a]);
        }else{
            deriv += k*T/(6*pi*pi)*matsubaraSum(kt2[a])*heaviside(kt2[a])*real(Gk[a]);
        }
    }

    return deriv;
}

double derivUpp(double **u, int i, int j, bool firstDir){
    double upp;

    if(firstDir){
        if(i==0){
            upp=(2*u[0][j]-5*u[1][j]+4*u[2][j]-u[3][j])/pow(dphi1,3);
        }else if(i==phi1Steps-1){
            upp=(2*u[i][j]-5*u[i-1][j]+4*u[i-2][j]-u[i-3][j])/pow(dphi1,3);
        }else{
            upp=(u[i+1][j]+u[i-1][j]-2*u[i][j])/pow(dphi1,2);
        }
    }else{
        if(j==0){
            upp=(2*u[i][0]-5*u[i][1]+4*u[i][2]-u[i][3])/pow(dphi2,3);
        }else if(j==phi2Steps-1){
            upp=(2*u[i][j]-5*u[i][j-1]+4*u[i][j-2]-u[i][j-3])/pow(dphi2,3);
        }else{
            upp=(u[i][j+1]+u[i][j-1]-2*u[i][j])/pow(dphi2,2);
        }
    }

    return upp;
}

void nextStep(double **u, double k, double **k1, double **k2, double **k3, double **k4){
    double upp1;
    double upp2;
    double phi1;
    double phi2;
    double** uEdit = new double*[phi1Steps];

    for(int i=0;i<phi1Steps;i++){
        uEdit[i] = new double[phi2Steps];
    }

    for(int i=0;i<phi1Steps;i++){
        phi1 = phi1Range[0]+i*dphi1;
        for(int j=0;j<phi2Steps;j++){
            phi2 = phi2Range[0]+j*dphi2;
            upp1 = derivUpp(u,i,j,true);
            upp2 = derivUpp(u,i,j,false);
            k1[i][j]=uk(k, phi1, phi2, upp1, upp2);
        }
    }

    for(int i=0;i<phi1Steps;i++){
        for(int j=0;j<phi2Steps;j++){
            uEdit[i][j]=u[i][j]+dk*k1[i][j]/2;
        }
    }

    for(int i=0;i<phi1Steps;i++){
        phi1 = phi1Range[0]+i*dphi1;
        for(int j=0;j<phi2Steps;j++){
            phi2 = phi2Range[0]+j*dphi2;
            upp1 = derivUpp(uEdit,i,j,true);
            upp2 = derivUpp(uEdit,i,j,false);
            k2[i][j]=uk(k+dk/2, phi1, phi2, upp1, upp2);
        }
    }

    for(int i=0;i<phi1Steps;i++){
        for(int j=0;j<phi2Steps;j++){
            uEdit[i][j]=u[i][j]+dk*k2[i][j]/2;
        }
    }

    for(int i=0;i<phi1Steps;i++){
        phi1 = phi1Range[0]+i*dphi1;
        for(int j=0;j<phi2Steps;j++){
            phi2 = phi2Range[0]+j*dphi2;
            upp1 = derivUpp(uEdit,i,j,true);
            upp2 = derivUpp(uEdit,i,j,false);
            k3[i][j]=uk(k+dk/2, phi1, phi2, upp1, upp2);
        }
    }

    for(int i=0;i<phi1Steps;i++){
        for(int j=0;j<phi2Steps;j++){
            uEdit[i][j]=u[i][j]+dk*k3[i][j];
        }
    }

    for(int i=0;i<phi1Steps;i++){
        phi1 = phi1Range[0]+i*dphi1;
        for(int j=0;j<phi2Steps;j++){
            phi2 = phi2Range[0]+j*dphi2;
            upp1 = derivUpp(uEdit,i,j,true);
            upp2 = derivUpp(uEdit,i,j,false);
            k4[i][j]=uk(k+dk, phi1, phi2, upp1, upp2);
        }
    }

    for(int i=0;i<phi1Steps;i++){
        for(int j=0;j<phi2Steps;j++){
            u[i][j]+= dk*(k1[i][j]+2*k2[i][j]+2*k3[i][j]+k4[i][j])/6;
        }
    }

    delete[] uEdit;

}

bool nextStepRK45(double **u, double k, double **k1, double **k2, double **k3, double **k4, double **k5, double **k6){
    double upp1;
    double upp2;
    double phi1;
    double phi2;
    double** uEdit = new double*[phi1Steps];

    for(int i=0;i<phi1Steps;i++){
        uEdit[i] = new double[phi2Steps];
    }

    for(int i=0;i<phi1Steps;i++){
        phi1 = phi1Range[0]+i*dphi1;
        for(int j=0;j<phi2Steps;j++){
            phi2 = phi2Range[0]+j*dphi2;
            upp1 = derivUpp(u,i,j,true);
            upp2 = derivUpp(u,i,j,false);
            k1[i][j]=uk(k+A[0]*dk, phi1, phi2, upp1, upp2);
        }
    }

    for(int i=0;i<phi1Steps;i++){
        for(int j=0;j<phi2Steps;j++){
            uEdit[i][j]=u[i][j]+B[1][0]*dk*k1[i][j];
        }
    }

    for(int i=0;i<phi1Steps;i++){
        phi1 = phi1Range[0]+i*dphi1;
        for(int j=0;j<phi2Steps;j++){
            phi2 = phi2Range[0]+j*dphi2;
            upp1 = derivUpp(uEdit,i,j,true);
            upp2 = derivUpp(uEdit,i,j,false);
            k2[i][j]=uk(k+A[1]*dk, phi1, phi2, upp1, upp2);
        }
    }

    for(int i=0;i<phi1Steps;i++){
        for(int j=0;j<phi2Steps;j++){
            uEdit[i][j]=u[i][j]+B[2][0]*dk*k1[i][j]+B[2][1]*dk*k2[i][j];
        }
    }

    for(int i=0;i<phi1Steps;i++){
        phi1 = phi1Range[0]+i*dphi1;
        for(int j=0;j<phi2Steps;j++){
            phi2 = phi2Range[0]+j*dphi2;
            upp1 = derivUpp(uEdit,i,j,true);
            upp2 = derivUpp(uEdit,i,j,false);
            k3[i][j]=uk(k+A[2]*dk, phi1, phi2, upp1, upp2);
        }
    }

    for(int i=0;i<phi1Steps;i++){
        for(int j=0;j<phi2Steps;j++){
            uEdit[i][j]=u[i][j]+B[3][0]*dk*k1[i][j]+B[3][1]*dk*k2[i][j]+B[3][2]*dk*k3[i][j];
        }
    }

    for(int i=0;i<phi1Steps;i++){
        phi1 = phi1Range[0]+i*dphi1;
        for(int j=0;j<phi2Steps;j++){
            phi2 = phi2Range[0]+j*dphi2;
            upp1 = derivUpp(uEdit,i,j,true);
            upp2 = derivUpp(uEdit,i,j,false);
            k4[i][j]=uk(k+A[3]*dk, phi1, phi2, upp1, upp2);
        }
    }

    for(int i=0;i<phi1Steps;i++){
        for(int j=0;j<phi2Steps;j++){
            uEdit[i][j]=u[i][j]+B[4][0]*dk*k1[i][j]+B[4][1]*dk*k2[i][j]+B[4][2]*dk*k3[i][j]+B[4][3]*dk*k4[i][j];
        }
    }

    for(int i=0;i<phi1Steps;i++){
        phi1 = phi1Range[0]+i*dphi1;
        for(int j=0;j<phi2Steps;j++){
            phi2 = phi2Range[0]+j*dphi2;
            upp1 = derivUpp(uEdit,i,j,true);
            upp2 = derivUpp(uEdit,i,j,false);
            k5[i][j]=uk(k+A[4]*dk, phi1, phi2, upp1, upp2);
        }
    }

    for(int i=0;i<phi1Steps;i++){
        for(int j=0;j<phi2Steps;j++){
            uEdit[i][j]=u[i][j]+B[5][0]*dk*k1[i][j]+B[5][1]*dk*k2[i][j]+B[5][2]*dk*k3[i][j]+B[5][3]*dk*k4[i][j]+B[5][4]*dk*k5[i][j];
        }
    }

    for(int i=0;i<phi1Steps;i++){
        phi1 = phi1Range[0]+i*dphi1;
        for(int j=0;j<phi2Steps;j++){
            phi2 = phi2Range[0]+j*dphi2;
            upp1 = derivUpp(uEdit,i,j,true);
            upp2 = derivUpp(uEdit,i,j,false);
            k6[i][j]=uk(k+A[5]*dk, phi1, phi2, upp1, upp2);
        }
    }

    double TE = abs(dk*(CT[0]*k1[0][0]+CT[1]*k2[0][0]+CT[2]*k3[0][0]+CT[3]*k4[0][0]+CT[4]*k5[0][0]+CT[5]*k6[0][0]));
    double epsilon = releps*abs(u[0][0])+abseps;

    for(int i=0;i<phi1Steps;i++){
        for(int j=0;j<phi2Steps;j++){
            double TEtemp = abs(dk*(CT[0]*k1[i][j]+CT[1]*k2[i][j]+CT[2]*k3[i][j]+CT[3]*k4[i][j]+CT[4]*k5[i][j]+CT[5]*k6[i][j]));
            double epstemp = releps*abs(u[i][j])+abseps;

            if(TEtemp > TE){
                TE = TEtemp;
            }

            if(epstemp<epsilon){
                epsilon = epstemp;
            }
        }
    }

    double dknew = sF*dk*pow(epsilon/TE, 1./5);

    if(abs(dknew) < abs(minKstep)){
        dknew = -minKstep;
        TE=sF*epsilon;
    }

    if(isnan(dknew)){
        dk=sF*dk;
        for(int i=0;i<phi1Steps;i++){
            delete[] uEdit[i];
        }
        delete[] uEdit;
        if(!err){
            err = true;
        }
        return false;
    }

    if(TE>epsilon){
        dk = dknew;
        for(int i=0;i<phi1Steps;i++){
            delete[] uEdit[i];
        }
        delete[] uEdit;
        /*if(abs(dk)<abs(minKstep)){
            err = true;
        }*/
        return false;
    }else{
        for(int i=0;i<phi1Steps;i++){
            for(int j=0;j<phi2Steps;j++){
                double du = dk*(CH[0]*k1[i][j]+CH[1]*k2[i][j]+CH[2]*k3[i][j]+CH[3]*k4[i][j]+CH[4]*k5[i][j]+CH[5]*k6[i][j]);
                if(isnan(du)){
                    err = true;
                }
                u[i][j]+= du;
            }
        }

        dk=dknew;
        for(int i=0;i<phi1Steps;i++){
            delete[] uEdit[i];
        }
        delete[] uEdit;
        return true;
    }



}

void nextStepAdaptive(double **u, double k, double **k1, double **k2, double **k3, double **k4, double **k5, double **k6){
    bool complete = nextStepRK45(u, k, k1, k2, k3, k4, k5, k6);
    cout<<"\r"<<k<<" "<<dk<<flush;
    if(!complete and !err){
        nextStepAdaptive(u, k, k1, k2, k3, k4, k5, k6);
    }
}

int main()
{
    string fileName;

    int zeroJ = (int)round(-phi2Range[0]/dphi2);

/*    cout<<"Phi steps: ";
    cin>>phi1Steps;

    phi2Steps=phi1Steps;

    cout<<"k steps: ";
    cin>>kSteps;*/

    cout<<"T: ";
    cin>>T;

    if(T==0.){
        zeroTemp = true;
    }

    cout<<"File name: ";
    cin>>fileName;

    dphi1 = abs(phi1Range[1]-phi1Range[0])/phi1Steps;
    dphi2 = abs(phi2Range[1]-phi2Range[0])/phi2Steps;
    dk = -kStart/kSteps;

    double kCurrent = kStart;

    double** u = new double*[phi1Steps];
    double** k1 = new double*[phi1Steps];
    double** k2 = new double*[phi1Steps];
    double** k3 = new double*[phi1Steps];
    double** k4 = new double*[phi1Steps];
    double** k5 = new double*[phi1Steps];
    double** k6 = new double*[phi1Steps];

    for(int i=0;i<phi1Steps;i++){
        u[i] = new double[phi2Steps];
        k1[i] = new double[phi2Steps];
        k2[i] = new double[phi2Steps];
        k3[i] = new double[phi2Steps];
        k4[i] = new double[phi2Steps];
        k5[i] = new double[phi2Steps];
        k6[i] = new double[phi2Steps];
    }

    for(int i=0;i<phi1Steps;i++){
        double phi1;
        phi1 = phi1Range[0]+i*dphi1;
        for(int j=0;j<phi2Steps;j++){
            double phi2;
            phi2 = phi2Range[0]+j*dphi2;
            u[i][j] = initV(phi1,phi2);
        }
    }

    ofstream initOutput("pot.csv");

    for(int i=0;i<phi1Steps;i++){
        double phi1;
        phi1 = phi1Range[0]+i*dphi1;
        initOutput << u[i][zeroJ] << "," << phi1 << endl;
    }

    initOutput.close();

    kSteps = 0;
    while(kCurrent > 0){
        kCurrent += dk;
        if(!err){
            if(kCurrent > 0){
                nextStepAdaptive(u, kCurrent, k1, k2, k3, k4, k5, k6);
            }else{
                nextStepAdaptive(u, 0, k1, k2, k3, k4, k5, k6);
                cout<<endl<<"Simulation Complete!";
            }
        }
        kSteps+=1;
    }

    ofstream output(fileName+".csv");

    for(int i=0;i<phi1Steps;i++){
        double phi1;
        phi1 = phi1Range[0]+i*dphi1;
        output << u[i][zeroJ] << "," << phi1 << endl;
    }

    output.close();
    cout<<endl;

    if(isnan(u[0][0])){
        cout<<endl<<"Error: Check phiSteps^2/kSteps < 240"<<endl;
    }

    for(int i=0;i<phi1Steps;i++){
        delete[] u[i];
        delete[] k1[i];
        delete[] k2[i];
        delete[] k3[i];
        delete[] k4[i];
        delete[] k5[i];
        delete[] k6[i];
    }
    
    delete[] u;
    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] k4;
    delete[] k5;
    delete[] k6;

}