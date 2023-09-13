#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <nlohmann/json.hpp>
#include "thread"
#include "progress.hpp"

using namespace std;

using json = nlohmann::json;

class QSEA{
    public:

        void run(int runID, json params, progress::progBar& bar, json& output)
        {   
            T = params["T"][runID].get<double>();
            N = params["N"].get<int>();
            double Gamma = params["Gamma"].get<double>();
            double vev = params["v"].get<double>();

            m2 = -2*pow(Gamma,4)/pow(vev,2);
            lambda = 12*pow(Gamma/vev,4);

            kStart = params["params"]["kStart"].get<double>();
            phiSteps = params["params"]["phiSteps"].get<int>();
            vector<double> phiLimits = params["params"]["phiRange"].get<vector<double>>();
            abseps = pow(10, params["params"]["logAbsErr"].get<double>());
            releps = pow(10, params["params"]["logRelErr"].get<double>());

            phiRange[0] = phiLimits[0];
            phiRange[1] = phiLimits[1];

            if(T==0){
                zeroTemp = true;
            }

            dphi = abs(phiRange[1]-phiRange[0])/phiSteps;
            dk = -kStart/kSteps;

            
            double kCurrent = kStart;

            double* u = new double[phiSteps];
            double* k1 = new double[phiSteps];
            double* k2 = new double[phiSteps];
            double* k3 = new double[phiSteps];
            double* k4 = new double[phiSteps];
            double* k5 = new double[phiSteps];
            double* k6 = new double[phiSteps];

            json fieldValues = json::array();
            json treeLevelPot = json::array();

            for(int i=0;i<phiSteps;i++){
                double phi;
                phi = phiRange[0]+i*dphi;
                u[i] = initV(phi);

                fieldValues.push_back(phi);
                treeLevelPot.push_back(u[i]);
            }

#pragma omp critical
            {
            output["phi"] = fieldValues;
            output["treeLevel"] = treeLevelPot;
            }

            kSteps = 0;
            while(kCurrent > 0){
                kCurrent += dk;
                if(kCurrent > 0){
                    nextStepAdaptive(u, kCurrent, k1, k2, k3, k4, k5, k6, bar);
                }else{
                    nextStepAdaptive(u, 0, k1, k2, k3, k4, k5, k6, bar);
                }
                kSteps+=1;
            }

            json finiteTempOutput;
            finiteTempOutput["pot"] = json::array();

            finiteTempOutput["T"] = T;

            for(int i=0;i<phiSteps;i++){
                finiteTempOutput["pot"].push_back(u[i]);
            }

#pragma omp critical
            {
            output["runs"].push_back(finiteTempOutput);
            }
            
            delete[] u;
            delete[] k1;
            delete[] k2;
            delete[] k3;
            delete[] k4;
            delete[] k5;
            delete[] k6;

        }

    private:
        double prevProg = 0;
        double m2 = -0.014;
        double lambda = 0.185;
        double alpha = 0;
        const double pi = 3.14159265358979323846;

        double phiRange[2] = {-15, 15};
        double kStart=10;

        int phiSteps = 10000;
        int kSteps = 1000;

        int kStop=kSteps;

        double T = 0;

        int N = 1;

        double dphi = abs(phiRange[1]-phiRange[0])/phiSteps;
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
        double abseps = pow(10, -8);

        /*"Safety Factor" determines the relative size of the next step compared to the maximum it's permitted to be.
        Used to avoid numerical errors. */
        double sF = 0.9;

        bool zeroTemp = false;
        bool err = false;

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

        double initV(double phi){
            return m2/2*phi*phi+alpha*pow(phi,3)/6+lambda*pow(phi,4)/24;
        }

        void nanerror(double quant, string message){
            if(isnan(quant)){
                if(!err){
                    cout<<endl<<message<<endl;
                }
                err = true;
            }
        }

        double sumV4(double phi){
            return (2+N)/3*lambda;
        }

        double uk(double k, double phi, double upp){
            double kt2 = k*k-m2-alpha*phi-lambda*phi*phi/2;
            double Et2 = kt2+upp;

            double deriv;

            if(zeroTemp){
                deriv = N*k*heaviside(kt2)/sumV4(phi)*(-Et2+sqrt(pow(Et2,2)+pow(kt2,2)*sumV4(phi)/(16*pi*pi)));
            }else{
                deriv = N*k*heaviside(kt2)/sumV4(phi)*(-Et2+sqrt(pow(Et2,2)+T*sumV4(phi)*matsubaraSum(kt2)/(3*pi*pi)));
            }
            
            return deriv;
        }

        double derivUpp(double *u, int i){
            double upp;

            if(i==0){
                upp=(2*u[0]-5*u[1]+4*u[2]-u[3])/pow(dphi,3);
            }else if(i==phiSteps-1){
                upp=(2*u[i]-5*u[i-1]+4*u[i-2]-u[i-3])/pow(dphi,3);
            }else{
                upp=(u[i+1]+u[i-1]-2*u[i])/pow(dphi,2);
            }

            return upp;
        }

        void nextStep(double *u, double k, double *k1, double *k2, double *k3, double *k4){
            //Old RK4 method. No longer used.
            double upp;
            double phi;
            double* uEdit = new double[phiSteps];

            for(int i=0;i<phiSteps;i++){
                phi = phiRange[0]+i*dphi;
                upp = derivUpp(u,i);
                k1[i]=uk(k, phi, upp);
            }

            for(int j=0;j<phiSteps;j++){
                uEdit[j]=u[j]+dk*k1[j]/2;
            }

            for(int i=0;i<phiSteps;i++){
                phi = phiRange[0]+i*dphi;
                upp = derivUpp(uEdit,i);
                k2[i]=uk(k+dk/2, phi, upp);
            }

            for(int j=0;j<phiSteps;j++){
                uEdit[j]=u[j]+dk*k2[j]/2;
            }

            for(int i=0;i<phiSteps;i++){
                phi = phiRange[0]+i*dphi;
                upp = derivUpp(uEdit,i);
                k3[i]=uk(k+dk/2, phi, upp);
            }

            for(int j=0;j<phiSteps;j++){
                uEdit[j]=u[j]+dk*k3[j];
            }

            for(int i=0;i<phiSteps;i++){
                phi = phiRange[0]+i*dphi;
                upp = derivUpp(uEdit,i);
                k4[i]=uk(k+dk, phi, upp);
            }

            for(int i=0;i<phiSteps;i++){
                u[i]+= dk*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;
                
            }

            delete[] uEdit;

        }

        bool nextStepRK45(double *u, double k, double *k1, double *k2, double *k3, double *k4, double *k5, double *k6){
            //New RK45 adaptive method that recalculates the step if the error is too large.
            double upp;
            double phi;
            double* uEdit = new double[phiSteps];

            for(int i=0;i<phiSteps;i++){
                phi = phiRange[0]+i*dphi;
                upp = derivUpp(u,i);
                k1[i]=uk(k+A[0]*dk, phi, upp);
            }

            for(int j=0;j<phiSteps;j++){
                uEdit[j]=u[j]+B[1][0]*dk*k1[j];
            }

            for(int i=0;i<phiSteps;i++){
                phi = phiRange[0]+i*dphi;
                upp = derivUpp(uEdit,i);
                k2[i]=uk(k+A[1]*dk, phi, upp);
            }

            for(int j=0;j<phiSteps;j++){
                uEdit[j]=u[j]+B[2][0]*dk*k1[j]+B[2][1]*dk*k2[j];
            }

            for(int i=0;i<phiSteps;i++){
                phi = phiRange[0]+i*dphi;
                upp = derivUpp(uEdit,i);
                k3[i]=uk(k+A[2]*dk, phi, upp);
            }

            for(int j=0;j<phiSteps;j++){
                uEdit[j]=u[j]+B[3][0]*dk*k1[j]+B[3][1]*dk*k2[j]+B[3][2]*dk*k3[j];
            }

            for(int i=0;i<phiSteps;i++){
                phi = phiRange[0]+i*dphi;
                upp = derivUpp(uEdit,i);
                k4[i]=uk(k+A[3]*dk, phi, upp);
            }

            for(int j=0;j<phiSteps;j++){
                uEdit[j]=u[j]+B[4][0]*dk*k1[j]+B[4][1]*dk*k2[j]+B[4][2]*dk*k3[j]+B[4][3]*dk*k4[j];
            }

            for(int i=0;i<phiSteps;i++){
                phi = phiRange[0]+i*dphi;
                upp = derivUpp(uEdit,i);
                k5[i]=uk(k+A[4]*dk, phi, upp);
            }

            for(int j=0;j<phiSteps;j++){
                uEdit[j]=u[j]+B[5][0]*dk*k1[j]+B[5][1]*dk*k2[j]+B[5][2]*dk*k3[j]+B[5][3]*dk*k4[j]+B[5][4]*dk*k5[j];
            }

            for(int i=0;i<phiSteps;i++){
                phi = phiRange[0]+i*dphi;
                upp = derivUpp(uEdit,i);
                k6[i]=uk(k+A[5]*dk, phi, upp);
            }

            double TE = abs(dk*(CT[0]*k1[0]+CT[1]*k2[0]+CT[2]*k3[0]+CT[3]*k4[0]+CT[4]*k5[0]+CT[5]*k6[0]));
            double epsilon = releps*abs(u[0])+abseps;

            for(int i=0;i<phiSteps;i++){
                double TEtemp = abs(dk*(CT[0]*k1[i]+CT[1]*k2[i]+CT[2]*k3[i]+CT[3]*k4[i]+CT[4]*k5[i]+CT[5]*k6[i]));
                double epstemp = releps*abs(u[i])+abseps;

                if(TEtemp > TE){
                    TE = TEtemp;
                }

                if(epstemp<epsilon){
                    epsilon = epstemp;
                }
            }

            double dknew = sF*dk*pow(epsilon/TE, 1./5);


            if(TE>epsilon){
                dk = dknew;
                delete[] uEdit;
                return false;
            }else{
                for(int i=0;i<phiSteps;i++){
                    u[i]+= dk*(CH[0]*k1[i]+CH[1]*k2[i]+CH[2]*k3[i]+CH[3]*k4[i]+CH[4]*k5[i]+CH[5]*k6[i]);
                }
                dk=dknew;
                delete[] uEdit;
                return true;
            }



        }

        void nextStepAdaptive(double *u, double k, double *k1, double *k2, double *k3, double *k4, double *k5, double *k6, progress::progBar& bar){
            //Runs the step and checks it's within the error tolerance.
            bool complete = nextStepRK45(u, k, k1, k2, k3, k4, k5, k6);
            
            if(complete){
#pragma omp critical
                {
                bar.add(100*((kStart-k)/kStart- prevProg));
                }
                prevProg = (kStart-k)/kStart;
            }else if(!complete){
                //If not, then it calls the step again.
                nextStepAdaptive(u, k, k1, k2, k3, k4, k5, k6, bar);
            }
        }
};

json readFile(std::string fileName, std::string errorType){
    std::ifstream file;
    file.open(fileName);

    if(file.fail()){
        throw std::runtime_error("error loading " + errorType + " file!");
    }

    json data = json::parse(file);

    file.close();

    return data;
}

int main(){
    json params = readFile("qsea-params.json", "input");
    json output;
    output["runs"] = json::array();

    string outputFileName = params["outputFileName"].get<string>();
    vector<double> temps = params["T"].get<vector<double>>();

    vector<QSEA> qseaVector;

    progress::progBar bar(100*6);

    thread progress;
    progress = bar.start();

    for(int i = 0; i < temps.size(); i++){
        QSEA qsea;
        qseaVector.push_back(qsea);
    }

#pragma omp parallel for shared(bar, qseaVector, temps, params, output)
    for(int i = 0; i < 6; i++){
        qseaVector[i].run(i, params, bar, output);

    }
    
    bar.finish();
    progress.join();


    ofstream outFile(outputFileName+".json");

    outFile << setw(4) << output << endl;

    outFile.close();

    return 0;
}