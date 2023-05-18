#include <iostream>
#include <vector>
#include <cmath>
#include "properties.h"
#include "MPcHX.h"
using namespace std;

// Air properties
double CalcRhoAir(double T, double P){
	//return 1.1774;
    T += 273.15;
    return P/(287*T);
}
double CalcCpAir(double T){
    return 1031.5 - 0.21*T + 4.143E-4*T*T;
	//return 1005.7;
}
double CalcMuAir(double T){
    T = T + 273.15;
    return (1.458E-6*pow(T, 1.5))/(T+110.40);
}
double CalcLambdaAir(double T){
	return 0.026;
}
void HXmultipipe(double T1in, double T2in, double &T1o, double &T2o, string fluid1, string fluid2, double mdot, double vair){
    cout << endl << endl;
    cout << "*********************************** HX MultiPipes counter flow ***********************************" << endl;
    cout << endl << endl;
    cout << "mdot: " << mdot << " vair: " << vair << endl<<endl;
    cout << "T1in: " << T1in << " mdot: " << mdot << endl;
    int N = 100; // Numero nodes
    double T1ini = T1in;
    double T2ini = T2in;
    double P1ini = 100000;
    double P2ini = 5E6;
    vector <double> T1(N, T1ini);
    vector <double> T1ant(N, T1ini);
    vector <double> T1w(N, T1ini);
    vector <double> T2(N, T2ini);
    vector <double> T2w(N, T2ini);
    vector <double> T2ant(N, T2ini);
    vector <double> P1(N, P1ini);
    vector <double> P2(N, P2ini);
    //vector <double> v1(N, 0.0);
    //vector <double> v2(N, 0.0);
    double Tw1 = T1ini;
    double Tw2 = T1ini;
    double Tw1min = Tw1;
    double Tw2min = Tw2;
    double mdot2 = mdot;
    double lambdas = 16.3;

    // Numerical parameters
    double t = 0.0;
    double tmax = 50000.20;
    double inct = 1;

    // Geometrical parameters
    double L = 1.5;
    double incz = L/(double)N;
    double D = 0.15;
    double Nt = 12;
    //double Di = 0.014;
    //double Do = 0.017;
    double Di = 0.025;
    double Do = 0.0282;
    double Ai = Nt * M_PI * Di * incz;
    double Ao = Nt * M_PI * Do * incz; 
    double S1 = Nt * M_PI * Di * Di / 4.0;
    double S2 = M_PI * D * D / 4.0 - Nt * M_PI * Do * Do / 4.0;
    double V1 = S1 * incz;
    double V2 = S2 * incz;
    double UO = 0.0;
    double Po = M_PI * Do;
    double D2h = 4*S2/(Po*Nt);

    // Calculate mdot1
    double mdot1 = CalcRhoAir(T1in, P1ini) * vair * S1;
    cout << "mdot1: " << mdot1 << endl;
    cout << "S1: " << S1 << " S2: " << S2 << " V1: " << V1 << " V2: " << V2 << " Ai: " << Ai << " Ao: " << Ao << endl;
    double caudal = vair * S1 * 3600;

    while(t<tmax)
    {
        bool itefin = false;
        while(!itefin)
        {
            itefin = true;
            for(int i = 0; i<N; i++)
            {
                // Initialize properties
                double cp1 = CalcCpAir(T1[i]);
                double rho1 = CalcRhoAir(T1[i], P1[i]);
                double mu1 = CalcMuAir(T1[i]);
                double mu1w = CalcMuAir(T1w[i]);
                double lambda1 = CalcLambdaAir(T1[i]);
                double cp2 = calcCp(T2[i], fluid1);
                double rho2 = calcRho(T2[i], fluid1);
                double mu2 = calcMu(T2[i], fluid1);
                double mu2w = calcMu(T2w[i], fluid1);
                double lambda2 = calcLambda(T2[i], fluid1);

                // Calculate alpha air
                double v1 = mdot1 / (rho1*S1);
                double v1i = v1;
                double Re1 = rho1 * v1 * Di / mu1;
                double Pr1 = mu1*cp1/lambda1;
                //double Nu1 = 0.023*pow(Re1, 0.8)*pow(Pr1, 0.4)*1;
                double Nu1 = 0.027*pow(Re1, 0.8)*pow(Pr1, 0.33)*pow(mu1/mu1w, 0.14);
                double alphai = Nu1*lambda1/Di;

                // Calculate alpha molten salt
                double v2 = mdot2 / (rho2*S2);
                double Re2 = rho2 * v2 * D2h / mu2;
                double Pr2 = mu2*cp2/lambda2;
                //double Nu2 = 0.023*pow(Re2, 0.8)*pow(Pr2, 0.4)*1;
                //double Nu2 = 1.86*pow(Re2, 1.0/3.0)*pow(Pr2, 1.0/3.0)*(pow(mu2/mu2w, 0.14))*pow((Do/L),(1.0/3.0));
                double Nu2 = 1.86*pow(Re2, 1.0/3.0)*pow(Pr2, 1.0/3.0)*(pow(mu2/mu2w, 0.14))*pow((D2h/L),(1.0/3.0));
                double alphao = Nu2*lambda2/D2h;

                // Calculate UO
                UO = pow( ((1.0/alphai)*Ao/Ai + (Ao/(incz*2*M_PI*lambdas))*log(Do/Di) + (1.0/alphao)), -1.0);
                //UO = pow( ((1.0/alphai) + (1.0/alphao)), -1.0);
                //cout << "Re1: " << Re1 << " Pr1: " << Pr1 <<" alphai: " << alphai << endl;
                //cout << "Re2: " << Re2 << " Pr2: " << Pr2 <<" alphao: " << alphao << " Nu2: " << Nu2 << endl;
                //cout << "UO: " << UO << endl;

                // Declarate some useful variables
                double T1actual = T1[i];
                double T2actual = T2[i];
                double T2n = 0.0;
                double T1s = 0.0;
                if(i == 0) T1s = T1in;
                else T1s = T1[i-1];
                if(i == N-1) T2n = T2in;
                else T2n = T2[i+1];
            
                // 1) Air side
                double Afi = rho1*cp1/inct*V1 + mdot1*cp1 + UO*Ao; 
                double Bfi = -mdot1*cp1;
                double Cfi = rho1*cp1/inct*V1*T1ant[i] + UO*Ao*T2[i];

                T1[i] = (Cfi - Bfi*T1s) / Afi;

                // 2) Molten Salt Side
                double Afo = rho2*cp2/inct*V2 + mdot2*cp2 + UO*Ao; 
                double Bfo = -mdot2*cp2;
                double Cfo = rho2*cp2/inct*V2*T2ant[i] + UO*Ao*T1[i];

                T2[i] = (Cfo - Bfo*T2n) / Afo;

                // Calculate Tw
                double Q = UO * Ao * (T1[i] - T2[i]);
                T1w[i] = T1[i] - Q / (Ai*alphai);
                T2w[i] = T2[i] + Q / (Ao*alphao);
                //cout << "Tw1: " << T1w[i] << " Tw2: " << T2w[i] << endl;

                // Calculate Plost
                double rug = 0.05/1000.0;
                double er1 = rug/Di;
                double f1 = 0.0625 / pow((log10(er1/3.7 + 5.74/(pow(Re1, 0.9)))),2.0);
                double v1mid = (v1+v1i) / 2.0;
                double tau1 = f1*rho1*v1mid*v1mid/2.0;
                double incp1 = ( (tau1*M_PI*Di*incz) + (mdot1/Nt*(v1-v1i)) )/ (M_PI*Di*Di/4.0);
                if (i > 0) P1[i] = P1[i-1] - incp1;
                else P1[i] = 100000;
                //cout << "i: " << P1[i] << endl;
                


                // Check error
                double error = 0.00001;
                if(fabs(T1actual - T1[i]) > error) itefin = false;
                if(fabs(T2actual - T2[i]) > error) itefin = false;

                // Print results
                //cout << "i: " << i << " T1: " << T1[i] << " T2: " << T2[i] << endl;
            }
        }
        T1ant = T1;
        T2ant = T2;
        t+=inct;
        //cout << "***************** t: " << t << " [s] " << endl;
        //cout << "T1o: " << T1[N-1] << " T2o: " << T2[0] << endl;
    }

    T1o = T1[N-1];
    T2o = T2[0];
    double P1o = P1[N-1];
    double P2o = P2[0];
    cout << "P1o: " << P1o << " P2o: " << P2[N-1] << " Plost: " << P1[N-1] - P1[0] <<endl;
    cout << "T1o: " << T1o << " T2o: " << T2o << " Tw2min: " << T2w[0] << " Tw1max: " << T1w[N-1] << endl;

    for(int i = 0; i<N; i++)
    {
        //cout << "T1: " << T1[i] << " Tw1: " << T1w[i] << " Tw2: " << T2w[i] << " T2: " << T2[i] << endl;
    }

    cout << "Caudal aire: " << caudal << " D2h: " << D2h << endl;
    cout << endl << endl;
    cout << "********************************* end HX MultiPipes counter flow *********************************" << endl;
    cout << endl << endl;
}

