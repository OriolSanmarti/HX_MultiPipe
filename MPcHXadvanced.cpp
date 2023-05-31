#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "properties.h"
#include "eHX.h"
#include "MPcHXadvanced.h"
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


void MPcHXadvanced(double T1in, double &T1o, double T2in, double &T2o, string fluid1, string fluid2){
        
    // Numerical parameters
    double inct = 10;
    double t = 0;
    double tmax = 50000;
    int N = 100;
    double error = 0.0001;

    // Geometrical parameters
    double L = 1.5;      // m
    double Ntt = 12;     // -
    double Di = 0.02245; // m 
    double Do = 0.02667; // m
    double De = 0.15;
    double incz = L / (double)N;
    double S1 = Ntt*M_PI*Di*Di/4.0;
    double S2 = M_PI*De*De/4.0 - Ntt*M_PI*Do*Do/4.0;
    double Ss =  Ntt * ( M_PI*Do*Do/4.0 - M_PI*Di*Di/4.0);
    double V1 = S1 * incz;
    double V2 = S2 * incz;
    double Vs = Ss * incz;
    double Pi = Ntt * M_PI * Di;
    double Po = Ntt * M_PI * Do;
    double Ai = Ntt * M_PI * Di * incz;
    double Ao = Ntt * M_PI * Do * incz;
    cout << " Ntt: " << Ntt << " Di: " << Di << " Do: " << Do << " De: " << De << " Di: " << Di << " Do: " << Di <<endl;
    cout << " S1: " << S1 << " S2: " << S2 << " Ss: " << Ss << " V1: " << V1 << " V2: " << V2 << " Vs: " << Vs << endl;
    double D2h = 4*S2/(Po);
    cout << " D2h: " << D2h << endl;


    // Operational parameters
    double P1in = 100000;
    double P2in = 100000;
    double v1 = 2.7;
    double mdot2 = 0.01;
    double mdot1 = CalcRhoAir(T1in, P1in)*S1*v1;
    double Caudalair = v1*S1*3600;
    cout << " Densitat aire: " << CalcRhoAir(T1in, P1in) << " Caudal aire: " << v1*S1*3600<<endl;
    
    // Initialize fields
    vector<double> T1(N, T1in);
    vector<double> T1ant(N, T1in);
    vector<double> T2(N, T2in);
    vector<double> T2ant(N, T2in);
    vector<double> Ts(N, T2in/2.0);
    vector<double> Tsant(N, T2in/2.0);
    vector<double> Tstdma(N, T2in/2.0);
    vector<double> alphai(N, 0.0);
    vector<double> alphao(N, 0.0);
    vector<double> P1(N, P1in);
    vector<double> P2(N, P2in);
    vector<double> v1a(N, v1);

    bool timefin = false;
    while(!timefin)
    {
        bool itefin = false;
        while(!itefin)
        {
            itefin = true;
            // Calc aire
            for(int i = 0; i<N; i++)
            {
                // Initialize properties
                double cp1 = CalcCpAir(T1[i]);
                double rho1 = CalcRhoAir(T1[i], P1[i]);
                double mu1 = CalcMuAir(T1[i]);
                double mu1w = CalcMuAir(Ts[i]);
                double lambda1 = CalcLambdaAir(T1[i]);

                // Calculate alpha air
                double v1i = v1;
                double v1 = mdot1 / (rho1*S1);
                v1a[i] = (mdot1 / (rho1*S1));
                double Re1 = rho1 * v1 * Di / mu1;
                double Pr1 = mu1*cp1/lambda1;
                //double Nu1 = 0.023*pow(Re1, 0.8)*pow(Pr1, 0.4)*1;
                double Nu1 = 0.027*pow(Re1, 0.8)*pow(Pr1, 0.33)*pow(mu1/mu1w, 0.14);
                //cout << " Re1: " << Re1 << " Pr1: " << Pr1 << " Nu1: " << Nu1 << endl;
                alphai[i] = Nu1*lambda1/Di;
                
                double Tn = T1in;
                if(i == 0) Tn = T1in;
                else Tn = T1[i-1];

                double Afi = rho1*cp1/inct*V1 + mdot1*cp1 + alphai[i]*Ai;
                double Bfi = -mdot1*cp1; 
                double Cfi = rho1*cp1/inct*V1*T1ant[i] + alphai[i]*Ai*Ts[i];
                T1[i] = (Cfi - Tn*Bfi) / (Afi);

                // Calculate Plost
                if(i > 0){
                    double rug = 0.05/1000.0;
                    double er1 = rug/Di;
                    double f1 = 0.0625 / pow((log10(er1/3.7 + 5.74/(pow(Re1, 0.9)))),2.0);
                    double v1mid = (v1a[i]);
                    double tau1 = f1*rho1*v1mid*v1mid/2.0;
                    double incp1 = ( (tau1*M_PI*Di*incz) + (mdot1/Ntt*(v1a[i]-v1a[i-1])) )/ (M_PI*Di*Di/4.0);
                    P1[i] = P1[i-1] - incp1;
                    //cout << "i : " << i << " v1mid: " << v1mid << " mdto1: " << mdot1 << " incp1: " << incp1 << " tau1: " << tau1 << " f1: " << f1 <<endl;
                }
            }
            
            // Calc MS
            for(int i = N-1; i>= 0; i--)
            {
                // Initialize properties
                double cp2 = calcCp(T2[i], fluid2);
                double rho2 = calcRho(T2[i], fluid2);
                double mu2 = calcMu(T2[i], fluid2);
                double mu2w = calcMu(Ts[i], fluid2);
                double lambda2 = calcLambda(T2[i], fluid2);

                // Calculate alpha molten salt
                double v2 = mdot2 / (rho2*S2);
                double Re2 = rho2 * v2 * D2h / mu2;
                double Pr2 = mu2*cp2/lambda2;
                double Gz2 = Re2*Pr2*D2h/L;
                //double Nu2 = 0.023*pow(Re2, 0.8)*pow(Pr2, 0.4)*1;
                double Nu2 = 3.66 + ((0.668*Re2*Pr2*(D2h/L)) / (1+0.40*(Re2*Pr2*pow((D2h/L),(2.0/3.0)))));
                //double Nu2 = 1.86*pow(Re2, 1.0/3.0)*pow(Pr2, 1.0/3.0)*(pow(mu2/mu2w, 0.14))*pow((D2h/L),(1.0/3.0));
                //double Nu2 = 3.66;
                cout << " Re2: " << Re2 << " Pr2: " << Pr2 << " Gz2: " << Gz2 << " Nu2: " << Nu2 << endl;
                alphao[i] = Nu2*lambda2/D2h;

                double Tsud = 0.0;
                if(i == N-1) Tsud = T2in;
                else Tsud = T2[i+1];
            
                double Afo = rho2*cp2/inct*V2 + mdot2*cp2 + alphao[i]*Ao;
                double Bfo = -mdot2*cp2;
                double Cfo = rho2*cp2/inct*V2*T2ant[i] + alphao[i]*Ao*Ts[i];
                T2[i] = (Cfo - Bfo*Tsud) / Afo;
            }

            // Calc Solid
            // Initialize properties
            double rhos = 7800;
            double cps = 500;
            double lambdas = 16.3;
        
            vector<double> Ptdma(N,0.0);
            vector<double> Rtdma(N,0.0);
            for(int i = 0; i<N; i++)
            {
                double ui = 0.0, ai = 0.0, li = 0.0, bi = 0.0;
                
                if(i == 0)
                {
                    ui = -lambdas*Ss/incz;
                    ai = rhos*cps/inct*Vs + lambdas*Ss/incz + alphai[i]*Ai + alphao[i]*Ao;
                    li = 0.0;
                    bi = rhos*cps/inct*Vs*Tsant[i] + alphai[i]*Ai*T1[i] + alphao[i]*Ao*T2[i];

                    Ptdma[i] = ui / ai;
                    Rtdma[i] = bi / ai;
                }
                else if(i == N-1)
                {
                    ui = 0.0;
                    ai = rhos*cps/inct*Vs + lambdas*Ss/incz + alphai[i]*Ai + alphao[i]*Ao;
                    li = -lambdas*Ss/incz;
                    bi = rhos*cps/inct*Vs*Tsant[i] + alphai[i]*Ai*T1[i] + alphao[i]*Ao*T2[i];

                    Rtdma[i] = (bi - li*Rtdma[i-1]) / (ai - li*Ptdma[i-1]);
                }
                else
                {
                    ui = -lambdas*Ss/incz;
                    ai = rhos*cps/inct*Vs + 2*lambdas*Ss/incz + alphai[i]*Ai + alphao[i]*Ao;
                    li = -lambdas*Ss/incz;
                    bi = rhos*cps/inct*Vs*Tsant[i] + alphai[i]*Ai*T1[i] + alphao[i]*Ao*T2[i];

                    Ptdma[i] = ui / (ai - li*Ptdma[i-1]);
                    Rtdma[i] = (bi - li*Rtdma[i-1]) / (ai - li*Ptdma[i-1]);
                }
            }
            for(int i = N-1; i>=0; i--)
            {   
                if(i == N-1) Ts[i] = Rtdma[i];
                else Ts[i] = Rtdma[i] - Ptdma[i]*Ts[i+1];

                if(fabs(Ts[i]-Tstdma[i]) > error) itefin = false;
            }
            
            Tstdma = Ts;
        }
        
        // Actualize fields
        T1ant = T1;
        T2ant = T2;
        Tsant = Ts;
        t += inct;

        // Check temporal convergence
        if(t>tmax) timefin = true;
        cout << "Temps: " << t << "\r"<<flush;
        //for(int i = 0; i<N; i++) cout << "z: " << incz*i << " T1: " << T1[i] << " T2: " << T2[i] << " Ts: " << Ts[i] << endl;
    }


    // PostProcess
    for(int i = 0; i<N; i++) cout << "z: " << incz*i << " T1: " << T1[i] << " Ts: " << Ts[i] << " T2: " << T2[i] << 
                            " alphai: " << alphai[i] << " alphao: " << alphao[i] <<endl;

    T2o = T2[0];
    T1o = T1[N-1];

    cout << "*********************** End cHX ********************" << endl;
    cout << "T1o: " << T1o << " T2o: " << T2o << " Twcritica: " << Ts[0] << " Caudal air: " << Caudalair << " P1o: " << P1[N-1] << " P1lost: " << P1[N-1] - P1[0] <<endl;
}

