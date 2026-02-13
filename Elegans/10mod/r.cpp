//=======================================================================================================================================================
#include <vector>
#include <cmath>
#include <complex>
#include <random>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;
//=======================================================================================================================================================
const int nmod = 10;
const int nosc1 = 68;
const int nosc2 = 33;
const int nosc3 = 76;
const int nosc4 = 11;
const int nosc5 = 13;
const int nosc6 = 7;
const int nosc7 = 6;
const int nosc8 = 8;
const int nosc9 = 8;
const int nosc10 = 18;
const int nosc = nosc1 + nosc2 + nosc3 + nosc4 + nosc5 + nosc6 + nosc7 + nosc8 + nosc9 + nosc10;

const double w1 = 0.0;
const double w2 = 1.0;
const double w3 = -3.5;
const double w4 = 2.0;
const double w5 = -2.0;
const double w6 = 0.6;
const double w7 = 0.4;
const double w8 = 1.2;
const double w9 = 1.5;
const double w10 = -1.5;

const double l_ini = 0.0;
const double l_fin = 10.0;

const double t1 = 0.0;
const double t2 = 500.0;
const double t3 = 1000.0;
const double dt = 0.01;

const double lin = 20.0;

double l;

vector<double> omega(nosc, 0.0);
vector<vector<double>> k(nmod, vector<double>(nmod));
vector<int> m1(nosc1);
vector<int> m2(nosc2);
vector<int> m3(nosc3);
vector<int> m4(nosc4);
vector<int> m5(nosc5);
vector<int> m6(nosc6);
vector<int> m7(nosc7);
vector<int> m8(nosc8);
vector<int> m9(nosc9);
vector<int> m10(nosc10);
vector<int> module(nosc);
vector<vector<int>> A(nosc, vector<int>(nosc, 0));
//=======================================================================================================================================================
void kuramoto(const vector<double>& y, vector<double>& dydt, double t){

    dydt = omega;

    for (int i = 0; i < nosc; i++){

        double soma = 0.0;

        int mod_i = module[i];

        for (int j = 0; j < nosc; j++){

            int mod_j = module[j];

            if (A[i][j] == 1){

                if (mod_i == mod_j){

                    soma += lin * sin(y[j] - y[i]) / k[mod_i-1][mod_j-1];

                } else {

                    soma += l * sin(y[j] - y[i]) / (2.0 * k[mod_i-1][mod_j-1]);

                }

            }

        }
        
        dydt[i] += soma;

    }

}
//=======================================================================================================================================================
void para_ord(const vector<double>& y, double& r){

    complex<double> I(0.0, 1.0);

    complex<double> z(0.0, 0.0);

    for (int i = 0; i < nosc; ++i) {
        complex<double> zi = exp(I * y[i]);
        z += zi;
    }

    z /= nosc;

    r = abs(z);
}
//=======================================================================================================================================================
struct OrderObserver {

    vector<double> &t_series;
    vector<double> &r_series;

    OrderObserver(vector<double>& t_,
                  vector<double>& r_)
        : t_series(t_), r_series(r_) {}

    void operator()(const vector<double>& y, double t)
    {
        double r;
        para_ord(y, r);

        t_series.push_back(t);
        r_series.push_back(r);
    }
};
//=======================================================================================================================================================
int main(int argc, char* argv[]){

    // Read the matrix
    ifstream file("A.txt");
    for (int i = 0; i < nosc; ++i) {
        for (int j = 0; j < nosc; ++j) {
            file >> A[i][j];
        }
    }
    file.close();

    // Read the normalizations
    file.open("K.txt");
    for (int i = 0; i < nmod; ++i) {
        for (int j = 0; j < nmod; ++j) {
            file >> k[i][j];
        }
    }
    file.close();

    // Read m1
    file.open("m1.txt");
    for (int i = 0; i < nosc1; ++i) {
        file >> m1[i];
    }
    file.close();

    // Read m2
    file.open("m2.txt");
    for (int i = 0; i < nosc2; ++i) {
        file >> m2[i];
    }
    file.close();

    // Read m3
    file.open("m3.txt");
    for (int i = 0; i < nosc3; ++i) {
        file >> m3[i];
    }
    file.close();

    // Read m4
    file.open("m4.txt");
    for (int i = 0; i < nosc4; ++i) {
        file >> m4[i];
    }
    file.close();

    // Read m5
    file.open("m5.txt");
    for (int i = 0; i < nosc5; ++i) {
        file >> m5[i];
    }
    file.close();

    // Read m6
    file.open("m6.txt");
    for (int i = 0; i < nosc6; ++i) {
        file >> m6[i];
    }
    file.close();

    // Read m7
    file.open("m7.txt");
    for (int i = 0; i < nosc7; ++i) {
        file >> m7[i];
    }
    file.close();

    // Read m8
    file.open("m8.txt");
    for (int i = 0; i < nosc8; ++i) {
        file >> m8[i];
    }
    file.close();

    // Read m9
    file.open("m9.txt");
    for (int i = 0; i < nosc9; ++i) {
        file >> m9[i];
    }
    file.close();

    // Read m10
    file.open("m10.txt");
    for (int i = 0; i < nosc10; ++i) {
        file >> m10[i];
    }
    file.close();

    // Read module
    file.open("module.txt");
    for (int i = 0; i < nosc; ++i) {
        file >> module[i];
    }
    file.close();

    int l_steps = stoi(argv[1]);
    int l_number = stoi(argv[2]);

    // lambda    
    l = l_ini + l_number * (l_fin - l_ini) / l_steps;

    // Cria arquivo rl.csv e rt.csv
    if (l_number == 1) filesystem::remove("rl.csv");
    ofstream rl_file("rl.csv", ios::app);
    
    // aleatorios
    random_device rd;
    default_random_engine generator(rd());
    uniform_real_distribution<double> uniform(0.0, 2.0 * M_PI); // Uniforme

    // vetor de posições
    vector<double> y(nosc);

    // integrador da parte de passo fixo
    runge_kutta4<vector<double>> fixed_stepper;

    // vetores do observador
    vector<double> t_obs, r_t;
    OrderObserver observer(t_obs, r_t);

    // posição inicial
    for (int i = 0; i < nosc; i++) y[i] = uniform(generator);

    // frequencia
    for (int i : m1) omega[i] += w1;
    for (int i : m2) omega[i] += w2;
    for (int i : m3) omega[i] += w3;
    for (int i : m4) omega[i] += w4;
    for (int i : m5) omega[i] += w5;
    for (int i : m6) omega[i] += w6;
    for (int i : m7) omega[i] += w7;
    for (int i : m8) omega[i] += w8;
    for (int i : m9) omega[i] += w9;
    for (int i : m10) omega[i] += w10;

    // integrador
    integrate_const(fixed_stepper, 
                    kuramoto,
                    y,
                    t1,
                    t2,
                    dt,
                    observer);

    // limpa vetores
    t_obs.clear();
    r_t.clear();

    // integrador
    integrate_const(fixed_stepper, 
                    kuramoto,
                    y,
                    t2,
                    t3,
                    dt,
                    observer);

    // media
    double mean_r = accumulate(r_t.begin(), r_t.end(), 0.0) / r_t.size();

    // variancia
    double var = 0.0;
    for (double x : r_t) {
        double d = x - mean_r;
        var += d * d;
    }

    // desvio padrão
    double sigma_r = sqrt(var / (r_t.size() - 1));

    // escreve valores no arquivo r_lambda.csv
    rl_file << l << "," << lin << "," << r_t.back() << "," << mean_r << "," << sigma_r << "\n";

    return 0;

}