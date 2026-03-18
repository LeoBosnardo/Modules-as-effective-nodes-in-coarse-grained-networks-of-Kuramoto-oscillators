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
const int nmod = 2;
const int nosc1 = 17;
const int nosc2 = 17;
const int nosc = nosc1 + nosc2;

const double q1 = 0.9;
const double q2 = 0.9;

const double sigma = 1.0;
const double w = 1.5;

const double l_ini = 0.0;
const double l_fin = 3.0;

const double t_intervalo = 50.0;
const double t_intervaloMaior = 100.0;
const double dt = 0.001;

double lin = 1.0;
const double lin_incremento = 0.1;
const double lin_max = 20.0;

double l;

vector<double> omega(nosc, 0.0);
vector<double> k(nmod*nmod);
vector<int> mrhi(nosc1);
vector<int> officer(nosc2);
vector<vector<int>> A(nosc, vector<int>(nosc, 0));
//=======================================================================================================================================================
void kuramoto(const vector<double>& y, vector<double>& dydt, double t){

    dydt = omega;

    // population mr hi
    for (int i : mrhi) {

        double soma = 0.0;

        for (int j : mrhi){
            if (A[i][j] == 1){
                soma += sin(y[j] - y[i]);
            }
        }

        dydt[i] += lin * soma / k[0];

        soma = 0.0;

        for (int j : officer){
            if (A[i][j] == 1){
                soma += sin(y[j] - y[i]);
            }
        }

        dydt[i] += l * soma / (2.0 * k[1]);
    }

    // population officer
    for (int i : officer) {

        double soma = 0.0;

        for (int j : mrhi){
            if (A[i][j] == 1){
                soma += sin(y[j] - y[i]);
            }
        }

        dydt[i] += l * soma / (2.0 * k[2]);

        soma = 0.0;

        for (int j : officer){
            if (A[i][j] == 1){
                soma += sin(y[j] - y[i]);
            }
        }

        dydt[i] += lin * soma / k[3];
    }
}
//=======================================================================================================================================================
void para_ord(const vector<double>& y, double& r1, double& r2, double& r){

    complex<double> I(0.0, 1.0);

    complex<double> z1(0.0, 0.0), z2(0.0, 0.0), z(0.0, 0.0);

    for (int i : mrhi) {
        complex<double> zi = exp(I * y[i]);
        z1 += zi;
        z += zi;
    }

    for (int i : officer) {
        complex<double> zi = exp(I * y[i]);
        z2 += zi;
        z += zi;
    }

    z1 /= nosc1;
    z2 /= nosc2;
    z /= nosc;

    r1 = abs(z1);
    r2 = abs(z2);
    r = abs(z);
}
//=======================================================================================================================================================
struct OrderObserver {

    vector<double> &t_series;
    vector<double> &r_series;
    vector<double> &r1_series;
    vector<double> &r2_series;

    OrderObserver(vector<double>& t_,
                  vector<double>& r_,
                  vector<double>& r1_,
                  vector<double>& r2_)
        : t_series(t_), r_series(r_), r1_series(r1_), r2_series(r2_) {}

    void operator()(const vector<double>& y, double t)
    {
        double r1, r2, r;
        para_ord(y, r1, r2, r);

        t_series.push_back(t);
        r_series.push_back(r);
        r1_series.push_back(r1);
        r2_series.push_back(r2);
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
    file.open("k.txt");
    for (int i = 0; i < nmod*nmod; ++i) {
        file >> k[i];
    }
    file.close();

    // Read Mr Hi
    file.open("mr_hi.txt");
    for (int i = 0; i < nosc1; ++i) {
        file >> mrhi[i];
    }
    file.close();

    // Read Officer
    file.open("officer.txt");
    for (int i = 0; i < nosc2; ++i) {
        file >> officer[i];
    }
    file.close();

    int run_number = stoi(argv[1]);
    int l_steps = stoi(argv[2]);
    int l_number = stoi(argv[3]);

    // lambda    
    l = l_ini + l_number * (l_fin - l_ini) / l_steps;

    // Cria pasta w*/p*/ e w*/p*/l*/
    ostringstream folder_stream_1;
    ostringstream folder_stream_2;
    folder_stream_1 << "w" << std::fixed << setprecision(2) << w << "/run_" << run_number << "/";
    folder_stream_2 << "w" << std::fixed << setprecision(2) << w << "/run_" << run_number << "/l" << setprecision(2) << l << "/";
    string folder_1 = folder_stream_1.str();
    string folder_2 = folder_stream_2.str();
    filesystem::create_directories(folder_1);
    filesystem::create_directories(folder_2);

    // Cria arquivo rl.csv e rt.csv
    if (l_number == 1) filesystem::remove(folder_1 + "rl.csv");
    ofstream rl_file(folder_1 + "rl.csv", ios::app);
    ofstream rt_file(folder_2 + "rt.csv");
    
    // aleatorios
    random_device rd;
    default_random_engine generator(rd());
    normal_distribution<double> norm_dist(0.0, sigma); // Normal(μ=0,σ)
    uniform_real_distribution<double> uniform(0.0, 2.0 * M_PI); // Uniforme

    // vetor de posições
    vector<double> y(nosc);

    // integrador da parte de passo fixo
    runge_kutta4<vector<double>> fixed_stepper;

    // integrador da parte de passo variado
    auto varying_stepper = make_controlled(
        1e-8,        // abs tol
        1e-8,        // rel tol
        runge_kutta_dopri5<vector<double>>() // integrador
    );

    // vetores do observador
    vector<double> t_obs, r_t, r1_t, r2_t;
    OrderObserver observer(t_obs, r_t, r1_t, r2_t);

    // posição inicial
    for (int i = 0; i < nosc; i++) y[i] = 0.0;//uniform(generator);

    // frequencia
    //double mean_freq = 0.0;
    //for (int i = 0; i < nosc; i++) {
    //    omega[i] = norm_dist(generator);
    //    mean_freq += omega[i];
    //}
    //mean_freq = mean_freq / nosc;
    //for (int i = 0; i < nosc; i++) omega[i] -= mean_freq;
    for (int i : mrhi) omega[i] += w;

    // menor valor inicial de r1 e r2
    double min_r1 = 0;
    double min_r2 = 0;

    // tempo inicial
    double t = 0.0;

    // loop dentro de cada intervalo de lin
    while ((min_r1 < q1 || min_r2 < q2) && (lin < lin_max)){

        // tempos
        double t1 = t;
        double t2 = t + t_intervalo;

        // limpa vetores com parametros de ordem do intervalo
        t_obs.clear();
        r_t.clear();
        r1_t.clear();
        r2_t.clear();

        // integrador
        integrate_const(
            fixed_stepper,
            kuramoto,
            y,
            t1,
            t2,
            dt,
            observer);

        if (run_number == 1) {
            for (int i = 0; i < r_t.size(); i++) {
                rt_file << t_obs[i] << "," << lin << "," << r1_t[i] << "," << r2_t[i] << "," << r_t[i] << "\n";
            }
        }

        // menores valores de r1 e r2 no intervalo
        min_r1 = *min_element(r1_t.begin(), r1_t.end());
        min_r2 = *min_element(r2_t.begin(), r2_t.end());

        // aumenta lin caso não tenha atingido sincronia interna suficiente
        t = t2;
        lin += lin_incremento;

        //cout << l << " " << lin << " " << t << " " << r1_t.back() << " " << r2_t.back() << "\n";

    }

    if (lin == lin_max){
        // escreve valores no arquivo r_lambda.csv
        rl_file << l << "lin did not converge" << "\n";
        return 0;
    }

    // volta no lin que convergiu
    lin -= 0.1;

    // intevalo de tempo maior
    double t1 = t;
    double t2 = t + t_intervaloMaior;

    // limpa vetores
    t_obs.clear();
    r_t.clear();
    r1_t.clear();
    r2_t.clear();

    // integrador
    integrate_const(fixed_stepper, 
                    kuramoto,
                    y,
                    t1,
                    t2,
                    dt,
                    observer);

    if (run_number == 1){
        for (int i = 0; i < r_t.size(); i++) {
            rt_file << t_obs[i] << "," << lin << "," << r1_t[i] << "," << r2_t[i] << "," << r_t[i] << "\n";
        }
    }

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
    rl_file << l << "," << lin << "," << r1_t.back() << "," << r2_t.back() << "," << r_t.back() << "," << mean_r << "," << sigma_r << "\n";

    return 0;

}