#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <functional>

#include "yaml-cpp/yaml.h"
#include "algebra.h"    
#include "matrizes.h"   
#include "CG.h"
#include "manufaturadas.cpp"

// Exportação de um campo escalar 2D em formato CSV (t, x, y, valor)
void escrever_campo_csv(std::ostream& out, const std::vector<double>& campo, int Nx, int Ny, double dx, double dy, double t) {
    out << std::scientific << std::setprecision(12);
    for (int i = 0; i < Nx; ++i) {
        double x = i * dx;
        for (int j = 0; j < Ny; ++j) {
            double y = j * dy;
            std::size_t k = static_cast<std::size_t>(i) * Ny + j;
            out << t << ',' << x << ',' << y << ',' << campo[k] << '\n';
        }
    }
}


int main() {
    YAML::Node config = YAML::LoadFile("params.yaml");

    // Parâmetros do modelo
    const double lambda_R = config["model_parameters"]["lambda_R"].as<double>();
    const double mu_A     = config["model_parameters"]["mu_A"].as<double>();
    const double mu_RM    = config["model_parameters"]["mu_RM"].as<double>();
    const double alpha_M  = config["model_parameters"]["alpha_M"].as<double>();
    const double v1       = config["model_parameters"]["v1"].as<double>();
    const double v2       = config["model_parameters"]["v2"].as<double>();
    const double lambda_M = config["model_parameters"]["lambda_M"].as<double>();
    const double rho      = config["model_parameters"]["rho"].as<double>();
    const double mu_MB    = config["model_parameters"]["mu_MB"].as<double>();
    const double alpha_B  = config["model_parameters"]["alpha_B"].as<double>();
    const double lambda_B = config["model_parameters"]["lambda_B"].as<double>();
    const double k_B      = config["model_parameters"]["k_B"].as<double>();

    // Malha e discretização
    const std::size_t Nx = config["mesh"]["Nx"].as<std::size_t>();
    const std::size_t Ny = config["mesh"]["Ny"].as<std::size_t>();
    const double Lx      = config["mesh"]["Lx"].as<double>();
    const double Ly      = config["mesh"]["Lx"].as<double>();
    const double dx      = Lx / (Nx - 1);
    const double dy      = Ly / (Ny - 1);

    // Configuração do CG
    const std::size_t max_iter =  config["CG"]["max_iter"].as<std::size_t>();
    const double tol           = config["CG"]["tol"].as<double>();

    // Tempo
    const std::size_t max_passos = config["time"]["max_steps"].as<std::size_t>();
    const double dt              = config["time"]["dt"].as<double>();

    //  Condições iniciais
    const double R0_val = config["initial_conditions"]["R0_val"].as<double>();
    const double M0_val = config["initial_conditions"]["M0_val"].as<double>();
    const double B0_val = config["initial_conditions"]["B0_val"].as<double>();
    const double dispM  = config["initial_conditions"]["dispM"].as<double>();   // largura da gaussiana de M
    const double dispB  = config["initial_conditions"]["dispB"].as<double>();   // largura da gaussiana de B
    
    // Arquivos de saída CSV
    std::string R_path     = config["output"]["R"].as<std::string>();
    std::string M_path     = config["output"]["M"].as<std::string>();
    std::string B_path     = config["output"]["B"].as<std::string>();
    std::string error_path = config["output"]["error"].as<std::string>();
    bool error_trigger     = config["output"]["error_trigger"].as<bool>();
    
    // Parâmetros das soluções exatas
    const double B_R = config["mms"]["B_R"].as<double>();
    const double B_M = config["mms"]["B_M"].as<double>();
    const double B_B = config["mms"]["B_B"].as<double>();
    const double PI = std::acos(-1.0);
    const double kx  = PI / Lx;   
    const double ky  = PI /Ly;

    // Função denominadora NSFD
    const double Lambda_max = std::max({mu_A, lambda_M, lambda_B});
    const double phi = (std::exp(Lambda_max * dt) - 1.0) / Lambda_max;

    std::vector<double> R, M, B, b_M;

    if(error_trigger)
      mms::solucao_exata(R, M, B, Nx, Ny, dx, dy, 0.0, B_R, B_M, B_B, kx, ky);
    else
      inicializar(R, M, B, Nx, Ny, dx, dy, R0_val, M0_val, B0_val, dispM, dispB);

    const std::size_t N = Nx * Ny;

    
    std::ofstream fR(R_path), fM(M_path), fB(B_path), fE(error_path);
    if (!fR || !fM || !fB || !fE) {
        std::cerr << "Erro ao abrir arquivos de saída.\n";
        return 1;
    }
    fR << "t,x,y,R\n";
    fM << "t,x,y,M\n";
    fB << "t,x,y,B\n";
    fE << "t,Erro Relativo R L2,Erro Relativo M L2,Erro Relativo B L2,Erro Relativo R L_inf,Erro Relativo M L_inf,Erro Relativo B L_inf\n";

    // Laço de evolução temporal
    double t = 0.0;

    // Estado inicial (t = 0)
    escrever_campo_csv(fR, R, Nx, Ny, dx, dy, t);
    escrever_campo_csv(fM, M, Nx, Ny, dx, dy, t);
    escrever_campo_csv(fB, B, Nx, Ny, dx, dy, t);
    std::cerr << "Passo 0 salvo.\n";

    for (int n = 0; n < max_passos; ++n) {
        t += dt;

        // Atualização explícita de R
        std::vector<double> R_new(N);
        std::vector<double> f_R(N); // vetor nulo
        
        if(error_trigger)
          mms::preencher_fonte_R(f_R, Nx, Ny, dx, dy, t+dt, B_R, B_M, kx, ky, lambda_R, mu_A, mu_RM);
        
        for (std::size_t k = 0; k < N; ++k) {
            double num = (1.0 + phi * lambda_R) * R[k] + phi * f_R[k];
            double denom = 1.0 + phi * mu_A + phi * mu_RM * M[k];
            R_new[k] = num / denom;
        }

        // Sistema linear para M (difusão implícita)
        // Monta a matriz pentadiagonal
        algebra::matriz<double> A_M = matriz_M(Nx, Ny, dx, dy, phi, alpha_M, lambda_M, rho, mu_MB, R, M, B);
        
        // Constroi o lado direito b_M (inclui advecção explícita upwind)
        b_M = geraFonteM(M, R_new, N, Nx, Ny, dx, dy, v1, v2, lambda_M, mu_RM, phi);
        
        if(error_trigger){
          std::vector<double> f_M;
          mms::preencher_fonte_M(f_M, Nx, Ny, dx, dy, t+dt, B_R, B_M, B_B, kx, ky, lambda_M, rho, mu_RM, mu_MB, alpha_M, v1, v2);
        for (std::size_t k = 0; k < N; ++k)
          b_M[k] += phi * f_M[k];
        }

        // Resolve o sistema para M usando CG (chute inicial = M atual)
        std::vector<double> M_new = CG(A_M, b_M, M, max_iter, tol);

        // Sistema linear para B (difusão implícita)
        auto A_B = matriz_B(Nx, Ny, dx, dy, phi, alpha_B, lambda_B, k_B, B);
        std::vector<double> b_B(N);
        for (std::size_t k = 0; k < N; ++k)
            b_B[k] = B[k] + phi * (lambda_B * B[k] + mu_MB * M_new[k] * B[k]);
            
        if(error_trigger){
          std::vector<double> f_B;
          mms::preencher_fonte_B(f_B, Nx, Ny, dx, dy, t+dt, B_M, B_B, kx, ky, lambda_B, k_B, mu_MB, alpha_B);
          for (std::size_t k = 0; k < N; ++k)
            b_B[k] += phi * f_B[k];
        }
        
        std::vector<double> B_new = CG(A_B, b_B, B, max_iter, tol);

        // Atualiza os vetores para o próximo passo
        R.swap(R_new);
        M.swap(M_new);
        B.swap(B_new);

        // Salva o estado
        escrever_campo_csv(fR, R, Nx, Ny, dx, dy, t);
        escrever_campo_csv(fM, M, Nx, Ny, dx, dy, t);
        escrever_campo_csv(fB, B, Nx, Ny, dx, dy, t);
        std::cerr << "Passo " << n + 1 << " salvo (t = " << t << ").\n";
        
        if (error_trigger) {
          // Vetores para a solução exata no instante t
          std::vector<double> R_ex, M_ex, B_ex;
          mms::solucao_exata(R_ex, M_ex, B_ex, Nx, Ny, dx, dy, t,
                             B_R, B_M, B_B, kx, ky);

          // Cálculo dos erros (norma L2 e/ou Linf)
          double eR_L2 = mms::erro_L2_relativo(R, R_ex, dx, dy);
          double eM_L2 = mms::erro_L2_relativo(M, M_ex, dx, dy);
          double eB_L2 = mms::erro_L2_relativo(B, B_ex, dx, dy);
          double eR_Linf = mms::erro_Linf_relativo(R, R_ex);
          double eM_Linf = mms::erro_Linf_relativo(M, M_ex);
          double eB_Linf = mms::erro_Linf_relativo(B, B_ex);

          // Escrever em fE
          fE << t
             << ',' << eR_L2 << ',' << eM_L2 << ',' << eB_L2
             << ',' << eR_Linf << ',' << eM_Linf << ',' << eB_Linf
             << '\n';
        }
        
    }

    // Fecha os arquivos
    fR.close();
    fM.close();
    fB.close();
    fE.close();
    std::cerr << "Simulação finalizada com sucesso.\n";
    return 0;
}
