#ifndef MANUFATURADAS_H
#define MANUFATURADAS_H

#include <vector>
#include <cmath>
#include <stdexcept>

namespace mms {

    // Soluções exatas manufaturadas (compatíveis com Neumann zero)
    inline double R_exata(double x, double y, double t,
                          double B_R, double kx, double ky) {
        return B_R * (2.0 + std::cos(kx * x) * std::cos(ky * y) * std::exp(-t));
    }

    inline double M_exata(double x, double y, double t,
                          double B_M, double kx, double ky) {
        return B_M * (2.0 + std::cos(kx * x) * std::cos(ky * y) * std::exp(-t));
    }

    inline double B_exata(double x, double y, double t,
                          double B_B, double kx, double ky) {
        return B_B * (2.0 + std::cos(kx * x) * std::cos(ky * y) * std::exp(-t));
    }

    // Derivadas temporais
    inline double dR_dt(double x, double y, double t,
                        double B_R, double kx, double ky) {
        return -B_R * std::cos(kx * x) * std::cos(ky * y) * std::exp(-t);
    }

    inline double dM_dt(double x, double y, double t,
                        double B_M, double kx, double ky) {
        return -B_M * std::cos(kx * x) * std::cos(ky * y) * std::exp(-t);
    }

    inline double dB_dt(double x, double y, double t,
                        double B_B, double kx, double ky) {
        return -B_B * std::cos(kx * x) * std::cos(ky * y) * std::exp(-t);
    }

    // Laplacianos
    inline double lap_R(double x, double y, double t,
                        double B_R, double kx, double ky) {
        return -(kx * kx + ky * ky) * B_R * std::cos(kx * x) * std::cos(ky * y) * std::exp(-t);
    }

    inline double lap_M(double x, double y, double t,
                        double B_M, double kx, double ky) {
        return -(kx * kx + ky * ky) * B_M * std::cos(kx * x) * std::cos(ky * y) * std::exp(-t);
    }

    inline double lap_B(double x, double y, double t,
                        double B_B, double kx, double ky) {
        return -(kx * kx + ky * ky) * B_B * std::cos(kx * x) * std::cos(ky * y) * std::exp(-t);
    }

    // Gradientes de M (para advecção)
    inline double dM_dx(double x, double y, double t,
                        double B_M, double kx, double ky) {
        return -kx * B_M * std::sin(kx * x) * std::cos(ky * y) * std::exp(-t);
    }

    inline double dM_dy(double x, double y, double t,
                        double B_M, double kx, double ky) {
        return -ky * B_M * std::cos(kx * x) * std::sin(ky * y) * std::exp(-t);
    }

    // Termos fonte (forçantes)
    inline double fonte_R(double x, double y, double t,
                          double B_R, double B_M, double kx, double ky,
                          double lambda_R, double mu_A, double mu_RM) {
        double R = R_exata(x, y, t, B_R, kx, ky);
        double M = M_exata(x, y, t, B_M, kx, ky);
        double dt = dR_dt(x, y, t, B_R, kx, ky);
        return dt - (lambda_R * R - mu_A * R - mu_RM * R * M);
    }

    inline double fonte_M(double x, double y, double t,
                          double B_R, double B_M, double B_B,
                          double kx, double ky,
                          double lambda_M, double rho, double mu_RM, double mu_MB,
                          double alpha_M, double v1, double v2) {
        double R = R_exata(x, y, t, B_R, kx, ky);
        double M = M_exata(x, y, t, B_M, kx, ky);
        double B = B_exata(x, y, t, B_B, kx, ky);
        double dt = dM_dt(x, y, t, B_M, kx, ky);
        double lap = lap_M(x, y, t, B_M, kx, ky);
        double grad_x = dM_dx(x, y, t, B_M, kx, ky);
        double grad_y = dM_dy(x, y, t, B_M, kx, ky);
        double adv = v1 * grad_x + v2 * grad_y;
        double reac = lambda_M * M * (1.0 - M / (rho * R)) - mu_MB * M * B + mu_RM * R * M;
        return dt - (alpha_M * lap - adv + reac);
    }

    inline double fonte_B(double x, double y, double t,
                          double B_M, double B_B,
                          double kx, double ky,
                          double lambda_B, double k_B, double mu_MB, double alpha_B) {
        double M = M_exata(x, y, t, B_M, kx, ky);
        double B = B_exata(x, y, t, B_B, kx, ky);
        double dt = dB_dt(x, y, t, B_B, kx, ky);
        double lap = lap_B(x, y, t, B_B, kx, ky);
        double reac = lambda_B * B * (1.0 - B / k_B) + mu_MB * M * B;
        return dt - (alpha_B * lap + reac);
    }

    // Preenchimento de vetores de fonte para a malha completa
    inline void preencher_fonte_R(std::vector<double>& f_R,
                                  int Nx, int Ny, double dx, double dy, double t,
                                  double B_R, double B_M, double kx, double ky,
                                  double lambda_R, double mu_A, double mu_RM) {
        f_R.resize(Nx * Ny);
        for (int i = 0; i < Nx; ++i) {
            double x = i * dx;
            for (int j = 0; j < Ny; ++j) {
                double y = j * dy;
                f_R[i * Ny + j] = fonte_R(x, y, t, B_R, B_M, kx, ky,
                                         lambda_R, mu_A, mu_RM);
            }
        }
    }

    inline void preencher_fonte_M(std::vector<double>& f_M,
                                  int Nx, int Ny, double dx, double dy, double t,
                                  double B_R, double B_M, double B_B,
                                  double kx, double ky,
                                  double lambda_M, double rho, double mu_RM, double mu_MB,
                                  double alpha_M, double v1, double v2) {
        f_M.resize(Nx * Ny);
        for (int i = 0; i < Nx; ++i) {
            double x = i * dx;
            for (int j = 0; j < Ny; ++j) {
                double y = j * dy;
                f_M[i * Ny + j] = fonte_M(x, y, t, B_R, B_M, B_B, kx, ky,
                                         lambda_M, rho, mu_RM, mu_MB,
                                         alpha_M, v1, v2);
            }
        }
    }

    inline void preencher_fonte_B(std::vector<double>& f_B,
                                  int Nx, int Ny, double dx, double dy, double t,
                                  double B_M, double B_B,
                                  double kx, double ky,
                                  double lambda_B, double k_B, double mu_MB, double alpha_B) {
        f_B.resize(Nx * Ny);
        for (int i = 0; i < Nx; ++i) {
            double x = i * dx;
            for (int j = 0; j < Ny; ++j) {
                double y = j * dy;
                f_B[i * Ny + j] = fonte_B(x, y, t, B_M, B_B, kx, ky,
                                         lambda_B, k_B, mu_MB, alpha_B);
            }
        }
    }

    // Preenchimento dos valores exatos na malha (em t)
    inline void solucao_exata(std::vector<double>& R_ex,
                              std::vector<double>& M_ex,
                              std::vector<double>& B_ex,
                              int Nx, int Ny, double dx, double dy, double t,
                              double B_R, double B_M, double B_B,
                              double kx, double ky) {
        R_ex.resize(Nx * Ny);
        M_ex.resize(Nx * Ny);
        B_ex.resize(Nx * Ny);
        for (int i = 0; i < Nx; ++i) {
            double x = i * dx;
            for (int j = 0; j < Ny; ++j) {
                double y = j * dy;
                std::size_t k = static_cast<std::size_t>(i) * Ny + j;
                R_ex[k] = R_exata(x, y, t, B_R, kx, ky);
                M_ex[k] = M_exata(x, y, t, B_M, kx, ky);
                B_ex[k] = B_exata(x, y, t, B_B, kx, ky);
            }
        }
    }

    // Cálculo de erros
    inline double erro_L2(const std::vector<double>& num,
                         const std::vector<double>& exata,
                         double dx, double dy) {
        if (num.size() != exata.size())
            throw std::invalid_argument("Vetores de tamanhos diferentes.");
        double soma = 0.0;
        for (std::size_t k = 0; k < num.size(); ++k) {
            double diff = num[k] - exata[k];
            soma += diff * diff;
        }
        return std::sqrt(soma * dx * dy);
    }

    inline double erro_Linf(const std::vector<double>& num,
                           const std::vector<double>& exata) {
        if (num.size() != exata.size())
            throw std::invalid_argument("Vetores de tamanhos diferentes.");
        double max_err = 0.0;
        for (std::size_t k = 0; k < num.size(); ++k) {
            double diff = std::abs(num[k] - exata[k]);
            if (diff > max_err) max_err = diff;
        }
        return max_err;
    }
    
    // Norma L2 de um vetor
    inline double norma_L2(const std::vector<double>& v, double dx, double dy) {
        double soma = 0.0;
        for (std::size_t k = 0; k < v.size(); ++k)
            soma += v[k] * v[k];
        return std::sqrt(soma * dx * dy);
    }
    
    // Norma L_inf de um vetor
    inline double norma_Linf(const std::vector<double>& v) {
        double max_val = 0.0;
        for (std::size_t k = 0; k < v.size(); ++k) {
            if (std::abs(v[k]) > max_val)
                max_val = std::abs(v[k]);
        }
        return max_val;
    }
    
    // Erro L2 relativo
    inline double erro_L2_relativo(const std::vector<double>& num,
                                  const std::vector<double>& exata,
                                  double dx, double dy) {
        double num_erro = erro_L2(num, exata, dx, dy);
        double den_norma = norma_L2(exata, dx, dy);
        return (den_norma > 0.0) ? (num_erro / den_norma) * 100.0 : 0.0;
    }
    
    // Erro L_inf relativo
    inline double erro_Linf_relativo(const std::vector<double>& num,
                                    const std::vector<double>& exata) {
        double num_erro = erro_Linf(num, exata);
        double den_norma = norma_Linf(exata);
        return (den_norma > 0.0) ? (num_erro / den_norma) * 100.0 : 0.0;
    }

} 

#endif
