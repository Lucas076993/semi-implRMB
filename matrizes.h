#include <vector>
#include <stdexcept>
#include <cmath>
#include "algebra.h"

// Gera o vetor fonte para o sistema linear de M
std::vector<double> geraFonteM(std::vector<double>& M, std::vector<double>& R_new, std::size_t N, std::size_t Nx, std::size_t Ny, double dx, double dy, double v1, double v2, double lambda_M, double mu_RM, double phi){
  std::vector<double> b_M(N);
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {
                std::size_t k = i * Ny + j;

                // Advecção upwind em x
                double adv_x = 0.0;
                if (v1 >= 0.0) {
                    if (i > 0)
                        adv_x = v1 * (M[k] - M[(i - 1) * Ny + j]) / dx;
                    else
                        adv_x = 0.0;      // fronteira oeste: fluxo nulo ⇒ dM/dx=0
                } else {
                    if (i < Nx - 1)
                        adv_x = v1 * (M[(i + 1) * Ny + j] - M[k]) / dx;
                    else
                        adv_x = 0.0;      // fronteira leste
                }

                // Advecção upwind em y
                double adv_y = 0.0;
                if (v2 >= 0.0) {
                    if (j > 0)
                        adv_y = v2 * (M[k] - M[i * Ny + (j - 1)]) / dy;
                    else
                        adv_y = 0.0;      // fronteira sul
                } else {
                    if (j < Ny - 1)
                        adv_y = v2 * (M[i * Ny + (j + 1)] - M[k]) / dy;
                    else
                        adv_y = 0.0;      // fronteira norte
                }
                double advec = adv_x + adv_y;

                // Lado direito
                b_M[k] = M[k] + phi * (-advec + lambda_M * M[k] + mu_RM * R_new[k] * M[k]);
            }
        }
        
        return b_M;
}

// Gera matriz de M
algebra::matriz<double> matriz_M(int Nx, int Ny, double dx, double dy, double phi, double alpha_M, double lambda_M, double rho, double mu_MB, const std::vector<double>& R_n, const std::vector<double>& M_n, const std::vector<double>& B_n) {
    const std::size_t N = static_cast<std::size_t>(Nx) * Ny;
    if (R_n.size() != N || M_n.size() != N || B_n.size() != N)
        throw std::invalid_argument("Tamanho dos vetores de estado incorreto.");

    algebra::matriz<double> A(Nx, Ny);

    const double coef_x = phi * alpha_M / (dx * dx);
    const double coef_y = phi * alpha_M / (dy * dy);

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            std::size_t k = static_cast<std::size_t>(i) * Ny + j;

            // Termo de reacção
            const double R_val = R_n[k];
            if (R_val <= 0.0) throw std::domain_error("R deve ser positivo.");
            double d = (lambda_M / (rho * R_val)) * M_n[k] + mu_MB * B_n[k];
            double diag_react = phi * d;

            // Contribuição da difusão para a diagonal
            double diff_diag = 0.0;
            double off_x = coef_x;
            double off_y = coef_y;

            // Vizinho leste (i+1) => supNy
            if (i < Nx - 1) {
                diff_diag += coef_x;
                A.supNy[k] = -coef_x;
            }
            // Vizinho oeste (i-1) => infNy
            if (i > 0) {
                diff_diag += coef_x;
                A.infNy[k] = -coef_x;
            }
            // Vizinho norte (j+1) => sup1
            if (j < Ny - 1) {
                diff_diag += coef_y;
                A.sup1[k] = -coef_y;
            }
            // Vizinho sul (j-1) => inf1
            if (j > 0) {
                diff_diag += coef_y;
                A.inf1[k] = -coef_y;
            }

            A.diag[k] = 1.0 + diff_diag + diag_react;
        }
    }
    return A;
}

algebra::matriz<double> matriz_B(int Nx, int Ny, double dx, double dy, double phi, double alpha_B, double lambda_B, double k_B, const std::vector<double>& B_n) {
    const std::size_t N = static_cast<std::size_t>(Nx) * Ny;
    if (B_n.size() != N)
        throw std::invalid_argument("Tamanho do vetor B_n incorreto.");

    algebra::matriz<double> A(Nx, Ny);

    const double coef_x = phi * alpha_B / (dx * dx);
    const double coef_y = phi * alpha_B / (dy * dy);

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            std::size_t k = static_cast<std::size_t>(i) * Ny + j;

            double reac = phi * (lambda_B / k_B) * B_n[k];
            double diff_diag = 0.0;

            if (i < Nx - 1) {
                diff_diag += coef_x;
                A.supNy[k] = -coef_x;
            }
            if (i > 0) {
                diff_diag += coef_x;
                A.infNy[k] = -coef_x;
            }
            if (j < Ny - 1) {
                diff_diag += coef_y;
                A.sup1[k] = -coef_y;
            }
            if (j > 0) {
                diff_diag += coef_y;
                A.inf1[k] = -coef_y;
            }

            A.diag[k] = 1.0 + diff_diag + reac;
        }
    }
    return A;
}

void inicializar(std::vector<double>& R, std::vector<double>& M, std::vector<double>& B, const int Nx, const int Ny, const double dx, const double dy, const double R0, const double M0, const double B0, const double dispM, const double dispB) {
    if (Nx <= 0 || Ny <= 0 || dx <= 0.0 || dy <= 0.0)
        throw std::invalid_argument("Dimensões da malha inválidas.");

    const std::size_t N = static_cast<std::size_t>(Nx) * Ny;
    R.assign(N, R0);          // uniforme
    M.resize(N);
    B.resize(N);

    // Centro da malha (considerando pontos de 0 a (Nx-1)*dx)
    const double cx = (Nx - 1) * dx / 2.0;
    const double cy = (Ny - 1) * dy / 2.0;

    for (int i = 0; i < Nx; ++i) {
        const double x = i * dx;
        const double dx2 = (x - cx) * (x - cx);
        for (int j = 0; j < Ny; ++j) {
            const double y = j * dy;
            const double dist2 = dx2 + (y - cy) * (y - cy);
            const std::size_t k = static_cast<std::size_t>(i) * Ny + j;
            M[k] = M0 * std::exp(-dispM * dist2);
            B[k] = B0 * std::exp(-dispB * dist2);
        }
    }
}
