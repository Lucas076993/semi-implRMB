#ifndef ALGEBRA_H
#define ALGEBRA_H

namespace algebra{
  // Estrutura de matriz para facilitar abstração e contas 
  // Uso de flatten matrix para contiguidade de alocação e otimização
  template <typename T> struct matriz{
      // Dimensões da malha
        std::size_t Nx;   // número de pontos em x
        std::size_t Ny;   // número de pontos em y
        std::size_t linhas; // Nx * Ny
        std::size_t colunas; // Nx * Ny

        // Armazenamento das 5 diagonais (todas com tamanho linhas)
        // Ordem: k = i * Ny + j
        // diag       : A(k, k)
        // sup1       : A(k, k+1)   (válido se j < Ny-1, caso contrário 0)
        // inf1       : A(k, k-1)   (válido se j > 0)
        // supNy      : A(k, k+Ny)  (válido se i < Nx-1)
        // infNy      : A(k, k-Ny)  (válido se i > 0)
        std::vector<T> diag;
        std::vector<T> sup1;
        std::vector<T> inf1;
        std::vector<T> supNy;
        std::vector<T> infNy;

        // Construtor padrão (matriz vazia)
        matriz() : Nx(0), Ny(0), linhas(0), colunas(0) {}

        // Construtor que aloca e inicializa com zeros
        matriz(std::size_t nx, std::size_t ny)
            : Nx(nx), Ny(ny), linhas(nx * ny), colunas(nx * ny),
              diag(linhas, T(0)),
              sup1(linhas, T(0)),
              inf1(linhas, T(0)),
              supNy(linhas, T(0)),
              infNy(linhas, T(0)) {}

        // Acesso a elementos (pode ser lento, útil para depuração)
        T operator()(std::size_t i, std::size_t j) const {
            if (i >= linhas || j >= colunas)
                throw std::out_of_range("Índice fora dos limites.");
            if (i == j) return diag[i];
            if (i + 1 == j && i % Ny != Ny - 1) return sup1[i];
            if (i == j + 1 && j % Ny != Ny - 1) return inf1[j]; // ou inf1[i]
            if (i + Ny == j) return supNy[i];
            if (i == j + Ny) return infNy[j]; // ou infNy[i]
            return T(0);
        }
    };

  // produto matriz vetor
  template <typename T> std::vector<T> operator*(const matriz<T> &A, const std::vector<T> &x){
    if (A.colunas != x.size())
            throw std::invalid_argument("Dimensões matriz/vetor incompatíveis.");

        std::vector<T> y(A.linhas, T(0));
        for (std::size_t k = 0; k < A.linhas; ++k) {
            T soma = A.diag[k] * x[k];
            // vizinho norte (k+1)
            if (k % A.Ny != A.Ny - 1)   // não é a última coluna
                soma += A.sup1[k] * x[k + 1];
            // vizinho sul (k-1)
            if (k % A.Ny != 0)          // não é a primeira coluna
                soma += A.inf1[k] * x[k - 1];
            // vizinho leste (k+Ny)
            if (k + A.Ny < A.linhas)
                soma += A.supNy[k] * x[k + A.Ny];
            // vizinho oeste (k-Ny)
            if (k >= A.Ny)
                soma += A.infNy[k] * x[k - A.Ny];
            y[k] = soma;
        }
        return y;
  }
  
  // produto interno canônico entre dois vetores
  template <typename T> T operator*(const std::vector<T> &y, const std::vector<T> &x){
    T soma = 0;
    
    for(std::size_t i = 0; i < x.size(); i++)
      soma += y[i] * x[i];
    
    return soma;
  }
  
  // produto escalar vetor
  template <typename T> std::vector<T> operator*(const T escalar, const std::vector<T> &x){
    std::vector<T> y(x.size());
    
    for(std::size_t i = 0; i < x.size(); i++)
      y[i] = escalar * x[i];
    return y;
  }
  
  
  // diferença entre vetores
  template <typename T> std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b){
    std::vector<T> c(b.size());    
    for(std::size_t i = 0; i < a.size(); i++)
      c[i] = a[i] - b[i];
      
    return c;
  }
  
  // soma entre vetores
  template <typename T> std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b){
    std::vector<T> c(b.size());    
    for(std::size_t i = 0; i < a.size(); i++)
      c[i] = a[i] + b[i];
      
    return c;
  }
}

#endif
