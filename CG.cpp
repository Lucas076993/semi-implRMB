#include <iostream>
#include <vector>

#include "CG.h"
#include "algebra.h"



std::vector<double> CG(const algebra::matriz<double> &A, const std::vector<double> &b, const std::vector<double> &x0, const std::size_t max_iter, const double tol){
  using algebra::operator+;
  using algebra::operator-;
  using algebra::operator*;
  
  // Verificação de dimensionalidade
  if(A.linhas != A.colunas || A.linhas != b.size() || A.colunas != x0.size()){
    std::cerr << "Dimensão inválida!!!\nRetornando valor inicial..." << std::endl;
    return x0;
  }
  // Declaração de variáveis
  std::vector<double> x(x0);
  std::vector<double> r(x0.size());
  std::vector<double> z(x0.size());
  double s, aux, aux1, m;
  
  
  // Cálculo inicial do resíduo: r = b - Ax
  r = b - A * x;
  std::vector<double> v(r);
  aux = r * r;
  
  for(std::size_t k = 0; k < max_iter; k++){
    z = A * v;
    s = aux / (v * z);
    
    x = x + s * v;
    r = r - s * z;
    
    aux1 = r * r;
    
    if(aux1 < tol){
      std::cerr << "Convergência no passo k = " << k << std::endl;
      break;
    }
    
    m = aux1 / aux;
    aux = aux1;
    v = r + m * v;
  }
  
  return x;
}
