#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <tuple>

// Tamaño de la malla (número de celdas)
const int Lx = 3;
const int Ly = 3;
const int numCells = Lx*Ly;

// Función para calcular el error cuadrático medio para una celda
double calculateError(double densityObserved, double densitySimulated) {
    double diff = densitySimulated - densityObserved;
    return diff * diff;
}

// Supongamos que esta función simula el proceso en cada celda para valores específicos de S y D
double simulateLatticeBoltzmannCell(int i, double S, double D) {
    // Código de simulación para una celda usando S y D
    // Retorna la densidad simulada para esa celda
    int ix = i/Lx;
    int iy = i%Lx;
    double rho = (ix+1)*(iy+1)*(S+D);
    return rho;
}

// Función para limitar los valores dentro de un intervalo
double clamp(double value, double min, double max) {
    return std::max(min, std::min(value, max));
}

// Función para ajustar S y D para cada celda
std::tuple<std::vector<double>, std::vector<double>>  optimizeParametersForCells(const std::vector<double>& densityObserved) {
    std::vector<double> S(numCells, 1.0); // Inicializa S para cada celda
    std::vector<double> D(numCells, 1.0); // Inicializa D para cada celda

    double learningRate = 0.0001;
    int maxIterations = 100000;
    double tolerance = 1e-6; // Criterio de parada basado en la tolerancia de error

    for (int iter = 0; iter < maxIterations; ++iter) {
        double totalError = 0.0;
        bool convergence = true;

        for (int cell = 0; cell < numCells; ++cell) {
            // Simula la densidad para la celda actual
            double densitySimulated = simulateLatticeBoltzmannCell(cell,S[cell], D[cell]);

            // Calcula el error para la celda actual
            double error = calculateError(densityObserved[cell], densitySimulated);
            totalError += error;

            // Calcula el gradiente de S y D para esta celda
            double delta = 0.0001;
            double S_gradient = (calculateError(densityObserved[cell], simulateLatticeBoltzmannCell(cell,S[cell] + delta, D[cell])) - error) / delta;
            double D_gradient = (calculateError(densityObserved[cell], simulateLatticeBoltzmannCell(cell,S[cell], D[cell] + delta)) - error) / delta;

            // Actualiza S y D para esta celda usando el gradiente
            double newS = S[cell] - learningRate * S_gradient;
            double newD = D[cell] - learningRate * D_gradient;

            // Clamping para mantener S y D dentro de sus intervalos permitidos
            S[cell] = clamp(newS, 1.0, 10.0);
            D[cell] = clamp(newD, 0.0, 4.0);

            // Si el cambio en S o D es menor que la tolerancia, marca
            if (std::abs(S_gradient) > tolerance || std::abs(D_gradient) > tolerance) {
                convergence = false;
            }

        }

        // Verifica el criterio de parada (error promedio)
        if (totalError / numCells < tolerance && convergence) {
            std::cout << "Convergencia alcanzada en iteración " << iter << std::endl;
            break;
        }

        // Opcional: Imprime el progreso
        if (iter % 100 == 0) {
            std::cout << "Iteración " << iter << ": Error promedio = " << totalError / numCells << std::endl;
        }
    }

    std::cout << "Optimización completada" << std::endl;
    std::cout << "Source:"<<"\n";
     for (double source : S) {
        std::cout << source << " ";
    }
    std::cout<<"\n";
    std::cout << "Diffusion:"<<"\n";
     for (double diff : D) {
        std::cout << diff << " ";
    }
    std::cout<<"\n";
    return std::make_tuple(S, D);
}

int main() {
  /*
 S=[1,10]
 D=[0,4]
 (0,1,3,4): S1= 8.2 , D1= 1.3
 (2,5,6,7,8): S2= 2.4 , D2= 3.9

 (ix+1)*(iy+1)*(S+D)
 */
  double S1 = 8.2;
  double S2 = 2.4;
  double D1 = 1.3;
  double D2 = 3.9;


  int rhoSize = Lx*Ly;
  std::vector<double> rhoObsVector(rhoSize);

  for(int i=0; i<rhoSize; i++){
    int ix = i/Lx;
    int iy = i%Lx;

    if(i==0 || i==1 || i==3 || i==4){
    rhoObsVector[i] = (ix+1)*(iy+1)*(S1+D1);
    }
    else{
    rhoObsVector[i] = (ix+1)*(iy+1)*(S2+D2);
    }
  }

    //std::vector<double> densityObserved(numCells, 1.0); // Ejemplo de densidades observadas

    std::vector<double> S(rhoSize);
    std::vector<double> D(rhoSize);
    std::vector<double> rhoSimulatedVector(rhoSize);
    std::tie(S, D) = optimizeParametersForCells(rhoObsVector);

    for(int i=0; i<rhoSize; i++){
    int ix = i/Lx;
    int iy = i%Lx;
    rhoSimulatedVector[i] = (ix+1)*(iy+1)*(S[i]+D[i]);
    }

    std::cout<<"Rho observado"<<"\n";
    for (double element : rhoObsVector) {
        std::cout << element << " ";
    }

    std::cout<<"\n Rho simulado"<<"\n";
    for (double element : rhoSimulatedVector) {
        std::cout << element << " ";
    }
               std::cout<<"\n";
    return 0;

}
