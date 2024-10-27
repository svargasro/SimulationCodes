#include "/home/sergio/Documentos/ceres-solver-2.2.0/include/ceres/ceres.h"
#include <vector>
#include <iostream>

// Convertir simulate_density en plantilla
template <typename T>
std::vector<T> simulate_density(T S) {
    std::vector<T> simulated_density = {S * T(1.0), S * T(0.9), S * T(0.8), S * T(1.1), S * T(1.2)};
    return simulated_density;
}

// Estructura que define la función de costo
struct DensityResidual {
    DensityResidual(const std::vector<double>& observed_density)
        : observed_density_(observed_density) {}

    template <typename T>
    bool operator()(const T* const S, T* residual) const {
        // Simular la densidad con el valor actual de S
        std::vector<T> simulated_density = simulate_density(*S);

        // Calcular el residual como la suma de diferencias cuadráticas
        T sum_residuals = T(0.0);
        for (size_t i = 0; i < observed_density_.size(); ++i) {
            T diff = simulated_density[i] - T(observed_density_[i]);
            sum_residuals += diff * diff;
        }
        residual[0] = sum_residuals;
        return true;
    }

private:
    const std::vector<double>& observed_density_;
};

int main() {
    // Datos de densidad observada para cada celda
    std::vector<double> density_observed = {1.0, 0.9, 0.8, 1.1, 1.2};

    int num_cells = density_observed.size();

    // Valores iniciales para S
    std::vector<double> S(num_cells, 1.0); // Valores iniciales para cada celda

    ceres::Problem problem;

    for (int i = 0; i < num_cells; ++i) {
        ceres::CostFunction* cost_function =
            new ceres::AutoDiffCostFunction<DensityResidual, 1, 1>(
                new DensityResidual(density_observed));

        // Cada celda tiene su propio parámetro S
        problem.AddResidualBlock(cost_function, nullptr, &S[i]);

        // Establecer límites para S
        problem.SetParameterLowerBound(&S[i], 0, 1.0);
        problem.SetParameterUpperBound(&S[i], 0, 10.0);
    }

    // Configuración del solucionador
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    // Imprimir los resultados
    std::cout << "Final Report:\n";
    std::cout << summary.BriefReport() << "\n";

    for (int i = 0; i < num_cells; ++i) {
        std::cout << "Optimized S[" << i << "] = " << S[i] << std::endl;
    }

    return 0;
}
