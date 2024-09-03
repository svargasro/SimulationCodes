import numpy as np
from scipy.optimize import minimize

# Tamaño de la malla (número de celdas)

Lx=3
Ly=3
num_cells = Lx*Ly

# Simulación de Lattice-Boltzmann para una celda específica
def simulate_lattice_boltzmann_cell(i, S, D):
    # Aquí se incluiría el código de simulación que retorna la densidad simulada
    ix = int(i/Lx);
    iy = i%Lx;
    rho = (ix+1)*(iy+1)*(S+D);
    return rho;

# Función de error cuadrático medio para una celda específica
def calculate_error(density_observed, density_simulated):
    diff = density_simulated - density_observed
    return diff**2

# Función objetivo para la optimización
def objective(params, density_observed):
    S = params[:num_cells]
    D = params[num_cells:]

    total_error = 0.0
    for cell in range(num_cells):
        density_simulated = simulate_lattice_boltzmann_cell(cell,S[cell], D[cell])
        total_error += calculate_error(density_observed[cell], density_simulated)

    return total_error

# Ajuste de parámetros S y D
def optimize_parameters_for_cells(density_observed):
    # Valores iniciales para S y D
    S_initial = np.full(num_cells, 1.0)
    D_initial = np.full(num_cells, 1.0)

    initial_params = np.concatenate([S_initial, D_initial])

    # Restricciones de los parámetros
    bounds = [(1.0, 10.0)] * num_cells + [(0.0, 4.0)] * num_cells

    # Optimización usando el método L-BFGS-B (que soporta restricciones de parámetros)
    result = minimize(objective, initial_params, args=(density_observed,),
                      method='L-BFGS-B', bounds=bounds)

    # Extraer los valores optimizados de S y D
    S_opt = result.x[:num_cells]
    D_opt = result.x[num_cells:]

    return S_opt, D_opt

# Ejemplo de densidades observadas

S1 = 8.2;
S2 = 2.4;
D1 = 1.3;
D2 = 3.9;

rhoSize = Lx*Ly;
rhoObsVector = np.zeros(rhoSize)

Sreal=np.zeros(rhoSize)
Dreal=np.zeros(rhoSize)

for i in range(0,rhoSize):
    ix=int(i/Lx)
    iy=i%Lx
    if (i==0 or i==1 or i==3 or i==4):
        rhoObsVector[i]= (ix+1)*(iy+1)*(S1+D1)
        Sreal[i]=S1
        Dreal[i]=D1
    else:
        rhoObsVector[i] = (ix+1)*(iy+1)*(S2+D2)
        Sreal[i]=S2
        Dreal[i]=D2

# Realizar la optimización
S_optimized, D_optimized = optimize_parameters_for_cells(rhoObsVector)
print("Valores optimizados de S:", S_optimized)
print("Valores optimizados de D:", D_optimized)
print("Valores de S reales:",Sreal)
print("Valores de D reales:",Dreal)
rhoSimulated = np.zeros(rhoSize)
for i in range(0,rhoSize):
    ix=int(i/Lx)
    iy=i%Lx
    rhoSimulated[i] = (ix+1)*(iy+1)*(S_optimized[i]+D_optimized[i])




print("rhoObservado")
print(rhoObsVector)
print("rhoSimulado")
print(rhoSimulated)
