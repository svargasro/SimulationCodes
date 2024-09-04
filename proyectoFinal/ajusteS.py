import numpy as np
from scipy.optimize import minimize

#Tamaño del sistema.
Lx=3
Ly=3
num_cells = Lx*Ly

#Función que hace el papel de la simulación de LatticeBoltzmann
def simulate_lattice_boltzmann_cell(i, S):
    rho = S*(i+1)
    return rho

#Cálculo del error cuadrático para una celda.
def calculate_error(density_observed, density_simulated):
    diff = density_simulated - density_observed
    return diff**2

# Función objetivo para la optimización
def objective(params, density_observed):
    S = params[:num_cells]
    total_error = 0.0
    for cell in range(num_cells): #Se suman los errores cuadráticos asociados a cada celda.
        density_simulated = simulate_lattice_boltzmann_cell(cell,S[cell])
        total_error += calculate_error(density_observed[cell], density_simulated)
    return total_error

#Ajuste de S.
def optimize_parameters_for_cells(density_observed):

    S_initial = np.full(num_cells, 1.0)   #Valores iniciales para S

    bounds = [(1.0, 10.0)] * num_cells    #Restricciones sobre S

    # Optimización usando el método L-BFGS-B (que soporta restricciones de parámetros)
    result = minimize(objective, S_initial, args=(density_observed), method='L-BFGS-B', bounds=bounds)

    # Extraer los valores optimizados de S
    S_opt = result.x[:num_cells]

    return S_opt



rhoSize = Lx*Ly;
rhoObsVector = np.zeros(rhoSize)


Sreal= np.random.uniform(1.0, 10.0, size=rhoSize) #Se crea un arreglo de S, que serán los valores reales.

#Se construye un vector de densidades observadas.
for i in range(0,rhoSize):
    rhoObsVector[i]= Sreal[i]*(i+1)

S_optimized = optimize_parameters_for_cells(rhoObsVector) #Se llama a la función de optimización.

rhoAdjusted = np.zeros(rhoSize)
for i in range(0,rhoSize):
    rhoAdjusted[i] = S_optimized[i]*(i+1)

print("\nValores reales:\n")
print("Valores de S reales:\n",Sreal)
print("rho_Obs:\n",rhoObsVector)
print("---------------------------------------------\n")
print("Valores ajustados:\n")
print("Valores optimizados de S:\n", S_optimized)
print("rho_Adjusted:\n",rhoAdjusted)
