# ========================================================
# SIMULACIÓN DE UNA CURVA ODMR DE UN CENTRO NV EN DIAMANTE
# ========================================================

import numpy as np                                  # librería 1: nos permite trabajar con vectores / cálculos numéricos
import matplotlib.pyplot as plt                     # librería 2: nos permite dibujar las gráficas

# DATOS INICIALES
D     = 2.87                                                    # Frecuencia de resonancia del centro NV (GHz)
N     = 1000                                                    # Número de puntos del barrido
B     = float(input('Introduce el campo magnético (mT): '))     # Campo magnético introducido por el usuario (mT)

if B < 0:                                                       # Pequeña comprobación
    print('El campo magnético no puede ser menor que cero. ')
    exit()

B     = B / 1000                                                # Pasamos de mT a T
gamma = 28                                                      # Relación giromagnética del electrón (GHz / T)

# Frecuencias de resonancia debido al efecto Zeeman
f1 = D - gamma*B
f2 = D + gamma*B

# El barrido se adapta automáticamente al campo magnético
margen = max(0.02, 2*gamma*B) # GHz

f_min = min(f1, f2) - margen 
f_max = max(f1, f2) + margen 

frec  = np.linspace(f_min, f_max, N)  # Vector de frecuencias del barrido


# FUNCIÓN CURVA ODMR
def lorentziana(f, f0, ancho, contr):              
    """
    Calcula una resonancia con forma lorentziana.

    Parámetros:
        · f  : vector de frecuencias
        · f0 : frecuencia central de la resonancia

        · ancho : anchura de la resonancia
        · contr : profundidad de la resonancia
    """

    return contr*ancho/((f - f0)**2 + ancho**2)

# CONSTRUCCIÓN CURVA ODMR

# Parámetros de la resonancia
ancho = 0.003 
contr = 0.080

fluor0 = 1 - lorentziana(frec, D, ancho, contr)      # Fluorescencia de la frecuencia; mínimo !!!!

# -----------------------------------------------------------------------------------------------------------------

# RESULTADOS

# Cálculo del campo magnético a partir de las resonancias
separacion = f2 - f1
B_calculado = separacion / (2 * gamma)

print("\n==============================")
print("          RESULTADOS")
print("==============================")

print(f"Campo introducido: {B*1000:.3f} mT")
print(f"Resonancia 1: {f1:.6f} GHz")
print(f"Resonancia 2: {f2:.6f} GHz")
print(f"Separación: {separacion*1000:.2f} MHz")
print(f"Campo calculado: {B_calculado*1000:.3f} mT")

# -----------------------------------------------------------------------------------------------------------------

# REPRESENTACIÓN CURVA ODMR SIN CAMPO
plt.figure(figsize=(8,5))

plt.plot(frec, fluor0)

plt.xlabel('Frecuencia de microondas (GHz)')
plt.ylabel('Flourescencia normalizada')
plt.title('Simulación de una curva ODMR sin campo magnético')

plt.axvline(D, linestyle='--', color = 'red', label=f'D = {D:.3f} GHz')

plt.grid(True)
plt.tight_layout()

plt.savefig('odmr_sin_campo_magnético.png')

# -----------------------------------------------------------------------------------------------------------------

# CONSTRUCCIÓN DE LA CURVA ODMR

# Fluorescencia con dos resonancias
flour = (1 - lorentziana(frec, f1, ancho, contr) - lorentziana(frec, f2, ancho, contr))

# REPRESENTACIÓN CURVA ODMR CON CAMPO
plt.figure(figsize=(8,5))

plt.plot(frec, flour)

plt.xlabel('Frecuencia de microondas (GHz)')
plt.ylabel('Flourescencia normalizada')
plt.title(f'Curva ODMR con campo magnético B = {B*1000:.1f} mT')

plt.axvline(f1, linestyle='--', color = 'orange', label=f'f1 = {f1:.3f} GHz')
plt.axvline(f2, linestyle='--', color = 'yellow', label=f'f2 = {f2:.3f} GHz')

plt.grid(True)
plt.tight_layout()

plt.savefig('odmr_con_campo_magnético.png')
plt.show()








