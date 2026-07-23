'''
Simulación de curvas ODMR de un centro NV en diamante.

Este programa complementa el manual:
"Introducción a los Sensores Cuánticos Basados en Centros NV en Diamante".

El código simula:

1. Una curva ODMR sin campo magnético, con una resonancia centrada
en la división de campo cero D = 2,87 GHz.

2. Una curva ODMR en presencia de un campo magnético externo,
mostrando el desdoblamiento Zeeman de las resonancias.

3. El cálculo del campo magnético a partir de la separación entre
las dos frecuencias de resonancia.

Autora: Sofía Núñez de Andrés
Año: 2026
Versión: 1.0
License: Educational use
'''

import numpy as np                                             # librería 1: nos permite trabajar con vectores / cálculos numéricos # type: ignore
import matplotlib.pyplot as plt                                # librería 2: nos permite dibujar las gráficas                       # type: ignore
from datetime import datetime                                  # librería 3: nos permite poner la fecha y la hora real              # type: ignore

print("\n")
print("\n===========================================================")
print(" SIMULADOR ODMR - CENTROS NV EN DIAMANTE")
print("===========================================================")

# ============================================================
# PARÁMETROS FÍSICOS DEL CENTRO NV
D         = 2.87                                               # División de campo cero (GHz)
gamma     = 28                                                 # Relación giromagnética (GHz/T)
# ============================================================

# ============================================================
# PARÁMETROS DE LA SIMULACIÓN
N            = 1000                                            # Número de puntos del barrido
ancho        = 0.003                                           # Anchura de la resonancia (GHz)
contraste    = 0.080                                           # Contraste óptico

rng = np.random.default_rng()                                  # Generador de números aleatorios
# ============================================================

# ============================================================
# ENTRADA DE PARÁMETROS


entrada_B    = input("\n-> Campos magnéticos (mT) [3]            : ")
B            = float(entrada_B) if entrada_B else 3.0

entrada_P    = input("-> Potencia de microondas (dBm) [0]      : ")
potencia_dBm = float(entrada_P) if entrada_P else 0.0

entrada_F    = input("-> Tasa de fotones (fotones/s) [1000000] : ")
tasa_fotones = float(entrada_F) if entrada_F else 1_000_000

entrada_R    = input("-> Nivel de ruido (0 - 0.01) [0.0008]    : ")
nivel_ruido  = float(entrada_R) if entrada_R else 0.0008

entrada_N    = input("-> Número de puntos [1000]               : ")
N            = float(entrada_N) if entrada_N else 1000

if B < 0:
    raise ValueError("El campo magnético debe ser positivo.")

if tasa_fotones <= 0:
    raise ValueError("La tasa de fotones debe ser positiva.")

if nivel_ruido < 0:
    raise ValueError("El nivel de ruido no puede ser negativo.")

if N < 100:
    raise ValueError("El número de puntos debe ser al menos 100.")

B = B / 1000
# ============================================================

# ============================================================
# FRECUENCIAS DE RESONANCIA
f1 = D - gamma * B
f2 = D + gamma * B
# ============================================================

# ============================================================
# BARRIDOS DE FRECUENCIA
# Barrido para la curva sin campo magnético
margen_sin_campo = 0.03 # GHz
frecuencias_sin_campo = np.linspace(D - margen_sin_campo, D + margen_sin_campo, N)

# Barrido para la curva con campo magnético
margen_con_campo = max(0.03, gamma * B + 0.02)
frecuencias_con_campo = np.linspace(f1 - margen_con_campo, f2 + margen_con_campo, N)
# ============================================================

# ============================================================
# FUNCIÓN CURVA ODMR
def lorentziana(f, f0, ancho, contraste):
    '''
    Calcula una resonancia lorentziana normalizada.

    Parámetros
    ----------
    · f         : array
                  Frecuencias.
    · f0        : float
                  Frecuencia central.
    · ancho     : float
                  Anchura de la resonancia.
    · contraste : float
                  Profundidad máxima de la resonancia.
    '''

    return contraste * (ancho**2) / ((f - f0)**2 + ancho**2)
# ============================================================

# ============================================================
# FUNCIÓN DE SENSIBILIDAD RELATIVA
def calcular_sensibilidad_relativa(ancho_resonancia,contraste_optico,tasa_fotones_detectados):
    """
    Calcula un indicador relativo de sensibilidad magnética.

    La expresión utilizada es:

    sensibilidad ∝ anchura / (contraste * sqrt(tasa de fotones))

    Un valor menor indica una mejor sensibilidad.

    Esta magnitud es relativa y no se expresa en nT/sqrt(Hz),
    porque sería necesaria una constante de calibración
    experimental para obtener una sensibilidad absoluta.
    """

    if ancho_resonancia <= 0:
        raise ValueError("La anchura de la resonancia debe ser positiva.")

    if contraste_optico <= 0:
        raise ValueError("El contraste óptico debe ser positivo.")

    if tasa_fotones_detectados <= 0:
        raise ValueError("La tasa de fotones debe ser positiva.")

    sensibilidad_relativa = (ancho_resonancia / (contraste_optico* np.sqrt(tasa_fotones_detectados)))

    return sensibilidad_relativa
# ============================================================

# ============================================================
# EFECTO DE LA POTENCIA DE MICROONDAS
def ajustar_parametros_microondas(potencia_dBm,ancho_base,contraste_base):

    if potencia_dBm < -30:
        potencia_dBm = -30
    if potencia_dBm > 20:
        potencia_dBm = 20

    # Contraste
    if potencia_dBm <= 0:
        contraste = contraste_base * (0.4 + (0.6 * (potencia_dBm + 30) / 30))
    else:
        contraste = contraste_base

    # Anchura
    if potencia_dBm <= 0:
        ancho = ancho_base
    else:
        ancho = ancho_base * (1 + potencia_dBm / 20)

    return ancho, contraste

ancho, contraste = ajustar_parametros_microondas(potencia_dBm, ancho, contraste)

# ============================================================
# CONSTRUCCIÓN DE LAS CURVAS ODMR
# Curva ODMR sin campo magnético:
# una única resonancia centrada en D
fluorescencia_sin_campo = (1 - lorentziana(frecuencias_sin_campo, D, ancho, contraste))

# -------------------------------------------------------
# Ruido experimental (sin campo)
# -------------------------------------------------------
fluorescencia_sin_campo_ruido = (fluorescencia_sin_campo + rng.normal(0, nivel_ruido, len(fluorescencia_sin_campo)))

# Curva ODMR con campo magnético:
# dos resonancias debidas al efecto Zeeman
fluorescencia_con_campo = (1 - lorentziana(frecuencias_con_campo, f1, ancho, contraste) - lorentziana(frecuencias_con_campo, f2, ancho, contraste))

# -------------------------------------------------------
# Ruido experimental (con campo)
# -------------------------------------------------------
fluorescencia_con_campo_ruido = (fluorescencia_con_campo + rng.normal(0, nivel_ruido, len(fluorescencia_con_campo)))
# ============================================================

# ============================================================
# LOCALIZACIÓN AUTOMÁTICA DE LAS RESONANCIAS
# Índices correspondientes a los mínimos experimentales
indice_f1 = np.argmin(fluorescencia_con_campo_ruido[:len(frecuencias_con_campo)//2])

indice_f2 = (np.argmin(fluorescencia_con_campo_ruido[len(frecuencias_con_campo)//2:]) + len(frecuencias_con_campo)//2)

# Frecuencias obtenidas experimentalmente
f1_exp = frecuencias_con_campo[indice_f1]
f2_exp = frecuencias_con_campo[indice_f2]
# ============================================================

# ============================================================
# CÁLCULO DE LA SENSIBILIDAD RELATIVA
sensibilidad_relativa = calcular_sensibilidad_relativa(ancho, contraste, tasa_fotones)
# ============================================================

# ============================================================
# RESULTADOS
# Separación obtenida a partir de la curva simulada
separacion = f2_exp - f1_exp

# Campo recuperado a partir de la simulación
B_calculado = separacion/(2*gamma)

# Error de la medida
error = abs(B_calculado - B)
error_relativo = 100*error/B if B != 0 else 0

print("\n===========================================================")
print(" PARÁMETROS DE LA SIMULACIÓN")

print(f"\n · Campo magnético             : {B*1000:.3f} mT")
print(f" · Potencia de microondas      : {potencia_dBm:.1f} dBm")
print(f" · Tasa de fotones             : {tasa_fotones:,.0f} fotones/s")
print(f" · Nivel de ruido              : {nivel_ruido:.4f}")
print(f" · Número de puntos            : {N}")

print("\n===========================================================")

print("\n===========================================================")
print(" RESULTADOS CON UN VALOR FIJO DE B")

print(f"\n · Campo introducido : {B * 1000:.3f} mT")

print(f"\n · Resonancia 1 (experimental) : {f1_exp:.6f} GHz")
print(f" · Resonancia 2 (experimental) : {f2_exp:.6f} GHz")
print(f" · Separación                  : {separacion * 1000:.2f} MHz")

print(f"\n · Campo calculado : {B_calculado * 1000:.3f} mT")

print(f"\n · Error absoluto : {error*1000:.4f} mT")
print(f" · Error relativo : {error_relativo:.3f} %")

print(f"\n · Anchura de resonancia : {ancho * 1000:.3f} MHz")
print(f" · Contraste óptico      : {contraste * 100:.2f} %")
print(f" · Tasa de fotones       : {tasa_fotones:,.0f} fotones/s")

print("\n · Indicador de sensibilidad relativa : " f"{sensibilidad_relativa:.3e}")
# ============================================================

# ============================================================
# REPRESENTACIÓN DE LA CURVA ODMR SIN CAMPO MAGNÉTICO
plt.figure(figsize=(9, 5.5))

plt.plot(frecuencias_sin_campo,    fluorescencia_sin_campo,       color="black",  linewidth=2.2,     label="Curva teórica")
plt.scatter(frecuencias_sin_campo, fluorescencia_sin_campo_ruido, s=8, alpha=0.5, color="royalblue", label="Datos simulados")

plt.axvline(D, linestyle="--", color="#007BFF", linewidth=1.8, alpha=0.7, label=f"D = {D:.3f} GHz")

plt.xlabel("Frecuencia de microondas (GHz)")
plt.ylabel("Fluorescencia normalizada")
plt.title(f"Simulación ODMR sin campo magnético", fontsize=14, fontweight="bold", pad=18)
plt.suptitle(f"Potencia = {potencia_dBm:.1f} dBm | " f"Fotones = {tasa_fotones:.1e} s⁻¹ | " f"Ruido = {nivel_ruido:.4f} | N = {N}", fontsize=10, fontstyle="italic", y=0.98)

plt.legend(loc="best")
plt.grid(True)
plt.tight_layout()

plt.savefig("odmr_sin_campo_magnetico.png", dpi=300)
# ============================================================

# ============================================================
# REPRESENTACIÓN DE LA CURVA ODMR CON CAMPO MAGNÉTICO
plt.figure(figsize=(9, 5.5))

plt.plot(frecuencias_con_campo,    fluorescencia_con_campo,       color="black",  linewidth=2.2,     label="Curva teórica")
plt.scatter(frecuencias_con_campo, fluorescencia_con_campo_ruido, s=8, alpha=0.5, color="royalblue", label="Datos simulados")

plt.axvline(f1, linestyle="--", color="#FF1493", linewidth=1.8, alpha=0.7, label=f"f₁ = {f1:.3f} GHz")
plt.axvline(f2, linestyle="--", color="#4B0082", linewidth=1.8, alpha=0.7, label=f"f₂ = {f2:.3f} GHz")

plt.xlabel("Frecuencia de microondas (GHz)")
plt.ylabel("Fluorescencia normalizada")
plt.title(f"Simulación ODMR con campo magnético", fontsize=14, fontweight="bold", pad=18)
plt.suptitle(f"B = {B*1000:.3f} mT | Potencia = {potencia_dBm:.1f} dBm | " f"Fotones = {tasa_fotones:.1e} s⁻¹ | " f"Ruido = {nivel_ruido:.4f} | N = {N}", fontsize=10, fontstyle="italic", y=0.98)

plt.legend(loc="best")
plt.grid(True)
plt.tight_layout()

plt.savefig("odmr_con_campo_magnetico.png", dpi=300)
# ============================================================

# ============================================================
# GUARDADO DE RESULTADOS EN UN ARCHIVO TXT

with open("resultados_ODMR.txt", "w", encoding="utf-8") as archivo:
    archivo.write("="*60 + "\n")
    archivo.write(" "*10 + "SIMULADOR ODMR - CENTROS NV EN DIAMANTE\n")
    archivo.write("="*60 + "\n\n")

    archivo.write("-> PARÁMETROS DE LA SIMULACIÓN\n")
    archivo.write("-"*60 + "\n")

    archivo.write(f"• Campo magnético        : {B*1000:.3f} mT\n")
    archivo.write(f"• Potencia de microondas : {potencia_dBm:.1f} dBm\n")
    archivo.write(f"• Tasa de fotones        : {tasa_fotones:,.0f} fotones/s\n")
    archivo.write(f"• Nivel de ruido         : {nivel_ruido:.4f}\n")
    archivo.write(f"• Número de puntos       : {N}\n")

    archivo.write("\n")

    archivo.write("-> PARÁMETROS ODMR\n")
    archivo.write("-"*60 + "\n")

    archivo.write(f"• Anchura                : {ancho*1000:.2f} MHz\n")
    archivo.write(f"• Contraste              : {contraste*100:.2f} %\n")

    archivo.write("\n")

    archivo.write("-> RESULTADOS\n")
    archivo.write("-"*60 + "\n")

    archivo.write(f"• Resonancia 1           : {f1_exp:.6f} GHz\n")
    archivo.write(f"• Resonancia 2           : {f2_exp:.6f} GHz\n")
    archivo.write(f"• Separación             : {separacion*1000:.2f} MHz\n")

    archivo.write("\n")

    archivo.write(f"• Campo calculado        : {B_calculado*1000:.3f} mT\n")
    archivo.write(f"• Error absoluto         : {error*1000:.4f} mT\n")
    archivo.write(f"• Error relativo         : {error_relativo:.3f} %\n")
    archivo.write(f"• Sensibilidad           : {sensibilidad_relativa:.3e} T/\u221AHz\n")

    archivo.write("\n")
    archivo.write("="*60 + "\n")
    archivo.write(" Generated by ODMR Simulator v1.0\n")
    archivo.write("="*60 + "\n")
    archivo.write(f" Fecha: {datetime.now().strftime('%d/%m/%Y %H:%M:%S')}\n")
    archivo.write(" Developed during the CINN research internship (2026)\n")
    archivo.write("-"*60 + "\n")

print("\n-> Los resultados se han guardado en 'resultados_ODMR.txt'")
print("===========================================================")
# ============================================================

# ============================================================
# COMPARACIÓN ESPECTROS PARA DISTINTOS VALORES DE B
# Cálculos + Representación Gráfica:
print("\n===========================================================")
print(" RESULTADOS CON DISTINTOS VALORES DE B")

entrada_C = input("\n-> Campos magnéticos (mT) [0,2.5,5,7.5,10] : ")      # El usuario introduce los campos

print("===========================================================")

if entrada_C:
    campos = [float(x) for x in entrada_C.split(",")]
else:
    campos = [0,2.5,5,7.5,10]

campo_maximo = max(abs(B_mT) for B_mT in campos) / 1000
margen       = 0.03     # GHz

frecuencia_min = D - gamma*campo_maximo - margen
frecuencia_max = D + gamma*campo_maximo + margen
frecuencias_comparacion = np.linspace(frecuencia_min, frecuencia_max, N)

plt.figure(figsize=(9,5.5))

for B_mT in campos:
    B_comp = B_mT / 1000      # pasar a teslas

    # Cálculo f1 y f2
    f1_comp = D - gamma * B_comp
    f2_comp = D + gamma * B_comp

    # Cálculo fluorescencia
    fluorescencia_comp = (1 - lorentziana(frecuencias_comparacion, f1_comp, ancho, contraste) - lorentziana(frecuencias_comparacion, f2_comp, ancho, contraste))

    # Dibujar la curva
    plt.plot(frecuencias_comparacion, fluorescencia_comp, linewidth=2.2, label=f"{B_mT:g} mT")

plt.xlabel("Frecuencia de microondas (GHz)")
plt.ylabel("Fluorescencia normalizada")

plt.title("Comparación de espectros ODMR para distintos campos magnéticos", fontsize=14, fontweight="bold", pad=18)
plt.suptitle("Espectros simulados para distintos campos magnéticos", fontsize=10, fontstyle="italic", y=0.98)

plt.legend(title="Campo magnético")
plt.grid(True)
plt.tight_layout()

plt.savefig("comparacion_campos_magneticos.png", dpi=300)
plt.show()
# ============================================================
# Fin del programa.









