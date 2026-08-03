"""
Simulación de espectros ODMR de un centro NV- en diamante.

Este programa complementa el manual:
"Introducción a los Sensores Cuánticos Basados en Centros NV en Diamante".

El programa permite:

1. Representar una curva ODMR sin campo magnético.
2. Representar una curva ODMR con un campo magnético externo.
3. Estudiar la variación de las frecuencias de resonancia con el campo.
4. Estudiar la variación de las resonancias con el ángulo entre el campo
magnético y el eje del centro NV.
5. Mostrar la evolución del espectro ODMR al aumentar el campo.
6. Guardar las gráficas y los resultados numéricos de la simulación.

Autora: Sofía Núñez de Andrés
Prácticas extracurriculares en el CINN, 2026
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

print("\n==============================================================")
print(" SIMULADOR ODMR - CENTROS NV EN DIAMANTE")
print("==============================================================")

# ============================================================
# PARÁMETROS FÍSICOS DEL CENTRO NV
D_GHz = 2.87 # División de campo cero (GHz)
gamma_e_GHz_T = 28.0 # Relación giromagnética del electrón (GHz/T)

# PARÁMETROS DE LA SIMULACIÓN
numero_puntos = 1000 # Número de puntos del barrido
ancho_fwhm_GHz = 0.003 # Anchura total a media altura (GHz)
contraste_por_resonancia = 0.080 # Contraste de cada resonancia ODMR
margen_frecuencia_GHz = 0.03 # Margen del barrido de frecuencia (GHz)

# PARÁMETROS DE DETECCIÓN
tasa_fotones_Hz = 250_000 # Fotones detectados por segundo
tiempo_integracion_s = 0.01 # Tiempo de integración por punto (s)
ruido_tecnico_relativo = 0.002 # Desviación típica relativa de ruido gaussiano

# PARÁMETROS PARA LA SENSIBILIDAD
factor_perfil_lorentziano = 0.77 # Perfil lorentziano según la convención utilizada

# SEMILLA DEL GENERADOR ALEATORIO
semilla_aleatoria = 7 # Permite reproducir la misma simulación
generador_aleatorio = np.random.default_rng(semilla_aleatoria)

carpeta_script = Path(__file__).resolve().parent
carpeta_resultados = carpeta_script / "results_odmr"
carpeta_resultados.mkdir(parents=True, exist_ok=True)
# ============================================================

# ============================================================
# MATRICES DE ESPÍN PARA S=1. Base: |+1>,|0>,|-1>
Sx = (1 / np.sqrt(2)) * np.array(
    [
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0]
    ], 
    dtype=complex
)
Sy = (1 / np.sqrt(2)) * np.array(
    [
        [0, -1j, 0],
        [1j, 0, -1j],
        [0, 1j, 0]
    ], 
    dtype=complex
)
Sz = np.array(
    [
        [1, 0, 0],
        [0, 0, 0],
        [0, 0, -1]
    ], 
    dtype=complex
)
# ============================================================

# ============================================================
# FUNCIONES
def calcular_resonancias(campo_x_T, campo_y_T, campo_z_T, perturbacion_E_GHz=0.0):
    """
    Construye el Hamiltoniano del estado fundamental del centro NV
    y calcula sus frecuencias de transición.

    Parameters
    ----------
    campo_x_T : float
    Componente x del campo magnético, en teslas.

    campo_y_T : float
    Componente y del campo magnético, en teslas.

    campo_z_T : float
    Componente z del campo magnético, en teslas.

    perturbacion_E_GHz : float, optional
    Perturbación transversal E, en GHz. Por defecto es cero.

    Returns
    -------
    frecuencia_inferior_GHz : float
    Menor frecuencia de transición desde el nivel fundamental,
    expresada en GHz.

    frecuencia_superior_GHz : float
    Mayor frecuencia de transición desde el nivel fundamental,
    expresada en GHz.

    energias_GHz : numpy.ndarray
    Autovalores ordenados del Hamiltoniano H/h, en GHz.
    """

    hamiltoniano_zfs = D_GHz * (Sz @ Sz)

    # Se omite el término -(2/3) D I porque solo produce un
    # desplazamiento común de todos los niveles de energía.

    hamiltoniano_transversal = perturbacion_E_GHz * ((Sx @ Sx) - (Sy @ Sy))

    hamiltoniano_zeeman = gamma_e_GHz_T * (campo_x_T * Sx + campo_y_T * Sy + campo_z_T * Sz)

    hamiltoniano_total = (hamiltoniano_zfs + hamiltoniano_transversal + hamiltoniano_zeeman)

    energias_GHz = np.linalg.eigvalsh(hamiltoniano_total)
    energias_GHz = np.sort(np.real(energias_GHz))

    frecuencia_inferior_GHz = energias_GHz[1] - energias_GHz[0]
    frecuencia_superior_GHz = energias_GHz[2] - energias_GHz[0]

    frecuencias_GHz = np.sort([frecuencia_inferior_GHz, frecuencia_superior_GHz])

    return frecuencias_GHz[0], frecuencias_GHz[1], energias_GHz


def calcular_perfil_lorentziano(frecuencias_GHz, frecuencia_central_GHz, ancho_fwhm_GHz=ancho_fwhm_GHz):
    """
    Calcula un perfil lorentziano normalizado.

    Parameters
    ----------
    frecuencias_GHz : numpy.ndarray
    Valores del barrido de frecuencia, en GHz.

    frecuencia_central_GHz : float
    Frecuencia central de la resonancia, en GHz.

    ancho_fwhm_GHz : float, optional
    Anchura total a media altura, FWHM, en GHz.

    Returns
    -------
    numpy.ndarray
    Perfil lorentziano con valor máximo igual a 1.
    """

    semiancho_GHz = ancho_fwhm_GHz / 2

    diferencia_frecuencia_GHz = (frecuencias_GHz - frecuencia_central_GHz)

    perfil = semiancho_GHz**2 / (diferencia_frecuencia_GHz**2 + semiancho_GHz**2)

    return perfil

def obtener_resonancias_unicas(resonancias_GHz, tolerancia_GHz):
    """
    Elimina resonancias coincidentes dentro de una tolerancia.

    Parameters
    ----------
    resonancias_GHz : iterable
    Frecuencias de resonancia, en GHz.

    tolerancia_GHz : float
    Diferencia máxima para considerar dos resonancias degeneradas.

    Returns
    -------
    list
    Lista ordenada de resonancias no duplicadas.
    """
    resonancias_ordenadas_GHz = sorted(resonancias_GHz)
    resonancias_unicas_GHz = []

    for frecuencia_GHz in resonancias_ordenadas_GHz:
        es_duplicada = any(np.isclose(frecuencia_GHz, frecuencia_anterior_GHz, atol=tolerancia_GHz, rtol=0.0) for frecuencia_anterior_GHz in resonancias_unicas_GHz)

        if not es_duplicada:
            resonancias_unicas_GHz.append(frecuencia_GHz)

    return resonancias_unicas_GHz

def construir_curva_odmr(frecuencias_GHz, resonancias_GHz, contraste=contraste_por_resonancia, ancho_fwhm_GHz=ancho_fwhm_GHz):
    """
    Construye una curva ODMR normalizada sin ruido.

    Parameters
    ----------
    frecuencias_GHz : numpy.ndarray
    Valores del barrido de frecuencia, en GHz.

    resonancias_GHz : iterable
    Frecuencias centrales de las transiciones ODMR, en GHz.

    contraste : float, optional
    Profundidad relativa de cada resonancia.

    ancho_fwhm_GHz : float, optional
    Anchura total a media altura de cada resonancia, en GHz.

    Returns
    -------
    numpy.ndarray
    Fluorescencia ODMR normalizada.
    """
    fluorescencia_normalizada = np.ones_like(frecuencias_GHz, dtype=float)

    tolerancia_degeneracion_GHz = ancho_fwhm_GHz / 100

    resonancias_unicas_GHz = obtener_resonancias_unicas(resonancias_GHz,tolerancia_degeneracion_GHz)

    for frecuencia_central_GHz in resonancias_unicas_GHz:
        perfil_lorentziano = calcular_perfil_lorentziano(frecuencias_GHz, frecuencia_central_GHz, ancho_fwhm_GHz)

        fluorescencia_normalizada -= (contraste * perfil_lorentziano)

    return fluorescencia_normalizada

def guardar_grafica(nombre_archivo, resolucion_dpi=300):
    """
    Guarda la figura activa dentro de la carpeta de resultados.

    Parameters
    ----------
    nombre_archivo : str
    Nombre del archivo de salida.

    resolucion_dpi : int, optional
    Resolución de la imagen en puntos por pulgada.
    Por defecto es 300 dpi.

    Returns
    -------
    pathlib.Path
    Ruta completa del archivo guardado.
    """

    nombre_archivo = str(nombre_archivo).strip()

    if not nombre_archivo:
        raise ValueError("El nombre del archivo de la gráfica no puede estar vacío.")

    ruta_archivo = Path(nombre_archivo)

    if ruta_archivo.suffix == "":
        ruta_archivo = ruta_archivo.with_suffix(".png")

    ruta_completa = carpeta_resultados / ruta_archivo.name

    plt.tight_layout()

    plt.savefig(ruta_completa, dpi=resolucion_dpi, bbox_inches="tight")

    plt.close()

    return ruta_completa

def calcular_componentes_campo(modulo_campo_T, angulo_rad):
    """
    Descompone el campo magnético en componentes cartesianas.

    Se supone que el campo está contenido en el plano XZ y que
    el eje Z coincide con el eje del centro NV.

    Parameters
    ----------
    modulo_campo_T : float
    Módulo del campo magnético, en teslas.

    angulo_rad : float
    Ángulo entre el campo magnético y el eje NV, en radianes.

    Returns
    -------
    campo_x_T : float
    Componente transversal del campo en el eje X.

    campo_y_T : float
    Componente del campo en el eje Y. En este modelo es cero.

    campo_z_T : float
    Componente longitudinal del campo, paralela al eje NV.
    """

    campo_x_T = modulo_campo_T * np.sin(angulo_rad)
    campo_y_T = 0.0
    campo_z_T = modulo_campo_T * np.cos(angulo_rad)

    return campo_x_T, campo_y_T, campo_z_T
# ============================================================

# ============================================================
# DATOS INTRODUCIDOS POR EL USUARIO
def pedir_float(mensaje, valor_por_defecto, valor_minimo=None,valor_maximo=None):
    """
    Solicita un número decimal al usuario y valida su intervalo.

    Parameters
    ----------
    mensaje : str
    Texto mostrado al solicitar el valor.

    valor_por_defecto : float
    Valor utilizado cuando el usuario pulsa Intro sin escribir nada.

    valor_minimo : float or None, optional
    Valor mínimo permitido. Si es None, no se establece un mínimo.

    valor_maximo : float or None, optional
    Valor máximo permitido. Si es None, no se establece un máximo.

    Returns
    -------
    float
    Valor numérico validado.
    """
    while True:
        entrada = input(mensaje).strip()

        if entrada == "":
            valor = float(valor_por_defecto)
        else:
            entrada = entrada.replace(",", ".")

            try:
                valor = float(entrada)
            except ValueError:
                print("Entrada no válida. Introduce un valor numérico.")
                continue

        if not np.isfinite(valor):
            print("El valor debe ser un número finito.")
            continue

        if valor_minimo is not None and valor < valor_minimo:
            print(f"El valor debe ser mayor o igual que {valor_minimo}.")
            continue

        if valor_maximo is not None and valor > valor_maximo:
            print(f"El valor debe ser menor o igual que {valor_maximo}.")
            continue

        return valor
# ============================================================

# ============================================================
# COMPONENTES DE LOS CAMPOS
print("\nIntroduce los parámetros de la simulación. ")
print("Pulsa Intro para utilizar el valor indicado entre corchetes. ")

campo_mT    = pedir_float("\n -> Campo magnético (mT) [3]: ", valor_por_defecto=3.0, valor_minimo=0.0)
campo_T     = campo_mT / 1000

angulo_grados = pedir_float("\n -> Ángulo respecto al eje NV en grados [0]: ", valor_por_defecto=0.0, valor_minimo=0.0, valor_maximo=90.0)
angulo_rad    = np.radians(angulo_grados)

perturbacion_E_MHz = pedir_float("\n -> Perturbación transversal E (MHz) [0]: ", valor_por_defecto=0.0, valor_minimo=0.0)
perturbacion_E_GHz = perturbacion_E_MHz / 1000

# Se supone que el campo está contenido en el plano XZ
campo_x_T, campo_y_T, campo_z_T = calcular_componentes_campo(campo_T, angulo_rad)
# ============================================================

# ============================================================
# CÁLCULO DE LAS RESONANCIAS
frecuencia_inferior_sin_campo_GHz, frecuencia_superior_sin_campo_GHz, energias_sin_campo_GHz = calcular_resonancias(campo_x_T=0.0, campo_y_T=0.0, campo_z_T=0.0, perturbacion_E_GHz=perturbacion_E_GHz)
frecuencia_inferior_GHz, frecuencia_superior_GHz, energias_con_campo_GHz = calcular_resonancias(campo_x_T=campo_x_T, campo_y_T=campo_y_T, campo_z_T=campo_z_T, perturbacion_E_GHz=perturbacion_E_GHz)

resonancias_sin_campo_GHz = [frecuencia_inferior_sin_campo_GHz, frecuencia_superior_sin_campo_GHz]
resonancias_con_campo_GHz = [frecuencia_inferior_GHz, frecuencia_superior_GHz]

todas_las_resonancias_GHz = (resonancias_sin_campo_GHz + resonancias_con_campo_GHz)

frecuencia_minima_GHz = min(todas_las_resonancias_GHz) - margen_frecuencia_GHz
frecuencia_maxima_GHz = max(todas_las_resonancias_GHz) + margen_frecuencia_GHz
if frecuencia_minima_GHz >= frecuencia_maxima_GHz:
    raise ValueError("El intervalo de frecuencias calculado no es válido. ")

frecuencias_GHz = np.linspace(frecuencia_minima_GHz, frecuencia_maxima_GHz, numero_puntos)

curva_sin_campo = construir_curva_odmr(frecuencias_GHz=frecuencias_GHz, resonancias_GHz=resonancias_sin_campo_GHz)
curva_con_campo = construir_curva_odmr(frecuencias_GHz=frecuencias_GHz, resonancias_GHz=resonancias_con_campo_GHz)
# ============================================================

# ============================================================
# GRÁFICA 1: CURVA ODMR SIN CAMPO MAGNÉTICO
figura_1, eje = plt.subplots(figsize=(8, 5))

eje.plot(frecuencias_GHz, curva_sin_campo, color="tab:blue", linewidth=2, label="Señal ODMR")

eje.axvline(D_GHz, linestyle="--", linewidth=1.5, color="tab:red", alpha=0.8, label=f"ZFS = {D_GHz:.2f} GHz")

eje.set_xlabel("Frecuencia de microondas (GHz)")
eje.set_ylabel("Fluorescencia normalizada")

eje.set_title("Espectro ODMR sin campo magnético")

eje.set_xlim(frecuencia_minima_GHz, frecuencia_maxima_GHz)
eje.set_ylim(0.88, 1.01)

eje.grid(True, alpha=0.3)
eje.legend(loc="best")

guardar_grafica("odmr_without_magnetic_field.png")
# ============================================================

# ============================================================
# GRÁFICA 2: CURVA ODMR CON CAMPO MAGNÉTICO
figura_2, eje = plt.subplots(figsize=(8, 5))

eje.plot(frecuencias_GHz, curva_con_campo, color="tab:green", linewidth=2, label="Señal ODMR")

resonancias_visibles_GHz = obtener_resonancias_unicas(resonancias_con_campo_GHz, tolerancia_GHz=ancho_fwhm_GHz / 100)

if len(resonancias_visibles_GHz) == 1:
    eje.axvline(resonancias_visibles_GHz[0], linestyle="--", linewidth=1.5, color="tab:red", alpha=0.8, label=f"Resonancia degenerada = {resonancias_visibles_GHz[0]:.6f} GHz")

else:
    eje.axvline(frecuencia_inferior_GHz, linestyle="--", linewidth=1.5, color="tab:red", alpha=0.8, label=f"Resonancia inferior = {frecuencia_inferior_GHz:.6f} GHz")
    eje.axvline(frecuencia_superior_GHz, linestyle="--", linewidth=1.5, color="tab:orange", alpha=0.8, label=f"Resonancia superior = {frecuencia_superior_GHz:.6f} GHz")

eje.set_xlabel("Frecuencia de microondas (GHz)")
eje.set_ylabel("Fluorescencia normalizada")

eje.set_title("Espectro ODMR con campo magnético\n" f"B = {campo_mT:.3f} mT, θ = {angulo_grados:.1f}° y E = {perturbacion_E_MHz:.3f} MHz")

eje.set_xlim(frecuencia_minima_GHz, frecuencia_maxima_GHz)
eje.set_ylim(0.88, 1.01)

eje.grid(True, alpha=0.3)
eje.legend(loc="best")

guardar_grafica("odmr_with_magnetic_field.png")
# ============================================================

# ============================================================
# GRÁFICA 3: FRECUENCIAS DE RESONANCIA FRENTE AL CAMPO
campo_maximo_mT = max(5.0, campo_mT * 2)
campos_mT = np.linspace(0.0, campo_maximo_mT, 201)

frecuencias_inferiores_campo_GHz = []
frecuencias_superiores_campo_GHz = []

for campo_actual_mT in campos_mT:
    campo_actual_T = campo_actual_mT / 1000

    campo_x_actual_T, campo_y_actual_T, campo_z_actual_T = calcular_componentes_campo(campo_actual_T,angulo_rad)

    frecuencia_inferior_actual_GHz, frecuencia_superior_actual_GHz, _ = calcular_resonancias(campo_x_T=campo_x_actual_T,campo_y_T=campo_y_actual_T,campo_z_T=campo_z_actual_T,perturbacion_E_GHz=perturbacion_E_GHz)

    frecuencias_inferiores_campo_GHz.append(frecuencia_inferior_actual_GHz)
    frecuencias_superiores_campo_GHz.append(frecuencia_superior_actual_GHz)

frecuencias_inferiores_campo_GHz = np.array(frecuencias_inferiores_campo_GHz)
frecuencias_superiores_campo_GHz = np.array(frecuencias_superiores_campo_GHz)

figura_3, eje = plt.subplots(figsize=(8, 5))

eje.plot(campos_mT,frecuencias_inferiores_campo_GHz,linewidth=2,color="tab:blue",label="Resonancia inferior")
eje.plot(campos_mT,frecuencias_superiores_campo_GHz,linewidth=2,linestyle="--",color="tab:orange",label="Resonancia superior")

eje.axvline(campo_mT, linestyle="--",linewidth=1.5,color="tab:red",alpha=0.8,label=f"Campo elegido = {campo_mT:.3f} mT")

eje.set_xlabel("Módulo del campo magnético (mT)")
eje.set_ylabel("Frecuencia de resonancia (GHz)")
eje.set_title("Dependencia de las resonancias con el campo magnético\n" f"θ = {angulo_grados:.1f}° y E = {perturbacion_E_MHz:.3f} MHz")

margen = 0.025
eje.set_xlim(0.0, campo_maximo_mT)
eje.set_ylim(min(frecuencias_inferiores_campo_GHz) - margen, max(frecuencias_superiores_campo_GHz) + margen)

eje.grid(True, alpha=0.3)
eje.legend(loc="best")

guardar_grafica("frequencies_vs_field.png")
# ============================================================

# ============================================================
# GRÁFICA 4: FRECUENCIAS DE RESONANCIA FRENTE AL ÁNGULO
angulos_grados = np.linspace(0.0, 90.0, 181)

frecuencias_inferiores_angulo_GHz = []
frecuencias_superiores_angulo_GHz = []

for angulo_actual_grados in angulos_grados:
    angulo_actual_rad = np.radians(angulo_actual_grados)

    campo_x_actual_T, campo_y_actual_T, campo_z_actual_T = calcular_componentes_campo(campo_T,angulo_actual_rad)

    frecuencia_inferior_actual_GHz, frecuencia_superior_actual_GHz, _ = calcular_resonancias(campo_x_T=campo_x_actual_T,campo_y_T=campo_y_actual_T,campo_z_T=campo_z_actual_T,perturbacion_E_GHz=perturbacion_E_GHz)

    frecuencias_inferiores_angulo_GHz.append(frecuencia_inferior_actual_GHz)
    frecuencias_superiores_angulo_GHz.append(frecuencia_superior_actual_GHz)

frecuencias_inferiores_angulo_GHz = np.array(frecuencias_inferiores_angulo_GHz)
frecuencias_superiores_angulo_GHz = np.array(frecuencias_superiores_angulo_GHz)

figura_4, eje = plt.subplots(figsize=(8, 5))

eje.plot(angulos_grados,frecuencias_inferiores_angulo_GHz,linewidth=2,color="tab:blue",label="Resonancia inferior")
eje.plot(angulos_grados,frecuencias_superiores_angulo_GHz,linewidth=2,color="tab:orange",label="Resonancia superior")

eje.axvline(angulo_grados,linestyle="--",linewidth=1.5,color="tab:red",alpha=0.8,label=f"Ángulo elegido = {angulo_grados:.1f}°")

eje.set_xlabel("Ángulo entre el campo magnético y el eje NV (grados)")
eje.set_ylabel("Frecuencia de resonancia (GHz)")
eje.set_title("Dependencia angular de las resonancias\n" f"B = {campo_mT:.3f} mT y E = {perturbacion_E_MHz:.3f} MHz")

eje.set_xlim(-2, 92)

eje.grid(True, alpha=0.3)
eje.legend(loc="best")

guardar_grafica("frequencies_vs_angle.png")
# ============================================================

# ============================================================
# GRÁFICA 5: EVOLUCIÓN DEL ESPECTRO ODMR CON EL CAMPO
campos_mostrados_mT = np.linspace(0.0, campo_maximo_mT, 6)
resonancias_evolucion_GHz = []

for campo_actual_mT in campos_mostrados_mT:
    campo_actual_T = campo_actual_mT / 1000

    campo_x_actual_T, campo_y_actual_T, campo_z_actual_T = calcular_componentes_campo(campo_actual_T,angulo_rad)

    frecuencia_inferior_actual_GHz, frecuencia_superior_actual_GHz, _ = calcular_resonancias(campo_x_T=campo_x_actual_T,campo_y_T=campo_y_actual_T,campo_z_T=campo_z_actual_T,perturbacion_E_GHz=perturbacion_E_GHz)

    resonancias_evolucion_GHz.append((campo_actual_mT,frecuencia_inferior_actual_GHz,frecuencia_superior_actual_GHz))

todas_las_resonancias_evolucion_GHz = []

for campo_actual_mT, frecuencia_inferior_actual_GHz, frecuencia_superior_actual_GHz in resonancias_evolucion_GHz:
    todas_las_resonancias_evolucion_GHz.extend([frecuencia_inferior_actual_GHz,frecuencia_superior_actual_GHz])

frecuencia_minima_evolucion_GHz = min(todas_las_resonancias_evolucion_GHz) - margen_frecuencia_GHz
frecuencia_maxima_evolucion_GHz = max(todas_las_resonancias_evolucion_GHz) + margen_frecuencia_GHz

frecuencias_evolucion_GHz = np.linspace(frecuencia_minima_evolucion_GHz,frecuencia_maxima_evolucion_GHz,numero_puntos)

figura_5, eje = plt.subplots(figsize=(8, 5))

for campo_actual_mT, frecuencia_inferior_actual_GHz, frecuencia_superior_actual_GHz in resonancias_evolucion_GHz:
    curva_actual = construir_curva_odmr(frecuencias_GHz=frecuencias_evolucion_GHz,resonancias_GHz=[frecuencia_inferior_actual_GHz,frecuencia_superior_actual_GHz,])

    eje.plot(frecuencias_evolucion_GHz,curva_actual,linewidth=2.2,label=f"B = {campo_actual_mT:.2f} mT")

eje.set_xlabel("Frecuencia de microondas (GHz)")
eje.set_ylabel("Fluorescencia normalizada")
eje.set_title("Evolución del espectro ODMR al aumentar el campo\n" f"θ = {angulo_grados:.1f}° y E = {perturbacion_E_MHz:.3f} MHz")

eje.set_xlim(frecuencia_minima_evolucion_GHz, frecuencia_maxima_evolucion_GHz)
eje.set_ylim(0.88, 1.01)

eje.grid(True, alpha=0.3)
eje.legend(loc="best")

guardar_grafica("evolution_odmr_with_magnetic_field.png")
# ============================================================

# ============================================================
# ESTIMACIÓN DEL CAMPO A PARTIR DEL DESDOBLAMIENTO
separacion_resonancias_GHz = abs(frecuencia_superior_GHz - frecuencia_inferior_GHz)

campo_longitudinal_estimado_T = separacion_resonancias_GHz / (2 * gamma_e_GHz_T)
campo_longitudinal_estimado_mT = campo_longitudinal_estimado_T * 1000

campo_longitudinal_real_mT = abs(campo_z_T) * 1000
error_estimacion_mT = campo_longitudinal_estimado_mT - campo_longitudinal_real_mT

if campo_longitudinal_real_mT > 0:
    error_estimacion_porcentual = 100 * error_estimacion_mT / campo_longitudinal_real_mT
else:
    error_estimacion_porcentual = 0.0
# ============================================================

# ============================================================
# ARCHIVO DE RESULTADOS
ruta_txt = carpeta_resultados / "results_odmr.txt"

linea = "=" * 72
sublinea = "-" * 72

with open(ruta_txt, "w", encoding="utf-8") as archivo:
    archivo.write(linea + "\n")
    archivo.write(" SIMULACIÓN ODMR DE UN CENTRO NV EN DIAMANTE\n")
    archivo.write(" Sofía Núñez de Andrés - CINN, 2026\n")
    archivo.write(linea + "\n\n")

    archivo.write("1. PARÁMETROS FÍSICOS DEL CENTRO NV\n")
    archivo.write(sublinea + "\n")
    archivo.write(f"División de campo cero, D     = {D_GHz:>12.6f} GHz\n")
    archivo.write(f"Relación giromagnética, gamma = {gamma_e_GHz_T:>12.6f} GHz/T\n")
    archivo.write(f"Perturbación transversal, E   = {perturbacion_E_GHz:>12.6f} GHz\n")
    archivo.write(f"Perturbación transversal, E   = {perturbacion_E_MHz:>12.6f} MHz\n\n")

    archivo.write("2. PARÁMETROS DEL ESPECTRO ODMR\n")
    archivo.write(sublinea + "\n")
    archivo.write(f"Número de puntos del barrido    = {numero_puntos:>12d}\n")
    archivo.write(f"Anchura FWHM de las resonancias = {ancho_fwhm_GHz:>12.6f} GHz\n")
    archivo.write(f"Contraste por resonancia        = {contraste_por_resonancia:>12.6f}\n")
    archivo.write(f"Margen del barrido              = {margen_frecuencia_GHz:>12.6f} GHz\n")
    archivo.write(f"Frecuencia mínima del barrido   = {frecuencia_minima_GHz:>12.9f} GHz\n")
    archivo.write(f"Frecuencia máxima del barrido   = {frecuencia_maxima_GHz:>12.9f} GHz\n\n")

    archivo.write("3. CAMPO MAGNÉTICO INTRODUCIDO\n")
    archivo.write(sublinea + "\n")
    archivo.write(f"Módulo del campo              = {campo_mT:>12.6f} mT\n")
    archivo.write(f"Ángulo respecto al eje NV     = {angulo_grados:>12.6f} grados\n")
    archivo.write(f"Componente Bx                 = {campo_x_T:>12.9f} T\n")
    archivo.write(f"Componente By                 = {campo_y_T:>12.9f} T\n")
    archivo.write(f"Componente Bz                 = {campo_z_T:>12.9f} T\n")
    archivo.write(f"Componente longitudinal, |Bz| = {campo_longitudinal_real_mT:>12.9f} mT\n\n")

    archivo.write("4. RESULTADOS SIN CAMPO MAGNÉTICO\n")
    archivo.write(sublinea + "\n")
    archivo.write("Autovalores de H/h = " + ", ".join(f"{energia:.9f}" for energia in energias_sin_campo_GHz) + " GHz\n")
    archivo.write(f"Frecuencia de resonancia inferior = {frecuencia_inferior_sin_campo_GHz:>12.9f} GHz\n")
    archivo.write(f"Frecuencia de resonancia superior = {frecuencia_superior_sin_campo_GHz:>12.9f} GHz\n\n")

    archivo.write("5. RESULTADOS CON CAMPO MAGNÉTICO\n")
    archivo.write(sublinea + "\n")
    archivo.write("Autovalores de H/h = " + ", ".join(f"{energia:.9f}" for energia in energias_con_campo_GHz) + " GHz\n")
    archivo.write(f"Frecuencia de resonancia inferior = {frecuencia_inferior_GHz:>12.9f} GHz\n")
    archivo.write(f"Frecuencia de resonancia superior = {frecuencia_superior_GHz:>12.9f} GHz\n")
    archivo.write(f"Separación entre resonancias      = {separacion_resonancias_GHz:>12.9f} GHz\n\n")

    archivo.write("6. ESTIMACIÓN DEL CAMPO MEDIANTE EL DESDOBLAMIENTO\n")
    archivo.write(sublinea + "\n")
    archivo.write(f"Componente longitudinal real, |Bz| = {campo_longitudinal_real_mT:>12.9f} mT\n")
    archivo.write(f"Componente longitudinal estimada   = {campo_longitudinal_estimado_mT:>12.9f} mT\n")
    archivo.write(f"Error de la estimación             = {error_estimacion_mT:>12.9f} mT\n")
    archivo.write(f"Error porcentual                   = {error_estimacion_porcentual:>12.6f} %\n")
    archivo.write("Nota: la expresión utilizada es una aproximación lineal y estima\n")
    archivo.write("principalmente la componente del campo paralela al eje NV.\n\n")

    archivo.write("7. ARCHIVOS GENERADOS\n")
    archivo.write(sublinea + "\n")
    archivo.write("1. odmr_without_magnetic_field.png\n")
    archivo.write("2. odmr_with_magnetic_field.png\n")
    archivo.write("3. frequencies_vs_field.png\n")
    archivo.write("4. frequencies_vs_angle.png\n")
    archivo.write("5. evolution_odmr_with_magnetic_field.png\n")
    archivo.write("6. results_odmr.txt\n\n")

    archivo.write(linea + "\n")
    archivo.write("Fin de la simulación\n")
    archivo.write(linea + "\n")
# ============================================================

# ============================================================
# RESULTADOS MOSTRADOS EN PANTALLA
print("\n==============================================================")
print(" RESULTADOS DE LA SIMULACIÓN ODMR")
print("==============================================================")

print("\n1. PARÁMETROS INTRODUCIDOS")
print(f"Campo magnético:            {campo_mT:.6f} mT")
print(f"Ángulo respecto al eje NV:  {angulo_grados:.6f} grados")
print(f"Perturbación transversal E: {perturbacion_E_MHz:.6f} MHz")

print("\n2. COMPONENTES DEL CAMPO")
print(f"Bx: {campo_x_T:.9f} T")
print(f"By: {campo_y_T:.9f} T")
print(f"Bz: {campo_z_T:.9f} T")
print(f"Componente longitudinal, |Bz|: {campo_longitudinal_real_mT:.9f} mT")

print("\n3. RESONANCIAS SIN CAMPO MAGNÉTICO")
print(f"Resonancia inferior: {frecuencia_inferior_sin_campo_GHz:.9f} GHz")
print(f"Resonancia superior: {frecuencia_superior_sin_campo_GHz:.9f} GHz")

print("\n4. RESONANCIAS CON CAMPO MAGNÉTICO")
print(f"Resonancia inferior:          {frecuencia_inferior_GHz:.9f} GHz")
print(f"Resonancia superior:          {frecuencia_superior_GHz:.9f} GHz")
print(f"Separación entre resonancias: {separacion_resonancias_GHz:.9f} GHz")

print("\n5. ESTIMACIÓN DEL CAMPO LONGITUDINAL")
print(f"Componente longitudinal real:     {campo_longitudinal_real_mT:.9f} mT")
print(f"Componente longitudinal estimada: {campo_longitudinal_estimado_mT:.9f} mT")
print(f"Error de la estimación:           {error_estimacion_mT:.9f} mT")
print(f"Error porcentual:                 {error_estimacion_porcentual:.6f} %")

print("\n6. ARCHIVOS GENERADOS")
print("odmr_without_magnetic_field.png")
print("odmr_with_magnetic_field.png")
print("frequencies_vs_field.png")
print("frequencies_vs_angle.png")
print("evolution_odmr_with_magnetic_field.png")
print("results_odmr.txt")

print("\n==============================================================")
print(f"Resultados guardados en: {carpeta_resultados.resolve()}")
print("\n · Simulación finalizada correctamente. ")
print("==============================================================")
# ============================================================
