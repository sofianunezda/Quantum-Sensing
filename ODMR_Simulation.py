"""
Simulación de espectros ODMR de un centro NV- en diamante.

Este programa complementa el manual:
"Introducción a los Sensores Cuánticos Basados en Centros NV en Diamante".

El programa permite:

1.  Calcular las resonancias ODMR a partir del Hamiltoniano de espín.
2.  Representar espectros sin campo y con campo magnético.
3.  Estudiar el efecto del módulo y la orientación del campo.
4.  Comparar las cuatro orientaciones cristalográficas de los centros NV.
5.  Analizar la dependencia con la temperatura.
6.  Simular el ensanchamiento producido por la potencia de microondas.
7.  Añadir ruido fotónico y ruido técnico.
8.  Estimar la sensibilidad magnética.
9.  Estudiar el efecto de la perturbación transversal E.
10. Guardar las gráficas y un informe de resultados.

Autora: Sofía Núñez de Andrés
Prácticas extracurriculares en el CINN, Julio- Agosto 2026
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

print("\n==========================================================================")
print(" SIMULADOR ODMR - CENTROS NV EN DIAMANTE")
print("==========================================================================")

# ==========================================================================
# CONSTANTES FÍSICAS DEL CENTRO NV
# ==========================================================================

D_GHz = 2.87
gamma_e_GHz_T = 28.0

# ==========================================================================
# DEPENDENCIA DE D CON LA TEMPERATURA
# ==========================================================================

temperatura_referencia_C = 25.0
coeficiente_temperatura_D_GHz_C = -74e-6

# ==========================================================================
# PARÁMETROS DEL ESPECTRO ODMR
# ==========================================================================

numero_puntos = 1000
ancho_fwhm_GHz = 0.003
contraste_por_resonancia = 0.080
margen_frecuencia_GHz = 0.03

# ==========================================================================
# POTENCIA DE MICROONDAS
# ==========================================================================

potencia_microondas_referencia = 1.0
potencias_microondas_comparacion = np.array([0.25, 0.5, 1.0, 2.0, 4.0, 9.0])

# ==========================================================================
# DETECCIÓN DE FOTONES Y RUIDO
# ==========================================================================

tasa_fotones_Hz = 2_500_000
tiempo_integracion_s = 0.01
ruido_tecnico_relativo = 0.002

semilla_aleatoria = 7
generador_aleatorio = np.random.default_rng(semilla_aleatoria)

# ==========================================================================
# SENSIBILIDAD MAGNÉTICA
# ==========================================================================

factor_perfil_lorentziano = 0.77

# ==========================================================================
# CARPETA DE RESULTADOS
# ==========================================================================

carpeta_script = Path(__file__).resolve().parent
carpeta_resultados = carpeta_script / "results_odmr"
carpeta_resultados.mkdir(parents=True, exist_ok=True)

# ==========================================================================
# MATRICES DE ESPÍN PARA S=1. 
# Base utilizada: |+1>, |0>, |-1>
# ==========================================================================

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

# ==========================================================================
# CUATRO ORIENTACIONES CRISTALOGRÁFICAS DE LOS CENTROS NV
# ==========================================================================

orientaciones_NV = {
    "NV [1, 1, 1]": np.array([1.0, 1.0, 1.0]) / np.sqrt(3),
    "NV [1, -1, -1]": np.array([1.0, -1.0, -1.0]) / np.sqrt(3),
    "NV [-1, 1, -1]": np.array([-1.0, 1.0, -1.0]) / np.sqrt(3),
    "NV [-1, -1, 1]": np.array([-1.0, -1.0, 1.0]) / np.sqrt(3),
}

# ==========================================================================
# FUNCIONES
# ==========================================================================

def calcular_resonancias(campo_x_T, 
                         campo_y_T, 
                         campo_z_T, 
                         perturbacion_E_GHz, 
                         D_actual_GHz
                        ):
    """
    Construye el Hamiltoniano del estado fundamental del centro NV
    y calcula sus frecuencias de transición.
    """

    hamiltoniano_zfs = D_actual_GHz * (Sz @ Sz)

    # Se omite el término -(2/3) D I porque solo produce un
    # desplazamiento común de todos los niveles de energía.

    hamiltoniano_transversal = perturbacion_E_GHz * ((Sx @ Sx) - (Sy @ Sy))

    hamiltoniano_zeeman = gamma_e_GHz_T * (campo_x_T * Sx + campo_y_T * Sy + campo_z_T * Sz)

    hamiltoniano_total = (hamiltoniano_zfs + hamiltoniano_transversal + hamiltoniano_zeeman)

    energias_GHz = np.linalg.eigvalsh(hamiltoniano_total)

    frecuencia_inferior_GHz = energias_GHz[1] - energias_GHz[0]
    frecuencia_superior_GHz = energias_GHz[2] - energias_GHz[0]

    return frecuencia_inferior_GHz, frecuencia_superior_GHz, energias_GHz


def calcular_perfil_lorentziano(frecuencias_GHz, 
                                frecuencia_central_GHz, 
                                ancho_fwhm_GHz
                                ):
    """
    Calcula un perfil lorentziano normalizado.

    El perfil vale 1 en la frecuencia central y disminuye 
    simétricamente al alejarnos de la resonancia.
    """

    semiancho_GHz = ancho_fwhm_GHz / 2

    diferencia_frecuencia_GHz = (frecuencias_GHz - frecuencia_central_GHz)

    perfil_lorentziano = semiancho_GHz**2 / (diferencia_frecuencia_GHz**2 + semiancho_GHz**2)

    return perfil_lorentziano

def calcular_ancho_fwhm_con_potencia(ancho_fwhm_referencia_GHz, 
                                     potencia_microondas_relativa
                                     ):
    """
    Calcula la anchura FWHM efectiva debida al ensanchamiento 
    producido por la potencia de microondas.
    """
    if potencia_microondas_relativa < 0:
        raise ValueError("La potencia relativa de microondas no puede ser negativa.")

    return ancho_fwhm_referencia_GHz * np.sqrt(potencia_microondas_relativa)

def calcular_sensibilidad_magnetica(ancho_fwhm_GHz, 
                                    contraste, 
                                    tasa_fotones_Hz, 
                                    gamma_e_GHz_T, 
                                    factor_perfil
                                    ):
    """
    Estima la sensibilidad magnética limitada por el ruido fotónico.
    """
    if ancho_fwhm_GHz <= 0:
        raise ValueError("La anchura FWHM debe ser mayor que cero.")

    if contraste <= 0:
        raise ValueError("El contraste debe ser mayor que cero.")

    if tasa_fotones_Hz <= 0:
        raise ValueError("La tasa de fotones debe ser mayor que cero.")

    sensibilidad_T_raiz_Hz = (factor_perfil * ancho_fwhm_GHz) / (gamma_e_GHz_T * contraste * np.sqrt(tasa_fotones_Hz))

    return sensibilidad_T_raiz_Hz

def obtener_resonancias_unicas(resonancias_GHz, 
                               tolerancia_GHz
                               ):
    """
    Elimina resonancias coincidentes dentro de una tolerancia.

    Esto evita representar dos veces una misma línea ODMR cuando
    dos transiciones tienen frecuencias prácticamente iguales.
    """
    resonancias_ordenadas_GHz = sorted(resonancias_GHz)
    resonancias_unicas_GHz = []

    for frecuencia_GHz in resonancias_ordenadas_GHz:
        es_duplicada = any(
            np.isclose(
                frecuencia_GHz, 
                frecuencia_anterior_GHz, 
                atol=tolerancia_GHz, 
                rtol=0.0
            ) 
            for frecuencia_anterior_GHz in resonancias_unicas_GHz
        )

        if not es_duplicada:
            resonancias_unicas_GHz.append(frecuencia_GHz)

    return resonancias_unicas_GHz

def construir_curva_odmr(frecuencias_GHz, 
                         resonancias_GHz, 
                         contraste, 
                         ancho_fwhm_utilizado_GHz
                         ):
    """
    Construye una curva ODMR ideal con fluorescencia normalizada.

    Cada resonancia produce una disminución lorentziana de la fluorescencia
    """
    fluorescencia_normalizada = np.ones_like(frecuencias_GHz, dtype=float)

    tolerancia_degeneracion_GHz = ancho_fwhm_utilizado_GHz / 100

    resonancias_unicas_GHz = obtener_resonancias_unicas(
        resonancias_GHz, 
        tolerancia_degeneracion_GHz
    )

    for frecuencia_central_GHz in resonancias_unicas_GHz:
        perfil_lorentziano = calcular_perfil_lorentziano(frecuencias_GHz, 
                                                         frecuencia_central_GHz, 
                                                         ancho_fwhm_utilizado_GHz
                                                         )

        fluorescencia_normalizada -= (contraste * perfil_lorentziano)

    return fluorescencia_normalizada

def añadir_ruido_experimental(curva_odmr, 
                              tasa_fotones_Hz, 
                              tiempo_integracion_s, 
                              ruido_tecnico_relativo
                              ):
    """
    Añade ruido fotónico de Poisson y ruido técnico gaussiano
    a una curva ODMR ideal.
    """
    numero_fotones_por_punto = tasa_fotones_Hz * tiempo_integracion_s

    if numero_fotones_por_punto <= 0:
        raise ValueError("El número de fotones por punto debe ser mayor que cero.")
    
    conteos_esperados = curva_odmr * numero_fotones_por_punto
    conteos_esperados = np.clip(conteos_esperados, 0.0, None)

    conteos_detectados = generador_aleatorio.poisson(conteos_esperados)

    curva_con_ruido_fotonico = conteos_detectados / numero_fotones_por_punto

    ruido_tecnico = generador_aleatorio.normal(loc=0.0, scale=ruido_tecnico_relativo, size=len(curva_odmr))

    curva_con_ruido = curva_con_ruido_fotonico + ruido_tecnico
    curva_con_ruido = np.clip(curva_con_ruido, 0.0, 1.02)

    return curva_con_ruido

def construir_espectro_odmr_completo(frecuencias_GHz, 
                                     resonancias_GHz,
                                     contraste, 
                                     ancho_fwhm_utilizado_GHz, 
                                     tasa_fotones_Hz, 
                                     tiempo_integracion_s, 
                                     ruido_tecnico_relativo
                                     ):
    """
    Construye una curva ODMR ideal y su versión con ruido experimental.
    """

    curva_ideal = construir_curva_odmr(frecuencias_GHz=frecuencias_GHz,
                                       resonancias_GHz=resonancias_GHz,
                                       contraste=contraste,
                                       ancho_fwhm_utilizado_GHz=ancho_fwhm_utilizado_GHz
                                       )
    curva_con_ruido = añadir_ruido_experimental(curva_odmr=curva_ideal,
                                                tasa_fotones_Hz=tasa_fotones_Hz,
                                                tiempo_integracion_s=tiempo_integracion_s,
                                                ruido_tecnico_relativo=ruido_tecnico_relativo
                                                )

    return curva_ideal, curva_con_ruido

def guardar_grafica(nombre_archivo, 
                    resolucion_dpi=300
                    ):
    """
    Guarda la figura activa dentro de la carpeta de resultados.

    Si el nombre no incluye extensión, se añade automáticamente
    la extensión PNG.
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

def calcular_componentes_campo(modulo_campo_T, 
                               angulo_rad
                               ):
    """
    Calcula las componentes cartesianas del campo magnético
    a partir de su módulo y del ángulo respecto al eje NV.

    Se supone que el campo está contenido en el plano XZ y que
    el eje Z coincide con el eje del centro NV.

        Bx = B·sin(\u03B8)
        By = 0.0
        Bz = B·cos(\u03B8)
    """

    campo_x_T = modulo_campo_T * np.sin(angulo_rad)
    campo_y_T = 0.0
    campo_z_T = modulo_campo_T * np.cos(angulo_rad)

    return campo_x_T, campo_y_T, campo_z_T

def calcular_vector_campo(modulo_campo_T, 
                          angulo_polar_rad, 
                          angulo_azimutal_rad
                          ):
    """
    Construye el vector del campo magnético en el sistema del cristal.

        Bx = B·sin(\u03B8)·cos(\u03C6)
        By = B·sin(\u03B8)·sin(\u03C6)
        Bz = B·cos(\u03B8)

    \u03B8 es el ángulo polar medido desde el eje Z, en radianes.
    \u03C6 es el ángulo azimutal medido desde el eje X en el plano XY, en radianes.
    """
    campo_x_T = modulo_campo_T * np.sin(angulo_polar_rad) * np.cos(angulo_azimutal_rad)
    campo_y_T = modulo_campo_T * np.sin(angulo_polar_rad) * np.sin(angulo_azimutal_rad)
    campo_z_T = modulo_campo_T * np.cos(angulo_polar_rad) 

    vector_campo_T = np.array([campo_x_T, campo_y_T, campo_z_T])

    return vector_campo_T

def construir_base_local_NV(eje_NV):
    """
    Construye una base cartesiana local cuyo eje z coincide con el eje NV.

    El eje Z local coincide con el eje del centro NV. Los ejes X e Y locales 
    son perpendiculares entre sí y al eje NV.
    """
    eje_z_local = eje_NV / np.linalg.norm(eje_NV)

    vector_referencia = np.array([0.0, 0.0, 1.0])

    # Evita un producto vectorial casi nulo.
    if abs(np.dot(eje_z_local, vector_referencia)) > 0.99:
        vector_referencia = np.array([1.0, 0.0, 0.0])

    eje_x_local = np.cross(vector_referencia, eje_z_local)
    eje_x_local /= np.linalg.norm(eje_x_local)

    eje_y_local = np.cross(eje_z_local, eje_x_local)
    eje_y_local /= np.linalg.norm(eje_y_local)

    return eje_x_local, eje_y_local, eje_z_local

def transformar_campo_a_base_NV(vector_campo_T, 
                                eje_NV
                                ):
    """
    Transforma el campo magnético global al 
    sistema de referencia local del centro NV.
    """
    eje_x_local, eje_y_local, eje_z_local = construir_base_local_NV(eje_NV)

    campo_x_local_T = float(np.dot(vector_campo_T, eje_x_local))
    campo_y_local_T = float(np.dot(vector_campo_T, eje_y_local))
    campo_z_local_T = float(np.dot(vector_campo_T, eje_z_local))

    return campo_x_local_T, campo_y_local_T, campo_z_local_T

def calcular_resonancias_cuatro_orientaciones(vector_campo_T, 
                                              perturbacion_E_GHz, 
                                              D_actual_GHz
                                              ):
    """
    Calcula las dos resonancias ODMR de cada una de las cuatro 
    orientaciones cristalográficas de los centros NV.
    """
    resultados_orientaciones = {}

    for nombre_NV, eje_NV in orientaciones_NV.items():
        campo_x_local_T, campo_y_local_T, campo_z_local_T = (
            transformar_campo_a_base_NV(
                vector_campo_T, 
                eje_NV
            )
        )

        frecuencia_inferior_GHz, frecuencia_superior_GHz, energias_GHz = calcular_resonancias(
            campo_x_T = campo_x_local_T,
            campo_y_T = campo_y_local_T,
            campo_z_T = campo_z_local_T,
            perturbacion_E_GHz = perturbacion_E_GHz,
            D_actual_GHz = D_actual_GHz
        )

        resultados_orientaciones[nombre_NV] = {
            "campo_x_local_T": campo_x_local_T,
            "campo_y_local_T": campo_y_local_T,
            "campo_z_local_T": campo_z_local_T,
            "frecuencia_inferior_GHz": frecuencia_inferior_GHz,
            "frecuencia_superior_GHz": frecuencia_superior_GHz,
            "energias_GHz": energias_GHz
        }

    return resultados_orientaciones

def construir_curva_odmr_cuatro_orientaciones(frecuencias_GHz, 
                                              resultados_orientaciones, 
                                              contraste_total, 
                                              ancho_fwhm_utilizado_GHz
                                              ):
    """
    Construye el espectro ODMR conjunto de cuatro familias NV.

    Se supone que las cuatro orientaciones están igualmente pobladas, 
    por lo que el contraste total se reparte entre ellas.
    """

    fluorescencia_normalizada = np.ones_like(frecuencias_GHz, dtype=float)
    numero_orientaciones = len(resultados_orientaciones)
    contraste_por_familia = contraste_total / numero_orientaciones
    tolerancia_degeneracion_GHz = ancho_fwhm_utilizado_GHz / 100

    for resultados_NV in resultados_orientaciones.values():
        resonancias_familia_GHz = [resultados_NV["frecuencia_inferior_GHz"], resultados_NV["frecuencia_superior_GHz"]]

        resonancias_unicas_GHz = obtener_resonancias_unicas(resonancias_familia_GHz, 
                                                            tolerancia_degeneracion_GHz
                                                            )
        
        for frecuencia_central_GHz in resonancias_familia_GHz:
            perfil_lorentziano = calcular_perfil_lorentziano(frecuencias_GHz, frecuencia_central_GHz, ancho_fwhm_utilizado_GHz)
            fluorescencia_normalizada -= contraste_por_familia * perfil_lorentziano

    return fluorescencia_normalizada

def calcular_D_con_temperatura(temperatura_C):
    """
    Calcula la división de campo cero D en función de la temperatura
    determinada.

    Se utiliza una aproximación lineal alrededor de la temperatura 
    de referencia.
    """

    cambio_temperatura_C = temperatura_C - temperatura_referencia_C

    desplazamiento_D_GHz = cambio_temperatura_C * coeficiente_temperatura_D_GHz_C

    D_temperatura_GHz = D_GHz + desplazamiento_D_GHz

    return D_temperatura_GHz

# ==========================================================================
# DATOS INTRODUCIDOS POR EL USUARIO
# ==========================================================================

def pedir_float(mensaje, 
                valor_por_defecto, 
                valor_minimo=None,
                valor_maximo=None
                ):
    """
    Solicita un número al usuario y comprueba que esté dentro 
    del intervalo permitido.

    Si el usuario pulsa Intro sin escribir nada, se utiliza el 
    valor predeterminado.
    """
    while True:
        entrada = input(mensaje).strip()

        if entrada == "":
            valor = float(valor_por_defecto)
        else:
            # Se acepta tanto la coma como el punto decimal.
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

# ==========================================================================
# COMPONENTES DE LOS CAMPOS
# ==========================================================================

print("\nIntroduce los parámetros de la simulación. ")
print("Pulsa Intro para utilizar el valor indicado entre corchetes. ")

# Campo magnético para el modelo de un único eje NV
campo_mT    = pedir_float("\n -> Campo magnético (mT) [3]: ", 
                          valor_por_defecto=3.0, 
                          valor_minimo=0.0)
campo_T     = campo_mT / 1000

angulo_grados = pedir_float("\n -> Ángulo respecto al eje NV (\u00B0) [0]: ", 
                            valor_por_defecto=0.0, 
                            valor_minimo=0.0, 
                            valor_maximo=90.0)
angulo_rad    = np.radians(angulo_grados)

# Orientación del campo en el sistema global del cristal 
angulo_polar_cristal_grados = pedir_float("\n -> Ángulo polar del campo respecto al eje Z del cristal (\u00B0) [0]: ", 
                                          valor_por_defecto=0.0, 
                                          valor_minimo=0.0, 
                                          valor_maximo=180.0)
angulo_polar_cristal_rad = np.radians(angulo_polar_cristal_grados)

angulo_azimutal_cristal_grados = pedir_float("\n -> Ángulo azimutal del campo en el plano XY (\u00B0) [0]: ", 
                                             valor_por_defecto=0.0, 
                                             valor_minimo=0.0, 
                                             valor_maximo=360.0)
angulo_azimutal_cristal_rad = np.radians(angulo_azimutal_cristal_grados)

# Perturbación transversal 
perturbacion_E_MHz = pedir_float("\n -> Perturbación transversal E (MHz) [0]: ", 
                                 valor_por_defecto=0.0, 
                                 valor_minimo=0.0)
perturbacion_E_GHz = perturbacion_E_MHz / 1000

# Temperatura
temperatura_C = pedir_float("\n -> Temperatura (\u00B0C) [25]: ", 
                            valor_por_defecto=25.0)

temperatura_minima_C = pedir_float("\n -> Temperatura mínima para el barrido (\u00B0C) [-50]: ", 
                                   valor_por_defecto=-50.0)
temperatura_maxima_C = pedir_float("\n -> Temperatura máxima para el barrido (\u00B0C) [150]: ", 
                                   valor_por_defecto=150.0)

while temperatura_maxima_C <= temperatura_minima_C:
    print("La temperatura máxima debe ser mayor que la temperatura mínima.")
    temperatura_maxima_C = pedir_float("\n -> Temperatura máxima para el barrido (\u00B0C) [150]: ", 
                                       valor_por_defecto=150.0)

# Potencia de microondas
potencia_microondas_relativa = pedir_float("\n -> Potencia relativa de microondas [1]: ", 
                                           valor_por_defecto=1.0, 
                                           valor_minimo=0.0)

# Detecciones de fotones y ruido
tasa_fotones_Hz        = pedir_float("\n -> Tasa de fotones detectados (fotones/s) [2500000]: ", 
                                     valor_por_defecto=2_500_000, 
                                     valor_minimo=1.0)
tiempo_integracion_s   = pedir_float("\n -> Tiempo de integración por punto (s) [0.01]: ", 
                                     valor_por_defecto=0.01, 
                                     valor_minimo=1e-6) 
ruido_tecnico_relativo = pedir_float("\n -> Ruido técnico relativo [0.002]: ", 
                                     valor_por_defecto=0.002, 
                                     valor_minimo=0.0) 

# ==========================================================================
# PARÁMETROS DERIVADOS DE LA SIMULACIÓN
# ==========================================================================

D_actual_GHz  = calcular_D_con_temperatura(temperatura_C)

ancho_fwhm_actual_GHz = calcular_ancho_fwhm_con_potencia(ancho_fwhm_referencia_GHz=ancho_fwhm_GHz, 
                                                         potencia_microondas_relativa=potencia_microondas_relativa)

factor_ensanchamiento = ancho_fwhm_actual_GHz / ancho_fwhm_GHz

numero_fotones_por_punto        = tasa_fotones_Hz * tiempo_integracion_s
ruido_fotonico_relativo         = 1 / np.sqrt(numero_fotones_por_punto)
ruido_total_relativo_aproximado = np.sqrt(ruido_fotonico_relativo**2 + ruido_tecnico_relativo**2)

# Campo en el sistema local de un único eje NV
campo_x_T, campo_y_T, campo_z_T = calcular_componentes_campo(campo_T, angulo_rad)

# Campo en el sistema global del cristal
vector_campo_cristal_T = calcular_vector_campo(modulo_campo_T=campo_T,
                                               angulo_polar_rad=angulo_polar_cristal_rad,
                                               angulo_azimutal_rad=angulo_azimutal_cristal_rad)

# ==========================================================================
# SENSIBILIDAD MAGNÉTICA
# ==========================================================================

sensibilidad_T_raiz_Hz  = calcular_sensibilidad_magnetica(ancho_fwhm_GHz=ancho_fwhm_actual_GHz, 
                                                          contraste=contraste_por_resonancia, 
                                                          tasa_fotones_Hz=tasa_fotones_Hz,
                                                          gamma_e_GHz_T=gamma_e_GHz_T, 
                                                          factor_perfil=factor_perfil_lorentziano)
sensibilidad_uT_raiz_Hz = sensibilidad_T_raiz_Hz * 1e6
sensibilidad_nT_raiz_Hz = sensibilidad_T_raiz_Hz * 1e9

campo_minimo_detectable_T  = sensibilidad_T_raiz_Hz / np.sqrt(tiempo_integracion_s)
campo_minimo_detectable_uT = campo_minimo_detectable_T * 1e6
campo_minimo_detectable_nT = campo_minimo_detectable_T * 1e9

# ==========================================================================
# RESONANCIAS DE UN ÚNICO EJE NV
# ==========================================================================

frecuencia_inferior_sin_campo_GHz, frecuencia_superior_sin_campo_GHz, energias_sin_campo_GHz = calcular_resonancias(
    campo_x_T=0.0, campo_y_T=0.0, campo_z_T=0.0, perturbacion_E_GHz=perturbacion_E_GHz, D_actual_GHz=D_actual_GHz)

frecuencia_inferior_GHz, frecuencia_superior_GHz, energias_con_campo_GHz = calcular_resonancias(
    campo_x_T=campo_x_T, campo_y_T=campo_y_T, campo_z_T=campo_z_T, perturbacion_E_GHz=perturbacion_E_GHz, D_actual_GHz=D_actual_GHz)

resonancias_sin_campo_GHz = [frecuencia_inferior_sin_campo_GHz, frecuencia_superior_sin_campo_GHz]
resonancias_con_campo_GHz = [frecuencia_inferior_GHz, frecuencia_superior_GHz]

todas_las_resonancias_GHz = (resonancias_sin_campo_GHz + resonancias_con_campo_GHz)

frecuencia_minima_GHz = min(todas_las_resonancias_GHz) - margen_frecuencia_GHz
frecuencia_maxima_GHz = max(todas_las_resonancias_GHz) + margen_frecuencia_GHz

if frecuencia_minima_GHz >= frecuencia_maxima_GHz:
    raise ValueError("El intervalo de frecuencias calculado no es válido.")

frecuencias_GHz = np.linspace(frecuencia_minima_GHz, frecuencia_maxima_GHz, numero_puntos)

# ==========================================================================
# ESPECTROS DE UN ÚNICO EJE NV
# ==========================================================================

curva_sin_campo,  curva_sin_campo_con_ruido = construir_espectro_odmr_completo(
    frecuencias_GHz=frecuencias_GHz, resonancias_GHz=resonancias_sin_campo_GHz, 
    contraste= contraste_por_resonancia, ancho_fwhm_utilizado_GHz=ancho_fwhm_actual_GHz, 
    tasa_fotones_Hz=tasa_fotones_Hz, tiempo_integracion_s=tiempo_integracion_s, 
    ruido_tecnico_relativo=ruido_tecnico_relativo)

curva_con_campo,  curva_con_campo_con_ruido = construir_espectro_odmr_completo(
    frecuencias_GHz=frecuencias_GHz, resonancias_GHz=resonancias_con_campo_GHz, 
    contraste= contraste_por_resonancia, ancho_fwhm_utilizado_GHz=ancho_fwhm_actual_GHz,
    tasa_fotones_Hz=tasa_fotones_Hz, tiempo_integracion_s=tiempo_integracion_s, 
    ruido_tecnico_relativo=ruido_tecnico_relativo)

# ==========================================================================
# CUATRO ORIENTACIONES CRISTALOGRÁFICAS NV
# ==========================================================================

resultados_cuatro_NV = calcular_resonancias_cuatro_orientaciones(
    vector_campo_T=vector_campo_cristal_T,
    perturbacion_E_GHz=perturbacion_E_GHz,
    D_actual_GHz=D_actual_GHz
    )

resonancias_cuatro_NV_GHz = []

for resultados_NV in resultados_cuatro_NV.values():
    resonancias_cuatro_NV_GHz.extend(
        [resultados_NV["frecuencia_inferior_GHz"],
         resultados_NV["frecuencia_superior_GHz"]]
    )

frecuencia_minima_cuatro_NV_GHz = min(resonancias_cuatro_NV_GHz) - margen_frecuencia_GHz
frecuencia_maxima_cuatro_NV_GHz = max(resonancias_cuatro_NV_GHz) + margen_frecuencia_GHz

frecuencias_cuatro_NV_GHz = np.linspace(frecuencia_minima_cuatro_NV_GHz, frecuencia_maxima_cuatro_NV_GHz, numero_puntos)

curva_cuatro_NV = construir_curva_odmr_cuatro_orientaciones(frecuencias_GHz=frecuencias_cuatro_NV_GHz, 
                                                            resultados_orientaciones=resultados_cuatro_NV,
                                                            contraste_total=contraste_por_resonancia,
                                                            ancho_fwhm_utilizado_GHz=ancho_fwhm_actual_GHz)

curva_cuatro_NV_con_ruido = añadir_ruido_experimental(curva_odmr=curva_cuatro_NV,
                                                      tasa_fotones_Hz=tasa_fotones_Hz,
                                                      tiempo_integracion_s=tiempo_integracion_s,
                                                      ruido_tecnico_relativo=ruido_tecnico_relativo)

# ==========================================================================
# DEPENDENCIA CON LA TEMPERATURA
# ==========================================================================

temperaturas_C = np.linspace(temperatura_minima_C, temperatura_maxima_C, 201)

valores_D_temperatura_GHz = []
frecuencias_inferiores_temperatura_GHz = []
frecuencias_superiores_temperatura_GHz = []

for temperatura_actual_C in temperaturas_C:
    D_temperatura_actual_GHz = calcular_D_con_temperatura(temperatura_actual_C)

    frecuencia_inferior_actual_GHz, frecuencia_superior_actual_GHz, _ = calcular_resonancias(
    campo_x_T=campo_x_T,
    campo_y_T=campo_y_T,
    campo_z_T=campo_z_T,
    perturbacion_E_GHz=perturbacion_E_GHz,
    D_actual_GHz=D_temperatura_actual_GHz,
    )

    valores_D_temperatura_GHz.append(D_temperatura_actual_GHz)
    frecuencias_inferiores_temperatura_GHz.append(frecuencia_inferior_actual_GHz)
    frecuencias_superiores_temperatura_GHz.append(frecuencia_superior_actual_GHz)

valores_D_temperatura_GHz = np.array(valores_D_temperatura_GHz)
frecuencias_inferiores_temperatura_GHz = np.array(frecuencias_inferiores_temperatura_GHz)
frecuencias_superiores_temperatura_GHz = np.array(frecuencias_superiores_temperatura_GHz)

desplazamiento_D_GHz = D_actual_GHz - D_GHz
desplazamiento_D_MHz = desplazamiento_D_GHz * 1000
desplazamiento_D_kHz = desplazamiento_D_GHz * 1_000_000

# ==========================================================================
# COMPARACIÓN DE POTENCIAS DE MICROONDAS
# ==========================================================================

anchos_fwhm_comparacion_GHz = []
curvas_potencia_microondas = []

for potencia_actual in potencias_microondas_comparacion:
    ancho_fwhm_potencia_actual_GHz = calcular_ancho_fwhm_con_potencia(
        ancho_fwhm_referencia_GHz=ancho_fwhm_GHz, 
        potencia_microondas_relativa=potencia_actual)

    curva_potencia_actual = construir_curva_odmr(
    frecuencias_GHz=frecuencias_GHz,
    resonancias_GHz=resonancias_con_campo_GHz,
    contraste=contraste_por_resonancia,
    ancho_fwhm_utilizado_GHz=ancho_fwhm_potencia_actual_GHz)

    anchos_fwhm_comparacion_GHz.append(ancho_fwhm_potencia_actual_GHz)
    curvas_potencia_microondas.append(curva_potencia_actual)

anchos_fwhm_comparacion_GHz = np.array(anchos_fwhm_comparacion_GHz)

# ==========================================================================
# SENSIBILIDAD FRENTE A LA TASA DE FOTONES
# ==========================================================================

tasas_fotones_sensibilidades_Hz = np.logspace(4, 8, 200)
sensibilidad_tasa_T_raiz_Hz   = []

for tasa_actual_Hz in tasas_fotones_sensibilidades_Hz:
    sensibilidad_actual_T_raiz_Hz = calcular_sensibilidad_magnetica(ancho_fwhm_GHz=ancho_fwhm_actual_GHz,
                                                                    contraste=contraste_por_resonancia,
                                                                    tasa_fotones_Hz=tasa_actual_Hz,
                                                                    gamma_e_GHz_T=gamma_e_GHz_T,
                                                                    factor_perfil=factor_perfil_lorentziano)
    sensibilidad_tasa_T_raiz_Hz.append(sensibilidad_actual_T_raiz_Hz)

sensibilidad_tasa_nT_raiz_Hz = np.array(sensibilidad_tasa_T_raiz_Hz) * 1e9

# ==========================================================================
# BARRIDO DE LA PERTURBACIÓN TRANSVERSAL 
# ==========================================================================

perturbaciones_E_MHz = np.linspace(0.0, 10.0, 200)

frecuencias_inferiores_E_GHz = []
frecuencias_superiores_E_GHz = []

for perturbacion_E_actual_MHz in perturbaciones_E_MHz:
    perturbacion_E_actual_GHz = perturbacion_E_actual_MHz / 1000
    
    frecuencia_inferior_actual_GHz, frecuencia_superior_actual_GHz, _ = calcular_resonancias(
        campo_x_T=0.0, 
        campo_y_T=0.0, 
        campo_z_T=0.0, 
        perturbacion_E_GHz=perturbacion_E_actual_GHz, 
        D_actual_GHz=D_actual_GHz)

    frecuencias_inferiores_E_GHz.append(frecuencia_inferior_actual_GHz)
    frecuencias_superiores_E_GHz.append(frecuencia_superior_actual_GHz)

frecuencias_inferiores_E_GHz = np.array(frecuencias_inferiores_E_GHz)
frecuencias_superiores_E_GHz = np.array(frecuencias_superiores_E_GHz)

# ==========================================================================
# GRÁFICA 1: CURVA ODMR SIN CAMPO MAGNÉTICO
# ==========================================================================

figura_1, eje = plt.subplots(figsize=(8, 5))

eje.plot(frecuencias_GHz, 
         curva_sin_campo, 
         color="tab:blue", 
         linewidth=2.0, 
         label="Espectro ideal"
         )
eje.plot(frecuencias_GHz, 
         curva_sin_campo_con_ruido, 
         color="tab:green", 
         linewidth=1.0, 
         alpha=0.65, 
         label="Espectro con ruido"
         )

eje.axvline(D_actual_GHz, 
            linestyle="--", 
            linewidth=1.5, 
            color="tab:red", 
            alpha=0.8, 
            label=f"D(T) = {D_actual_GHz:.6f} GHz"
            )

eje.set_xlabel("Frecuencia de microondas (GHz)")
eje.set_ylabel("Fluorescencia normalizada")

eje.set_title(f"Espectro ODMR sin campo magnético\n" 
              f"T = {temperatura_C:.1f} \u00B0C\n"
              f"E = {perturbacion_E_MHz:.3f} MHz"
              )

eje.set_xlim(frecuencia_minima_GHz, frecuencia_maxima_GHz)
eje.set_ylim(0.88, 1.01)

eje.grid(True, alpha=0.3)
eje.legend(loc="best")

guardar_grafica("odmr_without_magnetic_field.png")

# ==========================================================================
# GRÁFICA 2: CURVA ODMR CON CAMPO MAGNÉTICO
# ==========================================================================

figura_2, eje = plt.subplots(figsize=(8, 5))

eje.plot(frecuencias_GHz, 
         curva_con_campo, 
         color="tab:green", 
         linewidth=2.0, 
         label="Espectro ideal"
         )
eje.plot(frecuencias_GHz, 
         curva_con_campo_con_ruido, 
         color="tab:blue", 
         linewidth=1.0, 
         alpha=0.65, 
         label="Espectro con ruido"
         )

tolerancia_degeneracion_GHz = ancho_fwhm_actual_GHz / 100
resonancias_visibles_GHz    = obtener_resonancias_unicas(resonancias_con_campo_GHz, 
                                                         tolerancia_degeneracion_GHz
                                                         )

if len(resonancias_visibles_GHz) == 1:
    eje.axvline(resonancias_visibles_GHz[0], 
                linestyle="--", 
                linewidth=1.5, 
                color="tab:red", 
                alpha=0.8, 
                label=
                "Resonancia degenerada = "
                f"{resonancias_visibles_GHz[0]:.6f} GHz"
                )

else:
    eje.axvline(frecuencia_inferior_GHz, 
                linestyle="--", 
                linewidth=1.5, 
                color="tab:red", 
                alpha=0.8, 
                label=
                "Resonancia inferior = " 
                f"{frecuencia_inferior_GHz:.6f} GHz"
                )
    eje.axvline(frecuencia_superior_GHz, 
                linestyle="--", 
                linewidth=1.5, 
                color="tab:orange", 
                alpha=0.8, 
                label=
                "Resonancia superior = " 
                f"{frecuencia_superior_GHz:.6f} GHz"
                )

eje.set_xlabel("Frecuencia de microondas (GHz)")
eje.set_ylabel("Fluorescencia normalizada")

eje.set_title(f"Espectro ODMR con campo magnético\n" 
              f"B = {campo_mT:.3f} mT," 
              f"θ = {angulo_grados:.1f}\u00B0," 
              f"E = {perturbacion_E_MHz:.3f} MHz," 
              f"T = {temperatura_C:.1f} \u00B0C"
              )

eje.set_xlim(frecuencia_minima_GHz, frecuencia_maxima_GHz)
eje.set_ylim(0.88, 1.01)

eje.grid(True, alpha=0.3)
eje.legend(loc="best")

guardar_grafica("odmr_with_magnetic_field.png")

# ==========================================================================
# GRÁFICA 3: FRECUENCIAS DE RESONANCIA FRENTE AL CAMPO
# ==========================================================================

campo_maximo_mT = max(5.0, campo_mT * 2)
campos_mT = np.linspace(0.0, campo_maximo_mT, 201)

frecuencias_inferiores_campo_GHz = []
frecuencias_superiores_campo_GHz = []

for campo_actual_mT in campos_mT:
    campo_actual_T = campo_actual_mT / 1000

    campo_x_actual_T, campo_y_actual_T, campo_z_actual_T = calcular_componentes_campo(
        modulo_campo_T=campo_actual_T, 
        angulo_rad=angulo_rad
        )

    frecuencia_inferior_actual_GHz, frecuencia_superior_actual_GHz, _ = calcular_resonancias(
        campo_x_T=campo_x_actual_T,
        campo_y_T=campo_y_actual_T,
        campo_z_T=campo_z_actual_T,
        perturbacion_E_GHz=perturbacion_E_GHz,
        D_actual_GHz=D_actual_GHz
        )

    frecuencias_inferiores_campo_GHz.append(frecuencia_inferior_actual_GHz)
    frecuencias_superiores_campo_GHz.append(frecuencia_superior_actual_GHz)

frecuencias_inferiores_campo_GHz = np.array(frecuencias_inferiores_campo_GHz)
frecuencias_superiores_campo_GHz = np.array(frecuencias_superiores_campo_GHz)

figura_3, eje = plt.subplots(figsize=(8, 5))

eje.plot(campos_mT,
         frecuencias_inferiores_campo_GHz,
         linewidth=2.0,
         color="tab:blue",
         label="Resonancia inferior"
         )
eje.plot(campos_mT,
         frecuencias_superiores_campo_GHz,
         linewidth=2.0,
         linestyle="--",
         color="tab:orange",
         label="Resonancia superior"
         )

eje.axvline(campo_mT,
            linestyle="--",
            linewidth=1.5,
            color="tab:red",
            alpha=0.8,
            label=f"Campo elegido = {campo_mT:.3f} mT")

eje.set_xlabel("Módulo del campo magnético (mT)")
eje.set_ylabel("Frecuencia de resonancia (GHz)")

eje.set_title("Dependencia de las resonancias con el campo magnético\n" 
              f"θ = {angulo_grados:.1f}\u00B0," 
              f"E = {perturbacion_E_MHz:.3f} MHz,"
              f"T = {temperatura_C:.1f} \u00B0C")

margen_frecuencia_grafica_GHz = 0.025
eje.set_xlim(0.0, campo_maximo_mT)
eje.set_ylim(min(frecuencias_inferiores_campo_GHz) - margen_frecuencia_grafica_GHz, 
             max(frecuencias_superiores_campo_GHz) + margen_frecuencia_grafica_GHz)

eje.grid(True, alpha=0.3)
eje.legend(loc="best")

guardar_grafica("frequencies_vs_field.png")

# ==========================================================================
# GRÁFICA 4: FRECUENCIAS DE RESONANCIA FRENTE AL ÁNGULO
# ==========================================================================

angulos_grados = np.linspace(0.0, 90.0, 181)

frecuencias_inferiores_angulo_GHz = []
frecuencias_superiores_angulo_GHz = []

for angulo_actual_grados in angulos_grados:
    angulo_actual_rad = np.radians(angulo_actual_grados)

    campo_x_actual_T, campo_y_actual_T, campo_z_actual_T = calcular_componentes_campo(modulo_campo_T=campo_T,
                                                                                      angulo_rad=angulo_actual_rad
                                                                                      )

    frecuencia_inferior_actual_GHz, frecuencia_superior_actual_GHz, _ = calcular_resonancias(
        campo_x_T=campo_x_actual_T,
        campo_y_T=campo_y_actual_T,
        campo_z_T=campo_z_actual_T,
        perturbacion_E_GHz=perturbacion_E_GHz,
        D_actual_GHz=D_actual_GHz
        )

    frecuencias_inferiores_angulo_GHz.append(frecuencia_inferior_actual_GHz)
    frecuencias_superiores_angulo_GHz.append(frecuencia_superior_actual_GHz)

frecuencias_inferiores_angulo_GHz = np.array(frecuencias_inferiores_angulo_GHz)
frecuencias_superiores_angulo_GHz = np.array(frecuencias_superiores_angulo_GHz)

figura_4, eje = plt.subplots(figsize=(8, 5))

eje.plot(angulos_grados,
         frecuencias_inferiores_angulo_GHz,
         linewidth=2.0,
         color="tab:blue",
         label="Resonancia inferior"
         )
eje.plot(angulos_grados,
         frecuencias_superiores_angulo_GHz,
         linewidth=2.0,
         color="tab:orange",
         label="Resonancia superior"
         )

eje.axvline(angulo_grados,
            linestyle="--",
            linewidth=1.5,
            color="tab:red",
            alpha=0.8,
            label=f"Ángulo seleccionado = {angulo_grados:.1f}\u00B0"
            )

eje.set_xlabel("Ángulo entre el campo magnético y el eje NV (\u00B0)")
eje.set_ylabel("Frecuencia de resonancia (GHz)")
eje.set_title("Dependencia angular de las resonancias\n" 
              f"B = {campo_mT:.3f} mT," 
              f"E = {perturbacion_E_MHz:.3f} MHz,"
              f"T = {temperatura_C:.1f} \u00B0C")

eje.set_xlim(0.0, 90.0)

margen_frecuencia_grafica_GHz = 0.025
eje.set_ylim(min(frecuencias_inferiores_angulo_GHz) - margen_frecuencia_grafica_GHz, 
             max(frecuencias_superiores_angulo_GHz) + margen_frecuencia_grafica_GHz)

eje.grid(True, alpha=0.3)
eje.legend(loc="best")

guardar_grafica("frequencies_vs_angle.png")

# ==========================================================================
# GRÁFICA 5: EVOLUCIÓN DEL ESPECTRO ODMR CON EL CAMPO
# ==========================================================================

campos_mostrados_mT = np.linspace(0.0, campo_maximo_mT, 6)
resonancias_evolucion_GHz = []

for campo_actual_mT in campos_mostrados_mT:
    campo_actual_T = campo_actual_mT / 1000

    campo_x_actual_T, campo_y_actual_T, campo_z_actual_T = calcular_componentes_campo(
        modulo_campo_T=campo_actual_T,
        angulo_rad=angulo_rad
        )

    frecuencia_inferior_actual_GHz, frecuencia_superior_actual_GHz, _ = calcular_resonancias(
        campo_x_T=campo_x_actual_T,
        campo_y_T=campo_y_actual_T,
        campo_z_T=campo_z_actual_T,
        perturbacion_E_GHz=perturbacion_E_GHz,
        D_actual_GHz=D_actual_GHz
        )

    resonancias_evolucion_GHz.append(
        (campo_actual_mT,
         frecuencia_inferior_actual_GHz,
         frecuencia_superior_actual_GHz)
         )

todas_las_resonancias_evolucion_GHz = []

for campo_actual_mT, frecuencia_inferior_actual_GHz, frecuencia_superior_actual_GHz in resonancias_evolucion_GHz:
    todas_las_resonancias_evolucion_GHz.extend([frecuencia_inferior_actual_GHz,
                                                frecuencia_superior_actual_GHz]
                                                )

frecuencia_minima_evolucion_GHz = min(todas_las_resonancias_evolucion_GHz) - margen_frecuencia_GHz
frecuencia_maxima_evolucion_GHz = max(todas_las_resonancias_evolucion_GHz) + margen_frecuencia_GHz

frecuencias_evolucion_GHz = np.linspace(frecuencia_minima_evolucion_GHz, 
                                        frecuencia_maxima_evolucion_GHz,
                                        numero_puntos
                                        )

figura_5, eje = plt.subplots(figsize=(8, 5))

for campo_actual_mT, frecuencia_inferior_actual_GHz, frecuencia_superior_actual_GHz in resonancias_evolucion_GHz:
    curva_actual, curva_actual_con_ruido = construir_espectro_odmr_completo(
        frecuencias_GHz=frecuencias_evolucion_GHz,
        resonancias_GHz=[frecuencia_inferior_actual_GHz,
                         frecuencia_superior_actual_GHz],
        contraste=contraste_por_resonancia,
        ancho_fwhm_utilizado_GHz=ancho_fwhm_actual_GHz, 
        tasa_fotones_Hz=tasa_fotones_Hz, 
        tiempo_integracion_s=tiempo_integracion_s, 
        ruido_tecnico_relativo=ruido_tecnico_relativo
        )

    eje.plot(frecuencias_evolucion_GHz,
             curva_actual,
             linewidth=2.2,
             label=f"B = {campo_actual_mT:.2f} mT"
             )
    eje.plot(frecuencias_evolucion_GHz,
             curva_actual_con_ruido,
             linewidth=0.8, 
             alpha=0.35
             )

eje.set_xlabel("Frecuencia de microondas (GHz)")
eje.set_ylabel("Fluorescencia normalizada")
eje.set_title("Evolución del espectro ODMR al aumentar el campo\n" 
              f"θ = {angulo_grados:.1f}\u00B0," 
              f"E = {perturbacion_E_MHz:.3f} MHz,"
              f"T = {temperatura_C:.1f} \u00B0C")

eje.set_xlim(frecuencia_minima_evolucion_GHz, 
             frecuencia_maxima_evolucion_GHz
             )
eje.set_ylim(0.88, 1.01)

eje.grid(True, alpha=0.3)
eje.legend(loc="best")

guardar_grafica("evolution_odmr_with_magnetic_field.png")

# ==========================================================================
# GRÁFICA 6: ESPECTRO ODMR DE LAS CUATRO ORIENTACIONES NV
# ==========================================================================

figura_6, eje = plt.subplots(figsize=(8, 5))

eje.plot(frecuencias_cuatro_NV_GHz, 
         curva_cuatro_NV, 
         linewidth=2.0, 
         color="tab:blue", 
         label="Espectro conjunto"
         )
eje.plot(frecuencias_cuatro_NV_GHz, 
         curva_cuatro_NV_con_ruido, 
         linewidth=1.0, 
         alpha=0.65, 
         color="pink", 
         label="Espectro con ruido"
         )

colores_NV = ["tab:green", "tab:red", "tab:orange", "tab:purple"]

for color, (nombre_NV, resultados_NV) in zip(colores_NV, resultados_cuatro_NV.items()):
    eje.axvline(resultados_NV["frecuencia_inferior_GHz"], 
                linestyle="-", 
                linewidth=1.2, 
                color=color, 
                alpha=0.7
                ) # Resonancia inferior
    eje.axvline(resultados_NV["frecuencia_superior_GHz"], 
                linestyle="--", 
                linewidth=1.2, 
                color=color, 
                alpha=0.7, 
                label=nombre_NV
                ) # Resonancia superior

eje.set_xlabel("Frecuencia de microondas (GHz)")
eje.set_ylabel("Fluorescencia normalizada")
eje.set_title(f"Espectro ODMR de las cuatro orientaciones NV\n" 
              f"B = {campo_mT:.3f} mT,"
              f"θ = {angulo_polar_cristal_grados:.1f}\u00B0," 
              f"φ = {angulo_azimutal_cristal_grados:.1f}\u00B0,"
              f"T = {temperatura_C:.1f} \u00B0C"
              )

eje.set_xlim(frecuencia_minima_cuatro_NV_GHz, frecuencia_maxima_cuatro_NV_GHz)
eje.set_ylim(0.88, 1.01)

eje.grid(True, alpha=0.3)
eje.legend(loc="best")

guardar_grafica("odmr_four_nv_orientations.png")

# ==========================================================================
# GRÁFICA 7: FRECUENCIAS DE RESONANCIA FRENTE A LA TEMPERATURA
# ==========================================================================

figura_7, eje = plt.subplots(figsize=(8, 5))

eje.plot(temperaturas_C, 
         frecuencias_inferiores_temperatura_GHz, 
         linewidth=2.0, 
         color="tab:blue", 
         label="Resonancia inferior"
         )
eje.plot(temperaturas_C, 
         frecuencias_superiores_temperatura_GHz, 
         linewidth=2.0, 
         linestyle="--", 
         color="tab:orange", 
         label="Resonancia superior"
         )

eje.axvline(temperatura_C, 
            linestyle="--", 
            linewidth=1.5, 
            color="tab:red", 
            alpha=0.8, 
            label=f"Temperatura elegida = {temperatura_C:.1f} \u00B0C"
            )

eje.set_xlabel("Temperatura (\u00B0C)")
eje.set_ylabel("Frecuencia de resonancia (GHz)")
eje.set_title(f"Dependencia de las resonancias ODMR con la temperatura\n" 
              f"B = {campo_mT:.3f} mT," 
              f"θ = {angulo_grados:.1f}\u00B0," 
              f"E = {perturbacion_E_MHz:.3f} MHz")

eje.set_xlim(temperatura_minima_C, temperatura_maxima_C)
margen_temperatura_grafica_GHz = 0.005
eje.set_ylim(np.min(frecuencias_inferiores_temperatura_GHz) - margen_temperatura_grafica_GHz, 
             np.max(frecuencias_superiores_temperatura_GHz) + margen_temperatura_grafica_GHz)

eje.grid(True, alpha=0.3)
eje.legend(loc="best")

guardar_grafica("frequencies_vs_temperature.png")

# ==========================================================================
# GRÁFICA 8: ESPECTROS ODMR PARA DIFERENTES POTENCIAS
# ==========================================================================

figura_8, eje = plt.subplots(figsize=(8, 5))

for potencia_actual, ancho_actual_GHz, curva_actual in zip(potencias_microondas_comparacion, 
                                                           anchos_fwhm_comparacion_GHz, 
                                                           curvas_potencia_microondas):
    eje.plot(frecuencias_GHz, 
             curva_actual, 
             linewidth=1.8, 
             label=
             f"P = {potencia_actual:g}," 
             f"FWHM = {ancho_actual_GHz * 1000:.2f} MHz"
             )

eje.set_xlabel("Frecuencia de microondas (GHz)")
eje.set_ylabel("Fluorescencia normalizada")
eje.set_title(f"Ensanchamiento ODMR con la potencia de microondas\n" 
              f"B = {campo_mT:.3f} mT,"
              f"θ = {angulo_grados:.1f}\u00B0,"
              f"T = {temperatura_C:.1f} \u00B0C"
              )

eje.set_xlim(frecuencia_minima_GHz, frecuencia_maxima_GHz)
eje.set_ylim(0.88, 1.01)

eje.grid(True, alpha=0.3)
eje.legend(loc="best")

guardar_grafica("odmr_vs_microwave_power.png")

# ==========================================================================
# GRÁFICA 9: SENSIBILIDAD FRENTE A LA TASA DE FOTONES
# ==========================================================================

figura_9, eje = plt.subplots(figsize=(8, 5))

eje.plot(tasas_fotones_sensibilidades_Hz, 
         sensibilidad_tasa_nT_raiz_Hz, 
         linewidth=2.0, 
         color="tab:blue", 
         label="Sensibilidad magnética")

eje.scatter(tasa_fotones_Hz, 
            sensibilidad_nT_raiz_Hz, 
            color="tab:red", 
            s=55, 
            zorder=3, 
            label=
            "Valor seleccionado = "
            f"{sensibilidad_nT_raiz_Hz:.2f} nT/\u221AHz"
            )

eje.set_xscale("log")
eje.set_yscale("log")

eje.set_xlabel("Tasa de fotones detectados (fotones/s)")
eje.set_ylabel("Sensibilidad magnética (nT/\u221AHz)")
eje.set_title(f"Sensibilidad magnética limitada por ruido fotónico\n" 
              f"FWHM = {ancho_fwhm_actual_GHz * 1000:.2f} MHz," 
              f"C = {contraste_por_resonancia * 100:.1f} %")

eje.grid(True, which="both", alpha=0.3)
eje.legend(loc="best")

guardar_grafica("magnetic_sensitivity_vs_photon_rate.png")

# ==========================================================================
# GRÁFICA 10: DEPENDENCIA DE LAS RESONANCIAS CON LA PERTURBACIÓN TRANSVERSAL
# ==========================================================================

figura_10, eje = plt.subplots(figsize=(8, 5))

eje.plot(perturbaciones_E_MHz, 
         frecuencias_inferiores_E_GHz, 
         linewidth=2.0, 
         color="tab:blue", 
         label="Resonancia inferior"
         )
eje.plot(perturbaciones_E_MHz, 
         frecuencias_superiores_E_GHz, 
         linewidth=2.0, 
         color="tab:orange", 
         label="Resonancia superior"
         )

eje.axvline(perturbacion_E_MHz, 
            linestyle="--", 
            linewidth=1.5, 
            color="tab:red", 
            alpha=0.8, 
            label=f"E seleccionada = {perturbacion_E_MHz:.3f} MHz"
            )

eje.set_xlabel("Perturbación transversal E (MHz)")
eje.set_ylabel("Frecuencia de resonancia (GHz)")

eje.set_title("Dependencia de las resonancias ODMR con la perturbación transversal (sin campo magnético)\n"
              f"T = {temperatura_C:.1f} \u00B0C"
              )

eje.set_xlim(perturbaciones_E_MHz[0], 
             perturbaciones_E_MHz[-1]
             )

margen_frecuencia_E_GHz = 0.005
eje.set_ylim(np.min(frecuencias_inferiores_E_GHz) - margen_frecuencia_E_GHz, 
             np.max(frecuencias_superiores_E_GHz) + margen_frecuencia_E_GHz)

eje.grid(alpha=0.3)
eje.legend()

guardar_grafica("frequencies_vs_transverse_perturbation.png")

# ==========================================================================
# ESTIMACIÓN DEL CAMPO A PARTIR DEL DESDOBLAMIENTO
# ==========================================================================

separacion_resonancias_GHz = abs(frecuencia_superior_GHz - frecuencia_inferior_GHz)

campo_longitudinal_estimado_T  = separacion_resonancias_GHz / (2 * gamma_e_GHz_T)
campo_longitudinal_estimado_mT = campo_longitudinal_estimado_T * 1000

campo_longitudinal_real_mT = abs(campo_z_T) * 1000
error_estimacion_mT = campo_longitudinal_estimado_mT - campo_longitudinal_real_mT

if campo_longitudinal_real_mT > 0:
    error_estimacion_porcentual = 100 * error_estimacion_mT / campo_longitudinal_real_mT
else:
    error_estimacion_porcentual = 0.0

# ==========================================================================
# ARCHIVO DE RESULTADOS
# ==========================================================================

ruta_txt = carpeta_resultados / "results_odmr.txt"

linea = "=" * 78
sublinea = "-" * 78

tolerancia_degeneracion_GHz = ancho_fwhm_actual_GHz / 100

resonancias_distintas_cuatro_NV_GHz = obtener_resonancias_unicas(
    resonancias_cuatro_NV_GHz,
    tolerancia_degeneracion_GHz)

with open(ruta_txt, "w", encoding="utf-8") as archivo:
    # ======================================================================
    # CABECERA 
    # ======================================================================
    archivo.write(linea + "\n")
    archivo.write(" SIMULACIÓN ODMR DE CENTROS NV EN DIAMANTE\n")
    archivo.write(" Sofía Núñez de Andrés\n")
    archivo.write(" Prácticas externas: CINN (Julio - Agosto 2026)\n")
    
    archivo.write(linea + "\n\n")

    # ======================================================================
    # 1. PARÁMETROS FÍSICOS DEL CENTRO NV 
    # ======================================================================

    archivo.write("1. PARÁMETROS FÍSICOS DEL CENTRO NV\n")
    archivo.write(sublinea + "\n")
    archivo.write(f"División de campo cero de referencia = {D_GHz:.9f} GHz\n")
    archivo.write(f"Temperatura de referencia            = {temperatura_referencia_C:.3f} \u00B0C\n")
    archivo.write(f"Coeficiente térmico de D             = {coeficiente_temperatura_D_GHz_C:.9e} GHz/\u00B0C\n")
    archivo.write(f"Temperatura seleccionada             = {temperatura_C:.3f} \u00B0C\n")
    archivo.write(f"División de campo cero a T, D(T)     = {D_actual_GHz:.9f} GHz\n")
    archivo.write(f"Desplazamiento de D respecto a 25 \u00B0C = {desplazamiento_D_MHz:.6f} MHz\n")
    archivo.write(f"Relación giromagnética, \u03B3\u2091           = {gamma_e_GHz_T:.6f} GHz/T\n")
    archivo.write(f"Perturbación transversal, E          = {perturbacion_E_MHz:.6f} MHz\n")

    # ======================================================================
    # 2. PARÁMETROS DEL ESPECTRO ODMR 
    # ======================================================================

    archivo.write("2. PARÁMETROS DEL ESPECTRO ODMR\n")
    archivo.write(sublinea + "\n")
    archivo.write(f"Número de puntos del barrido    = {numero_puntos}\n")
    archivo.write(f"Anchura FWHM de las resonancias = {ancho_fwhm_GHz * 1000:.6f} MHz\n")
    archivo.write(f"Contraste por resonancia        = {contraste_por_resonancia * 100:.2f} %\n")
    archivo.write(f"Margen del barrido              = {margen_frecuencia_GHz:.6f} GHz\n")
    archivo.write(f"Frecuencia mínima del barrido   = {frecuencia_minima_GHz:.9f} GHz\n")
    archivo.write(f"Frecuencia máxima del barrido   = {frecuencia_maxima_GHz:.9f} GHz\n\n")

    # ======================================================================
    # 3. CAMPO MAGNÉTICO EN EL MODELO DE UN ÚNICO EJE NV 
    # ======================================================================

    archivo.write("3. CAMPO MAGNÉTICO EN EL MODELO DE UN ÚNICO EJE NV\n")
    archivo.write(sublinea + "\n")
    archivo.write(f"Módulo del campo              = {campo_mT:.6f} mT\n")
    archivo.write(f"Ángulo respecto al eje NV     = {angulo_grados:.3f}\u00B0\n")
    archivo.write(f"Componente Bx                 = {campo_x_T:.9f} T\n")
    archivo.write(f"Componente By                 = {campo_y_T:.9f} T\n")
    archivo.write(f"Componente Bz                 = {campo_z_T:.9f} T\n")
    archivo.write(f"Componente longitudinal, |Bz| = {campo_longitudinal_real_mT:.6f} mT\n\n")

    # ======================================================================
    # 4. RESULTADOS SIN CAMPO MAGNÉTICO
    # ======================================================================

    archivo.write("4. RESULTADOS SIN CAMPO MAGNÉTICO\n")
    archivo.write(sublinea + "\n")
    archivo.write("Autovalores de H/h = ")
    archivo.write(", ".join(f"{energia:.9f}" for energia in energias_sin_campo_GHz) + " GHz\n")
    archivo.write(f"Frecuencia de resonancia inferior      = {frecuencia_inferior_sin_campo_GHz:.9f} GHz\n")
    archivo.write(f"Frecuencia de resonancia superior      = {frecuencia_superior_sin_campo_GHz:.9f} GHz\n")
    archivo.write(f"Separación entre resonancias           = {(frecuencia_superior_sin_campo_GHz - 
                                                               frecuencia_inferior_sin_campo_GHz) * 1000:.6f} MHz\n")
    archivo.write(f"División de campo cero utilizada, D(T) = {D_actual_GHz:.9f} GHz\n\n")

    # ======================================================================
    # 5. RESULTADOS CON CAMPO MAGNÉTICO: UN ÚNICO EJE NV
    # ======================================================================

    archivo.write("5. RESULTADOS CON CAMPO MAGNÉTICO: UN ÚNICO EJE NV\n")
    archivo.write(sublinea + "\n")
    archivo.write("Autovalores de H/h = ")
    archivo.write(", ".join(f"{energia:.9f}" for energia in energias_con_campo_GHz) + " GHz\n")
    archivo.write(f"Frecuencia de resonancia inferior = {frecuencia_inferior_GHz:.9f} GHz\n")
    archivo.write(f"Frecuencia de resonancia superior = {frecuencia_superior_GHz:.9f} GHz\n")
    archivo.write(f"Separación entre resonancias      = {separacion_resonancias_GHz * 1000:.6f} MHz\n\n")

    # ======================================================================
    # 6. ESTIMACIÓN DEL CAMPO MEDIANTE EL DESDOBLAMIENTO
    # ======================================================================

    archivo.write("6. ESTIMACIÓN DEL CAMPO MEDIANTE EL DESDOBLAMIENTO\n")
    archivo.write(sublinea + "\n")
    archivo.write(f"Componente longitudinal real, |Bz| = {campo_longitudinal_real_mT:.6f} mT\n")
    archivo.write(f"Componente longitudinal estimada   = {campo_longitudinal_estimado_mT:.6f} mT\n")
    archivo.write(f"Error de la estimación             = {error_estimacion_mT:.6f} mT\n")
    archivo.write(f"Error porcentual                   = {error_estimacion_porcentual:.3f} %\n")
    archivo.write("Nota: la expresión utilizada es una aproximación lineal y estima\n")
    archivo.write("principalmente la componente del campo paralela al eje NV.\n\n")

    # ======================================================================
    # 7. CAMPO MAGNÉTICO EN EL SISTEMA DEL CRISTAL
    # ======================================================================

    archivo.write("7. CAMPO MAGNÉTICO EN EL SISTEMA DEL CRISTAL\n")
    archivo.write(sublinea + "\n")
    archivo.write(f"Módulo del campo     = {campo_mT:.6f} mT\n")
    archivo.write(f"Ángulo polar \u03B8       = {angulo_polar_cristal_grados:.3f}\u00B0\n")
    archivo.write(f"Ángulo azimutal \u03C6    = {angulo_azimutal_cristal_grados:.3f}\u00B0\n")
    archivo.write(f"Componente global Bx = {vector_campo_cristal_T[0]:.9f} T\n")
    archivo.write(f"Componente global By = {vector_campo_cristal_T[1]:.9f} T\n")
    archivo.write(f"Componente global Bz = {vector_campo_cristal_T[2]:.9f} T\n\n")

    # ======================================================================
    # 8. RESULTADOS DE LAS CUATRO ORIENTACIONES CRISTALOGRÁFICAS
    # ======================================================================

    archivo.write("8. RESULTADOS DE LAS CUATRO ORIENTACIONES CRISTALOGRÁFICAS\n")
    archivo.write(sublinea + "\n")
    archivo.write(f"Número total de transiciones calculadas = {len(resonancias_cuatro_NV_GHz)}\n")
    archivo.write(f"Número de resonancias distintas         = {len(resonancias_distintas_cuatro_NV_GHz)}\n")
    archivo.write(f"Frecuencia mínima del espectro conjunto = {frecuencia_minima_cuatro_NV_GHz:.9f} GHz\n")
    archivo.write(f"Frecuencia máxima del espectro conjunto = {frecuencia_maxima_cuatro_NV_GHz:.9f} GHz\n\n")

    for numero_NV, (nombre_NV, resultados_NV) in enumerate(resultados_cuatro_NV.items(), start=1):
        archivo.write(f"8.{numero_NV}. {nombre_NV}\n")
        archivo.write(f"Componente local Bx        = {resultados_NV['campo_x_local_T']:.9f} T\n")
        archivo.write(f"Componente local By        = {resultados_NV['campo_y_local_T']:.9f} T\n")
        archivo.write(f"Componente local Bz        = {resultados_NV['campo_z_local_T']:.9f} T\n")
        archivo.write(f"Proyección longitudinal Bz = {resultados_NV['campo_z_local_T'] * 1000:.6f} mT\n")
        archivo.write("Autovalores de H/h = ")
        archivo.write(", ".join(f"{energia:.9f}" for energia in resultados_NV["energias_GHz"]) + " GHz\n")
        archivo.write(f"Frecuencia de resonancia inferior = {resultados_NV['frecuencia_inferior_GHz']:.9f} GHz\n")
        archivo.write(f"Frecuencia de resonancia superior = {resultados_NV['frecuencia_superior_GHz']:.9f} GHz\n")
        archivo.write(f"Separación entre resonancias      = {abs(resultados_NV['frecuencia_superior_GHz'] - resultados_NV['frecuencia_inferior_GHz']) * 1000:.6f} MHz\n\n")

    # ======================================================================
    # 9. RESONANCIAS DISTINTAS DEL ESPECTRO CONJUNTO
    # ======================================================================

    archivo.write("9. RESONANCIAS DISTINTAS DEL ESPECTRO CONJUNTO\n")
    archivo.write(sublinea + "\n")

    for indice, frecuencia_GHz in enumerate(resonancias_distintas_cuatro_NV_GHz, start=1):
        archivo.write(f"Resonancia distinta {indice} = {frecuencia_GHz:.9f} GHz\n")

    archivo.write("\n")
    archivo.write("Nota: las cuatro familias NV pueden producir frecuencias coincidentes\n")
    archivo.write("cuando presentan proyecciones equivalentes del campo magnético.\n")
    archivo.write("En esos casos se calculan ocho transiciones, pero el número de líneas\n")
    archivo.write("distintas observables en el espectro puede ser inferior.\n\n")

    # ======================================================================
    # 10. DEPENDENCIA CON LA TEMPERATURA
    # ======================================================================

    archivo.write("10. DEPENDENCIA CON LA TEMPERATURA\n")
    archivo.write(sublinea + "\n")
    archivo.write(f"Temperatura seleccionada             = {temperatura_C:.3f} \u00B0C\n")
    archivo.write(f"Temperatura mínima del barrido       = {temperatura_minima_C:.3f} \u00B0C\n")
    archivo.write(f"Temperatura máxima del barrido       = {temperatura_maxima_C:.3f} \u00B0C\n")
    archivo.write(f"División de campo cero de referencia = {D_GHz:.9f} GHz\n")
    archivo.write(f"División de campo cero D(T)          = {D_actual_GHz:.9f} GHz\n")
    archivo.write(f"Desplazamiento térmico de D          = {desplazamiento_D_MHz:.9f} MHz\n")
    archivo.write(f"Coeficiente térmico utilizado        = {coeficiente_temperatura_D_GHz_C:.9e} GHz/\u00B0C\n")
    archivo.write(f"Resonancia inferior a T seleccionada = {frecuencia_inferior_GHz:.9f} GHz\n")
    archivo.write(f"Resonancia superior a T seleccionada = {frecuencia_superior_GHz:.9f} GHz\n\n")

    archivo.write("Nota: este módulo utiliza una aproximación lineal de D con la\n")
    archivo.write("temperatura alrededor de la temperatura de referencia.\n\n")

    # ======================================================================
    # 11. RUIDO EXPERIMENTAL
    # ======================================================================

    archivo.write("11. RUIDO EXPERIMENTAL\n")
    archivo.write(sublinea + "\n")
    archivo.write(f"Tasa de fotones detectados        = {tasa_fotones_Hz:.3f} fotones/s\n")
    archivo.write(f"Tiempo de integración por punto   = {tiempo_integracion_s:.6f} s\n")
    archivo.write(f"Número medio de fotones por punto = {numero_fotones_por_punto:.3f}\n")
    archivo.write("Ruido fotónico relativo           = " f"{ruido_fotonico_relativo:.6f} " f"({ruido_fotonico_relativo * 100:.3f} %)\n")
    archivo.write("Ruido técnico relativo            = " f"{ruido_tecnico_relativo:.6f} " f"({ruido_tecnico_relativo * 100:.3f} %)\n")
    archivo.write("Ruido total relativo aproximado   = " f"{ruido_total_relativo_aproximado:.6f} " f"({ruido_total_relativo_aproximado * 100:.3f} %)\n")
    archivo.write(f"Semilla aleatoria                 = {semilla_aleatoria}\n\n")
    archivo.write("Modelo utilizado:\n")
    archivo.write("- Ruido fotónico mediante distribución de Poisson.\n")
    archivo.write("- Ruido técnico mediante distribución normal.\n")
    archivo.write("- La semilla fija permite reproducir exactamente la simulación.\n\n")

    # ======================================================================
    # 12. POTENCIA DE MICROONDAS Y ENSANCHAMIENTO
    # ======================================================================

    archivo.write("12. POTENCIA DE MICROONDAS Y ENSANCHAMIENTO\n")
    archivo.write(sublinea + "\n")
    archivo.write(f"Potencia relativa seleccionada  = {potencia_microondas_relativa:.6f}\n")
    archivo.write(f"Potencia relativa de referencia = {potencia_microondas_referencia:.6f}\n")
    archivo.write(f"Anchura FWHM de referencia      = {ancho_fwhm_GHz * 1000:.6f} MHz\n")
    archivo.write(f"Anchura FWHM efectiva           = {ancho_fwhm_actual_GHz * 1000:.6f} MHz\n")
    archivo.write(f"Factor de ensanchamiento        = {factor_ensanchamiento:.6f}\n")
    archivo.write("Modelo utilizado: FWHM = FWHM\u2080 \u00B7 \u221A(P_relativa)\n\n")

    archivo.write("Comparación de potencias:\n")

    for potencia_actual, ancho_actual_GHz in zip(potencias_microondas_comparacion, anchos_fwhm_comparacion_GHz):
        archivo.write(f"Potencia relativa = {potencia_actual:.2f} | FWHM = {ancho_actual_GHz * 1000:.6f} MHz\n")

    archivo.write("\n")
    archivo.write("Nota: en este modelo educativo la potencia de microondas modifica\n")
    archivo.write("la anchura de las resonancias, pero no sus frecuencias centrales.\n\n")

    # ======================================================================
    # 13. SENSIBILIDAD MAGNÉTICA
    # ======================================================================

    archivo.write("13. SENSIBILIDAD MAGNÉTICA\n")
    archivo.write(sublinea + "\n")
    archivo.write(f"Factor del perfil Lorentziano      = {factor_perfil_lorentziano:.6f}\n")
    archivo.write(f"Anchura FWHM utilizada             = {ancho_fwhm_actual_GHz * 1000:.6f} MHz\n")
    archivo.write(f"Contraste utilizado                = {contraste_por_resonancia:.6f}\n")
    archivo.write(f"Tasa de fotones utilizada          = {tasa_fotones_Hz:.3f} fotones/s\n")
    archivo.write(f"Sensibilidad magnética             = {sensibilidad_nT_raiz_Hz:.6f} nT/\u221AHz\n")
    archivo.write(f"Tiempo de integración seleccionado = {tiempo_integracion_s:.6f} s\n")
    archivo.write(f"Campo mínimo detectable            = {campo_minimo_detectable_uT:.6f} \u00B5T\n\n")

    archivo.write("Modelo utilizado:\n")
    archivo.write("\u03B7_B = k \u00B7 FWHM / (\u03B3\u2091 \u00B7 C \u00B7 \u221AR)\n")
    archivo.write("\u03B4B  = \u03B7_B / \u221At\n\n")

    archivo.write("Nota: esta es una estimación limitada por ruido fotónico.\n")
    archivo.write("No incluye todas las fuentes de ruido ni imperfecciones de un\n")
    archivo.write("experimento real.\n\n")

    # ======================================================================
    # 14. EFECTO DE LA PERTURBACIÓN TRANSVERSAL
    # ======================================================================

    archivo.write("14. EFECTO DE LA PERTURBACIÓN TRANSVERSAL\n")
    archivo.write(sublinea + "\n")
    archivo.write(f"Perturbación transversal E        = {perturbacion_E_GHz * 1000:.6f} MHz\n")
    archivo.write(f"Frecuencia de resonancia inferior = {frecuencia_inferior_sin_campo_GHz:.9f} GHz\n")
    archivo.write(f"Frecuencia de resonancia superior = {frecuencia_superior_sin_campo_GHz:.9f} GHz\n")
    archivo.write(f"Separación entre resonancias      = {(frecuencia_superior_sin_campo_GHz - frecuencia_inferior_sin_campo_GHz) * 1000:.6f} MHz\n\n")

    archivo.write("Modelo físico:\n")
    archivo.write("La perturbación transversal E rompe la degeneración de los\n")
    archivo.write("estados m\u209B = \u00B11 incluso en ausencia de campo magnético.\n")
    archivo.write("Las frecuencias de resonancia se obtienen diagonalizando el Hamiltoniano\n")
    archivo.write("el Hamiltoniano para cada valor de E.\n")

    # ======================================================================
    # 15. ARCHIVOS GENERADOS
    # ======================================================================

    archivo.write("15. ARCHIVOS GENERADOS\n")
    archivo.write(sublinea + "\n")
    archivo.write("1.  odmr_without_magnetic_field.png\n")
    archivo.write("2.  odmr_with_magnetic_field.png\n")
    archivo.write("3.  frequencies_vs_field.png\n")
    archivo.write("4.  frequencies_vs_angle.png\n")
    archivo.write("5.  evolution_odmr_with_magnetic_field.png\n")
    archivo.write("6.  odmr_four_nv_orientations.png\n")
    archivo.write("7.  frequencies_vs_temperature.png\n")    
    archivo.write("8.  odmr_vs_microwave_power.png\n")
    archivo.write("9.  magnetic_sensitivity_vs_photon_rate.png\n")
    archivo.write("10. frequencies_vs_transverse_perturbation.png\n")
    archivo.write("11. results_odmr.txt\n\n")

    archivo.write(linea + "\n")
    archivo.write("Fin del informe.\n")
    archivo.write(linea + "\n")

# ==========================================================================
# RESULTADOS MOSTRADOS EN PANTALLA
# ==========================================================================

print("\n==========================================================================")
print(" RESULTADOS DE LA SIMULACIÓN ODMR")
print("==========================================================================")

print(f"\nResultados guardados en: {carpeta_resultados.name}")

print("\nPARÁMETROS UTILIZADOS:")
print(f"\n · Campo magnético = {campo_mT:.3f} mT")
print(f"\n · Ángulo respecto al eje NV = {angulo_grados:.3f}°")
print(f"\n · Temperatura = {temperatura_C:.3f} °C")
print(f"\n · Perturbación transversal E = {perturbacion_E_MHz:.3f} MHz")
print(f"\n · Potencia relativa = {potencia_microondas_relativa:.2f}")
print(f"\n · Tasa de fotones = {tasa_fotones_Hz:,.0f} fotones/s")
print(f"\n · Tiempo de integración = {tiempo_integracion_s:.3f} s")

print("\nARCHIVOS GENERADOS:")
print("\n · 10 gráficas")
print("\n · results_odmr.txt")

print("\nConsulte 'results_odmr.txt' para obtener el informe completo.")

print("\n==========================================================================")
print(" -> Simulación finalizada correctamente. ")
print("==========================================================================")