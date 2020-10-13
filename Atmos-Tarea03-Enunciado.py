# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Astrofísica Estelar – Tarea 03
#
# Usar la ecuación de Milne para calcular el espectro de flujo radiativo a diferentes profundidades en un modelo de una atmósfera ETL en equilibrio radiativo y con opacidad gris.

# ## Modelo de la atmósfera gris
#
# ### Propósito
#
# * Encontrar el espectro de flujo $F_\nu$ a diferentes profundidades en la atmósfera.
# * Esencialmente es reproducir la Figura 17.1 del libro de Hubeny & Mihalas (q.v.)
#
# ### Suposiciones
#
# * El hipótesis fundamental "gris" es que la opacidad es independiente de la frecuencia
#     - Entonces, la misma escala de profundidad óptica, $\tau$, aplica a todas frecuencias.
# * Suponemos equilibrio radiativo, de modo que la intensidad promedia integrada en frecuencia es igual a la función fuente: $ J = S $.
#     * Con un supuesto adicional de geometría plano-paralelo, esto significa que el flujo radiativo integrado en frecuencia, $ H $, es constante con la profundidad.
# * También asumimos *Equilibrio termodinámico local* y sin dispersión, de modo que la función fuente a cualquier frecuencia viene dada por la función de Planck a la temperatura local: $ S_\nu(\tau) = B_\nu[T(\tau)] $
#     - Donde $$B_\nu = \frac{2h \nu^3}{c^2} \frac{1}{e^{h\nu/kT} - 1}$$
#     - Además, la función Planck integrada en frecuencia es $$B = (\sigma/\pi) T^4$$
# * Además, usamos la *aproximación de Eddington*, $ J = 3 K $, como una relación de cierre para los momentos de la ecuación de transferencia radiativa integrada en frecuencia.

# ### Resultados del modelo
#
# * Flujo bolométrico es constante con profundidad: $H = \sigma T_\mathrm{eff}^4 \,/\, 4 \pi$
#     - La temperatura efectiva $T_\mathrm{eff}$ es el parámetro global único que caracteriza la atmósfera
# * Función fuente, función de Planck, e intensidad media (todas integrada en frecuencia): $S = B = J = 3 H (\tau + \frac{2}{3})$
#     - $\Rightarrow$ $B\big[\tau = \frac23\big] = B(T_\mathrm{eff}) = 4 H$

# ### Versión adimensional de las ecuaciones
#
# Para hacer los cálculos más limpios y generales, usamos la temperatura efectiva para normalizarlas:
#
# * Temperatura adimensional: $p(\tau) \equiv T(\tau) / T_\mathrm{eff}$
# * Frecuencia admensional: $\alpha = h\nu/k T_\mathrm{eff}$

# ### <font color=#c31> Ejercicio 1: Escribir una función para $p(\tau)$</font>
#
# Escriba una función:
# ```python
# def p(tau):
#     """T / T_eff en función de la profundidad óptica `tau`"""
#     return ????????
# ```
# y graficarla para el rango $\tau = [0, 4]$.  Checar que $p(2/3) = 1.0$ y que $p(0) \approx 0.841$.

# ### <font color=#c31> Ejercicio 2: Función Planck adimensional</font>
#
# La función Planck por unidad de frecuencia adimensional es $B_\alpha = 
# \frac{d\nu}{d\alpha} B_\nu$.  Muestre que se puede escribir:
# $$
# \frac{B_\alpha}{B(T_\mathrm{eff})} = \frac{C\, \alpha^3}{e^{\alpha / p(\tau)} - 1}
# $$
# donde el constante $C$ es un número puro (combinación de constantes físicos y matemáticos).  Encuentre el valor de $C$. 
#
# Una manera de hacerlo usando la librería `astropy.constants` es así:
#
# ```python
# import numpy as np
# from astropy.constants import c, k_B, h, sigma_sb
#
# bigC = ?????
# bigC_unit = bigC.unit.compose()
# bigC = bigC.value
# print('C =', bigC)
# print("Unidades de C:", bigC_unit)
# ```
#
# Note que los constantes de astropy son objetos especiales que tienen unidades físicas incorporados.  El renglón con `bigC.value` es para convertirla a un `float` normal.  Imprimimos también las unidades para checar que el resultado fue adimensional.

# ### <font color=#c31> Ejercicio 3: Graficar la función Planck para diferentes profundidades</font>
#
# Escriba una función para $B_\alpha / B(T_\mathrm{eff})$. Por ejemplo:
#
# ```python
# def planck(alpha, tau):
#     """Función de Planck normalizada (argumentos: frecuencia `alpha` y profundidad `tau`)"""
#     return bigC * ??????
# ```
#
# Grafique la función en cuatro diferentes maneras:
#
# #### <font color=#c31>(i) Graficar $B_\alpha / B(T_\mathrm{eff})$ contra frecuencia $\alpha$ para diferentes profundidades</font>
# Por ejemplo, $\tau$ = 0, 1, 2, 4, y 8. Usar un rango de frecuencia de $\alpha = [0, 12]$.
#
# * Checar que tanto el valor pico de $B_\alpha$ aumenta con $\tau$, así como el valor de la $\alpha$ donde el pico ocurre. Las curvas no deben de cruzarse.
#     
# #### <font color=#c31>(ii) Igual a la (i), pero graficando $B_\alpha / B(T(\tau))$.</font> 
# Es decir, normalizar por la función Planck integrada a la temperatura local de cada profundidad.
#
# - Ahora el area abajo de cada curva debe de ser igual. 
#     
# #### <font color=#c31>(iii) Graficar $B_\alpha / B(T_\mathrm{eff})$ contra $\tau$ para 3 frecencias particulares</font>
# Usar $\alpha$ = 1, 3, y 9 y un rango de $\tau = [0, 20]$.
#     
# - ¿por qué escoger estas 3 frecuencias?
# - Si $T_\mathrm{eff} = 6000\,\mathrm{K}$, a cuáles longitudes de onda corresponden?
#     
# #### <font color=#c31>(iv) Igual a la (iii), pero normalizando cada curva por el valor a $\tau = 2/3$.</font>
# Usar un rango más corto: $\tau = [0, 2]$.
#
# - Esto facilita la comparación del gradiente de la función fuente entre las 3 frecuencias
# - Explique cualitativamente qué está pasando en esta gráfica
#     
#

# ## La ecuación Milne para encontrar el flujo
#
# La ecuación general para el flujo es:
# $$
# H_\nu(\tau) = \frac12 \left[
# \int_{\tau}^{\infty} 
# S_\nu(t)\, E_2(t - \tau) \, dt
# \; -
# \int_0^{\tau}
# S_\nu(t)\, E_2(\tau - t) \, dt
# \right]
# $$
#
# La versión adimensional en ETL es
# $$
# H_{\alpha}(\tau) = \frac12 \left[
# \int_{\tau}^{\infty} 
# B_\alpha(t)\, E_2(t - \tau) \, dt
# \; -
# \int_0^{\tau}
# B_\alpha(t)\, E_2(\tau - t) \, dt
# \right]
# $$
# en donde $B_\alpha = 4 H C \alpha^3 / (e^{\alpha / p(\tau)} - 1)$ (ver arriba) y $E_2$ es la integral exponencial de orden dos. 
#
# Por lo tanto:
# $$
# \frac{H_{\alpha}(\tau)}{H} = 2 C \alpha^3
# \left[ \int_{\tau}^{\infty} 
# \frac{E_2(t - \tau) \, dt }{e^{\alpha / p(t)} - 1} 
# \; - 
# \int_0^{\tau} 
# \frac{E_2(\tau - t) \, dt }{e^{\alpha / p(t)} - 1}
# \right]
# = 
# \int_{\tau}^{\infty}
# 2 \frac{B_\alpha(t)}{B(T_\mathrm{eff})}\,
# E_2(t - \tau) \, dt
# \; -
# \int_0^{\tau}
# 2 \frac{B_\alpha(t)}{B(T_\mathrm{eff})}\,
# E_2(\tau - t) \, dt
# $$
#
#
#

# ### <font color=#c31> Ejercicio 4: Escribir una función para el integrando de la ecuación Milne</font>
#
# Se puede usar la función `scipy.special.expn` para la integral exponencial.  Además, el argumento de $E_2$ es siempre positivo en ambas integrales, entonces se puede usar el valor absoluto $|t - \tau|$, calculada con la función `abs`.
#
# ```python
# from scipy.special import expn
#
# def milne_integrand(t, alpha, tau):
#     """Dos veces B(t) por E_2(|t - tau|)"""
#     return ?????? 
# ```
#
# Note que el variable de integración va primero en esta función. Esto es necesario para que podemos usarlo con rutinas de integración numérica.

# ### Realizar las integrales
#
# * Escribimos dos funciones: una para la radiación en dirección descendente, otra para la radiación en dirección ascendente.  
# * Usamos la rutina `scipy.integrate.quad` que require al menos tres argumentos: la función para integrar, luego el límite inferior, luego el límite superior.  
#     - Se permite usar límites de integración finitas o infinitas. 
#     - Los argumentos auxiliares para el integrando están comunicados usando el argumento opcional `args`. 

# +
from scipy.integrate import quad

def downward(alpha, tau, integrand=milne_integrand):
    """Integrate the `integrand` between 0 and `tau` using quadpack"""
    result, error = quad(integrand, 0.0, tau, args=(alpha, tau))
    return result

def upward(alpha, tau, integrand=milne_integrand):
    """Integrate the `integrand` between `tau` and infinity using quadpack"""
    result, error = quad(integrand, tau, np.infty, args=(alpha, tau))
    return result


# -

# Ahora definimos una función `flux` que calcula el flujo neto como la diferencia entre el corriente ascendente y el descendente.
#
# Note que `quad` requiere que `alpha` y `tau` sean escalares.  Entonces, usamos la función `np.vectorize` como *decorador* de `flux` para que se pueda llamar con vectores (arreglos) para `alpha` y `tau`.  

@np.vectorize
def flux(alpha, tau):
    """Find net radiative flux as difference between 
    upward and downward streams"""
    return upward(alpha, tau) - downward(alpha, tau)


# ### <font color=#c31> Ejercicio 5: Graficar el espectro del flujo</font>
#
# * Utilice los mismos valores de $\tau$ y rango de $\alpha$ que usó en el 3(i) y 3(ii).  
# * Compare con la $B_\alpha / B(T(\tau))$ del 3(ii). 
#     - ¿Cómo es diferente el espectro del flujo a la función Planck a la misma profundidad?  ¿Por qué?
#

# ### <font color=#c31> [OPCIONAL] Ejercicio 6: Repitir para la intensidad promedia</font>
#
# * Escriba una función `schwarz_integrand` para el integrando en la ecuación de Schwarzschild para calcular $J_\alpha(\tau)$.
#     - Esta se utiliza en la función `meanJ` (proporcionado abajo) que evalua $J_\alpha(\tau)/H$
# * Grafique $J_\alpha(\tau)/H$ igual que en el ejercicio 5
# * Compare con la $B_\alpha / B(T_\mathrm{eff})$ del 3(i).  **OJO**: No es la misma comparación que en el ejercicio 5.
#     - ¿Qué pasa para $\tau > 2$?  ¿Por qué?
#

@np.vectorize
def meanJ(alpha, tau):
    """Mean intensity as sum of upward and downward streams"""
    return (upward(alpha, tau, integrand=schwarz_integrand) 
            + downward(alpha, tau, integrand=schwarz_integrand))


