import numpy as np
from numpy import exp, pi
from scipy.special import expn
from scipy.integrate import romberg, quad
import matplotlib.pyplot as plt

def p(tau):
    return (4.0 / (3.0*tau + 2.0))**0.25

CONSTANT = 0.308
CONSTANT = 4.0 * pi * 1.3806503e-16**4 / (
    6.62606876e-27**3 * 2.99792458e10**2 * 5.6703e-5)
def planck(alpha, tau):
    return (0.5*CONSTANT/pi) * alpha**3 * p(tau)**4 / (exp(alpha*p(tau)) - 1.0)

def planck_H(alpha, tau):
    return (2.0*CONSTANT/pi) * alpha**3 / (exp(alpha*p(tau)) - 1.0)

def milne_integrand(t, alpha, tau):
    return expn(2, abs(t - tau)) / (exp(alpha*p(tau)) - 1.0)

def downward(alpha, tau, integrand=milne_integrand):
    result, error = quad(integrand, 0.0, tau, args=(alpha, tau))
    return result

def upward(alpha, tau, integrand=milne_integrand):
    result, error = quad(integrand, tau, np.infty, args=(alpha, tau))
    return result

def flux(alpha, tau):
    result = upward(alpha, tau)
    if tau > 0.0:
        result -= downward(alpha, tau)
    result *= CONSTANT * alpha**3
    return result

def schwarz_integrand(t, alpha, tau):
    return expn(1, abs(t - tau)) / (exp(alpha*p(tau)) - 1.0)

def meanJ(alpha, tau):
    result = upward(alpha, tau, integrand=schwarz_integrand)
    if tau > 0.0:
        # Plus sign here, unlike for the flux 
        result += downward(alpha, tau, integrand=schwarz_integrand)
    result *= CONSTANT * alpha**3
    return result

def find_fluxes(alphas, tau):
    fluxes = []
    for alpha in alphas:
        fluxes.append(flux(alpha, tau))
    return np.array(fluxes)

def make_graph():
    alpha_pts = np.linspace(0.0, 12.0, 200)
    for tau, color in zip([0.0, 1.0, 2.0, 4.0, 8.0], "bgrcm"):
        T = 1./p(tau)
        flux_pts = find_fluxes(alpha_pts, tau)
        m = np.isfinite(flux_pts)
        flux_bolo = np.trapz(flux_pts[m], alpha_pts[m])
        plank_bolo = np.trapz(planck(alpha_pts[m], tau), alpha_pts[m])
        print('tau = {:.1f} T = {:.2f} F = {:.5f} B = {:.5f}'.format(tau, T, flux_bolo, plank_bolo))
        plt.plot(alpha_pts, flux_pts, "-" + color, 
                 label="H_alpha/H, tau = {}, T = {:.2f} T_ef".format(int(tau), T))
        plt.plot(alpha_pts, planck(alpha_pts, tau), "--" + color)
    plt.plot(alpha_pts, planck_H(alpha_pts, 2./3.), ":k",
             label="B_alpha(T_ef)/H")
    plt.legend(fontsize="x-small")
    plt.xlabel("alpha")
    plt.ylabel("Flux")
    plt.savefig("fluxes.pdf")

if __name__ == "__main__":
    make_graph()
