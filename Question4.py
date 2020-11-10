import sympy as sym

# LINERESIATION

# phi(F, x3,x4) -phi(F0, x30, x40) +(dphi/dF)((F0, x30, x40)
# define phi
# define psi
# define all involved symbolic functions

#constants

M, m, g, ell = sym.symbols('M, m, g, ell', real=True, postive=True)

#  System variables
x1, x2, x3, x4, F = sym.symbols('x1, x2, x3, x4, F')

# Define  Φ (phi)
phi = 4 * m * ell *x4**2 * sym.sin(x3) + 4 * F - 3 *m*g*sym.sin(x3) *sym.cos(x3)
phi /= 4 * (M+m) -3*m*sym.cos(x3)**2

# Define ψ  (psi)
psi =  -3 *(m * ell *x4**2 * sym.sin(x3)*sym.cos(x3) + F*sym.cos(x3) -(M+m)*g*sym.sin(x3))
psi /= (4*(M+m)-3*m*sym.cos(x3)**2)*ell


# determine the partial derivatives of Φ (phi), wrt to F, x3, x4

phi_deriv_F = phi.diff(F)
phi_deriv_x3 = phi.diff(x3)
phi_deriv_x4 = phi.diff(x4)

# determine the partial derivatives of ψ  (psi), wrt to F, x3, x4

psi_deriv_F = psi.diff(F)
psi_deriv_x3 = psi.diff(x3)
psi_deriv_x4 = psi.diff(x4)

# set the initial conditions

F0 = 0     # initial condition of F = 0

x30 = 0    # initial condition of x3 = 0

x40 = 0    # initial condition of x4 = 0

def evaluate_at_equilibrium(f):     # evaluting equilibrium function

    return f.subs([(F, F0), (x3, x30), (x4, x40)])

# Computes the derivatives of Φ (phi) at equilibrium point

phi_deriv_F_at_equlibrium = evaluate_at_equilibrium(phi.diff(F))    # Φ (phi) derivative of F at equilibrium

phi_deriv_x3_at_equlibrium = evaluate_at_equilibrium(phi.diff(x3))  # Φ (phi) derivative of x3 at equilibrium

phi_deriv_x4_at_equlibrium = evaluate_at_equilibrium(phi.diff(x4))  # Φ (phi) derivative of x4 at equilibrium

# Computes the derivatives of ψ  (psi) at equilibrium point

psi_deriv_F_at_equlibrium = evaluate_at_equilibrium(psi.diff(F))    # ψ  (psi) derivative of F at equilibrium

psi_deriv_x3_at_equlibrium = evaluate_at_equilibrium(psi.diff(x3))  # ψ  (psi) derivative of x3 at equilibrium

psi_deriv_x4_at_equlibrium = evaluate_at_equilibrium(psi.diff(x4))  # ψ  (psi) derivative of x4 at equilibrium

# x2' = aF - bx3
a = phi_deriv_F_at_equlibrium
b = -phi_deriv_x3_at_equlibrium

#x4' = -cF + dx3

c = -psi_deriv_F_at_equlibrium
d = psi_deriv_x3_at_equlibrium


M_value = 0.3
m_value = 0.1
ell_value = 0.35
g_value = 9.81

a_value = float(a.subs([(M, M_value), (m, m_value)]))                   # these need to be in float to specify type
b_value = float(b.subs([(M, M_value), (m, m_value), (g, g_value)]))
c_value = float(c.subs([(M, M_value), (m, m_value), (g, g_value), (ell, ell_value)]))
d_value = float(d.subs([(M, M_value), (m, m_value), (g, g_value), (ell, ell_value)]))


# ----- ONLY NUMERICAL VALUES --------

import control as ctrl                  # works for transfer functions
import matplotlib.pyplot as plt         #
import numpy as np
n_points = 500
t_final =0.2
t_span =np.linspace(0, t_final, n_points)                       # array of time instants

Force_acting_on_pendulum=np.sin(100*np.square(t_span))          # input signal, which  is the force acting on pendulum


#--------------------Trajectories of theta(t)------------------------
G_theta = ctrl.TransferFunction([-c_value], [1, 0, -d_value])   # numerator = -c
                                                                # denominator 1s^2+0s -d which is [1,0,-d] in vector

t_theta_out, y_theta_out, x_theta_out = ctrl.forced_response(G_theta, t_span,Force_acting_on_pendulum)

plt.plot(t_theta_out, y_theta_out)
plt.show()

#--------------------Trajectories of x(t)------------------------

G_x =ctrl.TransferFunction([a_value,0,a_value*d_value+b_value*c_value],[1,0,-d_value,0,0])

t_x_out, y_x_out, x_x_out = ctrl.forced_response(G_x, t_span,Force_acting_on_pendulum)

plt.plot(t_x_out, y_x_out)
plt.show()


