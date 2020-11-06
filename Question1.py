import sympy as sym

# LINERESIATION

# phi(F, x3,x4) -phi(F0, x30, x40) +(dphi/dF)((F0, x30, x40)
# define phi
# define psi
# define all involved symbolic functions

#constants

M, m, g, ell = sym.symbols('M, m, g, ell')

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


# prints the outputs of equation 3.3

print('Equation 3.3a Answer'),sym.pprint(phi_deriv_F_at_equlibrium)

print('Equation 3.3b Answer'),sym.pprint(phi_deriv_x3_at_equlibrium)

print('Equation 3.3c Answer'),sym.pprint(phi_deriv_x4_at_equlibrium)



print('Equation 3.3d Answer'),sym.pprint(psi_deriv_F_at_equlibrium)

print('Equation 3.3e Answer'),sym.pprint(psi_deriv_x3_at_equlibrium)

print('Equation 3.3f Answer'),sym.pprint(psi_deriv_x4_at_equlibrium)
