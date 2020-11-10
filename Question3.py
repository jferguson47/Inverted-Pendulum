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
a_value = a.subs([(M, M_value), (m, m_value)])
b_value = b.subs([(M, M_value), (m, m_value), (g, g_value)])
c_value = c.subs([(M, M_value), (m, m_value), (g, g_value), (ell, ell_value)])
d_value = d.subs([(M, M_value), (m, m_value), (g, g_value), (ell, ell_value)])

s, t = sym.symbols('s, t')
c, d = sym.symbols('c, d', real=True, positive=True)


# --------------------------------1. Impulse, Step and frequency Response of G_theta -----------------------------------
G_theta = - c / (s**2 - d)


# the value of A (amplitude) is set to 1

# 1.1. Impulse response for G_theta
# f_t = KICK (DIRAC PULSE) - impulse
Fs_impulse = 1                                 # Amplitude (A) = 1
X3s_impulse = G_theta * Fs_impulse             # the inverse of G_theta will be the the impulse response of the system
X3t_impulse = sym.inverse_laplace_transform(X3s_impulse, s, t)  # system output t domain
sym.pprint(X3t_impulse)    # put in report
print(sym.latex(X3t_impulse))
# 1.2. Step Response (PUSH) of G_theta: what is the system response when

# the input is a step pulse (Heaviside function)

# G = X / F

Fs_step = 1 / s                                             # system input Amplitude (1)/ s
X3s_step = G_theta * Fs_step                                # system output the s domain x=G*F
X3t_step = sym.inverse_laplace_transform(X3s_step, s, t)    # system output in the t domain
sym.pprint(X3t_step)    # put in report
print(sym.latex(X3t_step))

# 1.3. Frequency response (Shake) for G_theta: what happens to x3(t) when the input is
# sin(omega*t) for some omega?

# omega was given the value of 1    therefore the input u(t) is sin(t)
# the laplace transform of sin(t) is 1/*s**2 +1) which is what is used for F_S

# G=x/F

Fs_frequency = 1 / (s**2 + 1)
X3s_frequency = G_theta * Fs_frequency           # system output the s domain x=G*F
X3t_frequency = sym.inverse_laplace_transform(X3s_frequency, s, t)  # system output in the t domain
sym.pprint(X3t_frequency.simplify())    # put in report
#print(sym.latex(X3t_frequency))

# --------------------------------2. Impulse, Step and frequency Response of G_x ---------------------------------------

G_x=(a*s**2-a*d+b*c)/(s**4-d*s**2)

# the value of A (amplitude) is set to 1

# 2.1. Impulse response for G_x
# f_t = KICK (DIRAC PULSE) - impulse
Fsx_impulse = 1                                 # Amplitude (A) = 1
X3sx_impulse = G_x * Fsx_impulse                 # the inverse of G_x will be the the impulse response of the system
X3tx_impulse = sym.inverse_laplace_transform(X3sx_impulse, s, t)  # system output t domain
sym.pprint(X3tx_impulse)    # put in report
print(sym.latex(X3tx_impulse))

# 2.2. Step Response (PUSH) of G_x: what is the system response when

# the input is a step pulse (Heaviside function)

# G = X / F

Fsx_step = 1 / s                                            # system input Amplitude (1)/ s
X3sx_step = G_x * Fsx_step                                    # system output the s domain x=G*F
X3tx_step = sym.inverse_laplace_transform(X3sx_step, s, t)    # system output in the t domain
sym.pprint(X3tx_step)    # put in report
print(sym.latex(X3tx_step))

# 2.3. Frequency response (Shake) for G_x: what happens to x3(t) when the input is
# sin(omega*t) for some omega?

# omega was given the value of 1    therefore the input u(t) is sin(t)
# the laplace transform of sin(t) is 1/*s**2 +1) which is what is used for F_S

# G=x/F

Fsx_frequency = 1 / (s**2 + 1)
X3sx_frequency = G_x * Fsx_frequency           # system output the s domain x=G*F
X3tx_frequency = sym.inverse_laplace_transform(X3sx_frequency, s, t)  # system output in the t domain
sym.pprint(X3tx_frequency.simplify())    # put in report
print(sym.latex(X3t_frequency))
