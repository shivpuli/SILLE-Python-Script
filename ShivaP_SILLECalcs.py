# SILLE Calculations
import numpy as np
import math
import matplotlib.pyplot as plt
# Inputs to NASA CEA
p_chamber = 100 # psi
ethanol_k = 298 # temp (k) of ethanol (liquid), 100% proportion, no specified density or enthalpy
n2o = 220 # temp (k) of nitrous (gas), 100% proportion, no specified density or enthalpy
of = 2.5 # ox/fuel ratio is 2.5. Must be < 3
Pc_Pe = p_chamber/14.73 # pressure_chamber / pressure_exit 
# Solved parameters from NASA's CEA
mach = np.array([0, 1, 1.926]) # chamber,throat,exit mach numbers
gamma = np.array([1.2371, 1.2517, 1.2685]) # chamber,throat,exit heat capacity ratios 
p_array = np.array([6.8947, 3.8148, 1.0156]) # (bar) chamber,throat,exit pressures (bar = 100,000 Pa = 14.503 psi)
p_pa = p_array * 100000 # (Pa) converting bar to pascals for calcs
temp = np.array([2307.30, 2049.39, 1553.46]) # (K) chamber,throat,exit temperatrues
rho = np.array([7.6003e-1, 4.7379e-1, 1.6645e-1]) #  (kg/m^3) chamber,throat,exit densities
thrust = 222.411 # (N) thrust
Isp = np.array([0, 1006.3, 1695.7]) # chamber,throat,exit exhaust velocity (m/s)
dP = 0.2 * p_pa[0] # from the SILLE lecture
n2o_mm = 44.013 # molar mass of nitrous (g/mol)
ethanol_mm = 46.07 # molar mass of ethanol (g/mol)
Cd = 0.65 # discharge coefficient 
L_star = 1 # SILLE lecture
R_bar = 8.314 # universal gas constant (J/(mol*k))
Ac_ratio = 14 # kinda just vibed. huzel says > 3.5 
Cp = np.array([2.0352, 1.9315, 1.8474]) # chamber,throat,exit, specific heat (kJ/(kg*K))
rho_fuel = 789 # density of ethanol at room temp (kg/m^3)
conv_angle = 30 # converging section angle. 30-45, more means more chaotic flow but smaller chamber length. 
half_angle = 15 # chosen based on Huzel and Huang
n = 5 # number of holes for doublet
## Calculations
R_exit = 1000*(Cp[2] * (gamma[2]-1))/gamma[2] # gas constant of the mixture at exit
R_chamber = 1000 * Cp[0] * (gamma[0] - 1) / gamma[0] # gas constant in chamber
rho_ox = ((p_pa[0]+dP) * n2o_mm) / (R_bar * 298 * 1000) # solving for nitrous density, temp is 298
# Equation 6 from NASA Isentropic Flow
p_ratio = p_pa[2]/p_pa[0] # ratio of chamber vs throat pressure
mach_exit = math.sqrt((2 * (p_ratio**(-(gamma[1]-1)/(gamma[1])) - 1)) /  (gamma[1] - 1))
# Solving for speed of sounds
sos = math.sqrt(gamma[2] * R_exit * temp[2]) # speed of sound (m/s) using exit conditions
# Equation 1 for NASA. Used Isp instead, but both values were the same
v_exit = mach_exit * sos # velocity at exit (m/s) 
# Mass flow rate calculation. Mass flow has the units of kg/s so that's just F/V
mdot = thrust/Isp[2] # used isp instead of v_exit because both were the same 
mdot_ox = mdot * of / (of +1) # mdot for the oxidizer side
mdot_fuel = mdot / (of +1) # mdot for the fuel side
# Start Area calculations
A_ox = mdot_ox / (Cd * math.sqrt(2 * rho_ox * dP)) # area of the oxidizer (m^2)
A_fuel = mdot_fuel / (Cd * math.sqrt(2 * rho_fuel * dP)) # area of the fuel (m^2)
A_ox_hole = A_ox / n # dividing by n to find area of each hole (m^2)
A_fuel_hole = A_fuel / n # dividing by n to find area of each hole (m^2)
r_ox = math.sqrt(A_ox_hole/math.pi) # radius of each hole for ox (m)
r_fuel = math.sqrt(A_fuel_hole/math.pi) # radius of each for fuel (m)
print("---------- Injector Hole Sizing ---------- ")
print(f"Radius of ox holes are {39.3701 * r_ox:.6f}in") # sanity check values 
print(f"Radius of fuel holes are {39.3701 * r_fuel:.6f}in") # sanity check values
print()
# Area of throat
A_throat = (mdot * math.sqrt(temp[0])) / (p_pa[0] * math.sqrt(gamma[0]/R_chamber) * ((gamma[0] + 1)/2)**(-(gamma[0] + 1)/(2 * (gamma[0] - 1)))) # (m^2) we're using chamber values
# Equation 9 for NASA (Ratio of exit to throat)
A_et = ((gamma[0]+1)/2)**-((gamma[0]+1)/(2*(gamma[0]-1))) * ((1+(mach[2]**2*(gamma[0]-1))/2)**((gamma[0]+1)/(2*(gamma[0]-1))))/mach[2] # (m^2) equation 9 on isentropic flow
# Geometry for Engine
A_exit = A_et * A_throat # (m^2) solving for the nozzle area
r_exit = math.sqrt(A_exit/math.pi) # (m) radius of the nozzle 
r_throat = math.sqrt(A_throat/math.pi) # (m) radius of the throat 
L_nozzle = (r_exit - r_throat)/math.tan(math.radians(half_angle)) # (m) nozzle length 
A_chamber = Ac_ratio * A_throat # (m^2) Using the assumed ac_ratio of 10
r_chamber = math.sqrt(A_chamber/math.pi) # (m) radius of chamber
L_conv = (r_chamber - r_throat) / math.tan(math.radians(conv_angle)) # (m) length from r_c to r_t, 30 was chosen on vibes
V_chamber = A_throat * L_star # (m^3) volume of the chamber
V_conv = (math.pi * L_conv / 3) * (r_chamber**2 + r_chamber*r_throat + r_throat**2) # (m^3) converging part volume. volume of a conical frustum
L_cyl = (V_chamber - V_conv) / A_chamber # (m) length of rectangle part of chamber. Rearranging V_cyl / A_chamber, where V_cyl = V_chamber - V_converging
L_chamber = L_conv + L_cyl # (m) length of the chamber 
L_tot = L_chamber + L_nozzle
print("---------- Geometry of Engine ---------- ")
print(f"The total length is {L_tot * 39.3701:.4f}in") # sanity check to make sure it's < 6 in
print()
print(f"Chamber length is {39.3701 * L_chamber:.4f}in") # sanity check values 
print(f"Chamber radius is {39.3701 * r_chamber:.4f}in") # sanity check values
print()
print(f"Throat radius is {39.3701 * r_throat:.4f}in") # sanity check values
print()
print(f"Half angle is {15} degrees") # sanity check values
print(f"The convergence angle is {conv_angle} degrees") # to inform
print()
print(f"Nozzle length is {39.3701 * L_nozzle:.4f}in") # sanity check values 
print(f"Nozzle radius is {39.3701 * r_exit:.4f}in") # sanity check values 
print()
print(f"Expansion ratio is {A_et:.4f}") # sanity check values 
print(f"Contraction ratio is {Ac_ratio}") # sanity check values 
## RPA Style Plot
# Dimensions to inches
Rc = r_chamber * 39.3701
Rt = r_throat * 39.3701
Re = r_exit * 39.3701
Lcyl = L_cyl * 39.3701
Lconv = L_conv * 39.3701
Lc = L_chamber * 39.3701
Le = L_nozzle * 39.3701
Ltotal = L_tot * 39.3701
# x-locations
x0 = 0.0
x_cyl_end = Lcyl
x_throat = Lc
x_exit = Lc + Le
# blending
blend1 = min(0.15, 0.25 * Lconv)   # chamber-to-converging blend
blend2 = min(0.15, 0.25 * Lconv)   # converging-to-throat blend
# need bezier curves to gett curves
def quad_bezier(p0, p1, p2, n=30): 
    """Return x,y arrays for a quadratic Bezier curve."""
    t = np.linspace(0, 1, n)
    x = (1-t)**2*p0[0] + 2*(1-t)*t*p1[0] + t**2*p2[0]
    y = (1-t)**2*p0[1] + 2*(1-t)*t*p1[1] + t**2*p2[1]
    return x, y
# empty arrays
x_upper = []
y_upper = []
# chamber
x1 = np.linspace(x0, x_cyl_end - blend1, 20)
y1 = np.full_like(x1, Rc)
x_upper.extend(x1)
y_upper.extend(y1)
# rounded chamber converging corner
p0 = (x_cyl_end - blend1, Rc)
p1 = (x_cyl_end, Rc)
# approximate point on converging line
p2 = (x_cyl_end + blend1, Rc - (Rc - Rt) * (blend1 / Lconv))
xb, yb = quad_bezier(p0, p1, p2)
x_upper.extend(xb)
y_upper.extend(yb)
# converging line before throat blend
x3_start = x_cyl_end + blend1
x3_end = x_throat - blend2
x3 = np.linspace(x3_start, x3_end, 40)
y3 = Rc - (Rc - Rt) * ((x3 - x_cyl_end) / Lconv)
x_upper.extend(x3)
y_upper.extend(y3)
# rounded converging-to-throat corner
p0 = (x_throat - blend2, Rc - (Rc - Rt) * ((Lconv - blend2) / Lconv))
p1 = (x_throat, Rt)
p2 = (x_throat + blend2, Rt + (Re - Rt) * (blend2 / Le))
xb, yb = quad_bezier(p0, p1, p2)
x_upper.extend(xb)
y_upper.extend(yb)
# straight diverging nozzle line
x5 = np.linspace(x_throat + blend2, x_exit, 40)
y5 = Rt + (Re - Rt) * ((x5 - x_throat) / Le)
x_upper.extend(x5)
y_upper.extend(y5)
x_upper = np.array(x_upper)
y_upper = np.array(y_upper)
# symmetry is fun
x_lower = x_upper
y_lower = -y_upper
fig, ax = plt.subplots(figsize=(10, 4))
# contour
ax.plot(x_upper, y_upper, linewidth=2)
ax.plot(x_lower, y_lower, linewidth=2)
# centerline
ax.axhline(0, linestyle="--", linewidth=1)
# vertical lines
ax.axvline(x_cyl_end, linestyle=":", linewidth=1)
ax.axvline(x_throat, linestyle=":", linewidth=1)
ax.axvline(x_exit, linestyle=":", linewidth=1)
# Dimension helper function
def vertical_dimension(x, y0, y1, label, label_offset=0.08):
    ax.plot([x, x], [y0, y1], linewidth=1)
    ax.plot([x - 0.04, x + 0.04], [y0, y0], linewidth=1)
    ax.plot([x - 0.04, x + 0.04], [y1, y1], linewidth=1)
    ax.text(x - label_offset, (y0 + y1) / 2, label, ha="right", va="center")
# radius bars
vertical_dimension(x0 - 0.25, 0, Rc, "Rc")
vertical_dimension(x_throat - 0.18, 0, Rt, "Rt")
vertical_dimension(x_exit + 0.18, 0, Re, "Re", label_offset=-0.18)
# lengths
y_len_label = -1.25 * Rc
ax.annotate("", xy=(x0, y_len_label), xytext=(x_cyl_end, y_len_label),arrowprops=dict(arrowstyle="<->", linewidth=1))
ax.text((x0 + x_cyl_end) / 2, y_len_label - 0.08, "Lcyl", ha="center", va="top")
ax.annotate("", xy=(x_cyl_end, y_len_label), xytext=(x_throat, y_len_label),arrowprops=dict(arrowstyle="<->", linewidth=1))
ax.text((x_cyl_end + x_throat) / 2, y_len_label - 0.08, "Lconv", ha="center", va="top")
ax.annotate("", xy=(x_throat, y_len_label), xytext=(x_exit, y_len_label),arrowprops=dict(arrowstyle="<->", linewidth=1))
ax.text((x_throat + x_exit) / 2, y_len_label - 0.08, "Le", ha="center", va="top")
# tot chamber length
y_lc_label = 1.18 * Rc
ax.annotate("", xy=(x0, y_lc_label), xytext=(x_throat, y_lc_label),arrowprops=dict(arrowstyle="<->", linewidth=1))
ax.text((x0 + x_throat) / 2, y_lc_label + 0.08, "Lc", ha="center", va="bottom")
# angles
ax.text(x_cyl_end + 0.25 * Lconv, 0.75 * Rc, f"β = {conv_angle}°", ha="center")
ax.text(x_throat + 0.55 * Le, Rt + 0.12, f"θe = {half_angle}°", ha="center")
# axes + title
ax.set_title("RPA Style Engine Diagram ", pad=18)
ax.set_xlabel("Length (in)")
ax.set_ylabel("Radius (in)")
ax.set_aspect("equal", adjustable="box")
ax.grid(True, alpha=0.4)
# Dimensions
summary_text = (
    f"Rc = {Rc:.4f} in\n"
    f"Rt = {Rt:.4f} in\n"
    f"Re = {Re:.4f} in\n"
    f"Lcyl = {Lcyl:.4f} in\n"
    f"Lconv = {Lconv:.4f} in\n"
    f"Lc = {Lc:.4f} in\n"
    f"Le = {Le:.4f} in\n"
    f"Ltotal = {Ltotal:.4f} in\n"
    f"β = {conv_angle:.1f}°\n"
    f"θe = {half_angle:.1f}°\n"
    f"Ae/At = {A_et:.4f}\n"
    f"Contraction Ratio = {Ac_ratio:.1f}"
)

ax.text(
    1.03, 0.5, summary_text,
    transform=ax.transAxes,
    fontsize=9,
    va="center",
    ha="left",
    bbox=dict(boxstyle="round,pad=0.5", facecolor="white", edgecolor="black")
)
# space on the right
plt.subplots_adjust(right=0.72)
# margins
ax.set_xlim(x0 - 0.5, x_exit + 0.6)
ax.set_ylim(-1.55 * Rc, 1.45 * Rc)
plt.show()
plt.show()
