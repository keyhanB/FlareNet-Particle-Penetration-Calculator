from math import exp
from math import pi
from math import asin
from math import log
from math import cos
from math import sin
from math import fabs


def Area_with_Diameter (A):
    Area = pi * (A ** 2) / 4
    return Area


def LPM_to_M3perSec (A):
    M3persec = A / 60000
    return M3persec


def ActualFlowRateLPM (Q_Slpm, Pg_pa, T_0, P_0_pa, Tg_K):
    AFR = Q_Slpm / Pg_pa / T_0 * P_0_pa * Tg_K
    return AFR


def Inch_to_Meter (A):
    if A != 0:
        meter = A * 0.0254
    else:
        meter = 0
    return meter


def Deg_to_Radian (A):
    Rad = A * pi / 180
    return Rad


def Flow_Velocity (Flowrate, Area):
    Velocity = Flowrate / Area
    return Velocity


def Kn (mfp, dp):
    Kn = 2 * mfp / dp
    return Kn


def Cc (Kn):
    result = 1+Kn * (1.257+.4 * exp (-1.1 / Kn))
    return result


def visc (T, P):  # ASSUME INDEPENDENT OF P HERE, Sutherland Equation
    # Assume Air, https://www.cfd-online.com/Wiki/Sutherland%27s_law
    U0 = 1.716 * 10 ** (-5)
    T_ref = 273.15
    b = 1.458e-6
    S = 110.4
    Viscosity = U0 * ((T / T_ref) ** 1.5) * ((T_ref+S) / (T+S))
    return Viscosity


def mfp (visc, T, P):
    # For Air
    M = .029  # kg/mol/K
    R = 8.314  # J/mol/K
    l = 2 * visc / (P * (8 * M / pi / R / T) ** 0.5)
    return l


def KinViscosity (visc, rho):
    A = visc / rho
    return A


def Reynolds (rho, Velocity, Diameter, visc):
    A = rho * Velocity * Diameter / visc
    return A


def SpecificHeat (T):
    # For Air
    A = 287 * (3.653-0.001334 * T+0.000003291 * (T ** 2)-0.00000000191 * (T ** 3)+0.000000000000275 * (T ** 4))
    return A


def Prandtl (C, visc, K):
    A = C * visc / K
    return A


def ThermalCond (T):
    # For Air
    A = -0.0003933+0.00010184 * T-0.000000048574 * (T ** 2)+0.000000000015207 * (T ** 3)
    return A


def rho (T, P):
    result = P / 287.05 / T  # air density kg/m^3, R in J/kg/K T in K
    return result


def Mobility (dp, Cc, visc):
    B = Cc / (3 * visc * dp * pi)
    return B


def Particle_Diff_Coeff (Mobility, T):
    k = 1.3806e-23  # Boltzmann constant
    D = Mobility * k * T
    return D


def Zeta (Diffusion, Length, Q):
    A = pi * Diffusion * Length / Q
    return A


def Schmidt (visc, density, diffusion):
    A = visc / (density * diffusion)
    return A


def Effective_Density (K, Dm, dp):
    A = K * dp ** (Dm-3)
    return A


def Gravit_Depos (Length, V_ts, d_tube, U):
    A = Length * V_ts / (d_tube * U)
    return A


def K_Thermopho (Knud, Cc, Thermal_Cond_Fluid, Thermal_Cond_Particle):
    Cs = 1.17  # Talbot et al., 1980; Montassier et al., 1991
    Ct = 2.18
    Cm = 1.14
    K_ratio = Thermal_Cond_Fluid / Thermal_Cond_Particle
    K = (2 * Cs * (K_ratio+Ct * Knud) * Cc) / ((1+3 * Cm * Knud) * (1+2 * K_ratio+2 * Ct * Knud))
    return K


def Thermophoresis_Eff (Kth, Re, Pr, Re_t, T_wall, T_entrance, k, Cp, Q, rho_g, L, D):
    if T_entrance!=T_wall:
        if Re < Re_t:  # laminar
            # Jyh-Shyan Lin, Chuen-Jinn Tsai, 2003, Thermophoretic deposition efficiency in a cylindrical tube taking into account developing flow at the entrance region
            A = 1-0.783 * (Pr * Kth * (T_entrance-T_wall) / T_wall) ** (0.94)
            # Stratmann et al. (1994)
            # A= exp(-0.845*((Pr*Kth+0.025)/((T_wall/(T_entrance-T_wall))+0.28))**0.932)
            # Batchelor and Shen (1985)
            # A =1-(Pr*Kth*((T_entrance-T_wall)/T_entrance)*(1+(1-Pr*Kth)*((T_entrance-Tw)/T_entrance)))
            # Walker et al. (1979)
            # A= 1-((Pr*Kth/T_wall)*(T_entrance-T_wall))
        else:  # Turbulence
            f = (0.79 * log (Re)-1.64) ** (-2)  # smooth surface
            Nu = ((f / 8) * (Pr * (Re-1000))) / (1+12.7 * ((f / 8) ** (0.5)) * ((Pr ** (2 / 3))-1))  # fully developed turbulent flow, 0.5<Pr<2000, 3000<ReD<5*10e6 (Incropera and DeWitt, 1996)
            h = Nu * k / D
            A = ((T_wall+(T_entrance-T_wall) * (exp (-1 * ((pi * D * h * L) / (rho_g * Q * Cp))))) / T_entrance) ** (Pr * Kth)
    else:
        A=1
    return A


def Quality_Of_Sett_Eff (V_ts, theta, U):
    A = 1-(V_ts * sin (theta) / U)
    return A


def Settling_Eff (Zg, Re, Re_t, theta):
    if Re < Re_t:  # laminar flow expression Heyder and Gebhart (1977)
        ep = 0.75 * Zg * cos (theta)
        B = (1-(ep ** (2 / 3))) ** 0.5
        eta = 1-(2 / pi) * ((2 * ep * B)-((ep ** (1 / 3)) * B)+(asin (ep ** (1 / 3))))
    else:  # Turbulent flow expression Schwendiman et al. 1975
        eta = exp (-4 * Zg * cos (theta) / pi)
    return eta


def Terminal_Velocity (dp, rho_p, Cc, visc):
    A = ((dp ** 2) * rho_p * 9.81 * Cc / visc) / 18
    return A


def Vt_Turb_Depos (Stk, Re, U):
    A = U * (((6 * 10 ** (-4)) * (0.0395 * Stk * Re ** (3 / 4)) ** 2)+(2 * (10 ** (-8)) * Re)) / (5.03 * Re ** (1 / 8))
    return A


def Diffusion_Eff (zeta, Sch, Re, Re_t):
    if Re < Re_t:  # laminar flow expression from Holman (1972)
        Sh = 3.66+(0.2672 / (zeta+0.10079 * zeta ** (1 / 3)))  # Sherwood Number
        eta = exp (-1 * zeta * Sh)
    else:  # Turbulent flow expression Friedlander (1977)
        Sh = 0.0118 * Re ** (7 / 8) * Sch ** (1 / 3)  # Sherwood Number
        eta = exp (-1 * zeta * Sh)
    return eta


def Inertial_Turb_Depos_Eff (d_tube, Length, Vt_Turb, Q, Re, Re_t):
    if Re > Re_t:
        A = exp (-1 * ((pi * d_tube * Length * Vt_Turb) / Q))
    elif Re > 15600:  # not valid
        A = 1
    else:
        A = 1
    return A


def Inertial_Bend_Depos_Eff (Stk, Bend_Rad, Re, Re_t):
    if Re < Re_t: # check angle is in Rad or Deg ??? (two contradictory sources)
        A = (1+(Stk / 0.171) ** (0.452 * (Stk / 0.171)+2.242)) ** (-2 * (Bend_Rad * 180 / pi) / pi)
    else:
        A = exp (-2.823 * Stk * (Bend_Rad * 180 / pi))
    return A


def Inertial_Deposit_Contraction_Eff (Stk, Theta_Rad, A_o, A_i):
    # Validity: 0.001≤Stk(1−Ao/Ai)≤100 and 12≤θcont≤90,  Muyshondt et al. (1996)
    A = 1-(1 / (1+((2*Stk * (1-(A_o / A_i))) / (3.14 * exp (-0.0185 * (Theta_Rad * 180 / pi)))) ** (-1.24)))
    return A


def Stokes (V_ts, U, Length):
    A = (V_ts / 9.81) * U / Length
    return A


def Aspiration_Eff (Stk, U0, U, theta_Rad):
    # Liu, Zhang, and Kuehn (1989) correlation for isoaxial aspiration efficiency for 0.01 < Stk< 100 and 0.1 < U0/U < 10.
    theta = theta_Rad * 180 / pi
    Stk0 = (U0 / U) * Stk  # in correlation, Stokes is based on ambient velocity
    eta_asp = 1

    if theta == 0:  # isoaxial
        if (U0 / U) >= 1:  # sub-isokinetic sampling
            eta_asp = 1+(((U0 / U)-1) / (1+(0.418 / Stk0)))
        else:  # super-isokinetic sampling
            eta_asp = 1+(((U0 / U)-1) / (1+(0.506 * ((U0 / U) ** 0.5) / Stk0)))

    # Anisoaxial
    elif fabs (theta) <= 60:  # 0.02 < Stk <4 and 0.5 < U0/U < 2
        Stk0_prime = Stk0 * exp (0.022 * fabs (theta))  # Durham and Lundgren (1980)
        k = 2+(0.617 * (U / U0))
        eta_asp = 1+(((U0 / U) * cos (theta_Rad))-1) * ((1-((1+k * Stk0_prime) ** (-1))) / (1-((1+2.617 * Stk0_prime) ** (-1)))) * (1-((1+0.55 * Stk0_prime * exp (0.25 * Stk0_prime)) ** (-1)))

    elif fabs (theta) <= 90:  # 0.02 <Stk < 0.2 and 0.5 < U0/U < 2
        eta_asp = 1+(((U0 / U) * cos (theta_Rad))-1) * (3 * Stk0 ** ((U / U0) ** 0.5))
    return eta_asp


def InletTrans_Eff (Stk, U0, U, theta_Rad, Zg, Re):
    eta_inert = 1
    theta = theta_Rad * 180 / pi
    Stk0 = (U0 / U) * Stk  # in correlation, Stokes is based on ambient velocity

    # transmission efficiency for gravitational settling , Hangal and Willeke (1990b)
    K_theta = (Zg ** 0.5) * (Stk0 ** 0.5) * (Re ** (-0.25)) * cos (theta_Rad)
    eta_grav = exp (-4.7 * (K_theta ** (0.75)))

    # Hangal and Willeke (1990b) give the transmission efficiency for inertia
    if theta == 0:  # isoaxial
        if U0 >= U:  # sub-isokinetic sampling , 0.01 < Stk < 100 and 1 < U0/U < 10
            r = (U0 / U)-1
            eta_inert = (1+r / (1+2.66 * (Stk0 ** (-2 / 3)))) / (1+r / (1+(0.418 / Stk0)))
        elif U0 < U:  # super-isokinetic sampling ,0.02 < Stk < 4 and 0.25 < U0/U < 1.0
            IV = 0.09 * (Stk0 * (U-U0) / U0) ** 0.3
            eta_inert = exp (-75 * (IV ** 2))

    else:  # Anisoaxial
        if U0 <= U:
            IV = 0.09 * (Stk0 * cos (theta_Rad) * (U-U0) / U0) ** 0.3
        else:
            IV = 0
        alpha = 12 * ((1-(fabs (theta) / 90))-exp (-1 * fabs (theta)))
        if theta < 0:  # downward sampling
            IW = Stk0 * ((U0 / U) ** 0.5) * sin ((theta-alpha) * pi / 180) * sin (((theta-alpha) / 2) * pi / 180)
        if theta > 0:  # upward sampling
            IW = Stk0 * ((U0 / U) ** 0.5) * sin ((theta+alpha) * pi / 180) * sin (((theta+alpha) / 2) * pi / 180)
        eta_inert = exp (-75 * (IV+IW) ** 2)

    eta = eta_grav * eta_inert
    return eta
