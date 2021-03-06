# Programmed by Keyhan Babaee Under Prof. Steven Rogak supervision, https://github.com/KeyhanB
# Version 1.4
# Dec 22, 2018
import Functions as KN
import matplotlib.pyplot as plt
import numpy as np
import os
import xlrd
import xlsxwriter
from math import log
from stackedBarGraph import StackedBarGrapher
import logging

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', filename='Log.txt', level=logging.INFO, filemode='w')

# Important Note: all pressures are in Pascal and all temperatures are in Kelvin

if __name__ == "__main__":

    # Script Directory
    script_dir = os.path.dirname(os.path.realpath('__file__'))

    # Checking for directory existence and if there isn't any creation of the folder
    Excel_Directory = "Excel Output"
    Excel_Directory = os.path.join(script_dir, Excel_Directory)
    Graph_Directory = "Graph Output"
    Graph_Directory = os.path.join(script_dir, Graph_Directory)
    if not os.path.exists(Excel_Directory):
        os.makedirs(Excel_Directory)
    if not os.path.exists(Graph_Directory):
        os.makedirs(Graph_Directory)

    Input_Data_Excel = "Excel Input/Input.xlsx"
    Input_Data_Excel = os.path.join(script_dir, Input_Data_Excel)
    Output_Data_Excel = "Main Report.xlsx"
    Output_Data_Excel = os.path.join(Excel_Directory, Output_Data_Excel)

    # Reading Input file
    InputExcel_workbook = xlrd.open_workbook(Input_Data_Excel)
    InputExcel_Sheet1 = InputExcel_workbook.sheet_by_index(0)
    General_Input = []
    for row_index in range(0, InputExcel_Sheet1.nrows):
        General_Input.append([InputExcel_Sheet1.cell(row_index, col_index).value for col_index in range(InputExcel_Sheet1.ncols)])
    Input_Dimension_Row = len(General_Input)  # Number of rows of input data
    Input_Dimension_Col = len(General_Input[1])  # Number of columns of input data

    ####################################################### (Controlling Variables)
    T_0 = 294.26  # Kelvin for calculating actual flow rate
    P_0 = 101325  # Pa for calculating actual flow rate
    D_small = 10e-9  # smallest diameter of computation (10nm)
    D_large = 10e-6  # largest diameter of computation (10um)
    Nd = 199  # number of points (between above diameters) considered for the computation
    Reynolds_Lam_to_Turb = 3000  # transition between laminar to turbulence flow
    Particle_thermal_conduction = 0.2  # w/m/k
    eff_k = 0.418  # effective density variable for the result to be in kg/m^3
    eff_dm = 2.56  # effective density variable for the result to be in kg/m^3
    Dilution_Temp = 293  # Kelvin for calculating Dilution
    Computed_data_Column = 25 + 1  # Number of variables of data for each section needed to be calculated (change if you want to add more variables)
    Number_Particle_Var = 22 + 1  # Number of variables of data for each particle diameter needed to be calculated (change if you want to add more variables)
    Y_limit_lower = 0.5  # Graph limit _ Y axis
    Y_limit_upper = 1.05  # Graph limit _ Y axis
    Y_Axis_Unit = 0.1  # Graph unit _ Y axis
    Figure_DPI = 1200  # Graph dpi
    Polynomial_Deg = 6  # Graph fitting line polynomial degree
    Diffusion_Enable = 1  # Diffusion penetration plotting (0 for not plotting)
    Gravitational_Enable = 1  # Gravitational penetration plotting (0 for not plotting)
    Inertial_Enable = 1  # Turbulent inertial penetration plotting (0 for not plotting)
    Thermophoretic_Enable = 1  # Thermophoretic penetration plotting (0 for not plotting)
    Bend_Enable = 1  # Inertia Bend penetration plotting (0 for not plotting)
    Contr_Enable = 1  # Contraction penetration plotting (0 for not plotting)
    Asp_Enable = 1  # Aspiration penetration plotting (0 for not plotting)
    Inlet_Enable = 1  # Inlet transportation penetration plotting (0 for not plotting)
    Total_Enable = 1  # Total penetration plotting (0 for not plotting)
    Fit_Enable = 1  # Fitted polynomial plotting (0 for not plotting)

    #######################################################
    logging.info("FPPC Initiated")
    Flow_Properties = [[None for i in range(Computed_data_Column)] for j in range(Input_Dimension_Row)]

    for i in range(1, Input_Dimension_Row):
        if i == 1:
            i = i - 1
            Flow_Properties[i][1] = "Inner Diameter (meter)"
            Flow_Properties[i][2] = "Length (meter)"
            Flow_Properties[i][3] = "Bend Angle (Rad)"
            Flow_Properties[i][4] = "Theta Angle (Rad)"
            Flow_Properties[i][5] = "Actual Flow rate (lpm)"
            Flow_Properties[i][6] = "Actual Flow rate (m3/s)"
            Flow_Properties[i][7] = "Flow velocity rate (m/s)"
            Flow_Properties[i][8] = "Dynamic Viscosity (kg/m/s)"
            Flow_Properties[i][9] = "Mean free path (meter)"
            Flow_Properties[i][10] = "Density (kg/m^3)"
            Flow_Properties[i][11] = "Kinematic Viscosity (m^2/s)"
            Flow_Properties[i][12] = "Reynolds"
            Flow_Properties[i][13] = "Specific heat (J/kg/K)"
            Flow_Properties[i][14] = "Thermal Conductivity (W/m/K)"
            Flow_Properties[i][15] = "Prandtl"
            Flow_Properties[i][16] = "Tube Cross Section (m^2)"
            Flow_Properties[i][17] = "Inner Diameter (Second) (meter)"
            Flow_Properties[i][18] = "Half Theta (Second) (Rad)"
            Flow_Properties[i][19] = "Tube Cross Section (Second) (m^2)"
            Flow_Properties[i][20] = "Probe Theta (Rad)"
            Flow_Properties[i][21] = "Wall Overall Heat Coefficient (W/(m^2.K))"
            Flow_Properties[i][22] = "Tube Surface Area(m^2)"
            Flow_Properties[i][23] = "Heat loss (W)"
            Flow_Properties[i][24] = "Calculated Outlet Temperature (K)"
            Flow_Properties[i][25] = "Calculated Wall Temperature (K)"

            i = i + 1

        ####### Checking if user provided Inlet and gas temperature if not using outlet temperature of the last section
        if General_Input[i][11] == -1:  # T_Inlet_Calculator
            if General_Input[i - 1][12] != -1:
                General_Input[i][11] = General_Input[i - 1][12]
            else:
                logging.info(f"Invalid Temperature for the section number {i}")
                raise Exception("Invalid Temperature for calculation")

        if General_Input[i][9] == -1:  # T_Gas_Calculator
            if General_Input[i - 1][12] != -1:
                General_Input[i][9] = General_Input[i - 1][12]
            else:
                logging.info(f"Invalid Temperature for the section number {i}")
                raise Exception("Invalid Temperature for calculation")
        ####################
        ####### Checking for dilution
        if i != 1:
            if General_Input[i - 1][14] < General_Input[i][14]:  # We have dilution
                NewFlow = General_Input[i][14] - General_Input[i - 1][14]
                General_Input[i][11] = (NewFlow * Dilution_Temp + General_Input[i - 1][14] * General_Input[i][11]) / (NewFlow + General_Input[i - 1][14])
                General_Input[i][9] = General_Input[i][11]
        ####################

        Flow_Properties[i][0] = General_Input[i][4]  # ignore
        Flow_Properties[i][1] = KN.Inch_to_Meter(General_Input[i][5])  # ID (meter)
        Flow_Properties[i][2] = KN.Inch_to_Meter(General_Input[i][6])  # length (meter)
        Flow_Properties[i][3] = KN.Deg_to_Radian(General_Input[i][7])  # Bend (Rad)
        Flow_Properties[i][4] = KN.Deg_to_Radian(General_Input[i][8])  # Theta (Rad)
        Flow_Properties[i][5] = KN.ActualFlowRateLPM(Q_Slpm=General_Input[i][14], Pg_pa=General_Input[i][13], T_0=T_0, P_0_pa=P_0, Tg_K=General_Input[i][9])  # Actual flow rate (lpm)
        Flow_Properties[i][6] = KN.LPM_to_M3perSec(Flow_Properties[i][5])  # Actual flow rate (m^3/s)
        Flow_Properties[i][7] = KN.Flow_Velocity(Flowrate=Flow_Properties[i][6], Area=KN.Area_with_Diameter(Flow_Properties[i][1]))  # flow velocity derived from Actual flow rate (m/s)
        Flow_Properties[i][8] = KN.visc(T=General_Input[i][9], P=General_Input[i][13])  # dynamic Viscosity (kg/m/s)
        Flow_Properties[i][9] = KN.mfp(visc=Flow_Properties[i][8], T=General_Input[i][9], P=General_Input[i][13])  # mean free path (meter)
        Flow_Properties[i][10] = KN.rho(T=General_Input[i][9], P=General_Input[i][13])  # Density (kg/m^3)
        Flow_Properties[i][11] = KN.KinViscosity(visc=Flow_Properties[i][8], rho=Flow_Properties[i][10])  # Kinematic Viscosity (m2/s)
        Flow_Properties[i][12] = KN.Reynolds(rho=Flow_Properties[i][10], Velocity=Flow_Properties[i][7], Diameter=Flow_Properties[i][1], visc=Flow_Properties[i][8])  # Reynolds
        Flow_Properties[i][13] = KN.SpecificHeat(General_Input[i][9])  # Specific heat (J/kg/K)
        Flow_Properties[i][14] = KN.ThermalCond(General_Input[i][9])  # Thermal Conductivity (W/m/K)
        Flow_Properties[i][15] = KN.Prandtl(C=Flow_Properties[i][13], visc=Flow_Properties[i][8], K=Flow_Properties[i][14])  # Prandtl
        Flow_Properties[i][16] = KN.Area_with_Diameter(Flow_Properties[i][1])  # Area (m^2)
        Flow_Properties[i][17] = KN.Inch_to_Meter(General_Input[i][16])  # ID (Second) (meter)
        Flow_Properties[i][18] = KN.Deg_to_Radian(General_Input[i][17])  # Half_Theta (Second) (Rad)
        Flow_Properties[i][19] = KN.Area_with_Diameter(Flow_Properties[i][17])  # Area (Second) (m^2)
        Flow_Properties[i][20] = KN.Deg_to_Radian(General_Input[i][20])  # Probe Theta (Rad)
        Flow_Properties[i][21] = KN.Overall_Heat_Coefficient(Re=Flow_Properties[i][12], RSI=General_Input[i][22], ReT=Reynolds_Lam_to_Turb, Pr=Flow_Properties[i][15], k=Flow_Properties[i][14], D=Flow_Properties[i][1],
                                                             L=Flow_Properties[i][2])  # W/(m^2.K)
        Flow_Properties[i][22] = KN.Tube_Surface_Area(L=Flow_Properties[i][2], D=Flow_Properties[i][1])  # Area (m^2)
        Flow_Properties[i][23] = KN.Heat_Loss(U=Flow_Properties[i][21], A=Flow_Properties[i][22], T_Air_Duct=General_Input[i][11], T_ambient=General_Input[i][21])  # watt
        Flow_Properties[i][24] = KN.Outlet_Temperature(Q_Watt=Flow_Properties[i][23], Rho_Gas=Flow_Properties[i][10], Volumetric_Flow=Flow_Properties[i][6], Specific_Heat=Flow_Properties[i][13], T_Air_Inlet=General_Input[i][11])  # Kelvin
        Flow_Properties[i][25] = KN.Wall_Temperature_Calc(RSI=General_Input[i][22], A=Flow_Properties[i][22], Q_Watt=Flow_Properties[i][23], T_ambient=General_Input[i][21])  # Kelvin

        ####### Checking if user provided wall and oulet temperature
        if General_Input[i][10] == -1:  # T_Wall_Calculator
            General_Input[i][10] = Flow_Properties[i][25]
        if General_Input[i][12] == -1:  # T_Outlet_Calculator
            General_Input[i][12] = Flow_Properties[i][24]
        # TODO look for the better temperature calculator across the length of the pipe
        if General_Input[i][12] <= General_Input[i][21]:  # T_Ambient_Calculator
            General_Input[i][12] = General_Input[i][21]
        ####################

    logging.info(f"Reading and Calculating Variables for sections was successful {Input_Dimension_Row - 1}")
    #######  Efficiencies: Rows are for each section and columns are for each particle diameter
    eta_Diff = [[1 for i in range(Nd)] for j in range(Input_Dimension_Row + 1)]  # Diffusion penetration
    eta_Grav = [[1 for i in range(Nd)] for j in range(Input_Dimension_Row + 1)]  # Gravitational penetration
    eta_Inert = [[1 for i in range(Nd)] for j in range(Input_Dimension_Row + 1)]  # Turbulent inertial penetration
    eta_TP = [[1 for i in range(Nd)] for j in range(Input_Dimension_Row + 1)]  # Thermophoretic penetration
    eta_Bend = [[1 for i in range(Nd)] for j in range(Input_Dimension_Row + 1)]  # Inertia Bend penetration
    eta_Cont = [[1 for i in range(Nd)] for j in range(Input_Dimension_Row + 1)]  # Contraction penetration
    eta_Asp = [[1 for i in range(Nd)] for j in range(Input_Dimension_Row + 1)]  # Aspiration penetration
    eta_Inlet = [[1 for i in range(Nd)] for j in range(Input_Dimension_Row + 1)]  # Inlet transportation penetration
    e_Total = [[1 for i in range(Nd + 1)] for j in range(Input_Dimension_Row + 1)]  # Total penetration
    # +1 for the fitting variables

    los_Diff = [[1 for i in range(Nd)] for j in range(Input_Dimension_Row + 1)]  # Diffusion Loss
    los_Grav = [[1 for i in range(Nd)] for j in range(Input_Dimension_Row + 1)]  # Gravitational Loss
    los_Inert = [[1 for i in range(Nd)] for j in range(Input_Dimension_Row + 1)]  # Turbulent inertial Loss
    los_TP = [[1 for i in range(Nd)] for j in range(Input_Dimension_Row + 1)]  # Thermophoretic Loss
    los_Bend = [[1 for i in range(Nd)] for j in range(Input_Dimension_Row + 1)]  # Inertia Bend Loss
    los_Cont = [[1 for i in range(Nd)] for j in range(Input_Dimension_Row + 1)]  # Contraction Loss
    los_Asp = [[1 for i in range(Nd)] for j in range(Input_Dimension_Row + 1)]  # Aspiration Loss
    los_Inlet = [[1 for i in range(Nd)] for j in range(Input_Dimension_Row + 1)]  # Inlet transportation Loss
    los_Total = [[1 for i in range(Nd + 1)] for j in range(Input_Dimension_Row + 1)]  # Total Loss

    los_Norm_Sum = [1 for i in range(Nd)]  # TODO check for better names
    los_Norm_Sum1 = [1 for i in range(Nd)]  # TODO check for better names
    los_Diff_Norm = [1 for i in range(Nd)]  # Normalized Diffusion Loss
    los_Grav_Norm = [1 for i in range(Nd)]  # Normalized Gravitational Loss
    los_Inert_Norm = [1 for i in range(Nd)]  # Normalized Turbulent inertial Loss
    los_TP_Norm = [1 for i in range(Nd)]  # Normalized Thermophoretic Loss
    los_Bend_Norm = [1 for i in range(Nd)]  # Normalized Inertia Bend Loss
    los_Cont_Norm = [1 for i in range(Nd)]  # Normalized Contraction Loss
    los_Asp_Norm = [1 for i in range(Nd)]  # Normalized Aspiration Loss
    los_Inlet_Norm = [1 for i in range(Nd)]  # Normalized Inlet transportation Loss
    los_Total_Norm = [1 for i in range(Nd)]  # Normalized Total Loss

    Fitting_Coefficients = [[1 for i in range(Polynomial_Deg + 1)] for j in range(Input_Dimension_Row + 1)]  # Fitting line Coefficients for total penetration

    ######################
    Particle_Variables = [[None for i in range(Number_Particle_Var)] for j in range(Nd + 1)]
    ######################
    D_ratio = (D_large / D_small) ** (1 / (Nd - 1))
    diam = []  # diameters of particle in meter
    Diameter_Nano = []  # diameters of particle in nm
    Ln_Diam1 = []  # diameters of particle in natural logarithm
    Width_nm = []  # bin width in nm
    D_labels = []  # bin labels

    # Generating particle diameters
    Keyhan = 10 ** 9  # Conversion factor
    First = 0.1
    for i in range(Nd):
        d = D_small * D_ratio ** (i)
        diam.append(d)
        Particle_Variables[i + 1][0] = d
        Diameter_Nano.append(d * Keyhan)
        Ln_Diam1.append(log(d))
        Width_nm.append(log(Diameter_Nano[i]) - log(First))
        First = Diameter_Nano[i]
        if i % (int(Nd / 3)) == 0:
            D_labels.append(np.around(Diameter_Nano[i], decimals=1))
        else:
            D_labels.append("")
    Width_nm[0] = Width_nm[1]
    ######################

    # Naming particle variables
    i = -1
    Particle_Variables[i + 1][1] = "Knudsen Number"
    Particle_Variables[i + 1][2] = "Cunningham Correction Factor"
    Particle_Variables[i + 1][3] = "Mobility (s/kg)"
    Particle_Variables[i + 1][4] = "Particle Diffusion Coefficient (m^2/s)"
    Particle_Variables[i + 1][5] = "Zeta"
    Particle_Variables[i + 1][6] = "Schmidt Number"
    Particle_Variables[i + 1][7] = "Diffusion Penetration Efficiency"
    Particle_Variables[i + 1][8] = "Effective Density (kg/m^3)"
    Particle_Variables[i + 1][9] = "Terminal Velocity (m/s)"
    Particle_Variables[i + 1][10] = "Gravitational Deposition Parameter"
    Particle_Variables[i + 1][11] = "Gravitational Penetration Efficiency"
    Particle_Variables[i + 1][12] = "Validity of Gravitational Calculation"
    Particle_Variables[i + 1][13] = "Stokes Number"
    Particle_Variables[i + 1][14] = "Vt for Inertial Turbulent"
    Particle_Variables[i + 1][15] = "Inertial Turbulent Transportation Efficiency"
    Particle_Variables[i + 1][16] = "Kth for Thermophoresis"
    Particle_Variables[i + 1][17] = "Thermophoresis Transportation Efficiency"
    Particle_Variables[i + 1][18] = "Inertial Bend Transportation Efficiency"
    Particle_Variables[i + 1][19] = "Inertial Contraction Transportation Efficiency"
    Particle_Variables[i + 1][20] = "Inlet Aspiration Efficiency"
    Particle_Variables[i + 1][21] = "Inlet Transportation Efficiency"
    Particle_Variables[i + 1][22] = "Total Penetration Efficiency"
    ##########################################################################
    for Section in range(1, Input_Dimension_Row):
        for i in range(Nd):
            # Diffusional penetration for laminar and turbulent flow
            Particle_Variables[i + 1][1] = KN.Kn(mfp=Flow_Properties[Section][9], dp=diam[i])  # Knudsen number
            Particle_Variables[i + 1][2] = KN.Cc(Particle_Variables[i + 1][1])  # Cunningham correction factor
            Particle_Variables[i + 1][3] = KN.Mobility(dp=diam[i], Cc=Particle_Variables[i + 1][2], visc=Flow_Properties[Section][8])  # Mobility (s/kg)
            Particle_Variables[i + 1][4] = KN.Particle_Diff_Coeff(Mobility=Particle_Variables[i + 1][3], T=General_Input[Section][9])  # Particle diffusion coefficient (m2/s)
            Particle_Variables[i + 1][5] = KN.Zeta(Diffusion=Particle_Variables[i + 1][4], Length=Flow_Properties[Section][2], Q=Flow_Properties[Section][6])  # Zeta
            Particle_Variables[i + 1][6] = KN.Schmidt(visc=Flow_Properties[Section][8], density=Flow_Properties[Section][10], diffusion=Particle_Variables[i + 1][4])  # Schmidt Number
            Particle_Variables[i + 1][7] = KN.Diffusion_Eff(zeta=Particle_Variables[i + 1][5], Sch=Particle_Variables[i + 1][6], Re=Flow_Properties[Section][12], Re_t=Reynolds_Lam_to_Turb)  # Diffusion Penetration

            # Gravitational penetration
            Particle_Variables[i + 1][8] = KN.Effective_Density(eff_k, eff_dm, diam[i])  # Effective density (kg/m3)
            Particle_Variables[i + 1][9] = KN.Terminal_Velocity(dp=diam[i], rho_p=Particle_Variables[i + 1][8], Cc=Particle_Variables[i + 1][2], visc=Flow_Properties[Section][8])  # Terminal velocity m/s
            Particle_Variables[i + 1][10] = KN.Gravit_Depos(Length=Flow_Properties[Section][2], V_ts=Particle_Variables[i + 1][9], d_tube=Flow_Properties[Section][1], U=Flow_Properties[Section][7])  # Gravitational deposition parameter
            Particle_Variables[i + 1][11] = KN.Settling_Eff(Zg=Particle_Variables[i + 1][10], Re=Flow_Properties[Section][12], Re_t=Reynolds_Lam_to_Turb, theta=Flow_Properties[Section][4])  # Gravitational Penetration
            Particle_Variables[i + 1][12] = KN.Quality_Of_Sett_Eff(V_ts=Particle_Variables[i + 1][9], theta=Flow_Properties[Section][4], U=Flow_Properties[Section][7])  # If it reaches to 1 we have reliable results

            # TURBULENT INERTIAL penetration
            Particle_Variables[i + 1][13] = KN.Stokes(V_ts=Particle_Variables[i + 1][9], U=Flow_Properties[Section][7], Length=Flow_Properties[Section][1])  # Stokes number based on particle relaxation time, pipe velocity and radius
            Particle_Variables[i + 1][14] = KN.Vt_Turb_Depos(Stk=Particle_Variables[i + 1][13], Re=Flow_Properties[Section][12], U=Flow_Properties[Section][7])  # Vt for turbulence deposition
            Particle_Variables[i + 1][15] = KN.Inertial_Turb_Depos_Eff(d_tube=Flow_Properties[Section][1], Length=Flow_Properties[Section][2], Vt_Turb=Particle_Variables[i + 1][14], Q=Flow_Properties[Section][6], Re=Flow_Properties[Section][12],
                                                                       Re_t=Reynolds_Lam_to_Turb)  # Turbulence Penetration

            # THERMOPHORESIS penetration
            Particle_Variables[i + 1][16] = KN.K_Thermopho(Knud=Particle_Variables[i + 1][1], Cc=Particle_Variables[i + 1][2], Thermal_Cond_Fluid=Flow_Properties[Section][14], Thermal_Cond_Particle=Particle_thermal_conduction)
            Particle_Variables[i + 1][17] = KN.Thermophoresis_Eff(Kth=Particle_Variables[i + 1][16], Re=Flow_Properties[Section][12], Pr=Flow_Properties[Section][15], Re_t=Reynolds_Lam_to_Turb, T_wall=General_Input[Section][10],
                                                                  T_entrance=General_Input[Section][11], T_Exit=General_Input[Section][12],
                                                                  k=Flow_Properties[Section][14], Cp=Flow_Properties[Section][13], Q=Flow_Properties[Section][6], rho_g=Flow_Properties[Section][10], L=Flow_Properties[Section][2],
                                                                  D=Flow_Properties[Section][1])  # Thermophoresis Penetration
            # Inertial penetration bend
            Particle_Variables[i + 1][18] = KN.Inertial_Bend_Depos_Eff(Stk=Particle_Variables[i + 1][13], Bend_Rad=Flow_Properties[Section][3], Re=Flow_Properties[Section][12], Re_t=Reynolds_Lam_to_Turb)

            # Inertial deposition: contraction
            if Flow_Properties[Section][19] < Flow_Properties[Section][16]:
                Particle_Variables[i + 1][19] = KN.Inertial_Deposit_Contraction_Eff(Stk=Particle_Variables[i + 1][13], Theta_Rad=Flow_Properties[Section][18], A_o=Flow_Properties[Section][19], A_i=Flow_Properties[Section][16])
            else:
                Particle_Variables[i + 1][19] = 1

            # Check if this section is probe or not
            if General_Input[Section][18] == "Y":
                Particle_Variables[i + 1][7] = 1
                Particle_Variables[i + 1][11] = 1
                Particle_Variables[i + 1][15] = 1
                Particle_Variables[i + 1][17] = 1
                Particle_Variables[i + 1][18] = 1
                Particle_Variables[i + 1][19] = 1
                # INLET ASPIRATION penetration
                Particle_Variables[i + 1][20] = KN.Aspiration_Eff(Stk=Particle_Variables[i + 1][13], U0=General_Input[Section][19], U=Flow_Properties[Section][7], theta_Rad=Flow_Properties[Section][20])
                # INLET Transportation penetration
                Particle_Variables[i + 1][21] = KN.InletTrans_Eff(Stk=Particle_Variables[i + 1][13], U0=General_Input[Section][19], U=Flow_Properties[Section][7], theta_Rad=Flow_Properties[Section][20], Zg=Particle_Variables[i + 1][10],
                                                                  Re=Flow_Properties[Section][12])
            else:
                Particle_Variables[i + 1][20] = 1
                Particle_Variables[i + 1][21] = 1

            # Overall penetration
            Particle_Variables[i + 1][22] = Particle_Variables[i + 1][21] * Particle_Variables[i + 1][20] * Particle_Variables[i + 1][19] * Particle_Variables[i + 1][18] * Particle_Variables[i + 1][17] * Particle_Variables[i + 1][15] * \
                                            Particle_Variables[i + 1][11] * Particle_Variables[i + 1][7]

        # Saving Efficiencies
        for I in range(Nd):
            if Flow_Properties[Section][0] != 0:  # ignoring section
                eta_Diff[Section][I] = Flow_Properties[Section][0]
                eta_Grav[Section][I] = Flow_Properties[Section][0]
                eta_Inert[Section][I] = Flow_Properties[Section][0]
                eta_TP[Section][I] = Flow_Properties[Section][0]
                eta_Bend[Section][I] = Flow_Properties[Section][0]
                eta_Cont[Section][I] = Flow_Properties[Section][0]
                eta_Asp[Section][I] = Flow_Properties[Section][0]
                eta_Inlet[Section][I] = Flow_Properties[Section][0]
                e_Total[Section][I] = Flow_Properties[Section][0]
            else:
                eta_Diff[Section][I] = Particle_Variables[I + 1][7]
                eta_Grav[Section][I] = Particle_Variables[I + 1][11]
                eta_Inert[Section][I] = Particle_Variables[I + 1][15]
                eta_TP[Section][I] = Particle_Variables[I + 1][17]
                eta_Bend[Section][I] = Particle_Variables[I + 1][18]
                eta_Cont[Section][I] = Particle_Variables[I + 1][19]
                eta_Asp[Section][I] = Particle_Variables[I + 1][20]
                eta_Inlet[Section][I] = Particle_Variables[I + 1][21]
                e_Total[Section][I] = Particle_Variables[I + 1][22]
        logging.info(f"Calculating Particle variables were successful: S {Section}")
        ##########################################################################
        # Fitting polynomial on the total penetration for each section
        Z = np.polyfit(Ln_Diam1, e_Total[Section][0:Nd], Polynomial_Deg)
        Z1 = np.poly1d(Z)
        Fit_Expression = ""
        for Deg in range(Polynomial_Deg):
            K = Polynomial_Deg - Deg
            Fitting_Coefficients[Section][Deg] = Z[Deg]
            Fit_Expression += "{:.3E}".format(Z[Deg]) + " ." + "(" + "ln (dm)" + ")" + "$^{}$".format(K) + "+" + "\n"
        Fit_Expression += "{:.3E}".format(Z[Polynomial_Deg])
        Fit_Expression += "\n" + "(dm in meter)"
        Fitting_Coefficients[Section][Polynomial_Deg] = Z[Polynomial_Deg]
        ##########################################################################
        #  Saving Data for the Section
        Output_Data_Excel_Section = "Section_" + str(Section) + ".xlsx"
        Output_Data_Excel_Section = os.path.join(Excel_Directory, Output_Data_Excel_Section)
        workbook2 = xlsxwriter.Workbook(Output_Data_Excel_Section)
        P_Var = workbook2.add_worksheet("Particle Variables S_" + str(Section))
        for I in range(0, Nd + 1):
            for J in range(Number_Particle_Var):
                P_Var.write(I, J, Particle_Variables[I][J])
        for I in range(0, Polynomial_Deg + 1):
            J = Number_Particle_Var
            if I == 0:
                P_Var.write(I, J, "Fitting Coefficients")
            P_Var.write(I + 1, J, Fitting_Coefficients[Section][I])
        workbook2.close()
        ##########################################################################
        # Figures properties
        plt.figure(Section)
        plt.rcParams['mathtext.fontset'] = 'stix'
        plt.rcParams['font.family'] = 'STIXGeneral'
        if Diffusion_Enable == 1:
            plt.plot(Diameter_Nano, eta_Diff[Section][:], 'mD--', label='Diffusion', alpha=0.8, markersize=4)
        if Gravitational_Enable == 1:
            plt.plot(Diameter_Nano, eta_Grav[Section][:], 'gs:', label='Gravitational' + ', Re= ' + '{:.0f}'.format(Flow_Properties[Section][12]), alpha=0.7, markersize=4)
        if Inertial_Enable == 1:
            plt.plot(Diameter_Nano, eta_Inert[Section][:], 'co--', label='Inertial', alpha=0.6, markersize=3)
        if Thermophoretic_Enable == 1:
            plt.plot(Diameter_Nano, eta_TP[Section][:], 'k+:', label='Thermophoresis', alpha=0.5, markersize=3)
        if Bend_Enable == 1:
            plt.plot(Diameter_Nano, eta_Bend[Section][:], 'yx-', label='Bend', alpha=0.4, markersize=3)
        if Contr_Enable == 1:
            plt.plot(Diameter_Nano, eta_Cont[Section][:], 'm1-', label='Contraction', alpha=0.3, markersize=2)
        if Asp_Enable == 1:
            plt.plot(Diameter_Nano, eta_Asp[Section][:], 'cp-', label='Aspiration', alpha=0.2, markersize=1)
        if Inlet_Enable == 1:
            plt.plot(Diameter_Nano, eta_Inlet[Section][:], 'bs-', label='Probe Inlet', alpha=0.2, markersize=1)
        if Total_Enable == 1:
            plt.plot(Diameter_Nano, e_Total[Section][0:Nd], 'r*-', label='Total Penetration', alpha=1, markersize=2)
        if Fit_Enable == 1:
            plt.plot(Diameter_Nano, Z1(Ln_Diam1), 'b--', label=Fit_Expression, alpha=0.8, markersize=2)
        plt.title('Particle Penetration Coefficients in Section ' + str(Section))
        plt.xlabel('Mobility-equivalent Diameter (nm)')
        plt.xscale('log')
        plt.grid(True, which='minor', alpha=0.2)
        plt.ylabel('Penetration')
        plt.ylim((Y_limit_lower, Y_limit_upper))
        plt.yticks(np.arange(Y_limit_lower, Y_limit_upper, Y_Axis_Unit))
        plt.grid(True)
        plt.legend(title='Loss Mechanisms', bbox_to_anchor=(1.04, 0.5), loc='center left', fancybox=True, fontsize='small')

        # Saving figures
        name = "Section_" + str(Section) + ".jpg"
        Graph_Output_2 = os.path.join(Graph_Directory, name)
        plt.savefig(Graph_Output_2, format='jpg', dpi=Figure_DPI, bbox_inches='tight')
        plt.clf()
        plt.close()
    ##########################################################################
    # Calculating total penetration
    for Section in range(1, Input_Dimension_Row):
        for i in range(Nd):
            eta_Diff[Input_Dimension_Row][i] = eta_Diff[Section][i] * eta_Diff[Input_Dimension_Row][i]
            eta_Grav[Input_Dimension_Row][i] = eta_Grav[Section][i] * eta_Grav[Input_Dimension_Row][i]
            eta_Inert[Input_Dimension_Row][i] = eta_Inert[Section][i] * eta_Inert[Input_Dimension_Row][i]
            eta_TP[Input_Dimension_Row][i] = eta_TP[Section][i] * eta_TP[Input_Dimension_Row][i]
            eta_Bend[Input_Dimension_Row][i] = eta_Bend[Section][i] * eta_Bend[Input_Dimension_Row][i]
            eta_Cont[Input_Dimension_Row][i] = eta_Cont[Section][i] * eta_Cont[Input_Dimension_Row][i]
            eta_Asp[Input_Dimension_Row][i] = eta_Asp[Section][i] * eta_Asp[Input_Dimension_Row][i]
            eta_Inlet[Input_Dimension_Row][i] = eta_Inlet[Section][i] * eta_Inlet[Input_Dimension_Row][i]
            e_Total[Input_Dimension_Row][i] = e_Total[Section][i] * e_Total[Input_Dimension_Row][i]
    logging.info(f"Total Penetration calculated!")
    ##########################################################################
    # Fitting polynomial on the total penetration
    Z = np.polyfit(Ln_Diam1, e_Total[Input_Dimension_Row][0:Nd], Polynomial_Deg)
    Z1 = np.poly1d(Z)
    Fit_Expression = ""
    for Deg in range(Polynomial_Deg):
        K = Polynomial_Deg - Deg
        Fitting_Coefficients[Input_Dimension_Row][Deg] = Z[Deg]
        Fit_Expression += "{:.3E}".format(Z[Deg]) + " ." + "x" + "$^{}$".format(K) + "+" + "\n"
    Fit_Expression += "{:.3E}".format(Z[Polynomial_Deg])
    Fit_Expression += "\n" + "(x=ln(dm))"
    Fitting_Coefficients[Input_Dimension_Row][Polynomial_Deg] = Z[Polynomial_Deg]
    ##########################################################################
    for Section in range(1, Input_Dimension_Row + 1):
        for i in range(Nd):
            eta_Diff[Section][i] = 100 * eta_Diff[Section][i]
            eta_Grav[Section][i] = 100 * eta_Grav[Section][i]
            eta_Inert[Section][i] = 100 * eta_Inert[Section][i]
            eta_TP[Section][i] = 100 * eta_TP[Section][i]
            eta_Bend[Section][i] = 100 * eta_Bend[Section][i]
            eta_Cont[Section][i] = 100 * eta_Cont[Section][i]
            eta_Asp[Section][i] = 100 * eta_Asp[Section][i]
            eta_Inlet[Section][i] = 100 * eta_Inlet[Section][i]
            e_Total[Section][i] = 100 * e_Total[Section][i]

            los_Diff[Section][i] = (100 - eta_Diff[Section][i])
            los_Grav[Section][i] = (100 - eta_Grav[Section][i])
            los_Inert[Section][i] = (100 - eta_Inert[Section][i])
            los_TP[Section][i] = (100 - eta_TP[Section][i])
            los_Bend[Section][i] = (100 - eta_Bend[Section][i])
            los_Cont[Section][i] = (100 - eta_Cont[Section][i])
            los_Asp[Section][i] = (100 - eta_Asp[Section][i])
            los_Inlet[Section][i] = (100 - eta_Inlet[Section][i])
            los_Total[Section][i] = (100 - e_Total[Section][i])

    for i in range(Nd):
        los_Norm_Sum[i] = los_Diff[Input_Dimension_Row][i] + los_Grav[Input_Dimension_Row][i] + los_Inert[Input_Dimension_Row][i] + los_TP[Input_Dimension_Row][i] + los_Bend[Input_Dimension_Row][i] + los_Cont[Input_Dimension_Row][i] + \
                          los_Asp[Input_Dimension_Row][i] + los_Inlet[Input_Dimension_Row][i]
        los_Diff_Norm[i] = los_Diff[Input_Dimension_Row][i] / (los_Norm_Sum[i] / los_Total[Input_Dimension_Row][i])
        los_Grav_Norm[i] = los_Grav[Input_Dimension_Row][i] / (los_Norm_Sum[i] / los_Total[Input_Dimension_Row][i])
        los_Inert_Norm[i] = los_Inert[Input_Dimension_Row][i] / (los_Norm_Sum[i] / los_Total[Input_Dimension_Row][i])
        los_TP_Norm[i] = los_TP[Input_Dimension_Row][i] / (los_Norm_Sum[i] / los_Total[Input_Dimension_Row][i])
        los_Bend_Norm[i] = los_Bend[Input_Dimension_Row][i] / (los_Norm_Sum[i] / los_Total[Input_Dimension_Row][i])
        los_Cont_Norm[i] = los_Cont[Input_Dimension_Row][i] / (los_Norm_Sum[i] / los_Total[Input_Dimension_Row][i])
        los_Asp_Norm[i] = los_Asp[Input_Dimension_Row][i] / (los_Norm_Sum[i] / los_Total[Input_Dimension_Row][i])
        los_Inlet_Norm[i] = los_Inlet[Input_Dimension_Row][i] / (los_Norm_Sum[i] / los_Total[Input_Dimension_Row][i])

    ##########################################################################
    # Total Penetration Graph
    plt.figure(Input_Dimension_Row)
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'
    if Diffusion_Enable == 1:
        plt.plot(Diameter_Nano, eta_Diff[Input_Dimension_Row][:], 'mD--', label='Diffusion', alpha=0.8, markersize=4)
    if Gravitational_Enable == 1:
        plt.plot(Diameter_Nano, eta_Grav[Input_Dimension_Row][:], 'gs:', label='Gravitational', alpha=0.7, markersize=4)
    if Inertial_Enable == 1:
        plt.plot(Diameter_Nano, eta_Inert[Input_Dimension_Row][:], 'co--', label='Inertial', alpha=0.6, markersize=3)
    if Thermophoretic_Enable == 1:
        plt.plot(Diameter_Nano, eta_TP[Input_Dimension_Row][:], 'k+:', label='Thermophoresis', alpha=0.5, markersize=3)
    if Bend_Enable == 1:
        plt.plot(Diameter_Nano, eta_Bend[Input_Dimension_Row][:], 'yx-', label='Bend', alpha=0.4, markersize=3)
    if Contr_Enable == 1:
        plt.plot(Diameter_Nano, eta_Cont[Input_Dimension_Row][:], 'm1-', label='Contraction', alpha=0.3, markersize=2)
    if Asp_Enable == 1:
        plt.plot(Diameter_Nano, eta_Asp[Input_Dimension_Row][:], 'cp-', label='Aspiration', alpha=0.2, markersize=1)
    if Inlet_Enable == 1:
        plt.plot(Diameter_Nano, eta_Inlet[Input_Dimension_Row][:], 'bs-', label='Probe Inlet', alpha=0.2, markersize=1)
    if Total_Enable == 1:
        plt.plot(Diameter_Nano, e_Total[Input_Dimension_Row][0:Nd], 'r*-', label='Total Penetration', alpha=1, markersize=2)
    if Fit_Enable == 1:
        Z2 = Z1(Ln_Diam1) * 100
        plt.plot(Diameter_Nano, Z2, 'b--', label=Fit_Expression, alpha=0.8, markersize=2)
    plt.title('Particle Penetration Coefficients in ' + General_Input[Input_Dimension_Row - 1][1])
    plt.xlabel('Mobility-equivalent Diameter (nm)')
    plt.xscale('log')
    plt.grid(True, which='minor', alpha=0.2)
    plt.ylabel('Penetration(%)')
    plt.ylim((0, 105))
    plt.yticks(np.arange(0, 105, 10))
    plt.grid(True)
    plt.legend(title='Loss Mechanisms', bbox_to_anchor=(1.04, 0.5), loc='center left', fancybox=True, fontsize='small')

    # Saving Figures
    name = General_Input[Input_Dimension_Row - 1][1] + " - " + "Main Graph.jpg"
    Graph_Output_2 = os.path.join(Graph_Directory, name)
    plt.savefig(Graph_Output_2, format='jpg', dpi=Figure_DPI, bbox_inches='tight')
    plt.clf()
    plt.close()
    ##########################################################################
    # Bar plot
    SBG = StackedBarGrapher()
    Total_Losses = np.array([los_Diff_Norm[:], los_Grav_Norm[:], los_Inert_Norm[:], los_TP_Norm[:], los_Bend_Norm[:], los_Cont_Norm[:], los_Asp_Norm[:], los_Inlet_Norm[:]])
    Total_Losses = Total_Losses.transpose()
    Total_Losses = np.around(Total_Losses, decimals=1)
    Gender = ['Diffusion', 'Gravitational', 'Inertial', 'Thermophoresis', 'Bend', 'Contraction', 'Aspiration', 'Probe Inlet']
    # colors
    d_colors = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71", "#ffbf00", "#000000"]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    SBG.stackedBarPlot(ax, Total_Losses, d_colors, xLabels=D_labels, widths=Width_nm)
    plt.title('Particle loss Coefficients in ' + General_Input[Input_Dimension_Row - 1][1])
    plt.grid(True, which='major', alpha=0.1)
    plt.ylabel('Particle loss(%)')
    plt.xlabel('Mobility-equivalent Diameter (nm)')
    plt.grid(True)
    plt.legend(Gender, bbox_to_anchor=(1.04, 0.5), loc='center left', fancybox=True, fontsize='small')
    plt.tight_layout()
    name = General_Input[Input_Dimension_Row - 1][1] + " - " + "Main Graph_Loss.jpg"
    Graph_Output_2 = os.path.join(Graph_Directory, name)
    plt.savefig(Graph_Output_2, format='jpg', dpi=Figure_DPI, bbox_inches='tight')
    plt.clf()
    plt.close()
    ##########################################################################
    Shift = Computed_data_Column + Input_Dimension_Col + 2 - 1
    # Saving main report Data (Excel)
    Main_Excel_WorkBook = xlsxwriter.Workbook(Output_Data_Excel)
    # Diffusion Sheet
    ws_diff = Main_Excel_WorkBook.add_worksheet("Diffusion")
    for Section in range(0, Input_Dimension_Row):
        for i in range(Input_Dimension_Col):
            ws_diff.write(Section, i, General_Input[Section][i])
        for i in range(1, Computed_data_Column):
            ws_diff.write(Section, i + Input_Dimension_Col, Flow_Properties[Section][i])
        for i in range(Nd):
            if Section == 0:
                ws_diff.write(0, i + Shift, diam[i])
            ws_diff.write(Section + 1, i + Shift, eta_Diff[Section + 1][i])
    # Gravitational Sheet
    ws_grav = Main_Excel_WorkBook.add_worksheet("Gravity")
    for Section in range(0, Input_Dimension_Row):
        for i in range(Input_Dimension_Col):
            ws_grav.write(Section, i, General_Input[Section][i])
        for i in range(1, Computed_data_Column):
            ws_grav.write(Section, i + Input_Dimension_Col, Flow_Properties[Section][i])
        for i in range(Nd):
            if Section == 0:
                ws_grav.write(0, i + Shift, diam[i])
            ws_grav.write(Section + 1, i + Shift, eta_Grav[Section + 1][i])
    # Inertia Sheet
    ws_inert = Main_Excel_WorkBook.add_worksheet("Inertia")
    for Section in range(0, Input_Dimension_Row):
        for i in range(Input_Dimension_Col):
            ws_inert.write(Section, i, General_Input[Section][i])
        for i in range(1, Computed_data_Column):
            ws_inert.write(Section, i + Input_Dimension_Col, Flow_Properties[Section][i])
        for i in range(Nd):
            if Section == 0:
                ws_inert.write(0, i + Shift, diam[i])
            ws_inert.write(Section + 1, i + Shift, eta_Inert[Section + 1][i])
    # Thermophoresis Sheet
    ws_tp = Main_Excel_WorkBook.add_worksheet("Thermophoresis")
    for Section in range(0, Input_Dimension_Row):
        for i in range(Input_Dimension_Col):
            ws_tp.write(Section, i, General_Input[Section][i])
        for i in range(1, Computed_data_Column):
            ws_tp.write(Section, i + Input_Dimension_Col, Flow_Properties[Section][i])
        for i in range(Nd):
            if Section == 0:
                ws_tp.write(0, i + Shift, diam[i])
            ws_tp.write(Section + 1, i + Shift, eta_TP[Section + 1][i])
    # Bend Sheet
    ws_be = Main_Excel_WorkBook.add_worksheet("Bend")
    for Section in range(0, Input_Dimension_Row):
        for i in range(Input_Dimension_Col):
            ws_be.write(Section, i, General_Input[Section][i])
        for i in range(1, Computed_data_Column):
            ws_be.write(Section, i + Input_Dimension_Col, Flow_Properties[Section][i])
        for i in range(Nd):
            if Section == 0:
                ws_be.write(0, i + Shift, diam[i])
            ws_be.write(Section + 1, i + Shift, eta_Bend[Section + 1][i])
    # Contraction Sheet
    ws_cont = Main_Excel_WorkBook.add_worksheet("Contraction")
    for Section in range(0, Input_Dimension_Row):
        for i in range(Input_Dimension_Col):
            ws_cont.write(Section, i, General_Input[Section][i])
        for i in range(1, Computed_data_Column):
            ws_cont.write(Section, i + Input_Dimension_Col, Flow_Properties[Section][i])
        for i in range(Nd):
            if Section == 0:
                ws_cont.write(0, i + Shift, diam[i])
            ws_cont.write(Section + 1, i + Shift, eta_Cont[Section + 1][i])
    # Aspiration Sheet
    ws_asp = Main_Excel_WorkBook.add_worksheet("Aspiration")
    for Section in range(0, Input_Dimension_Row):
        for i in range(Input_Dimension_Col):
            ws_asp.write(Section, i, General_Input[Section][i])
        for i in range(1, Computed_data_Column):
            ws_asp.write(Section, i + Input_Dimension_Col, Flow_Properties[Section][i])
        for i in range(Nd):
            if Section == 0:
                ws_asp.write(0, i + Shift, diam[i])
            ws_asp.write(Section + 1, i + Shift, eta_Asp[Section + 1][i])
    # Inlet efficiency sheet
    ws_Inl = Main_Excel_WorkBook.add_worksheet("Inlet")
    for Section in range(0, Input_Dimension_Row):
        for i in range(Input_Dimension_Col):
            ws_Inl.write(Section, i, General_Input[Section][i])
        for i in range(1, Computed_data_Column):
            ws_Inl.write(Section, i + Input_Dimension_Col, Flow_Properties[Section][i])
        for i in range(Nd):
            if Section == 0:
                ws_Inl.write(0, i + Shift, diam[i])
            ws_Inl.write(Section + 1, i + Shift, eta_Inlet[Section + 1][i])
    # Total penetration Sheet
    ws_total = Main_Excel_WorkBook.add_worksheet("Total")
    for Section in range(0, Input_Dimension_Row):
        for i in range(Input_Dimension_Col):
            ws_total.write(Section, i, General_Input[Section][i])
        for i in range(1, Computed_data_Column):
            ws_total.write(Section, i + Input_Dimension_Col, Flow_Properties[Section][i])
        for i in range(Nd):
            if Section == 0:
                ws_total.write(0, i + Shift, diam[i])
            ws_total.write(Section + 1, i + Shift, e_Total[Section + 1][i])
        for I in range(0, Polynomial_Deg + 1):
            J = Nd + Shift
            if I == 0:
                ws_total.write(I, J, "Fitting Coefficients")
            ws_total.write(I + 1, J, Fitting_Coefficients[Input_Dimension_Row][I])

    # Diffusion loss Sheet
    ws_diff1 = Main_Excel_WorkBook.add_worksheet("Diffusion loss")
    for Section in range(0, Input_Dimension_Row):
        for i in range(Input_Dimension_Col):
            ws_diff1.write(Section, i, General_Input[Section][i])
        for i in range(1, Computed_data_Column):
            ws_diff1.write(Section, i + Input_Dimension_Col, Flow_Properties[Section][i])
        for i in range(Nd):
            if Section == 0:
                ws_diff1.write(0, i + Shift, diam[i])
            ws_diff1.write(Section + 1, i + Shift, los_Diff[Section + 1][i])
    # Gravitational loss Sheet
    ws_grav1 = Main_Excel_WorkBook.add_worksheet("Gravity loss")
    for Section in range(0, Input_Dimension_Row):
        for i in range(Input_Dimension_Col):
            ws_grav1.write(Section, i, General_Input[Section][i])
        for i in range(1, Computed_data_Column):
            ws_grav1.write(Section, i + Input_Dimension_Col, Flow_Properties[Section][i])
        for i in range(Nd):
            if Section == 0:
                ws_grav1.write(0, i + Shift, diam[i])
            ws_grav1.write(Section + 1, i + Shift, los_Grav[Section + 1][i])
    # Inertia loss Sheet
    ws_inert1 = Main_Excel_WorkBook.add_worksheet("Inertia loss")
    for Section in range(0, Input_Dimension_Row):
        for i in range(Input_Dimension_Col):
            ws_inert1.write(Section, i, General_Input[Section][i])
        for i in range(1, Computed_data_Column):
            ws_inert1.write(Section, i + Input_Dimension_Col, Flow_Properties[Section][i])
        for i in range(Nd):
            if Section == 0:
                ws_inert1.write(0, i + Shift, diam[i])
            ws_inert1.write(Section + 1, i + Shift, los_Inert[Section + 1][i])
    # Thermophoresis loss Sheet
    ws_tp1 = Main_Excel_WorkBook.add_worksheet("Thermophoresis loss")
    for Section in range(0, Input_Dimension_Row):
        for i in range(Input_Dimension_Col):
            ws_tp1.write(Section, i, General_Input[Section][i])
        for i in range(1, Computed_data_Column):
            ws_tp1.write(Section, i + Input_Dimension_Col, Flow_Properties[Section][i])
        for i in range(Nd):
            if Section == 0:
                ws_tp1.write(0, i + Shift, diam[i])
            ws_tp1.write(Section + 1, i + Shift, los_TP[Section + 1][i])
    # Bend loss Sheet
    ws_be1 = Main_Excel_WorkBook.add_worksheet("Bend loss")
    for Section in range(0, Input_Dimension_Row):
        for i in range(Input_Dimension_Col):
            ws_be1.write(Section, i, General_Input[Section][i])
        for i in range(1, Computed_data_Column):
            ws_be1.write(Section, i + Input_Dimension_Col, Flow_Properties[Section][i])
        for i in range(Nd):
            if Section == 0:
                ws_be1.write(0, i + Shift, diam[i])
            ws_be1.write(Section + 1, i + Shift, los_Bend[Section + 1][i])
    # Contraction loss Sheet
    ws_cont1 = Main_Excel_WorkBook.add_worksheet("Contraction loss")
    for Section in range(0, Input_Dimension_Row):
        for i in range(Input_Dimension_Col):
            ws_cont1.write(Section, i, General_Input[Section][i])
        for i in range(1, Computed_data_Column):
            ws_cont1.write(Section, i + Input_Dimension_Col, Flow_Properties[Section][i])
        for i in range(Nd):
            if Section == 0:
                ws_cont1.write(0, i + Shift, diam[i])
            ws_cont1.write(Section + 1, i + Shift, los_Cont[Section + 1][i])
    # Aspiration loss Sheet
    ws_asp1 = Main_Excel_WorkBook.add_worksheet("Aspiration loss")
    for Section in range(0, Input_Dimension_Row):
        for i in range(Input_Dimension_Col):
            ws_asp1.write(Section, i, General_Input[Section][i])
        for i in range(1, Computed_data_Column):
            ws_asp1.write(Section, i + Input_Dimension_Col, Flow_Properties[Section][i])
        for i in range(Nd):
            if Section == 0:
                ws_asp1.write(0, i + Shift, diam[i])
            ws_asp1.write(Section + 1, i + Shift, los_Asp[Section + 1][i])
    # Inlet loss sheet
    ws_Inl1 = Main_Excel_WorkBook.add_worksheet("Inlet loss")
    for Section in range(0, Input_Dimension_Row):
        for i in range(Input_Dimension_Col):
            ws_Inl1.write(Section, i, General_Input[Section][i])
        for i in range(1, Computed_data_Column):
            ws_Inl1.write(Section, i + Input_Dimension_Col, Flow_Properties[Section][i])
        for i in range(Nd):
            if Section == 0:
                ws_Inl1.write(0, i + Shift, diam[i])
            ws_Inl1.write(Section + 1, i + Shift, los_Inlet[Section + 1][i])
    # Total loss Sheet
    ws_total1 = Main_Excel_WorkBook.add_worksheet("Total loss")
    for Section in range(0, Input_Dimension_Row):
        for i in range(Input_Dimension_Col):
            ws_total1.write(Section, i, General_Input[Section][i])
        for i in range(1, Computed_data_Column):
            ws_total1.write(Section, i + Input_Dimension_Col, Flow_Properties[Section][i])
        for i in range(Nd):
            if Section == 0:
                ws_total1.write(0, i + Shift, diam[i])
            ws_total1.write(Section + 1, i + Shift, los_Total[Section + 1][i])

    Main_Excel_WorkBook.close()
