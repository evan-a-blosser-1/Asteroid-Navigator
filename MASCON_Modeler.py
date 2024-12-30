"""
Program: Asteroid_MASCON_Modeler.py
Description: Calculates the Center of mass for selected
                Models:
                        - MASCON I
                        - MASCON III
                        - MASCON VIII
            

MIT License

Copyright (c) [2024] [Evan Blosser]

"""
#############################
import os
import sys
import time
import threading
import trimesh
from sys import exit
# Mathmatical!!                     
import numpy as np               
import pandas as pd
# Plot 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from PySide6.QtWidgets import QVBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
########################################################
#################################### Personal Packages #    
sys.dont_write_bytecode = True
import OBJ_Pack as OBJ
#################################################################################
#################### Settings 
# Big G in km^3/kg.s^2
G = 6.67430e-20 
# # Scaling Settings
# rel_tol = 1e-09
# abs_tol = 1e-09
# Max_Iter = 100000
# # Mu calc settings 
# Max_Iter2 = 10000
# rel_tol2  = 1e-14
# abs_tol2 = 1e-14
# asteroid = '1950DA_Prograde'
# # 1950DA
# Vol_Accepted  = 0.797895   # km^3  Prograde   0.82 km^3 +/- 30% ( = 0.24599999999999997)  
#                            #       Retrograde 1.14 km^3 +/- 30%
# Mass_Accepted = 2.1e12     # 2.1 \times 10^{12} \pm 1.1~\mathrm{kg}
# Density       = 2.4e12     # kg/km^3    1.7\pm 0.7~\mathrm{g~cm}^{-3}
# ############################
def Precision_Delta(Perr_Diff):
        Delta = 0.01
        # print(f"Perr_Diff: {Perr_Diff}")
        if abs(Perr_Diff) > 1.0:
            Delta = 0.1
        if abs(Perr_Diff) < 1.0:
            Delta = 0.001
        if abs(Perr_Diff) < 0.5:
            Delta = 0.0001
        if abs(Perr_Diff) < 1e-2:
            Delta = 0.00001
        if abs(Perr_Diff) < 1e-3:
            Delta = 0.000001
        if abs(Perr_Diff) < 1e-4:
            Delta = 0.0000001
        if abs(Perr_Diff) < 1e-5:
            Delta = 0.00000001
        if abs(Perr_Diff) < 1e-6:
            Delta = 0.000000001
        if abs(Perr_Diff) < 1e-7:
            Delta = 0.0000000001
        if abs(Perr_Diff) < 1e-8:
            Delta = 0.00000000001
        if abs(Perr_Diff) < 1e-9:
            Delta = 0.000000000001
        if abs(Perr_Diff) < 1e-10:
            Delta = 0.0000000000001
        if abs(Perr_Diff) < 1e-11:
            Delta = 0.00000000000001
        if abs(Perr_Diff) < 1e-12:
            Delta = 0.000000000000001
        if abs(Perr_Diff) < 1e-13:
            Delta = 0.0000000000000001
        if abs(Perr_Diff) < 1e-14:
            Delta = 0.00000000000000001
        if abs(Perr_Diff) < 1e-15:
            Delta = 0.000000000000000001
        if abs(Perr_Diff) < 1e-16:
            Delta = 0.0000000000000000001
        if abs(Perr_Diff) < 1e-17:
            Delta = 0.00000000000000000001
        if abs(Perr_Diff) < 1e-18:
            Delta = 0.000000000000000000001
            
        return Delta
###################################
#################################
def MASCON_Modeler(Model_Save_Path,running_log_callback,log_callback,file_2_model,Spin_Rate, Density, Vol_Accepted, Mass_Accepted, 
                 rel_tol= 1e-09, abs_tol= 1e-09, Max_Iter = 100000,
                 rel_tol2= 1e-14,abs_tol2= 1e-14, Max_Iter2= 10000):
    # OBJ File                        
    Asteroid_Save = os.path.splitext(file_2_model)[0]
    
    #################################################################
    ################################################### Constants ###
    mu_const      = G*Mass_Accepted
    ###################################### Assign Vertices & Faces ##
    Verts_IN, Faces =  OBJ.OBJ_2_VertFace(file_2_model)
    Faces = Faces - 1
    ######################################################################
    ###################### trimesh Calculation ###########################
    mesh = trimesh.load(file_2_model)
    CM = mesh.center_mass

    polyhedron_volume = mesh.volume
    R_eff = (3 * polyhedron_volume / (4 * np.pi)) ** (1/3)
    # 
    tetra_count=np.shape(Faces)[0]
    # Use the output in Python
    log_callback(f"Poly CM: {CM}")
    #############################################################
    ############################################################
    #%% Scaling 
    running_log_callback("| Calculating Scaling Factor...")
    ####################################################################
    # Look at decimal function
    # x = np.zeros(3, dtpye="float64")
    # or np.array definition for more accurate value 
    #
    #
    Gamma_Coeff = 1
    for it_1 in range(Max_Iter):
        # Check % and decrease del_gamma?
        #
        Verts = Verts_IN*Gamma_Coeff
        #
        tetra_vol, total_vol = OBJ.Tetra_Volume(Verts, Faces, MASCON_Div=1,total_vol=0)
        ##############################
        Perr_Diff = (total_vol - Vol_Accepted)/Vol_Accepted

        Del_Gamma = Precision_Delta(Perr_Diff)    
        ##############################
        if total_vol > Vol_Accepted:
            Gamma_Coeff -= Del_Gamma
            # print(f"Gamma Coeff: {Gamma_Coeff}")
        ##############################
        elif total_vol < Vol_Accepted:
            Gamma_Coeff += Del_Gamma
            # print(f"Gamma Coeff: {Gamma_Coeff}")
        ##############################
        if np.isclose(total_vol,Vol_Accepted, rtol=rel_tol, atol=abs_tol):
            Converge_Message = f""" 
    |  Accepted Volume: {Vol_Accepted}
    | 
    | Calculated Volume: {total_vol}
    |
    | Gamma Coefficient: {Gamma_Coeff}
    |
    | Converged at Iteration: {it_1}        
            """
            log_callback(Converge_Message)    
            break
        ##############################
    else:
        Out_Message = f""" 
    |  Accepted Volume: {Vol_Accepted}
    | 
    | Calculated Volume: {total_vol}
    |
    | Gamma Coefficient: {Gamma_Coeff}
    |
    | Max Iterations Reached: {it_1 +1}        
        """
        log_callback(Out_Message) 
    #############################################



    #%% Volume Function Call
    #
    # Note: All tetra and Prism numbers are from outer layer 
    #       to inner from 1 to number 'n' respectively  
    #   
    #       i.e. 1 is outer, 2 is second layer, 3 is inner tetra in MIII
    ####################################################################
    ########################################################### MASCON I
    Full_Volume_Tetra_Array,Total_Volume_Out = OBJ.Tetra_Volume(Verts, Faces,MASCON_Div=1,total_vol=0)
    ####################################################################
    ######################################################### MASCON III
    ############# Tetrahedron Volumes ###############
    Volume_T_I   = OBJ.Tetra_Volume(Verts, Faces,MASCON_Div=1,total_vol=1)
    Volume_M3_T_II  = OBJ.Tetra_Volume(Verts, Faces,MASCON_Div=(2/3),total_vol=1)
    Volume_M3_T_III = OBJ.Tetra_Volume(Verts, Faces,MASCON_Div=(1/3),total_vol=1)
    ############## Prism Volumes ####################
    Vol_M3_Prism_I  =  Volume_T_I  - Volume_M3_T_II
    Vol_M3_Prism_II =  Volume_M3_T_II - Volume_M3_T_III
    # Check
    M1_Check = np.sum(Volume_T_I)
    M3_Check = np.sum(Volume_M3_T_III) + np.sum(Vol_M3_Prism_I) + np.sum(Vol_M3_Prism_II)
    ####################################################################
    ######################################################## MASCON VIII
    Volume_T_II   = OBJ.Tetra_Volume(Verts, Faces, MASCON_Div=(7/8),total_vol=1)
    Volume_T_III  = OBJ.Tetra_Volume(Verts, Faces, MASCON_Div=(6/8),total_vol=1)
    Volume_T_IV   = OBJ.Tetra_Volume(Verts, Faces, MASCON_Div=(5/8),total_vol=1)
    Volume_T_V    = OBJ.Tetra_Volume(Verts, Faces, MASCON_Div=(4/8),total_vol=1)
    Volume_T_VI   = OBJ.Tetra_Volume(Verts, Faces, MASCON_Div=(3/8),total_vol=1)
    Volume_T_VII  = OBJ.Tetra_Volume(Verts, Faces, MASCON_Div=(2/8),total_vol=1)
    Volume_T_VIII = OBJ.Tetra_Volume(Verts, Faces, MASCON_Div=(1/8),total_vol=1)
    #################################################
    ######################### Prism Volumes #########
    Vol_Prism_VII =  Volume_T_VII  - Volume_T_VIII
    Vol_Prism_VI  =  Volume_T_VI - Volume_T_VII 
    Vol_Prism_V   =  Volume_T_V  - Volume_T_VI
    Vol_Prism_IV  =  Volume_T_IV - Volume_T_V 
    Vol_Prism_III =  Volume_T_III  - Volume_T_IV
    Vol_Prism_II  =  Volume_T_II - Volume_T_III 
    Vol_Prism_I   =  Volume_T_I  - Volume_T_II
    # Check
    M8_Check = np.sum(Volume_T_VIII) + np.sum(Vol_Prism_VII) +\
                np.sum(Vol_Prism_VI) + np.sum(Vol_Prism_V) +\
                np.sum(Vol_Prism_IV) + np.sum(Vol_Prism_III) +\
                np.sum(Vol_Prism_II) + np.sum(Vol_Prism_I)
    ############################################################
    Volume_Check_Message = f"""
    {'-'*42}
    | Volume Check
    {'-'*42}
    |-----SUMS-----
    | MASCON I:    {M1_Check}
    | MASCON III:  {M3_Check}
    | MASCON VIII: {M8_Check}
    |
    | Total Volume: {Total_Volume_Out}
    {'-'*42}
    """
    running_log_callback(Volume_Check_Message)





    # Debug
        # if tt in [0,1,2,3,4,5]:
        #     print(f'Tetra: {Tetra_U[tt]}')
        #     print(U_vec)
        #     print(V_vec)
        #     print(W_vec)
        #     print(Vol_Full_Sum)
        #     print(A_Matrix)
    ########
    # x_coords_U = [tetra[0] for tetra in Tetra_U]
    # y_coords_U = [tetra[1] for tetra in Tetra_U]
    # z_coords_U = [tetra[2] for tetra in Tetra_U]
    # x_coords_V = [tetra[0] for tetra in Tetra_V]
    # y_coords_V = [tetra[1] for tetra in Tetra_V]
    # z_coords_V = [tetra[2] for tetra in Tetra_V]
    # x_coords_W = [tetra[0] for tetra in Tetra_W]
    # y_coords_W = [tetra[1] for tetra in Tetra_W]
    # z_coords_W = [tetra[2] for tetra in Tetra_W]
    # fig = plt.figure()
    # ax = fig.add_subplot(projection='3d')
    # ax.scatter(x_coords_U, y_coords_U, z_coords_U, c='r', marker='o', label='Tetra_U')
    # ax.scatter(x_coords_V, y_coords_V, z_coords_V, c='g', marker='^', label='Tetra_V')
    # ax.scatter(x_coords_W, y_coords_W, z_coords_W, c='b', marker='s', label='Tetra_W')
    # plt.show()
    ######################################################################################

    ###############################################
    #%% Grav Parameters 

    ####################################################################
    ####### Calculate Each Grav. Param. ################################
    #
    #
    Den = Density
    Del_Den = 0.01e12

    for it_2 in range(Max_Iter2):
        mu_i =  np.zeros(tetra_count, dtype="float64")
        Mass =  np.zeros(tetra_count, dtype="float64")
        mu_MIII = np.zeros((tetra_count, 3), dtype="float64")
        mu_MVIII = np.zeros((tetra_count, 8), dtype="float64")
        for i in range(tetra_count):
            ####################################################################
            ########################################################### MASCON I
            Mass[i] =     tetra_vol[i] * Den 
            mu_i[i] = G * tetra_vol[i] * Den
            ####################################################################
            ######################################################### MASCON III
            mu_MIII[i,0]  = G * Volume_M3_T_III[i] * Den
            mu_MIII[i,1]  = G * Vol_M3_Prism_II[i] * Den
            mu_MIII[i,2]  = G * Vol_M3_Prism_I[i]  * Den
            ####################################################################
            ######################################################## MASCON VIII
            mu_MVIII[i,0] = G * Volume_T_VIII[i] * Den
            mu_MVIII[i,1] = G * Vol_Prism_VII[i] * Den
            mu_MVIII[i,2] = G * Vol_Prism_VI[i]  * Den
            mu_MVIII[i,3] = G * Vol_Prism_V[i]   * Den
            mu_MVIII[i,4] = G * Vol_Prism_IV[i]  * Den
            mu_MVIII[i,5] = G * Vol_Prism_III[i] * Den
            mu_MVIII[i,6] = G * Vol_Prism_II[i]  * Den
            mu_MVIII[i,7] = G * Vol_Prism_I[i]   * Den
        ###############################################
        ######################################### Check
        mu = np.sum(mu_i)
        mu_3 = np.sum(mu_MIII)
        mu_8 = np.sum(mu_MVIII)    
        Mass_Check = np.sum(Mass)
        ##############################
        Perr_Diff_Mu = (mu - mu_const)/mu_const
        Del_Den = Precision_Delta(Perr_Diff)
        ##############################
        if mu > mu_const:
            Den -= Del_Den
            # print(f"Current Density (+): {Den}")
        ##############################
        elif mu < mu_const:
            Den += Del_Den
            # print(f"Current Density (-): {Den}")
        ##############################
        if np.isclose(mu, mu_const, rtol=rel_tol2, atol=abs_tol2):
            mu_converge_message = f"""
    {'-'*42}
    |Gravitational Param
    | MI:      {mu:.5e}
    | MIII:    {mu_3:.5e}
    | MVIII:   {mu_8:.5e}
    | True:    {mu_const:.5e}
    |Accepted Density: {Density:.5e}
    |The total Grav Parameter converged
    | in: {it_2} iterations
    | at:      {Den:.5e}
    | Density: {Density:.5e}
    {'-'*42}
    {'-'*42}
    |Accepted Volume: {Vol_Accepted}
    |Calculated Volume: {total_vol:.7f}
    {'-'*42}
    |Accepted Mass: {Mass_Accepted:.5e}
    |Calculated Mass: {Mass_Check:.5e}
    {'-'*42}
            """
            log_callback(mu_converge_message)
            break
    else:
        Vol_Mass_Den_Out_Message = f"""
    {'-'*42}
    | MAX Iterations Reached !
    {'-'*42}
    |Gravitational Param
    | MI:      {mu:.5e}
    | MIII:    {mu_3:.5e}
    | MVIII:   {mu_8:.5e}
    | True:    {mu_const:.5e}
    |The total Grav Parameter calculated
    | at: {Den:.5e}
    |Accepted Density: {Density:.5e}
    {'-'*42}
    {'-'*42}
    |Accepted Volume: {Vol_Accepted}
    |Calculated Volume: {total_vol:.7f}
    {'-'*42}
    |Accepted Mass: {Mass_Accepted:.5e}
    |Calculated Mass: {Mass_Check:.5e}
    {'-'*42}
        """ 
        log_callback(Vol_Mass_Den_Out_Message)

    #################################################


    #%% Tetra CM Calculation

    # Arrays 
    ####################################################################
    ########################################################### MASCON I
    Output_Array_M1 =  np.zeros((tetra_count, 3), dtype="float64")
    ####################################################################
    ######################################################### MASCON III
    CM_tetra_L1_M3 = np.zeros((tetra_count, 3), dtype="float64")
    CM_tetra_L2_M3 = np.zeros((tetra_count, 3), dtype="float64")
    CM_tetra_L3_M3 = np.zeros((tetra_count, 3), dtype="float64")
    Output_Array_MIII = np.zeros((tetra_count, 9), dtype="float64")
    ####################################################################
    ######################################################## MASCON VIII
    CM_tetra_M1 = np.empty((tetra_count,3))
    CM_tetra_M2 = np.empty((tetra_count,3))
    CM_tetra_M3 = np.empty((tetra_count,3))
    CM_tetra_M4 = np.empty((tetra_count,3))
    CM_tetra_M5 = np.empty((tetra_count,3))
    CM_tetra_M6 = np.empty((tetra_count,3))
    CM_tetra_M7 = np.empty((tetra_count,3))
    CM_tetra_M8 = np.empty((tetra_count,3))
    Output_Array_MVIII = np.zeros((tetra_count, 24), dtype="float64")
    #####################################################################
    ############## Center of Mass Calculations ##########################
    for it in range(0,tetra_count):
    ####################################################################
    ########################################################### MASCON I
        Center_mass_calc_x = (Verts[Faces[it,0],0] + Verts[Faces[it,1],0] + Verts[Faces[it,2],0] + CM[0])/4
        Center_mass_calc_y = (Verts[Faces[it,0],1] + Verts[Faces[it,1],1] + Verts[Faces[it,2],1] + CM[1])/4
        Center_mass_calc_z = (Verts[Faces[it,0],2] + Verts[Faces[it,1],2] + Verts[Faces[it,2],2] + CM[2])/4
        # Fill array
        Output_Array_M1[it] = (Center_mass_calc_x,Center_mass_calc_y,Center_mass_calc_z)
    ####################################################################
    ######################################################### MASCON III
        ############################################################################################### Tetra Center of Masses
        #######################################################################################################################
        ###################################################################################################### Center Layer ###
        CM_M3_calc_x3 = (Verts[Faces[it,0],0]*(1/3) + Verts[Faces[it,1],0]*(1/3) + Verts[Faces[it,2],0]*(1/3) + CM[0])/4
        CM_M3_calc_y3 = (Verts[Faces[it,0],1]*(1/3) + Verts[Faces[it,1],1]*(1/3) + Verts[Faces[it,2],1]*(1/3) + CM[1])/4
        CM_M3_calc_z3 = (Verts[Faces[it,0],2]*(1/3) + Verts[Faces[it,1],2]*(1/3) + Verts[Faces[it,2],2]*(1/3) + CM[2])/4
        # Fill array
        CM_tetra_L3_M3[it] =  (CM_M3_calc_x3,CM_M3_calc_y3,CM_M3_calc_z3)
        ######################################################################################################################
        #################################################################################################### Middle Layer ###
        # Tetrahedron CM
        CM_M3_Tetra_x2 = (Verts[Faces[it,0],0]*(2/3) + Verts[Faces[it,1],0]*(2/3) + Verts[Faces[it,2],0]*(2/3) + CM[0])/4
        CM_M3_Tetra_y2 = (Verts[Faces[it,0],1]*(2/3) + Verts[Faces[it,1],1]*(2/3) + Verts[Faces[it,2],1]*(2/3) + CM[1])/4
        CM_M3_Tetra_z2 = (Verts[Faces[it,0],2]*(2/3) + Verts[Faces[it,1],2]*(2/3) + Verts[Faces[it,2],2]*(2/3) + CM[2])/4
        # Fill array
        CM_tetra_L2_M3[it] =  (CM_M3_Tetra_x2,CM_M3_Tetra_y2,CM_M3_Tetra_z2)
        #####################################################################################################################
        ##################################################################################################### Outer Layer ###
        # Tetrahedron CM
        CM_M3_Tetra_x1 = (Verts[Faces[it,0],0] + Verts[Faces[it,1],0] + Verts[Faces[it,2],0] + CM[0])/4
        CM_M3_Tetra_y1 = (Verts[Faces[it,0],1] + Verts[Faces[it,1],1] + Verts[Faces[it,2],1] + CM[1])/4
        CM_M3_Tetra_z1 = (Verts[Faces[it,0],2] + Verts[Faces[it,1],2] + Verts[Faces[it,2],2] + CM[2])/4
        # Fill array
        CM_tetra_L1_M3[it] =  (CM_M3_Tetra_x1, CM_M3_Tetra_y1,CM_M3_Tetra_z1)
        ############################################################################################### Prism Center of Masses
        ######################################################################################################################
        #################################################################################################### Middle Layer ###
        # Prism CM
        CM_M3_Prism2_x2 = ( CM_tetra_L2_M3[it,0]*Volume_M3_T_II[it] - CM_tetra_L3_M3[it,0]*Volume_M3_T_III[it] )/Vol_M3_Prism_II[it]
        CM_M3_Prism2_y2 = ( CM_tetra_L2_M3[it,1]*Volume_M3_T_II[it] - CM_tetra_L3_M3[it,1]*Volume_M3_T_III[it] )/Vol_M3_Prism_II[it]
        CM_M3_Prism2_z2 = ( CM_tetra_L2_M3[it,2]*Volume_M3_T_II[it] - CM_tetra_L3_M3[it,2]*Volume_M3_T_III[it] )/Vol_M3_Prism_II[it]
        #####################################################################################################################
        ##################################################################################################### Outer Layer ###
        # Prism CM
        CM_M3_Prism1_x1 = ( CM_tetra_L1_M3[it,0]*Volume_T_I[it] - CM_tetra_L2_M3[it,0]*Volume_M3_T_II[it] )/Vol_M3_Prism_I[it]
        CM_M3_Prism1_y1 = ( CM_tetra_L1_M3[it,1]*Volume_T_I[it] - CM_tetra_L2_M3[it,1]*Volume_M3_T_II[it] )/Vol_M3_Prism_I[it]
        CM_M3_Prism1_z1 = ( CM_tetra_L1_M3[it,2]*Volume_T_I[it] - CM_tetra_L2_M3[it,2]*Volume_M3_T_II[it] )/Vol_M3_Prism_I[it]
        #################################################
        ### Fill Data Array 
        Output_Array_MIII[it] = (CM_M3_calc_x3,CM_M3_calc_y3,CM_M3_calc_z3,
                                CM_M3_Prism2_x2,CM_M3_Prism2_y2,CM_M3_Prism2_z2,
                                CM_M3_Prism1_x1,CM_M3_Prism1_y1,CM_M3_Prism1_z1) 
    #############################################################################################################################
    ######################################################## MASCON VIII ########################################################
        ###############################################################################################################
        # Tetrahedron VIII ############################################################################################
        CM_M8_Tetra_x = (Verts[Faces[it,0],0]*(1/8) + Verts[Faces[it,1],0]*(1/8) + Verts[Faces[it,2],0]*(1/8) + CM[0])/4
        CM_M8_Tetra_y = (Verts[Faces[it,0],1]*(1/8) + Verts[Faces[it,1],1]*(1/8) + Verts[Faces[it,2],1]*(1/8) + CM[1])/4
        CM_M8_Tetra_z = (Verts[Faces[it,0],2]*(1/8) + Verts[Faces[it,1],2]*(1/8) + Verts[Faces[it,2],2]*(1/8) + CM[2])/4
        # Fill array
        CM_tetra_M8[it] =  (CM_M8_Tetra_x,CM_M8_Tetra_y,CM_M8_Tetra_z)
        ###############################################################################################################
        # Tetrahedron VII #############################################################################################
        CM_Tetra_x7 = (Verts[Faces[it,0],0]*(2/8) + Verts[Faces[it,1],0]*(2/8) + Verts[Faces[it,2],0]*(2/8) + CM[0])/4
        CM_Tetra_y7 = (Verts[Faces[it,0],1]*(2/8) + Verts[Faces[it,1],1]*(2/8) + Verts[Faces[it,2],1]*(2/8) + CM[1])/4
        CM_Tetra_z7 = (Verts[Faces[it,0],2]*(2/8) + Verts[Faces[it,1],2]*(2/8) + Verts[Faces[it,2],2]*(2/8) + CM[2])/4
        # Fill array
        CM_tetra_M7[it] =  (CM_Tetra_x7,CM_Tetra_y7,CM_Tetra_z7)
        ###############################################################################################################
        # Tetrahedron VI  #############################################################################################
        CM_Tetra_x6 = (Verts[Faces[it,0],0]*(3/8) + Verts[Faces[it,1],0]*(3/8) + Verts[Faces[it,2],0]*(3/8) + CM[0])/4
        CM_Tetra_y6 = (Verts[Faces[it,0],1]*(3/8) + Verts[Faces[it,1],1]*(3/8) + Verts[Faces[it,2],1]*(3/8) + CM[1])/4
        CM_Tetra_z6 = (Verts[Faces[it,0],2]*(3/8) + Verts[Faces[it,1],2]*(3/8) + Verts[Faces[it,2],2]*(3/8) + CM[2])/4
        # Fill array
        CM_tetra_M6[it] =  (CM_Tetra_x6,CM_Tetra_y6,CM_Tetra_z6)
        ###############################################################################################################
        # Tetrahedron V  ##############################################################################################
        CM_Tetra_x5 = (Verts[Faces[it,0],0]*(4/8) + Verts[Faces[it,1],0]*(4/8) + Verts[Faces[it,2],0]*(4/8) + CM[0])/4
        CM_Tetra_y5 = (Verts[Faces[it,0],1]*(4/8) + Verts[Faces[it,1],1]*(4/8) + Verts[Faces[it,2],1]*(4/8) + CM[1])/4
        CM_Tetra_z5 = (Verts[Faces[it,0],2]*(4/8) + Verts[Faces[it,1],2]*(4/8) + Verts[Faces[it,2],2]*(4/8) + CM[2])/4
        # Fill array
        CM_tetra_M5[it] =  (CM_Tetra_x5,CM_Tetra_y5,CM_Tetra_z5)
        ###############################################################################################################
        # Tetrahedron IV  #############################################################################################
        CM_Tetra_x4 = (Verts[Faces[it,0],0]*(5/8) + Verts[Faces[it,1],0]*(5/8) + Verts[Faces[it,2],0]*(5/8) + CM[0])/4
        CM_Tetra_y4 = (Verts[Faces[it,0],1]*(5/8) + Verts[Faces[it,1],1]*(5/8) + Verts[Faces[it,2],1]*(5/8) + CM[1])/4
        CM_Tetra_z4 = (Verts[Faces[it,0],2]*(5/8) + Verts[Faces[it,1],2]*(5/8) + Verts[Faces[it,2],2]*(5/8) + CM[2])/4
        # Fill array
        CM_tetra_M4[it] =  (CM_Tetra_x4,CM_Tetra_y4,CM_Tetra_z4)
        ###############################################################################################################
        # Tetrahedron III #############################################################################################
        CM_Tetra_x3 = (Verts[Faces[it,0],0]*(6/8) + Verts[Faces[it,1],0]*(6/8) + Verts[Faces[it,2],0]*(6/8) + CM[0])/4
        CM_Tetra_y3 = (Verts[Faces[it,0],1]*(6/8) + Verts[Faces[it,1],1]*(6/8) + Verts[Faces[it,2],1]*(6/8) + CM[1])/4
        CM_Tetra_z3 = (Verts[Faces[it,0],2]*(6/8) + Verts[Faces[it,1],2]*(6/8) + Verts[Faces[it,2],2]*(6/8) + CM[2])/4
        # Fill array
        CM_tetra_M3[it] =  (CM_Tetra_x3,CM_Tetra_y3,CM_Tetra_z3)
        ###############################################################################################################
        # Tetrahedron it ##############################################################################################
        CM_Tetra_x2 = (Verts[Faces[it,0],0]*(7/8) + Verts[Faces[it,1],0]*(7/8) + Verts[Faces[it,2],0]*(7/8) + CM[0])/4
        CM_Tetra_y2 = (Verts[Faces[it,0],1]*(7/8) + Verts[Faces[it,1],1]*(7/8) + Verts[Faces[it,2],1]*(7/8) + CM[1])/4
        CM_Tetra_z2 = (Verts[Faces[it,0],2]*(7/8) + Verts[Faces[it,1],2]*(7/8) + Verts[Faces[it,2],2]*(7/8) + CM[2])/4
        # Fill array
        CM_tetra_M2[it] =  (CM_Tetra_x2,CM_Tetra_y2,CM_Tetra_z2)
        #############################################################################################
        # Total Tetrahedron #########################################################################
        CM_Tetra_x1 = (Verts[Faces[it,0],0] + Verts[Faces[it,1],0] + Verts[Faces[it,2],0] + CM[0])/4
        CM_Tetra_y1 = (Verts[Faces[it,0],1] + Verts[Faces[it,1],1] + Verts[Faces[it,2],1] + CM[1])/4
        CM_Tetra_z1 = (Verts[Faces[it,0],2] + Verts[Faces[it,1],2] + Verts[Faces[it,2],2] + CM[2])/4
        # Fill array
        CM_tetra_M1[it] =  (CM_Tetra_x1,CM_Tetra_y1,CM_Tetra_z1)
        #########################################################################################
        ################################# Prism Center of Masses ################################
        ######################################################################################################################
        # Prism CM VII Layer #################################################################################################
        CM_calc_x7 = ( CM_tetra_M7[it,0]*Volume_T_VII[it] - CM_tetra_M8[it,0]*Volume_T_VIII[it] )/Vol_Prism_VII[it]
        CM_calc_y7 = ( CM_tetra_M7[it,1]*Volume_T_VII[it] - CM_tetra_M8[it,1]*Volume_T_VIII[it] )/Vol_Prism_VII[it]
        CM_calc_z7 = ( CM_tetra_M7[it,2]*Volume_T_VII[it] - CM_tetra_M8[it,2]*Volume_T_VIII[it] )/Vol_Prism_VII[it]
        ######################################################################################################################
        # Prism CM VI Layer #################################################################################################
        CM_calc_x6 = ( CM_tetra_M6[it,0]*Volume_T_VI[it] - CM_tetra_M7[it,0]*Volume_T_VII[it] )/Vol_Prism_VI[it]
        CM_calc_y6 = ( CM_tetra_M6[it,1]*Volume_T_VI[it] - CM_tetra_M7[it,1]*Volume_T_VII[it] )/Vol_Prism_VI[it]
        CM_calc_z6 = ( CM_tetra_M6[it,2]*Volume_T_VI[it] - CM_tetra_M7[it,2]*Volume_T_VII[it] )/Vol_Prism_VI[it]
        ####################################################################################################################
        # Prism CM V Layer #################################################################################################
        CM_calc_x5 = ( CM_tetra_M5[it,0]*Volume_T_V[it] - CM_tetra_M6[it,0]*Volume_T_VI[it] )/Vol_Prism_V[it]
        CM_calc_y5 = ( CM_tetra_M5[it,1]*Volume_T_V[it] - CM_tetra_M6[it,1]*Volume_T_VI[it] )/Vol_Prism_V[it]
        CM_calc_z5 = ( CM_tetra_M5[it,2]*Volume_T_V[it] - CM_tetra_M6[it,2]*Volume_T_VI[it] )/Vol_Prism_V[it] 
        #####################################################################################################################
        # Prism CM IV Layer #################################################################################################
        CM_calc_x4 = ( CM_tetra_M4[it,0]*Volume_T_IV[it] - CM_tetra_M5[it,0]*Volume_T_V[it] )/Vol_Prism_IV[it]
        CM_calc_y4 = ( CM_tetra_M4[it,1]*Volume_T_IV[it] - CM_tetra_M5[it,1]*Volume_T_V[it] )/Vol_Prism_IV[it]
        CM_calc_z4 = ( CM_tetra_M4[it,2]*Volume_T_IV[it] - CM_tetra_M5[it,2]*Volume_T_V[it] )/Vol_Prism_IV[it]  
        ######################################################################################################################
        # Prism CM III Layer #################################################################################################
        CM_calc_x3 = ( CM_tetra_M3[it,0]*Volume_T_III[it] - CM_tetra_M4[it,0]*Volume_T_IV[it] )/Vol_Prism_III[it]
        CM_calc_y3 = ( CM_tetra_M3[it,1]*Volume_T_III[it] - CM_tetra_M4[it,1]*Volume_T_IV[it] )/Vol_Prism_III[it]
        CM_calc_z3 = ( CM_tetra_M3[it,2]*Volume_T_III[it] - CM_tetra_M4[it,2]*Volume_T_IV[it] )/Vol_Prism_III[it]   
        ######################################################################################################################
        # Prism CM it Layer  #################################################################################################
        CM_calc_x2 = ( CM_tetra_M2[it,0]*Volume_T_II[it] - CM_tetra_M3[it,0]*Volume_T_III[it] )/Vol_Prism_II[it]
        CM_calc_y2 = ( CM_tetra_M2[it,1]*Volume_T_II[it] - CM_tetra_M3[it,1]*Volume_T_III[it] )/Vol_Prism_II[it]
        CM_calc_z2 = ( CM_tetra_M2[it,2]*Volume_T_II[it] - CM_tetra_M3[it,2]*Volume_T_III[it] )/Vol_Prism_II[it]
        #####################################################################################################################
        # Prism CM Top Layer ################################################################################################
        CM_calc_x1 = ( CM_tetra_M1[it,0]*Volume_T_I[it] - CM_tetra_M2[it,0]*Volume_T_II[it] )/Vol_Prism_I[it]
        CM_calc_y1 = ( CM_tetra_M1[it,1]*Volume_T_I[it] - CM_tetra_M2[it,1]*Volume_T_II[it] )/Vol_Prism_I[it]
        CM_calc_z1 = ( CM_tetra_M1[it,2]*Volume_T_I[it] - CM_tetra_M2[it,2]*Volume_T_II[it] )/Vol_Prism_I[it]
    ################################################# 
    #################################################
        Output_Array_MVIII[it] = (CM_M8_Tetra_x,CM_M8_Tetra_y,CM_M8_Tetra_z,
                                    CM_calc_x7,CM_calc_y7,CM_calc_z7,
                                    CM_calc_x6,CM_calc_y6,CM_calc_z6,
                                    CM_calc_x5,CM_calc_y5,CM_calc_z5,
                                    CM_calc_x4,CM_calc_y4,CM_calc_z4,
                                    CM_calc_x3,CM_calc_y3,CM_calc_z3,
                                    CM_calc_x2,CM_calc_y2,CM_calc_z2,
                                    CM_calc_x1,CM_calc_y1,CM_calc_z1)  
    ###################################### fin ####################################
    ###############################################################################
    #%% File Save
    ####################################################################
    ########################################################### CM Save
    np.savetxt( Asteroid_Save +"_CM.in", Output_Array_M1,delimiter=' ')
    # np.savetxt( asteroid +"_M3CM.dat", Output_Array_MIII,delimiter=' ')
    # np.savetxt( asteroid +"_M8CM.dat", Output_Array_MVIII,delimiter=' ')

    #####################################################################################
    ######################### Gravitational Parameter Save ##############################
    np.savetxt( Asteroid_Save +"_mu.in", mu_i,delimiter=' ')
    # np.savetxt( asteroid +"_M3mu.in", mu_MIII,delimiter=' ')
    # np.savetxt( asteroid +"_M8mu.in", mu_MVIII,delimiter=' ')
    #####################################################################################
    ######################### Constants Save #############################################
    # write constants to file for feeding into calculations 
    #
    # Check if no spin rate is provided, keep space for it
    if not Spin_Rate or int(Spin_Rate) == 0:
        Spin = np.nan
    else:
        Spin = int(Spin_Rate)
        
    Const_File = pd.DataFrame({
        'Mass (kg)':[Mass_Check], 
        'Radius (Km)':[R_eff], 
        'Poly_x':CM[0], 
        'Poly_y':CM[1], 
        'Poly_z':CM[2],
        'Gamma':[Gamma_Coeff],
        'Spin (rot/hr)':[Spin],
    })
    Const_File.to_csv( Asteroid_Save +"_const.in", sep=' ',index=False)
    
    ####################################################################
    ######################################################### Volume 
    # Vol_DF_Out = pd.DataFrame({'Vol Tetra':Full_Volume_Tetra_Array[:],})
    # Vol_DF_Out_M3 = pd.DataFrame({
    # 'Vol P I'   :Vol_Prism_I[:],
    # 'Vol P II'  :Vol_Prism_II[:],
    # 'Vol T III' :Volume_T_III[:] 
    # })
    # Vol_DF_Out_M8 = pd.DataFrame({
    # 'Vol P I'   :Vol_Prism_I[:],
    # 'Vol P II'  :Vol_Prism_II[:] ,
    # 'Vol P III' :Vol_Prism_III[:] ,
    # 'Vol P IV'  :Vol_Prism_IV[:] ,
    # 'Vol P V'   :Vol_Prism_V[:] ,
    # 'Vol P VI'  :Vol_Prism_VI[:] ,
    # 'Vol P VII' :Vol_Prism_VII[:] ,
    # 'Vol T VIII':Volume_T_VIII[:] 
    # })
    ################################################################
    ########################################## Tetra/Prism Volume ##
    # Vol_DF_Out.to_csv(Aster_Vol_PATH + asteroid + '_VolM1.csv' ,sep=' ' ,index=False )
    # Vol_DF_Out_M3.to_csv(Aster_Vol_PATH_MIII + asteroid + '_VolM3.csv' ,sep=' ' ,index=False )
    # Vol_DF_Out_M8.to_csv(Aster_Vol_PATH_MVIII + asteroid + '_VolM8.csv' ,sep=' ' ,index=False )
    ###############################################################
    ######## Total Volume Save ###################
    # Append and save Constants file
    # Data_append = {'Volume Calc': Total_Volume_Out,
    #                'Volume Mirt': Volume_trimesh,
    #                'Scaling':     Gamma_Coeff,
    #                }
    # # Append data frame
    # Asteroid_Const = Asteroid_Const.assign(**Data_append)
    # print(Asteroid_Const.head())
    # # Save data frame
    # Asteroid_Const.to_csv(Aster_Const_PATH + asteroid  + "_const.in", sep=' ' ,index=False)
    ###################################################################
    Data_message = f"""
    {'-'*42}
    | Data ready, See respective directories
    {'-'*42}
    """
    running_log_callback(Data_message)
    ###########################
    ###############################################################
    ###############################################################
    return Output_Array_M1, Gamma_Coeff, Verts, Faces 
def plot_modeler_results(main_window, Output_Array_M1, Verts, Faces):
    ###############################################
    # verts already scaled by modeler and reused
    #  IFF new load of .obj is performed scaling 
    #  will be required
    Verts = Verts 
    Faces = Faces - 1
    # Plotting MASCON I
    xpM1 = Output_Array_M1[:, 0]
    ypM1 = Output_Array_M1[:, 1]
    zpM1 = Output_Array_M1[:, 2]
    # Adding the mesh to the MASCON I plot
    mesh1 = Poly3DCollection([Verts[ii] for ii in Faces],
                            edgecolor="yellow",
                            facecolors="white",
                            linewidth=0.015,
                            alpha=0.01)
    main_window.modeler_plot_canvas.set_canvas_configuration()
    main_window.modeler_plot_canvas.axis.add_collection3d(mesh1)
    main_window.modeler_plot_canvas.axis.scatter3D(xpM1, ypM1, zpM1, marker='.', color='#ff06b5')
    main_window.modeler_plot_canvas.axis.set_aspect('equal', 'box')
    main_window.modeler_plot_canvas.draw()