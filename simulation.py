# simulation.py
import numpy as np
import pandas as pd
import trimesh
import time
import sys  
from scipy.integrate import solve_ivp
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
##########################################
######################### Personal Modules
sys.dont_write_bytecode = True
import constants as C



    
def run_simulation(main_window, data_file, obj_file, mu_file,class_name,Const_file):
    asteroid = main_window.combo_folder.currentText()
    y0 = float(main_window.y0_input.text())
    H0 = float(main_window.h0_input.text())
    days = float(main_window.days_input.text())

    const = C.constants()
    A_Const = pd.read_csv(Const_file, sep=' ')
    Gamma = A_Const['Gamma'][0]
    Poly_CM = [A_Const['Poly_x'][0], 
               A_Const['Poly_y'][0], 
               A_Const['Poly_z'][0]]
    Spin = A_Const['Spin (rot/hr)'][0]
    omega = ((2*np.pi)/Spin) *(1/3600)
    CM = np.loadtxt(data_file, delimiter=' ', dtype=float)
    Terta_Count = len(CM)
    mu_I = np.loadtxt(mu_file, delimiter=' ')
    mu = np.sum(mu_I)
    ###
    mesh = trimesh.load_mesh(obj_file)
    mesh.apply_scale(Gamma)
    # Calculate the volume of the polyhedron
    polyhedron_volume = mesh.volume
    R_eff = (3 * polyhedron_volume / (4 * np.pi)) ** (1/3)
    ###
    def EOM_MASCON(t, a):
        x, y, z, vx, vy, vz = a
        r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        v_mag = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
        dxdt = vx
        dydt = vy
        dzdt = vz
        points = 1
        Ux = np.zeros(points, dtype="float64")
        Uy = np.zeros(points, dtype="float64")
        Uz = np.zeros(points, dtype="float64")
        for it in range(Terta_Count):
            R_x = CM[it, 0] - Poly_CM[0]
            R_y = CM[it, 1] - Poly_CM[1]
            R_z = CM[it, 2] - Poly_CM[2]
            r_x = a[0] - R_x
            r_y = a[1] - R_y
            r_z = a[2] - R_z
            vector = np.array([r_x, r_y, r_z])
            r_i = np.linalg.norm(vector)
            Ux += - (mu_I[it] * r_x) / (r_i ** 3)
            Uy += - (mu_I[it] * r_y) / (r_i ** 3)
            Uz += - (mu_I[it] * r_z) / (r_i ** 3)
        dvxdt = (omega ** 2) * x + 2 * omega * vy + Ux[0]
        dvydt = (omega ** 2) * y - 2 * omega * vx + Uy[0]
        dvzdt = Uz[0]
        dadt = [dxdt, dydt, dzdt, dvxdt, dvydt, dvzdt]
        return dadt

    def v_calc(Ham, omega, mu_I, CM, yp):
        U = np.zeros(1, dtype="float64")
        for it in range(len(CM)):
            x = 0 - CM[it, 0]
            y = yp - CM[it, 1]
            z = 0 - CM[it, 2]
            r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
            U += mu_I[it] / r
        psu = U[0]
        centri = (omega ** 2) * (y ** 2)
        arg = 2 * Ham + centri + 2 * psu
        if arg > 0:
            V = np.sqrt(arg)
        return V

    def Collision(t, a):
        global cond, critical
        mesh = trimesh.load_mesh(obj_file)
        mesh.apply_scale(Gamma)
        point_coords = np.array([a[0], a[1], a[2]])
        r_mag = np.sqrt(a[0] ** 2 + a[1] ** 2 + a[2] ** 2)
        is_inside = mesh.contains([point_coords])[0]
        if not is_inside:
            closest_point, distance, _ = mesh.nearest.on_surface([point_coords])
            tolerance = 0.001
            is_on_surface = distance[0] < tolerance
        else:
            is_on_surface = False
        cond = 0
        if is_inside:
            cond = 1
        elif is_on_surface:
            cond = 1
        if cond != 0:
            critical = 0
        else:
            critical = 1
        return critical

    Collision.direction = -1
    Collision.terminal = True

    def Escape(t, a):
        global cond, critical
        R_mag = np.sqrt(a[0] ** 2 + a[1] ** 2 + a[2] ** 2)
        cond = 0
        if R_mag >= 34.0:
            cond = 9
        if cond != 0:
            critical = 0
        else:
            critical = 1
        return critical

    Escape.direction = 1
    Escape.terminal = True

    def Plot_State(state, OBJ_File, gamma, Mesh_color, M_line=0.5, M_alpha=0.05):
        OBJ_Data = np.loadtxt(OBJ_File, delimiter=' ', dtype=str)
        vertices = np.array([line[1:].astype(float) for line in OBJ_Data if line[0] == 'v'])
        faces = np.array([line[1:].astype(int) for line in OBJ_Data if line[0] == 'f'])
        vertices = vertices * gamma
        faces = faces - 1
        mesh = Poly3DCollection([vertices[ii] for ii in faces],
                                edgecolor=Mesh_color,
                                facecolors="white",
                                linewidth=M_line,
                                alpha=M_alpha)
        X = state[0, :]
        Y = state[1, :]
        Z = state[2, :]
        main_window.simulation_plot_canvas.set_canvas_configuration()
        main_window.simulation_plot_canvas.axis.add_collection3d(mesh)
        main_window.simulation_plot_canvas.axis.plot(X, Y, Z, label='Trajectory', color='magenta')
        main_window.simulation_plot_canvas.axis.set_aspect('equal', 'box')
        main_window.simulation_plot_canvas.draw()

    T = 2 * np.pi * np.sqrt(R_eff ** 3 / mu)
    dt_min = T / T
    Start_Time = 0.0
    End_Time = days * const.day
    dN = round((End_Time - Start_Time) / dt_min)
    Time = np.linspace(start=Start_Time, stop=End_Time, num=dN)

    Calc_Start_Time = time.time()

    y = y0

    Ham = H0
    x_dot = v_calc(Ham, omega, mu_I, CM, y)
    a0 = [0.0, y, 0.0, x_dot, 0.0, 0.0]
    main_window.output.append(f"Calculating...")
    solution = solve_ivp(
        fun=EOM_MASCON,
        t_span=[Start_Time, End_Time],
        y0=a0,
        events=[Collision, Escape],
        method='DOP853',
        first_step=dt_min,
        rtol=1e-10,
        atol=1e-12,
        t_eval=Time,
    )
    state = solution.y
    sol_t = solution.t

    Calc_End_Time = time.time()
    Calculated_In = Calc_End_Time - Calc_Start_Time
    
    #####################################
    # User printout
    main_window.orbit_output.append(f"Calculation completed in {Calculated_In} seconds")
    if solution.status == 0:
        main_window.output.append(f"Calculation successful")
    if solution.status == 1:
        if solution.t_events[0].size > 0:
            main_window.output.append(f"Collision detected at {solution.t_events[0]}")
        if solution.t_events[1].size > 0:
            main_window.output.append(f"Escape detected at {solution.t_events[1]}")
            
    # main_window.orbit_output.append(f"Number of steps: {len(sol_t)}")
    # main_window.orbit_output.append(f"Final time: {sol_t[-1]}")
    main_window.orbit_output.append(f"Initial position: {state[:3,0]}")
    main_window.orbit_output.append(f"Initial velocity: {state[3:,0]}")
    main_window.orbit_output.append(f"Final position: {state[:3,-1]}")
    main_window.orbit_output.append(f"Final velocity: {state[3:,-1]}")
    # main_window.orbit_output.append(f"Final distance: {np.linalg.norm(state[:3,-1])}")
    

    Plot_State(state, obj_file, Gamma, 'red')