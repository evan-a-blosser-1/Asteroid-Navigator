"""

main_window.eqpt_output.append()


main_window.eqpt_logback.append()


"""
import sys
import time
import trimesh
import pandas as pd
import numpy as np 
from scipy.integrate import solve_ivp
from scipy.linalg import inv 
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
##########################################
######################### Personal Modules
sys.dont_write_bytecode = True
import constants as C
#################################################
def run_eqpt(main_window, data_file, obj_file, mu_file, file_name, Const_file, tol, max_iter, step_limit, gms):
    # load mesh
    mesh = trimesh.load_mesh(obj_file)
    # Scale the mesh
    A_Const = pd.read_csv(Const_file, sep=' ')
    Gamma = A_Const['Gamma'][0]
    mesh.apply_scale(Gamma)
    # Calculate the volume of the polyhedron
    polyhedron_volume = mesh.volume
    R_eff = (3 * polyhedron_volume / (4 * np.pi)) ** (1/3)
    # print(f"Volume: {polyhedron_volume} km^3")
    # print(f"Effective Radius: {R_eff} km")
    Spin = A_Const['Spin (rot/hr)'][0]
    omega = ((2*np.pi)/Spin) *(1/3600)
    ###########################################
    Guess_Lim = R_eff*1.5
    #
    CM_MI = np.loadtxt(data_file, delimiter=' ')
    mu_I = np.loadtxt(mu_file, delimiter=' ')
    ############################################
    def Hx_3D(vars,CM_MI,mu_I):
        x,y,z = vars
        Ux = np.zeros_like(x,dtype=np.float64)
        Uy = np.zeros_like(x,dtype=np.float64)
        Uz = np.zeros_like(x,dtype=np.float64)
        for i in range(len(CM_MI)):
            R_x = x - CM_MI[i,0]
            R_y = y - CM_MI[i,1]
            R_z = z - CM_MI[i,2]
            R = np.sqrt(R_x**2 + R_y**2 + R_z**2)
            ##############################
            Ux += - (mu_I[i]*R_x)/R**3
            Uy += - (mu_I[i]*R_y)/R**3
            Uz += - (mu_I[i]*R_z)/R**3
        Hx  =  x  + Ux
        Hy  =  y  + Uy
        Hz  =  z  + Uz
        Output = np.stack((Hx,Hy,Hz))
        return Output
    ###########################################
    def Hxx_3D(vars,CM_MI,mu_I):
        x,y,z = vars
        Uxx = np.zeros_like(x,dtype=np.float64)
        Uyy = np.zeros_like(x,dtype=np.float64)
        Uzz = np.zeros_like(x,dtype=np.float64)
        Uxy = np.zeros_like(x,dtype=np.float64)
        Uxz = np.zeros_like(x,dtype=np.float64)
        Uyx = np.zeros_like(x,dtype=np.float64)
        Uyz = np.zeros_like(x,dtype=np.float64)
        Uzx = np.zeros_like(x,dtype=np.float64)
        Uzy = np.zeros_like(x,dtype=np.float64)
        for i in range(len(CM_MI)):
            R_x = x - CM_MI[i,0]
            R_y = y - CM_MI[i,1]
            R_z = z - CM_MI[i,2]
            R = np.sqrt(R_x**2 + R_y**2 + R_z**2)
            ##############################
            Uxx += (3*mu_I[i]*R_x**2)/R**5
            Uyy += (3*mu_I[i]*R_y**2)/R**5
            Uzz += (3*mu_I[i]*R_z**2)/R**5
            ##############################  
            Uxy += (3*mu_I[i]*R_x*R_y)/R**5
            Uxz += (3*mu_I[i]*R_x*R_z)/R**5
            Uyx += (3*mu_I[i]*R_y*R_x)/R**5
            Uyz += (3*mu_I[i]*R_y*R_z)/R**5
            Uzx += (3*mu_I[i]*R_z*R_x)/R**5
            Uzy += (3*mu_I[i]*R_z*R_z)/R**5
        Top_Row = np.column_stack((Uxx,Uxy,Uxz))
        Mid_Row = np.column_stack((Uyx,Uyy,Uyz))
        Bot_Row = np.column_stack((Uzx,Uzy,Uzz))
        Output  = np.concatenate((Top_Row,Mid_Row,Bot_Row), axis =0)
        return Output
    ###########################################
    def New_Raph_3D(x0, tol, max_iter, step_limit):
        x = x0
        for i in range(max_iter):
            hx  = Hx_3D(x,CM_MI,mu_I)
            hxx = Hxx_3D(x,CM_MI,mu_I)
            # Debug
            # Check = np.dot(hxx,inv(hxx))
            # if hx[0] == 0.0 or hx[1] == 0.0:
            #     return x
            try:
                delta_x = np.matmul(-inv(hxx), hx)
                # delta_x = np.linalg.solve(hxx_reg, -hx) 
                
    
                    
                    
            except np.linalg.LinAlgError:
                main_window.eqpt_logback.append("Hessian is singular, stopping...")
                break
            
            
            # Limit the step size
            if np.linalg.norm(delta_x) > step_limit:
                delta_x = delta_x * (step_limit / np.linalg.norm(delta_x))
            
            
            x += delta_x
            if np.linalg.norm(delta_x) <tol:
                # main_window.eqpt_logback.append(f"| Converged to equilibrium point after {i+1} iterations")
                return x
            
        # main_window.eqpt_logback.append("| Did not converge within max iterations.")
        return x
    ####################################################
    ####################################################
    
    LB = -Guess_Lim
    UB =  Guess_Lim
    x_range = np.linspace(LB, UB, gms)
    y_range = np.linspace(LB, UB, gms)
    z_range = np.linspace(LB, UB, gms)
    X, Y, Z = np.meshgrid(x_range, y_range, z_range)
    # Flatten the meshgrid to create a list of initial guesses
    initial_guesses = [np.array([[x], [y], [z]]) for x, y, z in zip(X.ravel(), Y.ravel(), Z.ravel())]
    ####################################################
    Equi_Pts_3D = []
    for i, guess in enumerate(initial_guesses):
        #
        guess_copy = guess.copy()
        #
        # main_window.eqpt_logback.append(f"| Initial Guess:     x = {guess[0,0]:.6f}, y = {guess[1,0]:.6f}, z = {guess[2,0]:.6f}")
        Equilibrium = New_Raph_3D(guess_copy,tol=tol,max_iter=max_iter, step_limit=step_limit)
        #
        Equi_Pts_3D.append(Equilibrium)

        # main_window.eqpt_logback.append(f"| Equilibrium Reached at: x = {Equilibrium[0,0]:.6f}, y = {Equilibrium[1,0]:.6f}, z = {Equilibrium[2,0]:.6f}")
        
    main_window.eqpt_logback.append('Calculations complete, testing points for stabiltiy...')

    #################################################
    #################################################
    # Remove duplicate roots
    Equi_3D = np.unique(Equi_Pts_3D, axis=0)
    ########################################
    EQ3D_text = []
    for ii in range(len(Equi_3D)):
        EQ3D_text.append(f"E{ii}")
        ############################################
        # Guess the equilibrium point
        #vars = np.array([-0.653684,-0.423982,-0.609105],dtype=np.float64)
        vars = Equi_3D[ii]
        ###
        # Skip Eq. points inside asteroid
        vars_chk = Equi_3D[ii].ravel()[np.newaxis, :]
        # main_window.eqpt_logback.append(f"Equilibrium Point {vars_chk[0]} ")
        is_inside = mesh.contains([vars_chk[0]])
        if is_inside:
            main_window.eqpt_logback.append(f"Skipping point {vars_chk[0]} as its inside the asteroid")
            continue
        ########################################################
        ######## Testing Stability of Equilibrium Point ######## 
        ########################################################
        Output_A = Hxx_3D(vars,CM_MI,mu_I)
        # print(Output_A)
        U_xx = Hxx_3D(vars,CM_MI,mu_I)[0,0]
        U_xy = Hxx_3D(vars,CM_MI,mu_I)[0,1]
        U_xz = Hxx_3D(vars,CM_MI,mu_I)[0,2]
        U_yx = Hxx_3D(vars,CM_MI,mu_I)[1,0]
        U_yy = Hxx_3D(vars,CM_MI,mu_I)[1,1]
        U_yz = Hxx_3D(vars,CM_MI,mu_I)[1,2]
        U_zx = Hxx_3D(vars,CM_MI,mu_I)[2,0]
        U_zy = Hxx_3D(vars,CM_MI,mu_I)[2,1]
        U_zz = Hxx_3D(vars,CM_MI,mu_I)[2,2]
        # print(U_xx)
        # print(U_xy)
        ###################################
        # Characteristic coefficinats 
        Alpha = U_xx + U_yy + U_zz + 4*omega**2
        #
        Beta = U_xx*U_yy + U_yy*U_zz + U_zz*U_xx \
            - U_xy**2 - U_yz**2 - U_zx**2 \
                -U_zz*U_xy**2
        #
        Gamma = U_xx*U_yy*U_zz + 2*U_xy*U_yz*U_zx \
                    - U_xx*U_yz**2 - U_yy*U_zx**2 - U_zz*U_xy**2
        ###
        def Charac_Lam(Alpha, Beta, Gamma):
            
            # Coefficients of the polynomial Lambda**6 + Alpha*Lambda**4 + Beta*Lambda**2 + Gamma = 0
            coefficients = [1, 0, Alpha, 0, Beta, 0, Gamma]
            
            # Find the roots of the polynomial
            eigenvalues = np.roots(coefficients)
            
            return eigenvalues
        char_lam = Charac_Lam(Alpha, Beta, Gamma)
        main_window.eqpt_output.append(f"Eigenvalues: {char_lam}")

        ################################
        # Theorem 1 
        def is_positive_definite(matrix):
            """ Test eigenvalues of a matrix to 
                determine if it is positive definite.
                
            """
            eigenvalues = np.linalg.eigvals(matrix)
            return np.all(eigenvalues > 0)
        # Call the function
        pos_def = is_positive_definite(Output_A)
        main_window.eqpt_output.append(f"Positive Definite: {str(pos_def)}") 
        ###
        def Chol_positive_definite(matrix):
            """
            check for positive definite matrix
            using cholesky decomposition.
            
            This will only retunr true if the matrix is 
            positive definite 
            
            AND!!!
            
            symmetric, this if the eigenvalues retunr true,
            yet this returns false. Ten the matrix is 
            not symmetric, but positive definite. 


            """
            try:
                np.linalg.cholesky(matrix)
                return True
            except np.linalg.LinAlgError:
                return False
        # Call Cholesky function
        Chol_result = Chol_positive_definite(Output_A)
        main_window.eqpt_output.append(f"Symmetric: {str(Chol_result)}")  # Output: True
        
        ##########################################################
        ## Theorem 2
        #
        TH2_1 = U_xx + U_yy + U_zz  + 4*omega**2
        #
        TH2_2 = U_xx*U_yy + U_yy*U_zz + U_zz*U_xx + 4 *omega*U_zz
        TH2_2RS = U_xy**2 + U_yz**2 + U_xz**2 
        #
        TH2_3 = U_xx*U_yy*U_zz + 2*U_xy*U_yz*U_xz 
        TH2_3RS = U_xx*U_yz**2 + U_yy*U_xz**2 + U_zz*U_xy**2
        #
        TH2_4 = Alpha**2 + 18 *Alpha*Beta*Gamma 
        TH2_4RS = 4*Alpha**3*Gamma + 4 * Beta**3 + 27*Gamma**2
        #
        def Stability_Th1(TH2_1, TH2_2, TH2_2RS, TH2_3, TH2_3RS, TH2_4, TH2_4RS):
            if TH2_1 > 0:
                #main_window.eqpt_logback.append("Theorem 1: 1st condition met")
                if TH2_2 > TH2_2RS:
                    #main_window.eqpt_logback.append("Theorem 1: 2nd condition met")
                    if TH2_3 > TH2_3RS:
                        #main_window.eqpt_logback.append("Theorem 1: 3rd condition met")
                        if TH2_4 > TH2_4RS:
                            #main_window.eqpt_logback.append("Theorem 1: The equilibrium point is stable")
                            main_window.eqpt_output.append(f" Stable Eq. Point: x = {vars[0]}, y = {vars[1]}, z = {vars[2]}")
                            ############################################################################
        Stability_Th1(TH2_1, TH2_2, TH2_2RS, TH2_3, TH2_3RS, TH2_4, TH2_4RS)
        ############################################################################







def run_eqpt_simulation(main_window, data_file, obj_file, mu_file, file_name, Const_file, x_eqpt, y_eqpt, z_eqpt):
    #
    #
    days = float(main_window.eqpt_input_days.text())
    const = C.constants()
    A_Const = pd.read_csv(Const_file, sep=' ')
    Gamma = A_Const['Gamma'][0]
    CM_MI = np.loadtxt(data_file, delimiter=' ')
    mu_I = np.loadtxt(mu_file, delimiter=' ')
    mu   = np.sum(mu_I)
    ###
    Spin = A_Const['Spin (rot/hr)'][0]
    omega = ((2*np.pi)/Spin) *(1/3600)
    ###
    mesh = trimesh.load_mesh(obj_file)
    mesh.apply_scale(Gamma)
    # Calculate the volume of the polyhedron
    polyhedron_volume = mesh.volume
    R_eff = (3 * polyhedron_volume / (4 * np.pi)) ** (1/3)
    ############################################
    def Stability_EOM(t,a):
        x,y,z,vx,vy,vz = a
        dxdt = vx
        dydt = vy
        dzdt = vz
        #######################
        Uxx = np.zeros_like(x,dtype=np.float64)
        Uyy = np.zeros_like(x,dtype=np.float64)
        Uzz = np.zeros_like(x,dtype=np.float64)
        Uxy = np.zeros_like(x,dtype=np.float64)
        Uxz = np.zeros_like(x,dtype=np.float64)
        Uyx = np.zeros_like(x,dtype=np.float64)
        Uyz = np.zeros_like(x,dtype=np.float64)
        Uzx = np.zeros_like(x,dtype=np.float64)
        Uzy = np.zeros_like(x,dtype=np.float64)
        for i in range(len(CM_MI)):
            R_x = CM_MI[i,0] - x
            R_y = CM_MI[i,1] - y
            R_z = CM_MI[i,2] - z
            R = np.sqrt(R_x**2 + R_y**2 + R_z**2)
            ##############################
            Uxx += (3*mu_I[i]*R_x**2)/R**5
            Uyy += (3*mu_I[i]*R_y**2)/R**5
            Uzz += (3*mu_I[i]*R_z**2)/R**5
            ##############################  
            Uxy += (3*mu_I[i]*R_x*R_y)/R**5
            Uxz += (3*mu_I[i]*R_x*R_z)/R**5
            Uyx += (3*mu_I[i]*R_y*R_x)/R**5
            Uyz += (3*mu_I[i]*R_y*R_z)/R**5
            Uzx += (3*mu_I[i]*R_z*R_x)/R**5
            Uzy += (3*mu_I[i]*R_z*R_z)/R**5
        ########################
        # Stability Equation of Motion
        dvxdt =   2*omega*vy  - Uxx*x - Uxy*y - Uxz*z 
        dvydt = - 2*omega*vx  - Uyx*x - Uyy*y - Uyz*z
        dvzdt = - Uzx*x - Uzy*y - Uzz*z
        ###
        dadt = [dxdt, dydt, dzdt, dvxdt, dvydt, dvzdt]
        return dadt
        #################################
        #################################
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
    ###
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
    ###
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
    ###
    def Plot_State(state, OBJ_File, gamma, Mesh_color, M_line=0.5, M_alpha=0.05):
        # Load mesh using trimesh
        mesh = trimesh.load(OBJ_File)
        # Get vertices and faces
        vertices = mesh.vertices * gamma
        faces = mesh.faces - 1 
        ###
        mesh = Poly3DCollection([vertices[ii] for ii in faces],
                                edgecolor=Mesh_color,
                                facecolors="white",
                                linewidth=M_line,
                                alpha=M_alpha)
        X = state[0, :]
        Y = state[1, :]
        Z = state[2, :]
        main_window.eqpt_plot_canvas.set_canvas_configuration()
        main_window.eqpt_plot_canvas.axis.add_collection3d(mesh)
        main_window.eqpt_plot_canvas.axis.plot(X, Y, Z, label='Trajectory', color='magenta')
        main_window.eqpt_plot_canvas.axis.set_aspect('equal', 'box')
        main_window.eqpt_plot_canvas.draw()

    T = 2 * np.pi * np.sqrt(R_eff ** 3 / mu)
    dt_min = T / T
    Start_Time = 0.0
    End_Time = days * const.day
    dN = round((End_Time - Start_Time) / dt_min)
    Time = np.linspace(start=Start_Time, stop=End_Time, num=dN)

    Calc_Start_Time = time.time()
    ###
    
    #y = y0

    #Ham = H0
    
    # x_dot = v_calc(Ham, omega, mu_I, CM, y)
    a0 = [x_eqpt, y_eqpt, z_eqpt, 0.0, 0.0, 0.0]
    ###
    
    main_window.eqpt_logback.append(f"Testing Equilibrium Point...")
    solution = solve_ivp(
        fun=Stability_EOM,
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
    main_window.eqpt_logback.append(f"Calculation completed in {Calculated_In} seconds")
    if solution.status == 0:
        main_window.eqpt_logback.append(f"Calculation successful")
    if solution.status == 1:
        if solution.t_events[0].size > 0:
            main_window.eqpt_logback.append(f"Collision detected at {solution.t_events[0]}")
        if solution.t_events[1].size > 0:
            main_window.eqpt_logback.append(f"Escape detected at {solution.t_events[1]}")
            
    # main_window.eqpt_logback.append(f"Number of steps: {len(sol_t)}")
    # main_window.eqpt_logback.append(f"Final time: {sol_t[-1]}")
    main_window.eqpt_logback.append(f"Initial position: {state[:3,0]}")
    main_window.eqpt_logback.append(f"Initial velocity: {state[3:,0]}")
    main_window.eqpt_logback.append(f"Final position: {state[:3,-1]}")
    main_window.eqpt_logback.append(f"Final velocity: {state[3:,-1]}")
    # main_window.eqpt_logback.append(f"Final distance: {np.linalg.norm(state[:3,-1])}")
    

    Plot_State(state, obj_file, Gamma, 'red')