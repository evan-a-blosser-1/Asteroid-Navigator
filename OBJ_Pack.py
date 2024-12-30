"""
"""
import numpy as np
#%% volInt.c File Make

def volInt_Input_FileMake(OBJ_File, PATH_Name, asteroid):
    """Reads an OBJ file and converts the vertex and face data into arrays, then saves them to a file.
    
    Args:
        OBJ_File (str): Path to the .obj file
        Save_Path (str): Path to save the output file
    
    Returns:
        tuple: Two numpy arrays containing vertex and face data
    """
    # Read the OBJ file data
    with open(OBJ_File, 'r') as file:
        lines = file.readlines()
    
    # Extract vertex data
    vertices = np.array([list(map(float, line.strip().split()[1:])) for line in lines if line.startswith('v')])
    vertices = np.vstack(([0, 0, 0], vertices))
    
    # Extract face data
    faces = np.array([list(map(int, line.strip().split()[1:])) for line in lines if line.startswith('f')])
    faces = np.hstack((np.full((faces.shape[0], 1), 3), faces))
    
    
    
    # Save the data to a file
    Save_Path = PATH_Name + asteroid + ".in"    
    with open(Save_Path, "w") as Poly_Data_file:
        np.savetxt(Poly_Data_file, [len(vertices) - 1], fmt='%s', delimiter='\t')  # Number of vertices
        np.savetxt(Poly_Data_file, vertices, fmt='%.8f', delimiter='\t')  # Vertex data
        np.savetxt(Poly_Data_file, [len(faces)], fmt='%s', delimiter=' ')  # Number of faces
        np.savetxt(Poly_Data_file, faces, fmt='%s', delimiter=' ')  # Face data
    
    return vertices, faces
 

#%% OBJ file to Vertices & Faces 
###################################################################################
#################################### OBJ to Vertices & Faces ######################
def OBJ_2_VertFace(OBJ_File):
    """Reads an OBJ file and extracts the vertex and face data.
    
    Args:
        OBJ_File (str): Path to the .obj file
    
    Returns:
        tuple: Two numpy arrays containing vertex and face data
    """
    # Read the OBJ file data
    with open(OBJ_File, 'r') as file:
        lines = file.readlines()
    
    # Extract vertex data
    vertices = np.array([list(map(float, line.strip().split()[1:])) for line in lines if line.startswith('v')])
    
    # Extract face data
    faces = np.array([list(map(int, line.strip().split()[1:])) for line in lines if line.startswith('f')])
    
    return vertices, faces

#%% Tetrahedron Volume Calculation
###############################################
def Tetra_Volume(Verts,Faces,MASCON_Div,total_vol=0):
    """Tetrahedron Volume Calculation 
            - for polyhedron shape models. 
    Args:
        Verts (array): Polyhedron Vertices
        Faces (array): Polyhedron Face Data
        MASCON_Div (fracton): Divides the tetrahedron, set = 1 for total
        total_vol (bool): 0/1 for total volume calculations, or  just tetra calculations
    """
    Faces = Faces - 1
    tetra_count=np.shape(Faces)[0]
    ###################    
    # Make list
    Volume_Tetra_Array  = np.zeros(tetra_count, dtype="float64")

    Tetra_U = np.zeros((tetra_count, 3), dtype="float64")
    Tetra_V = np.zeros((tetra_count, 3), dtype="float64")
    Tetra_W = np.zeros((tetra_count, 3), dtype="float64")
    
    ############################### Analytical Method:
    for it in range(0,tetra_count):
    ##### Center of mass to vertex vectors
        U_vec = np.array([Verts[Faces[it,0],0]*MASCON_Div,
                          Verts[Faces[it,0],1]*MASCON_Div,
                          Verts[Faces[it,0],2]*MASCON_Div
                        ]) 
        V_vec = np.array([Verts[Faces[it,1],0]*MASCON_Div,
                          Verts[Faces[it,1],1]*MASCON_Div,
                          Verts[Faces[it,1],2]*MASCON_Div
                        ]) 
        W_vec = np.array([Verts[Faces[it,2],0]*MASCON_Div,
                          Verts[Faces[it,2],1]*MASCON_Div,
                          Verts[Faces[it,2],2]*MASCON_Div
                        ])
        Tetra_U[it] = U_vec
        Tetra_V[it] = V_vec 
        Tetra_W[it] = W_vec 
        ######################################################
        ############### Triple Scalar Product ################
        Vol_Full_Sum =  U_vec[0]*(V_vec[1]*W_vec[2] - V_vec[2]*W_vec[1]) -\
                        U_vec[1]*(V_vec[0]*W_vec[2] - V_vec[2]*W_vec[0]) +\
                        U_vec[2]*(V_vec[0]*W_vec[1] - V_vec[1]*W_vec[0]) 
        #
        Vol_tetra_full = (1/6) * abs((Vol_Full_Sum))
        Volume_Tetra_Array[it] = Vol_tetra_full
    ###############################################
    ######### Sum Tetra Volumes for Check #########
    Total_Volume_Out = np.sum(Volume_Tetra_Array)
    #######################################################
    ######### Output ######################################
    if total_vol == 0:
        Volume_Tetra_Array_Message = f"""
{'-'*42}
|
|  Total Volume Calculated as:
|    V = {Total_Volume_Out}
|
{'-'*42}
"""
        # print(Volume_Tetra_Array_Message )
        return Volume_Tetra_Array,Total_Volume_Out
    elif total_vol == 1:
        return Volume_Tetra_Array