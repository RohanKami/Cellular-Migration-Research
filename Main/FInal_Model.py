import gmsh
import math
import sys
import random
import numpy 
import copy
import time
import matplotlib.pyplot as plt
import matplotlib.cm as cm


# Initialize gmsh:
gmsh.initialize()

# Function used to generate the initial mesh for a cell
#N -> Node amount
#offset_x, offset_y -> Center point of generated cell
#Cell_Coordinates -> array used to store coordinates of each cell
#Radius -> Radius of cell
def Create_Cell_Mesh(N, lc, offset_x, offset_y, Cell_Coordinates, cell_number, radius):
    cell_points = {}
    Theta = 0
    dTheta = 2 * math.pi / N    
    Cell_Coordinates[0, cell_number, 0] = offset_x
    Cell_Coordinates[1, cell_number, 0] = offset_y
    for i in range(N):
        x = radius * math.cos(Theta) + offset_x
        y = radius * math.sin(Theta) + offset_y
        cell_points[str(i+1)] = gmsh.model.geo.add_point(x, y, 0, lc)
        Cell_Coordinates[0, cell_number, i+1] = x
        Cell_Coordinates[1, cell_number, i+1] = y
        Theta += dTheta
    return cell_points

# Function used to connect the points along a cell's mesh
#N -> Node amount
#cell_points -> G_Mesh Coordinates of points
def Connect_Mesh_Points(N, cell_points):
    line_dict = {}
    for i in range(N):
        if (i + 1 != N):
            line_dict[str(i+1)] = gmsh.model.geo.add_line(cell_points[str(i+1)], cell_points[str(i+2)])
        else:
            line_dict[str(i+1)] = gmsh.model.geo.add_line(cell_points[str(i+1)], cell_points[str(1)])
    return line_dict

# Function to generate random offsets ensuring no overlap for initial cell placement 
#used_offsets -> Coordinates of where cells have already been placed
#spacing -> minimum space we want between cells
#concentration_point_coordinates -> Coordinates of where chemo sources have been placed
#max_attempt -> ammount of attempts it does to place the point
import random
import math

import random
import math

def generate_random_offset(used_offsets, spacing, concentration_point_coordinates, R_Boundary, max_attempts=100):
    # Specify circular grid boundary (radius of the circle)
    Grid_Radius = R_Boundary  # Original radius of the grid
    effective_radius = Grid_Radius - 1  # Ensure points are at least 1 unit from the perimeter
    amount_concentration_signal = concentration_point_coordinates.shape[1]
    
    for _ in range(max_attempts):
        # Generate a random angle and distance within the effective radius
        angle = random.uniform(0, 2 * math.pi)
        distance = random.uniform(0, effective_radius)
        offset_x = distance * math.cos(angle)
        offset_y = distance * math.sin(angle)
        
        # Check against used_offsets
        if all(math.sqrt((offset_x - x)**2 + (offset_y - y)**2) >= spacing for x, y in used_offsets):
            # Check against concentration_point_coordinates
            if all(math.sqrt((offset_x - concentration_point_coordinates[0, i])**2 + 
                             (offset_y - concentration_point_coordinates[1, i])**2) >= 1.25 * spacing  # Increased chemo offset
                   for i in range(amount_concentration_signal)):
                used_offsets.append((offset_x, offset_y))
                return offset_x, offset_y
    
    raise ValueError("Failed to generate a valid random offset within the max_attempts limit.")



#Function used to solve for Mx
#node_amount -> Amount of nodes on each cell membrane or nucleus
#Membrane_cell_coordinates -> coordinates of points along membrane
#Nucleus_cell_cordinates -> Coordinates of points along nucleus
#cell_number -> Which cell are we solving Mx for
#rho_1 -> Density of nucleus
#rho_2 -> Density of membrane area
def solve_Mx(node_amount, membrane_cell_coordinates, nucleus_cell_coordinates, cell_number, rho_1, rho_2):
    sum1 = 0
    sum2 = 0
    for i in range(1, node_amount):
        sum1+= (membrane_cell_coordinates[0][cell_number][i] - membrane_cell_coordinates[0][cell_number][i + 1]) * ((membrane_cell_coordinates[1][cell_number][i] + 2*membrane_cell_coordinates[1][cell_number][i + 1]) * membrane_cell_coordinates[0][cell_number][i+1] + (2*membrane_cell_coordinates[1][cell_number][i] + membrane_cell_coordinates[1][cell_number][i+1]) * membrane_cell_coordinates[0][cell_number][i])
        sum2+= (nucleus_cell_coordinates[0][cell_number][i] - nucleus_cell_coordinates[0][cell_number][i + 1]) * ((nucleus_cell_coordinates[1][cell_number][i] + 2*nucleus_cell_coordinates[1][cell_number][i + 1]) * nucleus_cell_coordinates[0][cell_number][i+1] + (2*nucleus_cell_coordinates[1][cell_number][i] + nucleus_cell_coordinates[1][cell_number][i+1]) * nucleus_cell_coordinates[0][cell_number][i])
    
    sum1+= (membrane_cell_coordinates[0][cell_number][node_amount] - membrane_cell_coordinates[0][cell_number][1]) * ((membrane_cell_coordinates[1][cell_number][node_amount] + 2*membrane_cell_coordinates[1][cell_number][1]) * membrane_cell_coordinates[0][cell_number][1] + (2*membrane_cell_coordinates[1][cell_number][node_amount] + membrane_cell_coordinates[1][cell_number][1]) * membrane_cell_coordinates[0][cell_number][node_amount])
    sum2+= (nucleus_cell_coordinates[0][cell_number][node_amount] - nucleus_cell_coordinates[0][cell_number][1]) * ((nucleus_cell_coordinates[1][cell_number][node_amount] + 2*nucleus_cell_coordinates[1][cell_number][1]) * nucleus_cell_coordinates[0][cell_number][1] + (2*nucleus_cell_coordinates[1][cell_number][node_amount] + nucleus_cell_coordinates[1][cell_number][1]) * nucleus_cell_coordinates[0][cell_number][node_amount])
    Mx = rho_1/6 * sum1 + (rho_2-rho_1)/6 *  sum2
    return Mx

#Function used to solve for My
#node_amount -> Amount of nodes on each cell membrane or nucleus
#Membrane_cell_coordinates -> coordinates of points along membrane
#Nucleus_cell_cordinates -> Coordinates of points along nucleus
#cell_number -> Which cell are we solving Mx for
#rho_1 -> Density of cytoplasme
#rho_2 -> Density of nucleus
def solve_My(node_amount, membrane_cell_coordinates, nucleus_cell_coordinates, cell_number, rho_1, rho_2):
    sum1 = 0
    sum2 = 0
    for i in range(1, node_amount):
        sum1+= (membrane_cell_coordinates[0][cell_number][i] - membrane_cell_coordinates[0][cell_number][i+1]) * (membrane_cell_coordinates[1][cell_number][i]**2 + membrane_cell_coordinates[1][cell_number][i] * membrane_cell_coordinates[1][cell_number][i+1] + membrane_cell_coordinates[1][cell_number][i+1]**2)
        sum2+= (nucleus_cell_coordinates[0][cell_number][i] - nucleus_cell_coordinates[0][cell_number][i+1]) * (nucleus_cell_coordinates[1][cell_number][i]**2 + nucleus_cell_coordinates[1][cell_number][i] * nucleus_cell_coordinates[1][cell_number][i+1] + nucleus_cell_coordinates[1][cell_number][i+1]**2)

    sum1+= (membrane_cell_coordinates[0][cell_number][node_amount] - membrane_cell_coordinates[0][cell_number][1]) * (membrane_cell_coordinates[1][cell_number][node_amount]**2 + membrane_cell_coordinates[1][cell_number][node_amount] * membrane_cell_coordinates[1][cell_number][1] + membrane_cell_coordinates[1][cell_number][1]**2)
    sum2+= (nucleus_cell_coordinates[0][cell_number][node_amount] - nucleus_cell_coordinates[0][cell_number][1]) * (nucleus_cell_coordinates[1][cell_number][node_amount]**2 + nucleus_cell_coordinates[1][cell_number][node_amount] * nucleus_cell_coordinates[1][cell_number][1] + nucleus_cell_coordinates[1][cell_number][1]**2)
    My = rho_1/6 * sum1 + (rho_2-rho_1)/6 *  sum2
    return My

#Function used to solve for weighted area
#node_amount -> amount of nodes on cell membrane or nucleus boundary
#membrane_cell_coordinates -> array storing coordinates of the cell membrane
#nucleus_cell_coordinates -> array storing coordinates of nucleus membrane
#cell_number -> number of cell which we are calculating A for
#rho_1 -> density of cytoplasme
#rho_2 -> density of nucleus
def solve_A(node_amount, membrane_cell_coordinates, nucleus_cell_coordinates, cell_number, rho_1, rho_2):
    sum1 = 0
    sum2 = 0
    for i in range(1, node_amount):
        sum1+=(membrane_cell_coordinates[1][cell_number][i]+membrane_cell_coordinates[1][cell_number][i+1]) * (membrane_cell_coordinates[0][cell_number][i]-membrane_cell_coordinates[0][cell_number][i+1])
        sum2+=(nucleus_cell_coordinates[1][cell_number][i]+nucleus_cell_coordinates[1][cell_number][i+1]) * (nucleus_cell_coordinates[0][cell_number][i]-nucleus_cell_coordinates[0][cell_number][i+1])
   
    sum1+=(membrane_cell_coordinates[1][cell_number][node_amount]+membrane_cell_coordinates[1][cell_number][1]) * (membrane_cell_coordinates[0][cell_number][node_amount]-membrane_cell_coordinates[0][cell_number][1])    
    sum2+=(nucleus_cell_coordinates[1][cell_number][node_amount]+nucleus_cell_coordinates[1][cell_number][1]) * (nucleus_cell_coordinates[0][cell_number][node_amount]-nucleus_cell_coordinates[0][cell_number][1])
    A = rho_1/2 * sum1 + (rho_2-rho_1)/2 * sum2
    return A

#Function used to solve for center of mass of cel
#Mx -> Mx of cell
#My -> My of cell
#A -> A of cell
def solve_center_of_mass(Mx,My,A):
    x_c = Mx/A
    y_c = My/A
    return x_c, y_c

#Testing function to solve for concentration.
def solve_c(x, xs, y, ys, t, mu, teps, A):
    c = A/(mu*(t+teps))*math.exp(-((x-xs)**2+(y-ys)**2)/(4*mu*(t+teps)))
    return c

#Function to solve for concentration gradient of cell
#x -> x coordinates of node
#xs -> x coordinate of source
#y -> y coordinate of node
#ys -> y coordinate of surce
#mu -> diffusion constant
#teps -> Small value to avoid division by 0
#A -> Source intensity 
#Chemo_type -> decides if its chemo attractant or chemorepellant
def solve_grad_c(x, xs, y, ys, t, mu, teps, A, chemo_type):
    cx = -(A*(x-xs))/(2*mu**2*(t+teps)**2)*math.exp(-((x-xs)**2+(y-ys)**2)/(4*mu*(t+teps)))
    cy = -(A*(y-ys))/(2*mu**2*(t+teps)**2)*math.exp(-((x-xs)**2+(y-ys)**2)/(4*mu*(t+teps)))
    if (chemo_type == 0):
        return cx, cy
    else:
        return -cx, -cy



def run_gmesh():
    # Create the relevant Gmsh data structures from Gmsh model
    gmsh.model.geo.synchronize()

    # Generate mesh
    gmsh.model.mesh.generate()

    # Write mesh data to MSH file for Gmsh
    gmsh.write("GFG.msh")
    
#Function used to update cell position on gmesh
#Cell_nucleus -> coordinates of points on cell_nucleus
#cell_membrane -> coordinates of points on cell membrane
#Center_cell_points -> Coordinate of center of the cell
#Inner_Cell_Points -> ID of nucleus nodes
#Outer_Cell_Points -> ID of membrane nodes
#Outer_Cell_Lines -> ID of Lines connecting membrane
#Inner_Cell_Lines -> ID of lines connection nucleus boundary
def update_cell_coordinates(cell_nucleus, cell_membrane, Center_Cell_Points, Inner_Cell_Points, Outer_Cell_Points, Outer_Cell_Lines, Inner_Cell_Lines):
    for cell_number in range(1, Cell_Amount + 1):
        # Remove old lines
        if str(cell_number) in Outer_Cell_Lines:
            for line_id in Outer_Cell_Lines[str(cell_number)].values():
                gmsh.model.geo.remove([(1, line_id)])
        if str(cell_number) in Inner_Cell_Lines:
            for line_id in Inner_Cell_Lines[str(cell_number)].values():
                gmsh.model.geo.remove([(1, line_id)])
        
        # Remove old center points
        if str(cell_number) in Center_Cell_Points:
            gmsh.model.geo.remove([(0, Center_Cell_Points[str(cell_number)])])
        
        # Add new center points
        Center_Cell_Points[str(cell_number)] = gmsh.model.geo.add_point(
            cell_nucleus[0, cell_number, 0], cell_nucleus[1, cell_number, 0], 0, lc)
        
        # Remove old inner and outer points
        for node_number in range(1, Node_Amount + 1):
            if str(cell_number) in Inner_Cell_Points and str(node_number) in Inner_Cell_Points[str(cell_number)]:
                gmsh.model.geo.remove([(0, Inner_Cell_Points[str(cell_number)][str(node_number)])])
            if str(cell_number) in Outer_Cell_Points and str(node_number) in Outer_Cell_Points[str(cell_number)]:
                gmsh.model.geo.remove([(0, Outer_Cell_Points[str(cell_number)][str(node_number)])])
        
        # Add new inner and outer points
        for node_number in range(1, Node_Amount + 1):
            Inner_Cell_Points[str(cell_number)][str(node_number)] = gmsh.model.geo.add_point(
                cell_nucleus[0, cell_number, node_number], cell_nucleus[1, cell_number, node_number], 0, lc)
            Outer_Cell_Points[str(cell_number)][str(node_number)] = gmsh.model.geo.add_point(
                cell_membrane[0, cell_number, node_number], cell_membrane[1, cell_number, node_number], 0, lc)
        
        # Connect new points with lines
        Outer_Cell_Lines[str(cell_number)] = Connect_Mesh_Points(Node_Amount, Outer_Cell_Points[str(cell_number)])
        Inner_Cell_Lines[str(cell_number)] = Connect_Mesh_Points(Node_Amount, Inner_Cell_Points[str(cell_number)])
        
    
def write_mesh_to_vtk(filename="time.vtk"):
    try:
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(2)  # Change to (3) if it's a 3D mesh
        
        # Write the mesh to a VTK file
        gmsh.write(filename)
        print(f"Mesh has been written to {filename}")
    
    except Exception as e:
        # Print the error found for gmesh
        print(f"An error occurred: {e}")
    
#Rework this to be a stochastic term added to ODE
def random_move(Cell_Membrane_Coordinates, Cell_Nucleus_Coordinates, Node_Amount, cell_number, step_size=0.01):
    # Move the center of the cell randomly
    dx = numpy.random.uniform(-step_size, step_size)
    dy = numpy.random.uniform(-step_size, step_size)
    
    # Update the center point of the cell
    Cell_Nucleus_Coordinates[0, cell_number, 0] += dx
    Cell_Nucleus_Coordinates[1, cell_number, 0] += dy
    
    # Update all the nodes of the cell accordingly using vectorized operations
    Cell_Membrane_Coordinates[0, cell_number, 1:Node_Amount + 1] += dx
    Cell_Membrane_Coordinates[1, cell_number, 1:Node_Amount + 1] += dy
    Cell_Nucleus_Coordinates[0, cell_number, 1:Node_Amount + 1] += dx
    Cell_Nucleus_Coordinates[1, cell_number, 1:Node_Amount + 1] += dy
    
#Function used to calculate the leonard jones forces on each membrane node
#cell_membrane_coords -> coordinates of points on membrane
#epsilon -> strenght of leonard-jones interaction
#sigma -> range of leonard-jones interaction
#dt -> time step size
#num_cells -> number of cells in simulaltion
#num_nodes -> number of nodes on membrane boundary.
def calculate_lj_forces_vectorized(cell_membrane_coords, epsilon, sigma, dt, num_cells, num_nodes):


    for cell_number in range(1, num_cells):
        for other_cell_number in range(cell_number + 1, num_cells+1):
            # Extract coordinates for the two cells
            cell_coords = cell_membrane_coords[:, cell_number, 1:num_nodes + 1]
            other_cell_coords = cell_membrane_coords[:, other_cell_number, 1:num_nodes + 1]

            # Calculate pairwise distances
            dx = cell_coords[0][:, numpy.newaxis] - other_cell_coords[0]
            dy = cell_coords[1][:, numpy.newaxis] - other_cell_coords[1]
            distances = numpy.sqrt(dx**2 + dy**2)

            # Avoid division by zero
            mask = distances > 0
            inv_distances = numpy.zeros_like(distances)
            inv_distances[mask] = 1.0 / distances[mask]

            # Calculate Lennard-Jones forces
            inv_distances_6 = inv_distances**6
            inv_distances_12 = inv_distances_6**2
            forces = 4 * epsilon * (    ((sigma / distances)**12) - (sigma / distances)**6)

            # Calculate force components
            fx = forces * dx
            fy = forces  * dy

            # Update coordinates
            cell_membrane_coords[0, cell_number, 1:num_nodes + 1] += numpy.sum(fx, axis=1) * dt
            cell_membrane_coords[1, cell_number, 1:num_nodes + 1] += numpy.sum(fy, axis=1) * dt
            cell_membrane_coords[0, other_cell_number, 1:num_nodes + 1] -= numpy.sum(fx, axis=0) * dt
            cell_membrane_coords[1, other_cell_number, 1:num_nodes + 1] -= numpy.sum(fy, axis=0) * dt
            

#Calculate mean curvature of cell
#Cell_Membrane_Coordinates -> Coordinates of each node on cell membrane
#Cell_Nucleus_Coordinates -> Coordinates of each node on cell nucleus
#cell_number -> Amount of cells in simulation
#Node_Amount -> Amount of nodes on cell membrane
def calculate_mean_curvature(Cell_Membrane_Coordinates, Cell_Nucleus_Coordinates, cell_number, Node_Amount):
    # Calculate the mean curvature for each node in the cell membrane
    mean_curvature = []
    for i in range(1, Node_Amount + 1):
        # Get the coordinates of the current, previous, and next points
        x_i = Cell_Membrane_Coordinates[0][cell_number][i]
        y_i = Cell_Membrane_Coordinates[1][cell_number][i]
        
        x_prev = Cell_Membrane_Coordinates[0][cell_number][i - 1] if i > 1 else Cell_Membrane_Coordinates[0][cell_number][Node_Amount]
        y_prev = Cell_Membrane_Coordinates[1][cell_number][i - 1] if i > 1 else Cell_Membrane_Coordinates[1][cell_number][Node_Amount]
        
        x_next = Cell_Membrane_Coordinates[0][cell_number][i + 1] if i < Node_Amount else Cell_Membrane_Coordinates[0][cell_number][1]
        y_next = Cell_Membrane_Coordinates[1][cell_number][i + 1] if i < Node_Amount else Cell_Membrane_Coordinates[1][cell_number][1]

        # Calculate vectors
        vec1 = numpy.array([x_prev - x_i, y_prev - y_i])
        vec2 = numpy.array([x_next - x_i, y_next - y_i])

        # Calculate the angle between vectors
        angle = math.acos(numpy.dot(vec1, vec2) / (numpy.linalg.norm(vec1) * numpy.linalg.norm(vec2)))
        
        # Curvature is inversely proportional to the radius of the circumscribing circle
        if numpy.linalg.norm(vec1) != 0 and numpy.linalg.norm(vec2) != 0:
            curvature = 2 * numpy.sin(angle) / (numpy.linalg.norm(vec1) + numpy.linalg.norm(vec2))
        else:
            curvature = 0
        
        mean_curvature.append(curvature)
    
    return numpy.mean(mean_curvature)
    
#N -> Number of points on boundary

def Create_Boundary_Mesh(N, lc, offset_x, offset_y, R_Boundary):
    circle_points = {}
    Theta = 0
    dTheta = 2 * math.pi / N  # Angle step for each point on the circle

    for i in range(N):
        x = R_Boundary * math.cos(Theta) + offset_x
        y = R_Boundary * math.sin(Theta) + offset_y
        circle_points[str(i+1)] = gmsh.model.geo.add_point(x, y, 0, lc)
        Theta += dTheta

    return circle_points
# Circle points:
lc = 1e-2
Node_Amount = 500 #Amout of nodes in simulation
Cell_Amount = 10 #Amount of cells in simulation
rho_1 = 1.2 #Density of cytoplasme
rho_2 = 1.4 #Density of Nucleus
mu = 2 #Diffusion constant of chemical signal
teps = 1*10**(-6) #Small number to avoid division by 0


A = 17.5 #15 #Source intensity Controlling size of chemical signal
Beta = 4 #2.5 #Magnitude of cell reponse to external stimuli
alpha =  50 #10 #Internal Cell Membranes Relaxation Coefficient
Delta = 30 #30 #External Cell Membranes Relaxation Coefficient 
delta_n = 300 #External Nucleus Relaxation Coefficient
alpha_n = 100 #Internal Nucleus Relaxation Coefficient 

#Chemoatt
grad_cx = 0
grad_cy = 0
dcx = 0
dcy =0

#0.008
#1/3
epsilon = 0.0005 # Strength of Lennard Jones-interaction 0.0005 
sigma = 2/3 # Range of Lennard-Jones Interaction 1/2

#Boundary
x_center = 0
y_center = 0 
R_Boundary = 8 #Radius of boundary
k_boundary = 200 #Spring constant of boundary

# Parameters for Brownian motion
sigma_x = 0.2  # Intensity for x-direction
sigma_y = 0.2  # Intensity for y-direction


# Create a dictionary to store points of all cells
Center_Cell_Points = {}
Outer_Cell_Points = {}
Inner_Cell_Points = {}
Outer_Cell_Lines = {}
Inner_Cell_Lines = {}
Center_to_Inner_Lines = {}
Inner_to_Outer_Lines = {}
Mean_Curvature = {cell_number: [] for cell_number in range(1, Cell_Amount + 1)}
Cell_Area= {cell_number: [] for cell_number in range(1, Cell_Amount + 1)}
time_values = []



# Parameters for cell spacing to avoid overlap
cell_spacing = 3  # Minimum distance between cell centers to avoid overlap
used_offsets = []
Cell_Membrane_Coordinates = numpy.zeros((2, Cell_Amount + 1, Node_Amount + 1)) #x dimension is coordinates y dimension is cell # and z dimension is node #
Original_Cell_Membrane_Coordinates = numpy.zeros((2, Cell_Amount + 1, Node_Amount + 1))
Cell_Nucleus_Coordinates = numpy.zeros((2, Cell_Amount + 1, Node_Amount + 1))
Original_Cell_Nucleus_Coordinates = numpy.zeros((2, Cell_Amount + 1, Node_Amount + 1))

centroid_x = {cell_number: [] for cell_number in range(1, Cell_Amount + 1)}
centroid_y = {cell_number: [] for cell_number in range(1, Cell_Amount + 1)}




Amount_Concentration_Signal = 1
Chemical_Signal_Points = numpy.zeros(Amount_Concentration_Signal)
Concentration_Point_Coordinates = numpy.zeros((5, Amount_Concentration_Signal)) #x,y,time_on,time_off which concentration point, 1 = repel / 0 = attract
Concentration_Point_flags = numpy.zeros(Amount_Concentration_Signal)

#List where you want concentration points and on/off time of each
Concentration_Point_Coordinates[0][0] = 6
Concentration_Point_Coordinates[1][0] = 0
Concentration_Point_Coordinates[2][0] = 0
Concentration_Point_Coordinates[3][0] = 40
Concentration_Point_Coordinates[4][0] = 0


# Concentration_Point_Coordinates[0][1] = 35
# Concentration_Point_Coordinates[1][1] = 13
# Concentration_Point_Coordinates[2][1] = 0
# Concentration_Point_Coordinates[3][1] = 40
# Concentration_Point_Coordinates[4][1] = 0

# Concentration_Point_Coordinates[0][2] = 35
# Concentration_Point_Coordinates[1][2] = 11
# Concentration_Point_Coordinates[2][2] = 0
# Concentration_Point_Coordinates[3][2] = 40
# Concentration_Point_Coordinates[4][2] = 0

# Concentration_Point_Coordinates[0][3] = 35
# Concentration_Point_Coordinates[1][3] = 9
# Concentration_Point_Coordinates[2][3] = 0
# Concentration_Point_Coordinates[3][3] = 40
# Concentration_Point_Coordinates[4][3] = 0

# Concentration_Point_Coordinates[0][4] = 35
# Concentration_Point_Coordinates[1][4] = 7
# Concentration_Point_Coordinates[2][4] = 0
# Concentration_Point_Coordinates[3][4] = 40
# Concentration_Point_Coordinates[4][4] = 0

# Concentration_Point_Coordinates[0][5] = 35
# Concentration_Point_Coordinates[1][5] = 5
# Concentration_Point_Coordinates[2][5] = 0
# Concentration_Point_Coordinates[3][5] = 40
# Concentration_Point_Coordinates[4][5] = 0

# Concentration_Point_Coordinates[0][6] = 35
# Concentration_Point_Coordinates[1][6] = 3
# Concentration_Point_Coordinates[2][6] = 0
# Concentration_Point_Coordinates[3][6] = 40
# Concentration_Point_Coordinates[4][6] = 0

# Concentration_Point_Coordinates[0][7] = 35
# Concentration_Point_Coordinates[1][7] = 1
# Concentration_Point_Coordinates[2][7] = 0
# Concentration_Point_Coordinates[3][7] = 40
# Concentration_Point_Coordinates[4][7] = 0

# Concentration_Point_Coordinates[0][8] = 35
# Concentration_Point_Coordinates[1][8] = -1
# Concentration_Point_Coordinates[2][8] = 0
# Concentration_Point_Coordinates[3][8] = 40
# Concentration_Point_Coordinates[4][8] = 0

# Concentration_Point_Coordinates[0][9] = 35
# Concentration_Point_Coordinates[1][9] = -3
# Concentration_Point_Coordinates[2][9] = 0
# Concentration_Point_Coordinates[3][9] = 40
# Concentration_Point_Coordinates[4][9] = 0

# Concentration_Point_Coordinates[0][10] = 35
# Concentration_Point_Coordinates[1][10] = -5
# Concentration_Point_Coordinates[2][10] = 0
# Concentration_Point_Coordinates[3][10] = 40
# Concentration_Point_Coordinates[4][10] = 0

# Concentration_Point_Coordinates[0][11] = 35
# Concentration_Point_Coordinates[1][11] = -5
# Concentration_Point_Coordinates[2][11] = 0
# Concentration_Point_Coordinates[3][11] = 40
# Concentration_Point_Coordinates[4][11] = 0

# Concentration_Point_Coordinates[0][12] = 35
# Concentration_Point_Coordinates[1][12] = -7
# Concentration_Point_Coordinates[2][12] = 0
# Concentration_Point_Coordinates[3][12] = 40
# Concentration_Point_Coordinates[4][12] = 0

# Concentration_Point_Coordinates[0][13] = 35
# Concentration_Point_Coordinates[1][13] = -9
# Concentration_Point_Coordinates[2][13] = 0
# Concentration_Point_Coordinates[3][13] = 40
# Concentration_Point_Coordinates[4][13] = 0

# Concentration_Point_Coordinates[0][14] = 35
# Concentration_Point_Coordinates[1][14] = -11
# Concentration_Point_Coordinates[2][14] = 0
# Concentration_Point_Coordinates[3][14] = 40
# Concentration_Point_Coordinates[4][14] = 0

# Concentration_Point_Coordinates[0][15] = 35
# Concentration_Point_Coordinates[1][15] = -13
# Concentration_Point_Coordinates[2][15] = 0
# Concentration_Point_Coordinates[3][15] = 40
# Concentration_Point_Coordinates[4][15] = 0

# Concentration_Point_Coordinates[0][16] = 35
# Concentration_Point_Coordinates[1][16] = -15
# Concentration_Point_Coordinates[2][16] = 0
# Concentration_Point_Coordinates[3][16] = 40
# Concentration_Point_Coordinates[4][16] = 0

#Generate Boundary
circle_points = Create_Boundary_Mesh(500, lc, x_center, y_center, R_Boundary)

# Generate mesh for each cell with a random offset
for i in range(1, Cell_Amount + 1):
    offset_x, offset_y = generate_random_offset(used_offsets, cell_spacing, Concentration_Point_Coordinates, R_Boundary)
    Center_Cell_Points[str(i)] = gmsh.model.geo.add_point(offset_x, offset_y, 0, lc)
    Outer_Cell_Points[str(i)] = Create_Cell_Mesh(Node_Amount, lc, offset_x, offset_y, Cell_Membrane_Coordinates, i, 1)
    Inner_Cell_Points[str(i)] = Create_Cell_Mesh(Node_Amount, lc, offset_x, offset_y, Cell_Nucleus_Coordinates, i, 0.25)
    

# Connect mesh points for each cell
for i in range(1, Cell_Amount + 1):
    Outer_Cell_Lines[str(i)] = Connect_Mesh_Points(Node_Amount, Outer_Cell_Points[str(i)])
    Inner_Cell_Lines[str(i)] = Connect_Mesh_Points(Node_Amount, Inner_Cell_Points[str(i)])
    

Original_Cell_Membrane_Coordinates = copy.deepcopy(Cell_Membrane_Coordinates)
Original_Cell_Nucleus_Coordinates = copy.deepcopy(Cell_Nucleus_Coordinates)


dt = 0.002 # 0.002
t = 0
T_Max = 20 #20
PV_TIME = 0

#run_gmesh()
#create_axes()


#Outer loops checks current time
while t < T_Max:
    
    #How often we export a vtk file
    if int(t/dt) % 100 == 0:
        for i in range(0, Amount_Concentration_Signal):

            if t >= Concentration_Point_Coordinates[2][i] and t <= Concentration_Point_Coordinates[3][i] and (Concentration_Point_flags[i] == 0):
                    Chemical_Signal_Points[i] = gmsh.model.geo.add_point(Concentration_Point_Coordinates[0][i], Concentration_Point_Coordinates[1][i], 0, lc)
                    Concentration_Point_flags[i] = 1
            elif (Concentration_Point_flags[i] == 1) and (t < Concentration_Point_Coordinates[2][i] or t > Concentration_Point_Coordinates[3][i]):
                gmsh.model.geo.remove([(0, Chemical_Signal_Points[i])])
                Chemical_Signal_Points[i] = 0
                Concentration_Point_flags[i] = 0
        
        update_cell_coordinates(Cell_Nucleus_Coordinates, Cell_Membrane_Coordinates, Center_Cell_Points, Inner_Cell_Points, Outer_Cell_Points, Outer_Cell_Lines, Inner_Cell_Lines)
        write_mesh_to_vtk(f"time_{int(PV_TIME)}.vtk")
        PV_TIME += 1
    
    
    calculate_lj_forces_vectorized(Cell_Membrane_Coordinates, epsilon, sigma, dt, Cell_Amount, Node_Amount)


    #Iterate for each cell
    for cell_number in range(1, Cell_Amount + 1):
        original_c_x = Original_Cell_Nucleus_Coordinates[0][cell_number][0]
        original_c_y = Original_Cell_Nucleus_Coordinates[1][cell_number][0]
        Mx = solve_Mx(Node_Amount, Cell_Membrane_Coordinates, Cell_Nucleus_Coordinates, cell_number, rho_1, rho_2)
        My = solve_My(Node_Amount, Cell_Membrane_Coordinates, Cell_Nucleus_Coordinates, cell_number, rho_1, rho_2)
        A1 = solve_A(Node_Amount, Cell_Membrane_Coordinates, Cell_Nucleus_Coordinates, cell_number, rho_1, rho_2)
        
        Cell_Area[cell_number].append(A1)
    
    
        Cell_Nucleus_Coordinates[0][cell_number][0], Cell_Nucleus_Coordinates[1][cell_number][0] = solve_center_of_mass(Mx, My, A1)
        c_x = Cell_Nucleus_Coordinates[0][cell_number][0]
        c_y = Cell_Nucleus_Coordinates[1][cell_number][0]
        
        centroid_x[cell_number].append(c_x)
        centroid_y[cell_number].append(c_y)
        
        F_Brown_x = sigma_x * numpy.random.normal()
        F_Brown_y = sigma_y * numpy.random.normal()
        
        #Iterate over each node of each cell
        for i in range(1, Node_Amount + 1):
            #Iterate for each concentration signal
            
            #Calculate boundary conditions
            d = numpy.sqrt((Cell_Membrane_Coordinates[0][cell_number][i] - x_center)**2 + (Cell_Membrane_Coordinates[1][cell_number][i]-y_center)**2)
            
            if(d > R_Boundary):
                F_xBound = -k_boundary * (d-R_Boundary) * (Cell_Membrane_Coordinates[0][cell_number][i] - x_center)/d
                F_yBound = -k_boundary * (d-R_Boundary) * (Cell_Membrane_Coordinates[1][cell_number][i] - y_center)/d
            else:
                F_xBound = 0
                F_yBound = 0
               
            # Calculate total chemoattractant force
            for j in range(0, Amount_Concentration_Signal):
                if t >= Concentration_Point_Coordinates[2][j] and t <= Concentration_Point_Coordinates[3][j]:
                    grad_cx, grad_cy = solve_grad_c(Cell_Membrane_Coordinates[0][cell_number][i], Concentration_Point_Coordinates[0][j], Cell_Membrane_Coordinates[1][cell_number][i], Concentration_Point_Coordinates[1][j], t - Concentration_Point_Coordinates[2][j], mu, teps, A, Concentration_Point_Coordinates[4][j])                    
                    dcx = dcx + grad_cx
                    dcy = dcy + grad_cy

                
            # Read the cell membrane node coordinates
            m_x = Cell_Membrane_Coordinates[0][cell_number][i]
            m_y = Cell_Membrane_Coordinates[1][cell_number][i]
            n_x = Cell_Nucleus_Coordinates[0][cell_number][i]
            n_y = Cell_Nucleus_Coordinates[1][cell_number][i]
            
            ##Read the original lcell membrane node coordinates
            original_m_x = Original_Cell_Membrane_Coordinates[0][cell_number][i]
            original_m_y = Original_Cell_Membrane_Coordinates[1][cell_number][i]
            original_n_x = Original_Cell_Nucleus_Coordinates[0][cell_number][i]
            original_n_y = Original_Cell_Nucleus_Coordinates[1][cell_number][i]
            
            #If i == 1 then we remember previous node will be node N
            if i == 1:
                prev_x = Cell_Membrane_Coordinates[0][cell_number][Node_Amount]
                prev_y = Cell_Membrane_Coordinates[1][cell_number][Node_Amount]
                original_prev_x = Original_Cell_Membrane_Coordinates[0][cell_number][Node_Amount]
                original_prev_y = Original_Cell_Membrane_Coordinates[1][cell_number][Node_Amount]
            else:
                prev_x = Cell_Membrane_Coordinates[0][cell_number][i-1]
                prev_y = Cell_Membrane_Coordinates[1][cell_number][i-1]
                original_prev_x = Original_Cell_Membrane_Coordinates[0][cell_number][i-1]
                original_prev_y = Original_Cell_Membrane_Coordinates[1][cell_number][i-1]
            #If i == node_amount then we remember next node will be node 1.
            if i == Node_Amount:
                next_x = Cell_Membrane_Coordinates[0][cell_number][1]
                next_y = Cell_Membrane_Coordinates[1][cell_number][1]
                original_next_x = Original_Cell_Membrane_Coordinates[0][cell_number][1]
                original_next_y = Original_Cell_Membrane_Coordinates[1][cell_number][1]
            else:
                next_x = Cell_Membrane_Coordinates[0][cell_number][i+1]
                next_y = Cell_Membrane_Coordinates[1][cell_number][i+1]
                original_next_x = Original_Cell_Membrane_Coordinates[0][cell_number][i+1]
                original_next_y = Original_Cell_Membrane_Coordinates[1][cell_number][i+1]

            #Explicit time step on node summing all forces.
            Cell_Membrane_Coordinates[0][cell_number][i] += (Beta * dcx * dt +
                                                            alpha * (n_x - m_x - (original_n_x - original_m_x)) * dt + #
                                                            Delta * (original_m_x - original_next_x - (original_next_x - original_m_x) - (original_prev_x - original_m_x)  + (original_m_x - original_prev_x)) * dt +
                                                            F_xBound * dt +
                                                            F_Brown_x * numpy.sqrt(dt) ) 

            Cell_Membrane_Coordinates[1][cell_number][i] += (Beta * dcy * dt +
                                                            alpha * (n_y - m_y - (original_n_y - original_m_y)) * dt +
                                                            Delta * (original_m_y - original_next_y - (original_next_y - original_m_y) - (original_prev_y - original_m_y)  + (original_m_y - original_prev_y)) * dt +
                                                            F_yBound * dt +
                                                            F_Brown_y * numpy.sqrt(dt) )
            
            # Update nucleus coordinates
                #If i == 1 then we remember previous node will be node N
            if i == 1:
                prev_n_x = Cell_Nucleus_Coordinates[0][cell_number][Node_Amount]
                prev_n_y = Cell_Nucleus_Coordinates[1][cell_number][Node_Amount]
                original_prev_n_x = Original_Cell_Nucleus_Coordinates[0][cell_number][Node_Amount]
                original_prev_n_y = Original_Cell_Nucleus_Coordinates[1][cell_number][Node_Amount]
            else:
                prev_n_x = Cell_Nucleus_Coordinates[0][cell_number][i-1]
                prev_n_y = Cell_Nucleus_Coordinates[1][cell_number][i-1]
                original_prev_n_x = Original_Cell_Nucleus_Coordinates[0][cell_number][i-1]
                original_prev_n_y = Original_Cell_Nucleus_Coordinates[1][cell_number][i-1]

            #If i == node_amount then we remember next node will be node 1.
            if i == Node_Amount:
                next_n_x = Cell_Nucleus_Coordinates[0][cell_number][1]
                next_n_y = Cell_Nucleus_Coordinates[1][cell_number][1]
                original_next_n_x = Original_Cell_Nucleus_Coordinates[0][cell_number][1]
                original_next_n_y = Original_Cell_Nucleus_Coordinates[1][cell_number][1]
            else:
                next_n_x = Cell_Nucleus_Coordinates[0][cell_number][i+1]
                next_n_y = Cell_Nucleus_Coordinates[1][cell_number][i+1]
                original_next_n_x = Original_Cell_Nucleus_Coordinates[0][cell_number][i+1]
                original_next_n_y = Original_Cell_Nucleus_Coordinates[1][cell_number][i+1]
                
            
            #Get original values of difference between nodes
            original_hat_y_x = (original_n_x - original_next_n_x)
            original_hat_y_y = (original_n_y - original_next_n_y)

            original_hat_x_x = (original_n_x - original_c_x)
            original_hat_x_y = (original_n_y - original_c_y)

            #Explicit time step
            Cell_Nucleus_Coordinates[0][cell_number][i] += (alpha_n * (c_x - n_x + original_hat_x_x) * dt -
                                                            alpha * (n_x - m_x - (original_n_x - original_m_x)) * dt +
                                                            delta_n * (original_hat_y_x - (original_next_n_x - original_n_x) - (original_prev_n_x - original_n_x) + (original_n_x - original_prev_n_x)) * dt)

            Cell_Nucleus_Coordinates[1][cell_number][i] += (alpha_n * (c_y - n_y + original_hat_x_y) * dt -
                                                            alpha * (n_y - m_y - (original_n_y - original_m_y)) * dt +
                                                            delta_n * (original_hat_y_y - (original_next_n_y - original_n_y) - (original_prev_n_y - original_n_y) + (original_n_y - original_prev_n_y)) * dt)
            
            dcx = 0
            dcy = 0
        #Append for plotting meancurv
        MeanCurv = calculate_mean_curvature(Cell_Membrane_Coordinates, Cell_Nucleus_Coordinates, cell_number, Node_Amount)
        Mean_Curvature[cell_number].append(MeanCurv)
    time_values.append(t)
    t += dt
    
    """
# Set the size and DPI of the figure for a two-column format
plt.figure(figsize=(3.5, 2.5), dpi=300)

# Loop through each cell and plot their mean curvature
for cell_number in range(1, Cell_Amount + 1):
    plt.plot(time_values, Mean_Curvature[cell_number], label=f'Cell {cell_number}')

# Add labels with LaTeX formatting, title, legend, and grid
plt.xlabel('Time, $t$ (h)', fontsize=10)
plt.ylabel('Mean Curvature, $k$', fontsize=10)
plt.title('Mean Curvature $k$ as a Function of Time for All Cells', fontsize=12)
plt.legend(loc='best', fontsize=8)
plt.grid(True)

# Improve the style of the grid and axes
plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.minorticks_on()
plt.tick_params(axis='both', which='major', labelsize=8)

# Show the plot
plt.tight_layout()
plt.show()
    
plt.figure(figsize=(3.5, 2.5), dpi=300)
for cell_number in range(1, Cell_Amount + 1):
    plt.plot(time_values, Cell_Area[cell_number], label=f'Cell {cell_number}')

# Add labels with LaTeX formatting, title, legend, and grid
plt.xlabel('Time, $t$ (h)', fontsize=10)
plt.ylabel('Area $A$', fontsize=10)
plt.title('Area $A$ as a Function of Time for All Cells', fontsize=12)
plt.legend(loc='best', fontsize=8)
plt.grid(True)

# Improve the style of the grid and axes
plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.minorticks_on()
plt.tick_params(axis='both', which='major', labelsize=8)

# Show the plot
plt.tight_layout()
plt.show()
    
    
    
    """
    
    
# Import necessary libraries
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm

# Define the boundary radius
R_Boundary = 10  # Replace 10 with your desired radius

# Set the size and DPI of the figure for a two-column format
fig, ax = plt.subplots(figsize=(3.5, 2.5), dpi=300)  # Explicitly define the figure and axes

# Create a colormap for the time gradient
colors = cm.viridis(np.linspace(0, 1, len(time_values)))

# Loop through each cell and plot their centroid movement
for cell_number in range(1, Cell_Amount + 1):
    x_values = centroid_x[cell_number]
    y_values = centroid_y[cell_number]
    
    # Plot the movement of the centroid using a color gradient for time
    for i in range(1, len(time_values)):
        ax.plot(x_values[i-1:i+1], y_values[i-1:i+1], '-', color=colors[i])

# Add a circular boundary
circle = plt.Circle((0, 0), R_Boundary, color='black', fill=False, linestyle='--', linewidth=0.8)
ax.add_patch(circle)

# Set axis limits to match the boundary
ax.set_xlim(-R_Boundary, R_Boundary)
ax.set_ylim(-R_Boundary, R_Boundary)

# Maintain equal aspect ratio for proper circular shape
ax.set_aspect('equal', adjustable='box')

# Add labels with LaTeX formatting, title, and grid
ax.set_xlabel('Centroid X Coordinate, $x$', fontsize=10)
ax.set_ylabel('Centroid Y Coordinate, $y$', fontsize=10)
ax.set_title('Movement of Cell Centroids Over Time', fontsize=12)
ax.grid(True)

# Improve the style of the grid and axes
ax.grid(which='both', linestyle='--', linewidth=0.5)
ax.minorticks_on()
ax.tick_params(axis='both', which='major', labelsize=8)

# Add a color bar to represent time
norm = plt.Normalize(vmin=min(time_values), vmax=max(time_values))
sm = plt.cm.ScalarMappable(cmap=cm.viridis, norm=norm)
sm.set_array([])  # ScalarMappable requires a set array
cbar = fig.colorbar(sm, ax=ax)  # Explicitly link the colorbar to the axes
cbar.set_label('Time', fontsize=10)

# Adjust layout and show the plot
plt.tight_layout()
plt.show()
