import gmsh
import sys
import math
import random
import numpy

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
#cell_points ->
def Connect_Mesh_Points(N, cell_points):
    line_dict = {}
    for i in range(N):
        if (i + 1 != N):
            line_dict[str(i+1)] = gmsh.model.geo.add_line(cell_points[str(i+1)], cell_points[str(i+2)])
        else:
            line_dict[str(i+1)] = gmsh.model.geo.add_line(cell_points[str(i+1)], cell_points[str(1)])
    return line_dict

def create_center_to_inner_lines(N, center_point_id, inner_points):
    lines = {}
    for i in range(1, N + 1):
        lines[str(i)] = gmsh.model.geo.add_line(center_point_id, inner_points[str(i)])
    return lines

def create_inner_to_outer_radial_lines(N, inner_points, outer_points):
    lines = {}
    for i in range(1, N + 1):
        lines[str(i)] = gmsh.model.geo.add_line(inner_points[str(i)], outer_points[str(i)])
    return lines

# Function to generate random offsets ensuring no overlap for initial cell placement 
def generate_random_offset(used_offsets, spacing, max_attempts=100):
    Size_Of_Grid = 25
    for _ in range(max_attempts):
        offset_x = random.uniform(-Size_Of_Grid, Size_Of_Grid)
        offset_y = random.uniform(-Size_Of_Grid, Size_Of_Grid)
        if all(math.sqrt((offset_x - x)**2 + (offset_y - y)**2) >= spacing for x, y in used_offsets):
            used_offsets.append((offset_x, offset_y))
            return offset_x, offset_y
    raise ValueError("Could not find a non-overlapping position after multiple attempts")

#This function will take in a specified cell number and node number alongside the nodes updated x and y positions and then redraw the cell
def update_inner_cell_shape(cell_number, node_number, Cell_Points, Cell_Lines, Center_ID, Cell_Radial_Lines, new_x, new_y):

    for line_id in Cell_Lines[str(cell_number)].values():
        gmsh.model.geo.remove([(1, line_id)])  
    for line_id in Cell_Radial_Lines[str(cell_number)].values():
        gmsh.model.geo.remove([(1, line_id)])  

    old_point_id = Cell_Points[str(cell_number)][str(node_number)]
    gmsh.model.geo.remove([(0, old_point_id)])
    
    # Add a new point with updated coordinates
    new_point_id = gmsh.model.geo.add_point(new_x, new_y, 0, lc)
    # Update the Cell_Points dictionary with the new point ID
    Cell_Points[str(cell_number)][str(node_number)] = new_point_id
    
    Cell_Lines[str(cell_number)] = Connect_Mesh_Points(Node_Amount, Cell_Points[str(cell_number)])
    Cell_Radial_Lines[str(cell_number)] = create_center_to_inner_lines(Node_Amount, Center_ID, Cell_Points[str(cell_number)])
    

# This function will take in a specified cell number and node number alongside the nodes updated x and y positions and then redraw the outer cell
def update_outer_cell_shape(cell_number, node_number, Cell_Points, Cell_Lines, Inner_Cell_Points, Inner_Cell_Radial_Lines, new_x, new_y):
    for line_id in Cell_Lines[str(cell_number)].values():
        gmsh.model.geo.remove([(1, line_id)])  
    for line_id in Inner_Cell_Radial_Lines[str(cell_number)].values():
        gmsh.model.geo.remove([(1, line_id)])  
        
    old_point_id = Cell_Points[str(cell_number)][str(node_number)]
    gmsh.model.geo.remove([(0, old_point_id)])
    
    # Add a new point with updated coordinates
    new_point_id = gmsh.model.geo.add_point(new_x, new_y, 0, lc)
    # Update the Cell_Points dictionary with the new point ID
    Cell_Points[str(cell_number)][str(node_number)] = new_point_id
    

    Cell_Lines[str(cell_number)] = Connect_Mesh_Points(Node_Amount, Cell_Points[str(cell_number)])
    Inner_Cell_Radial_Lines[str(cell_number)] = create_inner_to_outer_radial_lines(Node_Amount, Inner_Cell_Points[str(cell_number)], Cell_Points[str(cell_number)])
    
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

def solve_center_of_mass(Mx,My,A):
    x_c = Mx/A
    y_c = My/A
    return x_c, y_c

def solve_c(x, xs, y, ys, t, mu, teps, A):
    c = A/(mu*(t+teps))*math.exp(-((x-xs)**2+(y-ys)**2)/(4*mu*(t+teps)))
    return c

def solve_grad_c(x, xs, y, ys, t, mu, teps, A):
    cx = -(A*(x-xs))/(2*mu**2*(t+teps)**2)*math.exp(-((x-xs)**2+(y-ys)**2)/(4*mu*(t+teps)))
    cy = -(A*(y-ys))/(2*mu**2*(t+teps)**2)*math.exp(-((x-xs)**2+(y-ys)**2)/(4*mu*(t+teps)))
    return cx, cy

def run_gmesh():
    # Create the relevant Gmsh data structures from Gmsh model
    gmsh.model.geo.synchronize()

    # Generate mesh
    gmsh.model.mesh.generate()

    # Write mesh data
    gmsh.write("GFG.msh")

    # Creates graphical user interface
    if 'close' not in sys.argv:
        gmsh.fltk.run()

    # Finalize the Gmsh API
    gmsh.finalize()
    

    
def test_center_of_mass(node_amount, membrane_cell_coordinates, nucleus_cell_coordinates, cell_number, rho_1, rho_2):
    Mx = solve_Mx(Node_Amount, Cell_Membrane_Coordinates, Cell_Nucleus_Coordinates, 1, rho_1, rho_2)
    My = solve_My(Node_Amount, Cell_Membrane_Coordinates, Cell_Nucleus_Coordinates, 1, rho_1, rho_2)
    A1 = solve_A(Node_Amount, Cell_Membrane_Coordinates, Cell_Nucleus_Coordinates, 1, rho_1, rho_2)
    xc, yc = solve_center_of_mass(Mx,My, A1)
    print(Cell_Nucleus_Coordinates[0][1][0])
    print(Cell_Nucleus_Coordinates[1][1][0])
    print(xc)
    print(yc)
    
# Function to update the coordinates of existing points
def update_cell_coordinates(cell_nucleus, cell_membrane, Cell_Points):
    for cell_number, nucleus_coords in enumerate(cell_nucleus, start=1):
        for node_number, (x, y) in enumerate(nucleus_coords, start=1):
            point_id = Cell_Points[str(cell_number)][str(node_number)]
            gmsh.model.geo.mesh.movePoint(point_id, x, y, 0)
    
    for cell_number, membrane_coords in enumerate(cell_membrane, start=1):
        for node_number, (x, y) in enumerate(membrane_coords, start=1):
            point_id = Cell_Points[str(cell_number)][str(node_number)]
            gmsh.model.geo.mesh.movePoint(point_id, x, y, 0)
            


# Circle points:
lc = 1e-2
Node_Amount = 5
Cell_Amount = 1
rho_1 = 1.2
rho_2 = 1.4
mu = 2
teps = 1*10**(-6)
xs = 1
ys = 0

A = 15
Beta = 2.5
alpha = 0.5
Delta = 30
delta_n = 30*200
alpha_n = 100
n = 10000

# Create a dictionary to store points of all cells
Center_Cell_Points = {}
Outer_Cell_Points = {}
Inner_Cell_Points = {}
Outer_Cell_Lines = {}
Inner_Cell_Lines = {}
Center_to_Inner_Lines = {}
Inner_to_Outer_Lines = {}



# Parameters for cell spacing to avoid overlap
cell_spacing = 3  # Minimum distance between cell centers to avoid overlap
used_offsets = []
Cell_Membrane_Coordinates = numpy.zeros((2, Cell_Amount + 1, Node_Amount + 1)) #x dimension is coordinates y dimension is cell # and z dimension is node #
Cell_Nucleus_Coordinates = numpy.zeros((2, Cell_Amount + 1, Node_Amount + 1))

# Generate mesh for each cell with a random offset
for i in range(1, Cell_Amount + 1):
    offset_x, offset_y = generate_random_offset(used_offsets, cell_spacing)
    Center_Cell_Points[str(i)] = gmsh.model.geo.add_point(offset_x, offset_y, 0, lc)
    Outer_Cell_Points[str(i)] = Create_Cell_Mesh(Node_Amount, lc, offset_x, offset_y, Cell_Membrane_Coordinates, i, 1)
    Inner_Cell_Points[str(i)] = Create_Cell_Mesh(Node_Amount, lc, offset_x, offset_y, Cell_Nucleus_Coordinates, i, 0.25)
   # Center_to_Inner_Lines[str(i)] = create_center_to_inner_lines(Node_Amount, Center_Cell_Points[str(i)], Inner_Cell_Points[str(i)])
   # Inner_to_Outer_Lines[str(i)] = create_inner_to_outer_radial_lines(Node_Amount, Inner_Cell_Points[str(i)], Outer_Cell_Points[str(i)])
    

# Connect mesh points for each cell
for i in range(1, Cell_Amount + 1):
    Outer_Cell_Lines[str(i)] = Connect_Mesh_Points(Node_Amount, Outer_Cell_Points[str(i)])
    Inner_Cell_Lines[str(i)] = Connect_Mesh_Points(Node_Amount, Inner_Cell_Points[str(i)])
    

cell_number = 1

#dt = 0.01
dt = 2
t = 0
T_Max = 20

#The center is blowing up probably due to coordinates themselves blowing up? stops working at dt = 1.
while t < T_Max:
    for cell_number in range(1, Cell_Amount + 1):
        
        Mx = solve_Mx(Node_Amount, Cell_Membrane_Coordinates, Cell_Nucleus_Coordinates, cell_number, rho_1, rho_2)
        My = solve_My(Node_Amount, Cell_Membrane_Coordinates, Cell_Nucleus_Coordinates, cell_number, rho_1, rho_2)
        A1 = solve_A(Node_Amount, Cell_Membrane_Coordinates, Cell_Nucleus_Coordinates, cell_number, rho_1, rho_2)
        Cell_Nucleus_Coordinates[0][cell_number][0], Cell_Nucleus_Coordinates[1][cell_number][0] = solve_center_of_mass(Mx,My,A1)
        
        c_x = Cell_Nucleus_Coordinates[0][cell_number][0]
        c_y = Cell_Nucleus_Coordinates[1][cell_number][0]
        test_center_of_mass(Node_Amount, Cell_Membrane_Coordinates, Cell_Nucleus_Coordinates, cell_number, rho_1, rho_2)
        for i in range(1, Node_Amount + 1):    
            
            c = solve_c(Cell_Membrane_Coordinates[0][cell_number][i], xs, Cell_Membrane_Coordinates[1][cell_number][i], ys, t, mu, teps, A)
            dcx, dcy = solve_grad_c(Cell_Membrane_Coordinates[0][cell_number][i], xs, Cell_Membrane_Coordinates[1][cell_number][i], ys, t, mu, teps, A)
            
            # Update cell membrane coordinates
            m_x = Cell_Membrane_Coordinates[0][cell_number][i]
            m_y = Cell_Membrane_Coordinates[1][cell_number][i]
            n_x = Cell_Nucleus_Coordinates[0][cell_number][i]
            n_y = Cell_Nucleus_Coordinates[1][cell_number][i]
            
            if i == 1:
                prev_x = Cell_Membrane_Coordinates[0][cell_number][Node_Amount]
                prev_y = Cell_Membrane_Coordinates[1][cell_number][Node_Amount]
            else:
                prev_x = Cell_Membrane_Coordinates[0][cell_number][i-1]
                prev_y = Cell_Membrane_Coordinates[1][cell_number][i-1]

            if i == Node_Amount:
                next_x = Cell_Membrane_Coordinates[0][cell_number][1]
                next_y = Cell_Membrane_Coordinates[1][cell_number][1]
            else:
                next_x = Cell_Membrane_Coordinates[0][cell_number][i+1]
                next_y = Cell_Membrane_Coordinates[1][cell_number][i+1]

            Cell_Membrane_Coordinates[0][cell_number][i] += Beta * dcx * dt + alpha * (n_x - m_x) * dt + Delta * (m_x - prev_x + m_x - next_x) * dt
            Cell_Membrane_Coordinates[1][cell_number][i] += Beta * dcy * dt + alpha * (n_y - m_y) * dt + Delta * (m_y - prev_y + m_y - next_y) * dt

            # Update nucleus coordinates
            if i == 1:
                prev_n_x = Cell_Nucleus_Coordinates[0][cell_number][Node_Amount]
                prev_n_y = Cell_Nucleus_Coordinates[1][cell_number][Node_Amount]
            else:
                prev_n_x = Cell_Nucleus_Coordinates[0][cell_number][i-1]
                prev_n_y = Cell_Nucleus_Coordinates[1][cell_number][i-1]

            if i == Node_Amount:
                next_n_x = Cell_Nucleus_Coordinates[0][cell_number][1]
                next_n_y = Cell_Nucleus_Coordinates[1][cell_number][1]
            else:
                next_n_x = Cell_Nucleus_Coordinates[0][cell_number][i+1]
                next_n_y = Cell_Nucleus_Coordinates[1][cell_number][i+1]

            hat_y_x = n_x - prev_n_x
            hat_y_y = n_y - prev_n_y

            Cell_Nucleus_Coordinates[0][cell_number][i] += alpha_n * (c_x - n_x - (n_x - c_x )) * dt - alpha * (n_x - m_x) * dt + delta_n * (hat_y_x - (n_x - next_n_x)) * dt
            Cell_Nucleus_Coordinates[1][cell_number][i] += alpha_n * (c_y - n_y - (n_y - c_y)) * dt - alpha * (n_y - m_y) * dt + delta_n * (hat_y_y - (n_y - next_n_y)) * dt

    t += dt
    
# gmsh.model.geo.synchronize() 


# # Define new coordinates for the point
# new_x = 1.0
# new_y = 2.0
# #update_inner_cell_shape(1,1, Inner_Cell_Points, Inner_Cell_Lines, Center_Cell_Points[str(1)], Center_to_Inner_Lines, new_x, new_y)
# #update_outer_cell_shape(1, 1, Outer_Cell_Points, Outer_Cell_Lines, Inner_Cell_Points, Inner_to_Outer_Lines, 1.5, 1.5)



run_gmesh()

