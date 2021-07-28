# ------------------------------------------------------------------------------
#
#  Gmsh Python 
#
# ------------------------------------------------------------------------------

# The Python API is entirely defined in the `gmsh.py' module (which contains the
# full documentation of all the functions in the API):
import gmsh
import sys
import numpy as np

# Before using any functions in the Python API, Gmsh must be initialized:
gmsh.initialize()

# Add a new model 
gmsh.model.add("gmsh_mesh")

# The Python API provides direct access to each supported geometry (CAD)
# kernel. The built-in kernel is used in this first tutorial: the corresponding
# API functions have the `gmsh.model.geo' prefix.

# Variables
# -------------
# Building dimensions
B_BUILDING = 32
H_BUILDING = 160

# Wind tunnel bounding box dimensions
L_WT_1  =  8 * H_BUILDING
L_WT_2  = 16 * H_BUILDING
B_WT    = 10 * H_BUILDING
H_WT    =  6 * H_BUILDING

# Mesh Elements
# "Radial" from the Building inner perim.
n_Ri = 20 # 20
# Progression in radial direction inner perim.
p_Ri = 1.1
# "Radial" from the Building outter perim.
n_Ro = 20 # 20
# Progression in radial direction outter perim.
p_Ro = 1.05
# "Tangential" along the Building
n_T = 20 # 20
# Progression in tangential direction
p_T = 1.1

# ABL-Extrusion
# Number of cells
n_abl = 20 # 20 
# Total expansion ratio delta_n / delta_1
r_abl = 1.15 #10 

# Extrusion above building
# Number of cells
n_Ra = 20 # 20
# Total expansion ratio delta_n / delta_1
r_Ra = 1.1 #10 


# Add vertices
# -------------
# gmsh.model.geo.addPoint():
# - the first 3 arguments are the point coordinates (x, y, z)
# - the next (optional) argument is the target mesh size close to the point
# - the last (optional) argument is the point tag (a stricly positive integer
#   that uniquely identifies the point)
# Center Point
p00 = gmsh.model.geo.addPoint(0                         , 0                         , 0 )

# block 01 vertices
p01 = gmsh.model.geo.addPoint(B_BUILDING/2	            , B_BUILDING/2		        , 0 )
p02 = gmsh.model.geo.addPoint(0			                , B_BUILDING/2		        , 0 )
p03 = gmsh.model.geo.addPoint(0			                , H_BUILDING 		        , 0 )
p04 = gmsh.model.geo.addPoint(np.sqrt(2)/2*H_BUILDING   , np.sqrt(2)/2*H_BUILDING   , 0 )

# block 02 new vertices
p09 = gmsh.model.geo.addPoint(0                         , B_WT/2                    , 0 )
p10 = gmsh.model.geo.addPoint(B_WT/2                    , B_WT/2                    , 0 )


# Add Lines
# -------------
# Straight lines =  {Start, End}
# Circle         =  {Start, Center, End}
# block 01 lines
line01  = gmsh.model.geo.addLine(p01, p02) 
line02  = gmsh.model.geo.addLine(p02, p03)
circle03= gmsh.model.geo.addCircleArc(p03, p00, p04) 
line04  = gmsh.model.geo.addLine(p01, p04)

# block 02 lines
line09  = gmsh.model.geo.addLine(p03, p09)
line10  = gmsh.model.geo.addLine(p09, p10)
line11  = gmsh.model.geo.addLine(p04, p10)

# Curve loops
# -------------
loop01 = gmsh.model.geo.addCurveLoop([line01, line02, circle03, -line04])
loop02 = gmsh.model.geo.addCurveLoop([-circle03, line09, line10, -line11])

# Surfaces
# -------------
surface01 = gmsh.model.geo.addPlaneSurface([-loop01])
surface02 = gmsh.model.geo.addPlaneSurface([-loop02])

# Before they can be meshed (and, more generally, before they can be used by API
# functions outside of the built-in CAD kernel functions), the CAD entities must
# be synchronized with the Gmsh model.
gmsh.model.geo.synchronize()

# Meshing constraints
# -------------

# We put nodes following a geometric progression on the curves
# Transfinite lines on inner circle
gmsh.model.geo.mesh.setTransfiniteCurve(line01 , n_T, meshType="Progession", coef=1.0)
gmsh.model.geo.mesh.setTransfiniteCurve(line04, n_Ri, meshType="Progession", coef=p_Ri)
gmsh.model.geo.mesh.setTransfiniteCurve(circle03, n_T, meshType="Progession", coef=1.0)
gmsh.model.geo.mesh.setTransfiniteCurve(line02 , n_Ri, meshType="Progession", coef=p_Ri)

# Transfinite lines on outter perimeter
gmsh.model.geo.mesh.setTransfiniteCurve(line10, n_T, meshType="Progession", coef=p_T)
gmsh.model.geo.mesh.setTransfiniteCurve(line09, n_Ro, meshType="Progession", coef=p_Ro)
gmsh.model.geo.mesh.setTransfiniteCurve(line11, n_Ro, meshType="Progession", coef=p_Ro)

gmsh.model.geo.synchronize()

points  = gmsh.model.get_entities(0)      # Get Points
curves  = gmsh.model.get_entities(1)      # Get Curves
surfaces= gmsh.model.get_entities(2)      # Get Surfaces

for surface in surfaces:
    # The `setTransfiniteSurface()' meshing constraint uses a transfinite
    # interpolation algorithm in the parametric plane of the surface to connect the
    # nodes on the boundary using a structured grid.
    gmsh.model.geo.mesh.setTransfiniteSurface(surface[1], "Left")

    # To create quadrangles instead of triangles, one can use the `setRecombine'
    # constraint:
    gmsh.model.geo.mesh.setRecombine(dim=2, tag=surface[1])

# Set physical groups
# -------------

# Mirror all entities
# -------------
# Important: Copy meshing method (unstructured or transfinite) when duplicating 
# geometrical entities with built-in geometry kernel
gmsh.option.setNumber("Geometry.CopyMeshingMethod", 1)

def symmetrizeAndCopy(model, a, b, c, d):
    # Important: Copy meshing method (unstructured or transfinite) when duplicating 
    # geometrical entities with built-in geometry kernel
    gmsh.option.setNumber("Geometry.CopyMeshingMethod", 1)
    
    # Select all entities
    entities = model.get_entities()

    # Return copy 
    copied_entities = model.geo.copy(entities)
    # For ABL-Extrusion, try the following out
    # for copy in copied_entities:
    #     if copy[0] == 2:
    #         gmsh.model.geo.mesh.setReverse(copy[0], copy[1])
    
    # Reverse surface normals
    model.geo.symmetrize(copied_entities, a, b, c, d)

    # Synchronize the model
    model.geo.synchronize()

# Mirror along XY-Axis
symmetrizeAndCopy(gmsh.model, -1, 1, 0, 0)

# Mirror along Y-Axis
symmetrizeAndCopy(gmsh.model, -1, 0, 0, 0)

# Create Inlet Volume
# -------------
# block 01 vertices
pi01 = gmsh.model.geo.addPoint(-L_WT_1                   , B_WT/2		            , 0 )
pi02 = gmsh.model.geo.addPoint(-B_WT/2                   , B_WT/2		            , 0 )
pi03 = gmsh.model.geo.addPoint(-B_WT/2                   , 0	    	            , 0 )
pi04 = gmsh.model.geo.addPoint(-L_WT_1                   , 0 	    	            , 0 )

linei01  = gmsh.model.geo.addLine(pi01, pi02) 
linei02  = gmsh.model.geo.addLine(pi02, pi03)
linei03  = gmsh.model.geo.addLine(pi04, pi03)
linei04  = gmsh.model.geo.addLine(pi04, pi01)

gmsh.model.geo.mesh.setTransfiniteCurve(linei02, n_T, meshType="Progession", coef=p_T)
gmsh.model.geo.mesh.setTransfiniteCurve(linei04, n_T, meshType="Progession", coef=p_T)
gmsh.model.geo.mesh.setTransfiniteCurve(linei01, int(n_Ro/2), meshType="Progession", coef=p_Ro)
gmsh.model.geo.mesh.setTransfiniteCurve(linei03, int(n_Ro/2), meshType="Progession", coef=p_Ro)

# Curve loops
# -------------
loopi01 = gmsh.model.geo.addCurveLoop([linei01, linei02, -linei03, linei04])

# Surfaces
# -------------
surfaceo01 = gmsh.model.geo.addPlaneSurface([-loopi01])

# Create Outlet Volume
# -------------
# block 01 vertices
po01 = gmsh.model.geo.addPoint( L_WT_2                   , B_WT/2		            , 0 )
po02 = gmsh.model.geo.addPoint( B_WT/2                   , B_WT/2		            , 0 )
po03 = gmsh.model.geo.addPoint( B_WT/2                   , 0	    	            , 0 )
po04 = gmsh.model.geo.addPoint( L_WT_2                   , 0 	    	            , 0 )

lineo01  = gmsh.model.geo.addLine(po02, po01) 
lineo02  = gmsh.model.geo.addLine(po02, po03)
lineo03  = gmsh.model.geo.addLine(po03, po04)
lineo04  = gmsh.model.geo.addLine(po04, po01)

gmsh.model.geo.mesh.setTransfiniteCurve(lineo02, n_T, meshType="Progession", coef=p_T)
gmsh.model.geo.mesh.setTransfiniteCurve(lineo04, n_T, meshType="Progession", coef=p_T)
gmsh.model.geo.mesh.setTransfiniteCurve(lineo01, n_Ro, meshType="Progession", coef=p_Ro)
gmsh.model.geo.mesh.setTransfiniteCurve(lineo03, n_Ro, meshType="Progession", coef=p_Ro)

# Curve loops
# -------------
loopo01 = gmsh.model.geo.addCurveLoop([-lineo01, lineo02, lineo03, lineo04])

# Surfaces
# -------------
surfaceo01 = gmsh.model.geo.addPlaneSurface([-loopo01])

gmsh.model.geo.synchronize()

surfaces= gmsh.model.get_entities(2)      # Get Surfaces

for surface in surfaces:
    # The `setTransfiniteSurface()' meshing constraint uses a transfinite
    # interpolation algorithm in the parametric plane of the surface to connect the
    # nodes on the boundary using a structured grid.
    gmsh.model.geo.mesh.setTransfiniteSurface(surface[1], "Left")

    # To create quadrangles instead of triangles, one can use the `setRecombine'
    # constraint:
    gmsh.model.geo.mesh.setRecombine(dim=2, tag=surface[1])

gmsh.model.geo.synchronize()

# Mirror along X-Axis
symmetrizeAndCopy(gmsh.model, 0, 1, 0, 0) 

# Create the building surface
building_curves = gmsh.model.getEntitiesInBoundingBox(-B_BUILDING/2, -B_BUILDING/2, 0, +B_BUILDING/2, +B_BUILDING/2, 0, dim=1)
building_curves_tags = []
for building_curve in building_curves:
    building_curves_tags.append(building_curve[1])

loop = gmsh.model.geo.addCurveLoop(building_curves_tags, reorient=True)
surface03 = gmsh.model.geo.addPlaneSurface([loop])

gmsh.model.geo.synchronize()

def createExtrusion(n, r, H):
    """
    # L       = H_BUILDING        # Total Length
    # n       = n_abl             # Number of cells / layers
    # r       = r_abl             # Total expansion ratio delta_n / delta_1
    """
    # List w/ extrusion
    heights   = np.zeros(n)     # Height of each extrusion layer
    cells     = np.zeros(n)     # Number of cells per layer

    # Calc expansion ratio
    beta    = r**(1/(n-1))      # Cell to cell expansion ratio

    # Solve for width of start cell
    # delta_1 = (1-beta) / (1-beta ** n) * H              # Block Mesh Grading
    delta_1 = (r - 1) / (r ** n - 1) * H                # GMSH Transfinite Grading

    # Starting values
    cells[0]    = 1           # Layers per cell   
    heights[0]  = delta_1   # First cell Height

    for i in range(1, n_abl):
        cells[i]    = 1
        delta_i     = delta_1 * beta ** i
        # heights[i]  = heights[i-1] + delta_i            # Block Mesh Grading
        heights[i]  = heights[i-1] + delta_1 * r **i    # GMSH Transfinite Grading

    # Normalize heights to 1
    heights = heights / heights[-1]

    return heights, cells

# Heights & cells for ABL Extrusion from 0 - H_BUILDING
heights_abl, cells_abl = createExtrusion(n_abl, r_abl, H_BUILDING)

# Height & cells for Extrusion above Building from H_BUILDING - H_WT
heights_above, cells_above = createExtrusion(n_Ra, r_Ra, H_WT - H_BUILDING)

# Get all Surfaces
surfaces = gmsh.model.get_entities(2)     
 
# Extrude the surfaces
# Important note: Do not use the gmsh.model.geo.extrudeBoundaryLayer() command here
# it will create the volumes along the mesh normals
# gmsh.model.geo.extrudeBoundaryLayer(surfaces,cells,heights, recombine=True)
gmsh.model.geo.extrude(surfaces, 0, 0, H_BUILDING, cells_abl, heights_abl, recombine=True)
gmsh.model.geo.synchronize()

# Get top surfaces
surfaces = gmsh.model.getEntitiesInBoundingBox(-L_WT_1, -B_WT/2, H_BUILDING, L_WT_2, B_WT/2, H_BUILDING, dim=2)

# Extrude the surfaces
gmsh.model.geo.extrude(surfaces, 0, 0, H_WT-H_BUILDING, cells_above, heights_above, recombine=True)
gmsh.model.geo.synchronize()

# Remove volume inside building
# ---------
# Vol = gmsh.model.getEntitiesInBoundingBox(-B_BUILDING/2, -B_BUILDING/2, 0, B_BUILDING/2, B_BUILDING/2, H_BUILDING, dim=3)

# Procedure using symmetry
# -------
# gmsh.option.setNumber("Geometry.Normals", 100)


def transfiniteSurface(z_coord):
    # Surface
    surface  = gmsh.model.getEntitiesInBoundingBox(-B_BUILDING/2, -B_BUILDING/2, z_coord,  B_BUILDING/2,  B_BUILDING/2, z_coord, dim=2)[0][1]

    # Corner points
    corner01 = gmsh.model.getEntitiesInBoundingBox(-B_BUILDING/2,  B_BUILDING/2, z_coord, -B_BUILDING/2,  B_BUILDING/2, z_coord, dim=0)[0][1]
    corner02 = gmsh.model.getEntitiesInBoundingBox( B_BUILDING/2,  B_BUILDING/2, z_coord,  B_BUILDING/2,  B_BUILDING/2, z_coord, dim=0)[0][1]
    corner03 = gmsh.model.getEntitiesInBoundingBox( B_BUILDING/2, -B_BUILDING/2, z_coord,  B_BUILDING/2, -B_BUILDING/2, z_coord, dim=0)[0][1]
    corner04 = gmsh.model.getEntitiesInBoundingBox(-B_BUILDING/2, -B_BUILDING/2, z_coord, -B_BUILDING/2, -B_BUILDING/2, z_coord, dim=0)[0][1]

    # The `setTransfiniteSurface()' meshing constraint uses a transfinite
    # interpolation algorithm in the parametric plane of the surface to connect the
    # nodes on the boundary using a structured grid.
    gmsh.model.geo.mesh.setTransfiniteSurface(surface, "Left", [corner01, corner02, corner03, corner04] )

    # To create quadrangles instead of triangles, one can use the `setRecombine'
    # constraint:
    gmsh.model.geo.mesh.setRecombine(dim=2, tag=surface)
    gmsh.model.geo.synchronize()


transfiniteSurface(0)
transfiniteSurface(H_BUILDING)
transfiniteSurface(H_WT)


# Physical Groups
# -------------

# At this level, Gmsh knows everything to display the rectangular surface 1 and
# to mesh it. An optional step is needed if we want to group elementary
# geometrical entities into more meaningful groups, e.g. to define some
# mathematical ("domain", "boundary"), functional ("left wing", "fuselage") or
# material ("steel", "carbon") properties.

# Such groups are called "Physical Groups" in Gmsh. By default, if physical
# groups are defined, Gmsh will export in output files only mesh elements that
# belong to at least one physical group. (To force Gmsh to save all elements,
# whether they belong to physical groups or not, set the `Mesh.SaveAll' option
# to 1.) Physical groups are also identified by tags, i.e. stricly positive
# integers, that should be unique per dimension (0D, 1D, 2D or 3D). Physical
# groups can also be given names.

# Here we define a physical curve that groups the left, bottom and right curves
# in a single group (with prescribed tag 5); and a physical surface with name
# "My surface" (with an automatic tag) containing the geometrical surface 1:

# Inlet
# a = gmsh.model.getBoundingBox([(1,-1)],[-L_WT_1, -B_WT, 0,  -L_WT_1, -B_WT, 0])

# I suggest you define your volume to be a "physical volume". Once any 
# physical items are present in the geometry, only those physical items 
# are exported in the mesh. I.e. if you don't mark the boundary as 
# physical as well, the triangles of the boundary will not be exported 
# separately.


surf_inlet   = gmsh.model.getEntitiesInBoundingBox(-L_WT_1, -B_WT/2, 0,     -L_WT_1,  B_WT/2, H_WT,  dim=2)
surf_outlet  = gmsh.model.getEntitiesInBoundingBox( L_WT_2, -B_WT/2, 0,      L_WT_2,  B_WT/2, H_WT,  dim=2)
surf_back    = gmsh.model.getEntitiesInBoundingBox(-L_WT_1,  B_WT/2, 0,      L_WT_2,  B_WT/2, H_WT,  dim=2)
surf_front   = gmsh.model.getEntitiesInBoundingBox(-L_WT_1, -B_WT/2, 0,      L_WT_2, -B_WT/2, H_WT,  dim=2)
surf_top     = gmsh.model.getEntitiesInBoundingBox(-L_WT_1, -B_WT/2, H_WT,   L_WT_2,  B_WT/2, H_WT,  dim=2)
surf_bottom  = gmsh.model.getEntitiesInBoundingBox(-L_WT_1, -B_WT/2, 0,      L_WT_2,  B_WT/2, 0,     dim=2)
surf_building= gmsh.model.getEntitiesInBoundingBox(-B_BUILDING/2, -B_BUILDING/2, 0, B_BUILDING/2, B_BUILDING/2, H_BUILDING, dim=2)
surf_building.remove((2, surface03))

vol_building = gmsh.model.getEntitiesInBoundingBox(-B_BUILDING/2, -B_BUILDING/2, 0, B_BUILDING/2, B_BUILDING/2, H_BUILDING, dim=3)
vol_internal = gmsh.model.getEntitiesInBoundingBox(-L_WT_1, -B_WT/2, 0, L_WT_2, B_WT/2, H_WT, dim=3)

vol_internal.remove(vol_building[0])


def setPhysicalName(model, dim_tags, dim, name):
    # for dim_tag in dim_tags:
    tags = list(zip(*dim_tags))[1]
    ps = gmsh.model.addPhysicalGroup(dim, tags)
    model.setPhysicalName(dim, ps, name)
    model.geo.synchronize()
        
setPhysicalName(gmsh.model, surf_inlet, 2, "inlet")    
setPhysicalName(gmsh.model, surf_outlet, 2, "outlet")    
setPhysicalName(gmsh.model, surf_back, 2, "back")    
setPhysicalName(gmsh.model, surf_front, 2, "front")    
setPhysicalName(gmsh.model, surf_bottom, 2, "bottom") 
setPhysicalName(gmsh.model, surf_top, 2, "top")    
setPhysicalName(gmsh.model, surf_building, 2, "building") 

setPhysicalName(gmsh.model, vol_internal, 3, "internal")

# Meshing
# -------------
# To create curvilinear elements mesh with second order elements;
# meaning adding an additional interpolation point per line
# MEMO: OF DOES NOT SUPPORT HIGHER ORDER SCHEMES
# gmsh.option.setNumber("Mesh.SecondOrderLinear", 0)
# gmsh.option.setNumber("Mesh.ElementOrder", 2)
# gmsh.option.setNumber("Mesh.HighOrderOptimize", 2)
# gmsh.option.setNumber("Mesh.HighOrderFastCurvingNewAlgo", 1)


# In order to compare with first order elements uncomment the following
gmsh.option.setNumber("Mesh.SecondOrderLinear", 0)
gmsh.option.setNumber("Mesh.ElementOrder", 1)

# Increase the number of sub-edges for a nicer display of the geometry:
gmsh.option.setNumber("Mesh.NumSubEdges", gmsh.option.get_number("Geometry.NumSubEdges"))

# Mesh Algorithm set to 
# 3D mesh algorithm (1: Delaunay, 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT)
# gmsh.option.setNumber("Mesh.Algorithm3D", 1)

# We can then generate a 3D mesh...
gmsh.model.mesh.generate(dim=3)

# We can then generate a 2D mesh...
# gmsh.model.mesh.generate(dim=2)

# ... and save it to disk
gmsh.write("gmsh.msh")

# Remember that by default, if physical groups are defined, Gmsh will export in
# the output mesh file only those elements that belong to at least one physical
# group. To force Gmsh to save all elements, you can use
#
# gmsh.option.setNumber("Mesh.SaveAll", 1)

# By default, Gmsh saves meshes in the latest version of the Gmsh mesh file
# format (the `MSH' format). You can save meshes in other mesh formats by
# specifying a filename with a different extension. For example
#d
#   gmsh.write("t1.unv")
#
# will save the mesh in the UNV format. You can also save the mesh in older
# versions of the MSH format: simply set
#
#   gmsh.option.setNumber("Mesh.MshFileVersion", x)
#
# for any version number `x'. As an alternative, you can also not specify the
# format explicitly, and just choose a filename with the `.msh2' or `.msh4'
# extension.

# To visualize the model we can run the graphical user interface with
# `gmsh.fltk.run()'. Here we run it only if the "-nopopup" is not provided in
# the command line arguments:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

# Note that starting with Gmsh 3.0, models can be built using other geometry
# kernels than the default "built-in" kernel. To use the OpenCASCADE CAD kernel
# instead of the built-in kernel, you should use the functions with the
# `gmsh.model.geo' prefix.
#
# Different CAD kernels have different features. With OpenCASCADE, instead of
# defining the surface by successively defining 4 points, 4 curves and 1 curve
# loop, one can define the rectangular surface directly with
#
# gmsh.model.geo.addRectangle(.2, 0, 0, .1, .3)
#
# After synchronization with the Gmsh model with
#
# gmsh.model.geo.synchronize()
#
# the underlying curves and points could be accessed with
# gmsh.model.getBoundary().
#
# See e.g. `t16.py', `t18.py', `t19.py' or `t20.py' for complete examples based
# on OpenCASCADE, and `demos/api' for more.

# This should be called when you are done using the Gmsh Python API:
gmsh.finalize()