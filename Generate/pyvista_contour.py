import pyvista as pv
from pyvista import examples


simulationResultFile = "../Data/1/simulation/result/case_t0001.vtu"
mesh = pv.read(simulationResultFile)

data = mesh
cpos = [-2, 5, 3]

# vol = data.threshold_percent(30, invert=1)
# vol.plot(show_edges=True, cpos=cpos)
# surf = data.extract_geometry()
# smooth = surf.smooth (n_iter = 10)
# smooth.plot()

# mesh = smooth
mesh = mesh.extract_geometry()
# contours = mesh.contour()
# exit()
# Extract the edges above a 12 degree feature angle
# edges = mesh.extract_feature_edges(200)
edges = mesh.extract_feature_edges(feature_angle=10, boundary_edges=True,
    non_manifold_edges=True, feature_edges=True, manifold_edges=True,
    inplace=False)
# edges = mesh.extract_feature_edges(boundary_edges=False,
#                            feature_edges=True,
#                            manifold_edges=False)

# Render the edge lines ontop of the original mesh
p = pv.Plotter()
# p.add_mesh(mesh.outline(), color='k')
p.add_mesh(mesh, color=True)
p.add_mesh(edges, color="red", line_width=5)
# p.add_mesh(contours, color="white", line_width=5)
# Define a camera position that will zoom to her eye
# p.camera_position = [(96.0, -197.0, 45.0), (7.0, -109.0, 22.0), (0, 0, 1)]
p.show()
