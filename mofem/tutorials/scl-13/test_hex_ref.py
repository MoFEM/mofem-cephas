import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

# Define the vertices of the hexahedron (cube) based on the corrected labeling:
#   v8 ------ v7
#    |\        |\
#    | \       | \
#    |  v5 ------ v6
#    |  |      |  |
#   v4 -|----- v3 |
#     \ |       \ |
#      v1 ------ v2

# Define the vertices of the hexahedron (cube) based on the updated labeling
hex_vertices_corrected = np.array([
    [0, 0, 0],  # v1
    [1, 0, 0],  # v2
    [1, 1, 0],  # v3
    [0, 1, 0],  # v4
    [1, 0, 1],  # v5
    [0, 0, 1],  # v6
    [0, 1, 1],  # v7
    [1, 1, 1]   # v8
])

# Define tetrahedron nodal connectivity using 0-based indexing, aligned with the corrected labeling
tetrahedrons_cubit_corrected = [
    [0, 2, 3, 6],  # Tet 1
    [0, 2, 6, 7],  # Tet 2
    [0, 1, 2, 7],  # Tet 3
    [0, 1, 7, 4],  # Tet 4
    [0, 4, 7, 6],  # Tet 5
    [0, 4, 6, 5]   # Tet 6
]

# Function to plot a tetrahedron
def plot_tetrahedron(ax, vertices, color):
    faces = [[vertices[0], vertices[1], vertices[2]],
             [vertices[0], vertices[1], vertices[3]],
             [vertices[0], vertices[2], vertices[3]],
             [vertices[1], vertices[2], vertices[3]]]
    poly3d = Poly3DCollection(faces, facecolors=color, linewidths=1, edgecolors='r', alpha=1.)
    ax.add_collection3d(poly3d)

# Function to plot the division of the hexahedron
def plot_hex_division(tetrahedrons, title):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Colors for the tetrahedrons
    colors = ['cyan', 'magenta', 'yellow', 'blue', 'green', 'orange'][:len(tetrahedrons)]

    # Plot each tetrahedron
    for i, tet in enumerate(tetrahedrons):
        plot_tetrahedron(ax, hex_vertices_corrected[tet], colors[i])

    # Set the limits and labels for the plot
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    ax.set_zlim([0, 1])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.title(title)
    plt.show()

# Plot the correct 6-tetrahedron split using the updated nodal connectivity
plot_hex_division(tetrahedrons_cubit_corrected, "Hexahedron Split into 6 Tetrahedrons (Corrected Nodal Connectivity)")