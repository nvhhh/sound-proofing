"""Build a printable two-leaf micro-perforated absorber coupon (single material)."""
import numpy as np
import trimesh
from shapely.geometry import Polygon, Point
from shapely.ops import unary_union

SIZE = 120.0          # coupon width/depth, mm
T_PLATE = 0.8         # perforated leaf thickness, mm
HOLE_D = 0.5          # hole diameter, mm
PITCH1, PITCH2 = 3.0, 6.0
GAP1, GAP2 = 25.0, 25.0
T_BACK = 2.0          # solid back wall, mm
RIB_T = 1.2           # rib wall thickness, mm
RIB_PITCH = 20.0      # compartment size, mm


def holed_plate(z0, pitch, circle_res=8):
    base = Polygon([(0, 0), (SIZE, 0), (SIZE, SIZE), (0, SIZE)])
    n = int(SIZE // pitch)
    off = (SIZE - (n - 1) * pitch) / 2
    holes = [Point(off + i * pitch, off + j * pitch).buffer(HOLE_D / 2, quad_segs=circle_res)
             for i in range(n) for j in range(n)]
    poly = base.difference(unary_union(holes))
    mesh = trimesh.creation.extrude_polygon(poly, T_PLATE)
    mesh.apply_translation([0, 0, z0])
    return mesh


def rib_layer(z0, height):
    parts = []
    n = int(SIZE // RIB_PITCH)
    for i in range(n + 1):
        x = min(i * RIB_PITCH, SIZE - RIB_T)
        parts.append(trimesh.creation.box([RIB_T, SIZE, height],
                     trimesh.transformations.translation_matrix([x + RIB_T / 2, SIZE / 2, z0 + height / 2])))
        parts.append(trimesh.creation.box([SIZE, RIB_T, height],
                     trimesh.transformations.translation_matrix([SIZE / 2, x + RIB_T / 2, z0 + height / 2])))
    return parts


z = 0.0
parts = [holed_plate(z, PITCH1)]
z += T_PLATE
parts += rib_layer(z, GAP1)
z += GAP1
parts.append(holed_plate(z, PITCH2))
z += T_PLATE
parts += rib_layer(z, GAP2)
z += GAP2
back = trimesh.creation.box([SIZE, SIZE, T_BACK],
                            trimesh.transformations.translation_matrix([SIZE / 2, SIZE / 2, z + T_BACK / 2]))
parts.append(back)

mesh = trimesh.util.concatenate(parts)
mesh.merge_vertices()
mesh.fix_normals()
print('watertight (before boolean):', mesh.is_watertight, 'faces', len(mesh.faces))
mesh.export('mpp_absorber.stl')
print('bbox mm:', np.round(mesh.extents, 2), 'volume mm^3:', round(mesh.volume, 1))

# Fuse overlapping bodies into one watertight solid.
from trimesh.boolean import union
solid = union(parts, engine='manifold')
solid.fix_normals()
print('watertight (after union):', solid.is_watertight, 'faces', len(solid.faces))
solid.export('mpp_absorber.stl')
print('bbox mm:', np.round(solid.extents, 2), 'volume mm^3:', round(solid.volume, 1))
