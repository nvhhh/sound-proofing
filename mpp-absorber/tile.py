"""Parametric, stackable micro-perforated absorber tile.

One tile is a square cell of the absorber. Tiles butt edge to edge in a grid;
3 mm dowels in the side sockets keep them aligned. Each tile mounts by four
countersunk M3 screws at the corner bosses, or by tape on the flat back.

Usage:
    python3 tile.py                      # default 6 x 6 cells, screws and tape
    python3 tile.py --cells 11 --mount tape
    python3 tile.py --bed 220            # largest tile that fits a 220 mm bed

All dimensions in mm. The hole grid sits between the rib lines and the cell size
is a multiple of both hole pitches, so no hole is blocked by a rib.
"""
import argparse
from dataclasses import dataclass, replace

import numpy as np
import trimesh
from manifold3d import Manifold


@dataclass(frozen=True)
class Params:
    cells: int = 6               # rib cells per side; tile width = cells * cell
    cell: float = 18.0           # rib cell size; must be a multiple of both pitches
    pitch_front: float = 3.0
    pitch_rear: float = 6.0
    hole_d: float = 0.5
    t_leaf: float = 0.8          # perforated leaf thickness
    gap_front: float = 25.0      # cavity behind the front leaf
    gap_rear: float = 25.0       # cavity behind the rear leaf
    rib_t: float = 1.2
    back: float = 5.0            # solid back; thick enough for dowel sockets
    boss: float = 8.0            # corner boss size (screw and tile corner)
    screw_d: float = 3.4         # M3 clearance
    csk_d: float = 6.5           # M3 countersink head diameter
    dowel_d: float = 3.2         # socket for a 3 mm dowel
    dowel_depth: float = 6.0
    mount: str = "both"          # "screw", "tape" or "both"

    @property
    def width(self):
        return self.cells * self.cell

    @property
    def depth(self):
        return 2 * self.t_leaf + self.gap_front + self.gap_rear + self.back


def check(p):
    for pitch in (p.pitch_front, p.pitch_rear):
        if p.cell % pitch:
            raise ValueError(f"cell {p.cell} must be a multiple of pitch {pitch}")
    if p.mount not in ("screw", "tape", "both"):
        raise ValueError("mount must be screw, tape or both")
    if p.boss + p.rib_t > p.cell:
        raise ValueError("boss overlaps the first rib cell")


def box(size, origin):
    return Manifold.cube(size, center=False).translate(origin)


def cyl(r, h, origin, axis="z"):
    m = Manifold.cylinder(h, r, r, 24, center=False)
    if axis == "x":
        m = m.rotate([0, 90, 0])
    elif axis == "y":
        m = m.rotate([90, 0, 0])
    return m.translate(origin)


def hole_centres(width, pitch, boss, cell, rib_t, hole_r):
    """Hole centres on the half-pitch grid, minus boss corners and any blocked by a rib."""
    n = int(round(width / pitch))
    c = (np.arange(n) + 0.5) * pitch
    ribs = np.arange(0, width, cell)
    # a hole is blocked if it touches a rib (rib centreline +/- rib_t/2)
    nearest = np.min(np.abs(c[:, None] - ribs[None, :]), axis=1)
    if np.any(nearest < hole_r + rib_t / 2):
        raise ValueError("a hole would be blocked by a rib")
    pts = [(x, y) for x in c for y in c
           if not ((x < boss or x > width - boss) and (y < boss or y > width - boss))]
    return pts


def leaf(p, pitch, z0, thickness):
    """Solid plate with the perforation, returned as a Manifold."""
    w = p.width
    plate = box([w, w, thickness], [0, 0, z0])
    pts = hole_centres(w, pitch, p.boss, p.cell, p.rib_t, p.hole_d / 2)
    holes = Manifold.compose([cyl(p.hole_d / 2, thickness + 2, [x, y, z0 - 1]) for x, y in pts])
    return plate - holes


def rib_walls(p, z0, height):
    w = p.width
    parts = []
    for x in np.arange(0, w, p.cell):
        x = min(x, w - p.rib_t)
        parts.append(box([p.rib_t, w, height], [x, 0, z0]))
        parts.append(box([w, p.rib_t, height], [0, x, z0]))
    return parts


def build(p):
    check(p)
    w, d = p.width, p.depth
    z = 0.0
    parts = []
    parts.append(leaf(p, p.pitch_front, z, p.t_leaf))
    z += p.t_leaf
    parts += rib_walls(p, z, p.gap_front)
    z += p.gap_front
    parts.append(leaf(p, p.pitch_rear, z, p.t_leaf))
    z += p.t_leaf
    parts += rib_walls(p, z, p.gap_rear)
    z += p.gap_rear
    back_z = z
    parts.append(box([w, w, p.back], [0, 0, back_z]))

    # Corner bosses run the full depth: they carry the screws and the tile corner.
    for x in (0, w - p.boss):
        for y in (0, w - p.boss):
            parts.append(box([p.boss, p.boss, d], [x, y, 0]))

    tile = parts[0]
    for part in parts[1:]:
        tile = tile + part

    cutters = []
    # Dowel sockets along the two side faces, inside the solid back.
    zc = back_z + p.back / 2
    for y in (w / 4, 3 * w / 4):
        cutters.append(cyl(p.dowel_d / 2, p.dowel_depth, [0, y, zc], axis="x"))
        cutters.append(cyl(p.dowel_d / 2, p.dowel_depth, [w, y, zc], axis="x").translate([-p.dowel_depth, 0, 0]))
    # Countersunk screw holes from the back face, one per corner boss.
    if p.mount in ("screw", "both"):
        m = p.boss / 2
        for x in (m, w - m):
            for y in (m, w - m):
                cutters.append(cyl(p.screw_d / 2, d + 2, [x, y, -1]))
                csk_h = (p.csk_d - p.screw_d) / 2           # 90 degree countersink
                cutters.append(Manifold.cylinder(csk_h, p.screw_d / 2, p.csk_d / 2, 24, center=False)
                               .translate([x, y, d - csk_h]))
                cutters.append(cyl(p.csk_d / 2, 3, [x, y, d]))
    tile = tile - Manifold.compose(cutters)
    return tile


def to_trimesh(m):
    mesh = m.to_mesh()
    return trimesh.Trimesh(vertices=np.asarray(mesh.vert_properties)[:, :3],
                           faces=np.asarray(mesh.tri_verts), process=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=int)
    ap.add_argument("--bed", type=float, help="printer bed side, mm; picks the largest tile that fits")
    ap.add_argument("--mount", choices=["screw", "tape", "both"])
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    p = Params()
    if a.bed:
        p = replace(p, cells=int((a.bed - 10) // p.cell))
    if a.cells:
        p = replace(p, cells=a.cells)
    if a.mount:
        p = replace(p, mount=a.mount)
    mesh = to_trimesh(build(p))
    out = a.out or f"mpp_tile_{int(p.width)}mm_{p.mount}.stl"
    mesh.export(out)
    print(f"cells {p.cells}  width {p.width:.0f} mm  depth {p.depth:.1f} mm  mount {p.mount}")
    print("bbox mm", np.round(mesh.extents, 2), " watertight", mesh.is_watertight,
          " volume mm^3", round(mesh.volume), " wrote", out)


if __name__ == "__main__":
    main()
