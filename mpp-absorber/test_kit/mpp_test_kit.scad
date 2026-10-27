// MPP test kit: impedance-tube sample, rigid reference, hole test strip.
// Units: mm. Choose `part`, press F6 (CGAL) or F5, then export STL.
// Not yet rendered in OpenSCAD in this repo: check it once in your install.
//
// Sample: fits a 29 mm ID impedance tube (ISO 10534-2 small tube). The stack follows the
// design in geometry/build_data: front leaf, cavity 25 mm, rear leaf, cavity 30 mm, back.
// Note: the small sample has no rib walls, so it is a single open cavity per gap.
// The rib cell effect is therefore not in the sample; test the tile separately if needed.

part = "sample";          // "sample", "rigid" or "holes"

tube_d    = 29;          // inner diameter of the impedance tube
clearance = 0.2;         // sample is printed slightly smaller than the bore
hole_d    = 0.5;         // perforation diameter
pitch     = 2.5;         // hole pitch, must match the design (front leaf)
t_leaf    = 1.0;         // leaf thickness
gap_front = 25;          // cavity behind the front leaf
gap_rear  = 30;          // cavity behind the rear leaf
back      = 5;           // solid back
wall      = 2.4;         // holder wall: about 6 perimeters at a 0.4 mm nozzle, for a stiff ring

$fn = 64;

r_s = (tube_d - clearance) / 2;   // outer radius of the sample

// Hole centres on the half-pitch grid, inside radius r, keeping a hole edge inside the disc.
module holes2d(r, p, d) {
    n = ceil(r / p) + 1;
    for (i = [-n : n - 1], j = [-n : n - 1])
        let (x = (i + 0.5) * p, y = (j + 0.5) * p)
            if (sqrt(x * x + y * y) + d / 2 < r)
                translate([x, y]) circle(d = d);
}

module leaf(r, t) {
    linear_extrude(height = t)
        difference() {
            circle(r = r);
            holes2d(r, pitch, hole_d);
        }
}

// Holder ring: outer radius r_s, inner radius r_s - wall.
module ring(h) {
    difference() {
        cylinder(r = r_s, h = h);
        translate([0, 0, -1]) cylinder(r = r_s - wall, h = h + 2);
    }
}

module sample() {
    z1 = t_leaf;                    // top of cavity 1
    z2 = z1 + gap_front;            // rear leaf bottom
    z3 = z2 + t_leaf;               // top of cavity 2
    z4 = z3 + gap_rear;             // back wall bottom
    leaf(r_s, t_leaf);
    translate([0, 0, z1]) ring(gap_front);
    translate([0, 0, z2]) leaf(r_s, t_leaf);
    translate([0, 0, z3]) ring(gap_rear);
    translate([0, 0, z4]) cylinder(r = r_s, h = back);
}

// Solid reference disc for the rigid-wall check (absorption should be near zero).
module rigid() {
    cylinder(r = r_s, h = 10);
}

// Hole test strip: five diameters on the design pitch, to measure printed hole size.
module holes_strip() {
    dias = [0.40, 0.45, 0.50, 0.55, 0.60];
    difference() {
        translate([-2, -6, 0]) cube([5 * pitch + 4, 12, t_leaf]);
        for (k = [0 : 4])
            translate([k * pitch + 0.5 * pitch - 2, 0, -1])
                cylinder(d = dias[k], h = t_leaf + 2);
    }
}

if (part == "sample") sample();
if (part == "rigid") rigid();
if (part == "holes") holes_strip();
