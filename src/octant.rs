type Point = [f64; 3];

#[derive(PartialEq, Eq, PartialOrd, Ord, Copy, Clone, Debug, Default, Hash)]
pub(crate) struct OctantId(pub usize);

#[derive(Clone, Debug, Default)]
/// A specific node in octree
pub(crate) struct Octant {
    /// Tree related attributes. For root octant, the parent is None.
    pub parent: Option<OctantId>,
    pub children: Vec<OctantId>,

    /// The actual data which will be stored within the tree.
    pub center: Point,
    /// The extent of octant (in radius).
    pub extent: f64,
    /// Child point indices in point cloud.
    pub ipoints: Vec<usize>,
    /// The ranking within sibling octants.
    pub ranking: usize,
}

impl Octant {
    /// Construct an Octant cube with extent (in radius).
    ///
    /// # Panic
    ///
    /// * Panics if extent is negative.
    ///
    pub fn new(extent: f64) -> Self {
        assert!(extent.is_sign_positive());
        Octant {
            extent: extent,
            ..Default::default()
        }
    }

    /// Construct a root octant from 3D points
    pub fn from_points(points: &[Point]) -> Self {
        use vecfx::*;

        // define the boundary in XYZ directions
        let xs: Vec<_> = points.iter().map(|[x, _, _]| *x).collect();
        let ys: Vec<_> = points.iter().map(|[_, y, _]| *y).collect();
        let zs: Vec<_> = points.iter().map(|[_, _, z]| *z).collect();

        let (xmin, ymin, zmin) = (xs.min(), ys.min(), zs.min());
        let (xmax, ymax, zmax) = (xs.max(), ys.max(), zs.max());

        let (wx, wy, wz) = (xmax - xmin, ymax - ymin, zmax - zmin);
        let width = [wx, wy, wz].max();

        // Construct the root octant containg all points
        let mut o = Octant::new(0.5 * width);
        o.center = [(xmax + xmin) / 2., (ymax + ymin) / 2., (zmax + zmin) / 2.];

        o.ipoints = (0..points.len()).collect();
        o
    }

    /// test if two octants are neighboring
    pub fn neighboring(&self, other: &Octant) -> bool {
        let e = other.extent + self.extent;

        for i in 0..3 {
            let v = (other.center[i] - self.center[i]).abs() - e;
            if v > 0.001 {
                return false;
            }
        }

        true
    }
}

#[test]
fn test_octant() {
    // test from points
    // coordinates from CH4 molecule
    let points = vec![
        [1.21333641, 1.25660414, 1.37365150],
        [1.56999084, 0.24779414, 1.37365150],
        [1.57000925, 1.76100233, 2.24730300],
        [1.57000925, 1.76100233, 0.50000000],
        [0.14333641, 1.25661732, 1.37365150],
    ];

    let octant = Octant::from_points(&points);
    assert_eq!(octant.ipoints.len(), 5);
    assert_relative_eq!(octant.extent, (2.24730300 - 0.5) / 2.0, epsilon = 1e-4);
    let x = (1.57000925 + 0.14333641) / 2.0;
    let y = (1.76100233 + 0.24779414) / 2.0;
    let z = (2.24730300 + 0.50000000) / 2.0;
    assert_relative_eq!(octant.center[0], x, epsilon = 1e-4);
    assert_relative_eq!(octant.center[1], y, epsilon = 1e-4);
    assert_relative_eq!(octant.center[2], z, epsilon = 1e-4);

    // test neighboring
    let mut octant1 = Octant::new(45.051);
    octant1.center = [44.798034, 18.452375500000002, 121.71243299999999];
    let mut octant2 = Octant::new(22.525501);
    octant2.center = [67.32353549999999, -4.073125999999995, 144.2379345];
    assert!(octant1.neighboring(&octant2));

    let mut octant3 = Octant::new(22.525501);
    octant3.center = [22.272532500000004, -4.073125999999995, 144.2379345];
    assert!(octant2.neighboring(&octant3));

    let mut octant4 = Octant::new(2.8156876874999996);
    octant4.center = [13.825469437500008, -1.2574383124999953, 113.26536993749998];
    assert!(!octant3.neighboring(&octant4));
}
