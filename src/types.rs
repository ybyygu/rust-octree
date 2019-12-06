// set up aliases for convenience
pub type Point = [f64; 3];
pub type Points = Vec<Point>;

#[derive(PartialEq, Eq, PartialOrd, Ord, Copy, Clone, Debug, Default, Hash)]
pub struct OctantId(pub usize);

#[derive(Clone, Debug, Default)]
/// A specific node of octree
pub struct Octant {
    /// tree related attributes
    /// for root octant, the parent is None
    pub parent: Option<OctantId>,
    pub children: Vec<OctantId>,

    /// The actual data which will be stored within the tree
    pub center: Point,
    /// The extent of octant (in radius)
    pub extent: f64,
    /// child points indicated with indices to point cloud
    pub ipoints: Vec<usize>,
    /// the ranking in sibling octant
    pub ranking: usize,
}

impl Octant {
    /// Construct an octant with defined extent (in radius)
    /// Panic if extent is negative
    pub fn new(extent: f64) -> Self {
        assert!(extent.is_sign_positive());
        Octant {
            extent: extent,
            ..Default::default()
        }
    }

    /// Construct root octant from a set of 3D points
    pub fn from_points(points: &[Point]) -> Self {
        // define the boundary in XYZ directions
        let mut p_min = points[0];
        let mut p_max = points[0];
        for p in points {
            for i in 0..3 {
                if p[i] > p_max[i] {
                    p_max[i] = p[i];
                } else if p[i] < p_min[i] {
                    p_min[i] = p[i];
                }
            }
        }

        // Construct the root octant containg all points
        let mut distance = 0.0f64;
        for i in 0..3 {
            distance = distance.max(p_max[i] - p_min[i]);
        }
        let mut octant = Octant::new(0.5 * distance);

        octant.center = [
            (p_max[0] + p_min[0]) / 2.,
            (p_max[1] + p_min[1]) / 2.,
            (p_max[2] + p_min[2]) / 2.,
        ];
        octant.ipoints = (0..points.len()).collect();

        octant
    }
}

impl Octant {
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
fn test_octree_octant_neighboring() {
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

#[derive(Debug)]
pub struct Query {
    pub center: Point,
    pub radius: f64,
}

/// Four possible relations of a query ball with an octant
#[derive(Debug, PartialEq)]
pub enum QORelation {
    /// the query ball has no common space with the octant
    Disjoint,
    /// the query ball is partially overlapping with the octant
    Overlaps,
    /// the query ball completely contains the octant
    Contains,
    /// the query ball is completely within the octant
    Within,
}

impl Query {
    pub fn new(r: f64) -> Self {
        assert!(r.is_sign_positive(), "radius has to be positive: {}", r);
        Query {
            center: [0.0; 3],
            radius: r,
        }
    }

    /// calculate the relation of query ball with octant
    pub fn relation(&self, octant: &Octant) -> QORelation {
        let extent = octant.extent;
        let radius = self.radius;

        let x = (self.center[0] - octant.center[0]).abs();
        let y = (self.center[1] - octant.center[1]).abs();
        let z = (self.center[2] - octant.center[2]).abs();

        // 1. cheap case: xyz > e+r
        let max_dist = extent + radius;
        if (x > max_dist || y > max_dist || z > max_dist) {
            return QORelation::Disjoint;
        }

        // 2. overlaps or not
        if (x < extent || y < extent || z < extent) {
            // expected to be common: e >= r
            // expected to be rare  : e < r
            if extent >= radius {
                // 2.1 Within
                // cheap case: xyz < e-r < e+r
                let min_dist = extent - radius;
                if (x <= min_dist && y <= min_dist && z <= min_dist) {
                    return QORelation::Within;
                }
            } else {
                if (x <= extent && y <= extent && z <= extent) {
                    // distance to the farthest corner point
                    let r_sqr = radius * radius;
                    let e = extent;
                    let d_sqr = (x + e) * (x + e) + (y + e) * (y + e) + (z + e) * (z + e);
                    // 2.2 Contains
                    if d_sqr <= r_sqr {
                        return QORelation::Contains;
                    }
                }
            }
            // cheap case: e < x < e+r || e < y < e+r || z < e < e+r
            return QORelation::Overlaps;
        }

        // 3. corner case: Disjoint or Overlaps?
        // FIXME: can we just assume "Overlaps" to improve efficiency?
        // expensive case: e < xyz < e+r
        // distance to the nearest corner point
        let r_sqr = radius * radius;
        let e = extent;
        let d_sqr = (x - e) * (x - e) + (y - e) * (y - e) + (z - e) * (z - e);
        if d_sqr > r_sqr {
            return QORelation::Disjoint;
        }

        QORelation::Overlaps
    }
}

#[test]
fn test_octree_query_relations() {
    let octant = Octant::new(2.5);
    let mut query = Query::new(1.4);
    let r = query.relation(&octant);
    assert_eq!(r, QORelation::Within);

    query.radius = 4.4; // 2.5*sqrt(3)
    let r = query.relation(&octant);
    assert_eq!(r, QORelation::Contains);

    let octant = Octant::new(2.5);
    let mut query = Query::new(0.4);
    query.center = [2.7, 2.7, 2.7];
    let r = query.relation(&octant);
    assert_eq!(r, QORelation::Overlaps);

    query.center = [2.7, -2.7, -2.7];
    let r = query.relation(&octant);
    assert_eq!(r, QORelation::Overlaps);

    query.center = [2.8, 2.8, 2.8];
    let r = query.relation(&octant);
    assert_eq!(r, QORelation::Disjoint);

    let query = Query {
        center: [31.079695, 10.200508, 146.169464],
        radius: 3.0,
    };

    let mut octant = Octant::new(22.525501499999997);
    octant.center = [67.32353549999999, 40.977877, 144.2379345];
    let r = query.relation(&octant);
    assert_eq!(r, QORelation::Disjoint);
}

pub const XYZ_TXT: &str = " N                  0.49180679   -7.01280337   -3.37298245
 H                  1.49136679   -7.04246937   -3.37298245
 C                 -0.19514721   -5.73699137   -3.37298245
 H                 -0.81998021   -5.66018837   -4.26280545
 C                 -1.08177021   -5.59086937   -2.14084145
 C                  0.79533179   -4.58138037   -3.37298245
 H                 -0.46899721   -5.65651737   -1.24178645
 H                 -1.58492621   -4.62430837   -2.16719845
 H                 -1.82600521   -6.38719137   -2.13160945
 O                  2.03225779   -4.81286537   -3.37298245
 H                  0.43991988   -3.57213195   -3.37298245
 H                 -0.03366507   -7.86361434   -3.37298245 ";

#[test]
fn test_octree_init() {
    use super::*;

    let points = get_positions_from_xyz_stream(&XYZ_TXT).unwrap();
    let octant = Octant::from_points(&points);
    assert_relative_eq!(octant.center[0], 0.103126, epsilon=1e-4);
    assert_relative_eq!(octant.center[1], -5.717873145, epsilon=1e-4);
    assert_relative_eq!(octant.center[2], -2.75229595, epsilon=1e-4);
    assert_relative_eq!(octant.extent, 2.145741195, epsilon=1e-4);
}
