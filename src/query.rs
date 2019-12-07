type Point = [f64; 3];

use crate::octant::*;

#[derive(Debug)]
pub(crate) struct Query {
    pub center: Point,
    pub radius: f64,
}

/// Four possible relations of a query ball with an octant
#[derive(Debug, PartialEq)]
pub(crate) enum QORelation {
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
        if x > max_dist || y > max_dist || z > max_dist {
            return QORelation::Disjoint;
        }

        // 2. overlaps or not
        if x < extent || y < extent || z < extent {
            // expected to be common: e >= r
            // expected to be rare  : e < r
            if extent >= radius {
                // 2.1 Within
                // cheap case: xyz < e-r < e+r
                let min_dist = extent - radius;
                if x <= min_dist && y <= min_dist && z <= min_dist {
                    return QORelation::Within;
                }
            } else {
                if x <= extent && y <= extent && z <= extent {
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

#[test]
fn test_octree_init() {
    const XYZ_TXT: &str = " N                  0.49180679   -7.01280337   -3.37298245
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

    fn read_points(txt: &str) -> Vec<[f64; 3]> {
        let mut positions = Vec::new();
        for line in txt.lines() {
            let attrs: Vec<_> = line.split_whitespace().collect();
            let (_symbol, position) = attrs.split_first().expect("empty line");
            assert_eq!(position.len(), 3,);
            let p: Vec<f64> = position.iter().map(|x| x.parse().unwrap()).collect();
            positions.push([p[0], p[1], p[2]]);
        }
    
        positions
    }
    let points = read_points(&XYZ_TXT);
    let octant = Octant::from_points(&points);
    assert_relative_eq!(octant.center[0], 0.103126, epsilon = 1e-4);
    assert_relative_eq!(octant.center[1], -5.717873145, epsilon = 1e-4);
    assert_relative_eq!(octant.center[2], -2.75229595, epsilon = 1e-4);
    assert_relative_eq!(octant.extent, 2.145741195, epsilon = 1e-4);
}
