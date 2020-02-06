use std::ops::{Index, IndexMut};

use log::*;

use crate::octant::*;
use crate::query::*;

type Point = [f64; 3];

#[derive(Clone, Debug)]
pub struct Octree {
    /// reference points in 3D space
    pub points: Vec<Point>,

    /// private data storing all octants in octree
    octants: Vec<Octant>,
    /// root octant index to Octree.octans
    root: OctantId,
    // /// for quick access octant containing certain point
    // mapping_octants: HashMap<usize, usize>,
}

impl Octree {
    /// Construct octree from points in 3D space
    pub fn new(points: impl IntoIterator<Item = Point>) -> Self {
        let points: Vec<_> = points.into_iter().collect();
        let octant = Octant::from_points(&points);
        let octants = vec![octant];
        let root = OctantId(0);

        Octree {
            points,
            octants: octants,
            root: root,
            // mapping_octants: HashMap::new(),
        }
    }

    fn root(&self) -> OctantId {
        OctantId(0)
    }

    /// Count octants in octree.
    pub fn count(&self) -> usize {
        self.octants.len()
    }

    /// Returns true if octree has no octant, false otherwise
    pub fn is_empty(&self) -> bool {
        self.octants.is_empty()
    }
}

impl Index<OctantId> for Octree {
    type Output = Octant;

    fn index(&self, node: OctantId) -> &Octant {
        &self.octants[node.0]
    }
}

impl IndexMut<OctantId> for Octree {
    fn index_mut(&mut self, node: OctantId) -> &mut Octant {
        &mut self.octants[node.0]
    }
}

impl Octree {
    /// Add octant as orphan node in tree, return OctantId for further operation
    fn new_node(&mut self, octant: Octant) -> OctantId {
        let next_index = self.octants.len();
        self.octants.push(octant);

        OctantId(next_index)
    }

    /// Append a new child octant to parent node
    fn append_child(&mut self, parent_node: OctantId, mut octant: Octant) -> OctantId {
        octant.parent = Some(parent_node);
        octant.ranking = self[parent_node].children.len();
        let n = self.new_node(octant);

        // get parent octant, update children attributes
        let parent_octant = &mut self[parent_node];
        parent_octant.children.push(n);

        n
    }

    /// split parent octant into 8 child octants, and append them into octree
    fn split_octant(&mut self, parent_node: OctantId) {
        let child_octants = octree_create_child_octants(&self[parent_node], &self.points);

        for octant in child_octants {
            self.append_child(parent_node, octant);
        }
    }
}

/// octant: octree node data
/// points: reference points in 3D space
fn octree_create_child_octants(octant: &Octant, points: &[Point]) -> Vec<Octant> {
    let extent = octant.extent as f64 / 2f64;

    // initialize 8 child octants
    // 1. update center
    // Note: rayon's par_iter is slow
    let mut child_octants: Vec<_> = (0..8usize)
        .into_iter()
        .map(|i| {
            let mut o = Octant::new(extent);
            let factors = get_octant_cell_factor(i);
            // j = 0, 1, 2 => x, y, z
            for j in 0..3 {
                o.center[j] += extent * factors[j] + octant.center[j]
            }
            o
        })
        .collect();

    // 2. update point indices
    if octant.ipoints.len() > 1 {
        let (x0, y0, z0) = (octant.center[0], octant.center[1], octant.center[2]);
        // 1. scan xyz
        for &i in octant.ipoints.iter() {
            let p = points[i];
            let (x, y, z) = (p[0] - x0, p[1] - y0, p[2] - z0);
            let index = get_octant_cell_index(x, y, z);
            child_octants[index].ipoints.push(i);
        }
    }

    child_octants
}

// zyx: +++ => 0
// zyx: ++- => 1
// zyx: --- => 7
// morton encode
#[inline]
fn get_octant_cell_index(x: f64, y: f64, z: f64) -> usize {
    // create lookup table, which could be faster
    match (
        z.is_sign_positive(),
        y.is_sign_positive(),
        x.is_sign_positive(),
    ) {
        (true, true, true) => 0,
        (true, true, false) => 1,
        (true, false, true) => 2,
        (true, false, false) => 3,
        (false, true, true) => 4,
        (false, true, false) => 5,
        (false, false, true) => 6,
        (false, false, false) => 7,
    }

    // another way: using bit shift
    // let bits = [z.is_sign_negative(), y.is_sign_negative(), x.is_sign_negative()];
    // bits.iter().fold(0, |acc, &b| acc*2 + b as usize)
}

#[test]
fn test_octree_cell_index() {
    let index = get_octant_cell_index(1.0, 1.0, 1.0);
    assert_eq!(index, 0);

    let index = get_octant_cell_index(-1.0, -1.0, -1.0);
    assert_eq!(index, 7);

    let index = get_octant_cell_index(-1.0, 1.0, 1.0);
    assert_eq!(index, 1);

    let index = get_octant_cell_index(-1.0, -1.0, 1.0);
    assert_eq!(index, 3);
}

// useful for calculate center of child octant
// morton decode
#[inline]
fn get_octant_cell_factor(index: usize) -> Point {
    debug_assert!(index < 8);
    [
        match (index & 0b001) == 0 {
            true => 1.0,
            false => -1.0,
        },
        match ((index & 0b010) >> 1) == 0 {
            true => 1.0,
            false => -1.0,
        },
        match ((index & 0b100) >> 2) == 0 {
            true => 1.0,
            false => -1.0,
        },
    ]
}

#[test]
fn test_octree_factor() {
    let x = get_octant_cell_factor(0);
    assert_eq!(1.0, x[0]);
    assert_eq!(1.0, x[1]);
    assert_eq!(1.0, x[2]);

    let x = get_octant_cell_factor(7);
    assert_eq!(-1.0, x[0]);
    assert_eq!(-1.0, x[1]);
    assert_eq!(-1.0, x[2]);

    let x = get_octant_cell_factor(2);
    assert_eq!(1.0, x[0]);
    assert_eq!(-1.0, x[1]);
    assert_eq!(1.0, x[2]);
}

impl Octree {
    /// Build octree by recursively dividing child octants
    ///
    /// * Parameters
    ///
    /// - bucket_size: the max number of points each octant holds before
    /// stopping recursively dividing.
    pub fn build(&mut self, bucket_size: usize) {
        debug_assert!(bucket_size > 0, "invalid bucket_size: {}!", bucket_size);

        let root = self.root();
        let npoints = self.points.len();
        if npoints > bucket_size {
            let mut need_split = vec![root];
            for depth in 0.. {
                // 1. split into child octants
                let mut remained = vec![];
                for &parent_node in need_split.iter() {
                    self.split_octant(parent_node);
                    for &child_node in &self[parent_node].children {
                        let octant = &self[child_node];
                        let n = octant.ipoints.len();
                        if n > bucket_size {
                            remained.push(child_node);
                        }
                    }
                }

                // 2. drill down to process child octants
                need_split.clear();
                need_split.extend(remained);

                // 3. loop control
                if need_split.is_empty() {
                    debug!("octree built after {:?} cycles.", depth);
                    break;
                }
            }
        }

        // cache octants
        // create mapping of point => octant
        // for (i, ref octant) in self.octants.iter().enumerate() {
        //     for &j in octant.ipoints.iter() {
        //         self.mapping_octants.insert(j, i);
        //     }
        // }
    }
}

// FIXME: useful or not?
/// calculate max allowed depth according min octant extent
fn max_depth(max_extent: f64, min_extent: f64) -> usize {
    assert!(
        min_extent < max_extent,
        "invalid parameters: {} {}",
        max_extent,
        min_extent
    );
    assert!(min_extent.is_sign_positive());
    let max_depth = ((max_extent / min_extent).ln() / 2f64.ln()).floor() as usize;
    max_depth
}

#[test]
fn test_octree_struct() {
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
            assert_eq!(position.len(), 3, "{:?}", position);
            let p: Vec<f64> = position.iter().map(|x| x.parse().unwrap()).collect();
            positions.push([p[0], p[1], p[2]]);
        }
    
        positions
    }

    let points = read_points(&XYZ_TXT);
    let mut octree = Octree::new(points);

    // test octree
    let root = octree.root();
    let octant = Octant::new(1.2);
    let child1 = octree.append_child(root, octant);
    assert_eq!(&octree[child1].parent, &Some(root));

    let octant = Octant::new(1.2);
    let child2 = octree.append_child(root, octant);
    let octant = Octant::new(1.3);
    let child3 = octree.append_child(child1, octant);

    let root_octant = &octree[root];
    assert!(root_octant.children.contains(&child1));
    assert!(root_octant.children.contains(&child2));
    assert_eq!(octree.count(), 4);

    let octant1 = &octree[child1];
    assert_eq!(octant1.extent, 1.2);
    assert!(octant1.children.contains(&child3));
}

#[test]
fn test_octree_split_children() {
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
            assert_eq!(position.len(), 3, "{:?}", position);
            let p: Vec<f64> = position.iter().map(|x| x.parse().unwrap()).collect();
            positions.push([p[0], p[1], p[2]]);
        }
    
        positions
    }

    let points = read_points(&XYZ_TXT);
    let mut octree = Octree::new(points);
    let root = octree.root();
    octree.split_octant(root);

    // reborrow as immutable
    let octree = octree;
    // root octant
    let octant = &octree[root];

    let children = &octant.children;
    let child = &octree[children[0]];
    let x = child.center[0] - octant.center[0];
    let y = child.center[1] - octant.center[1];
    let z = child.center[2] - octant.center[2];
    assert_relative_eq!(x, y, epsilon = 1e-4);
    assert_relative_eq!(x, z, epsilon = 1e-4);

    assert_eq!(octant.extent, child.extent * 2.0);

    let child = &octree[children[1]];
    let x = child.center[0] - octant.center[0];
    let y = child.center[1] - octant.center[1];
    let z = child.center[2] - octant.center[2];
    assert_relative_eq!(x, child.extent * -1., epsilon = 1e-4);
    assert_relative_eq!(y, child.extent * 1., epsilon = 1e-4);
    assert_relative_eq!(z, child.extent * 1., epsilon = 1e-4);

    let child7 = &octree[children[7]];
    assert!(child7.ipoints.contains(&2));
    assert_eq!(child7.parent, Some(root));
}

impl Octree {
    /// Search nearby points within radius of center.
    ///
    /// Parameters
    /// ----------
    /// - p: the coordinates of the point searching for neighbors.
    /// - radius: the cutoff radius for neighbors.
    ///
    /// Return
    /// ------
    /// indices of nearby points and distances
    pub fn search(&self, p: Point, radius: f64) -> impl Iterator<Item = (usize, f64)> + '_ {
        let mut query = Query::new(radius);
        query.center = p;

        // step 1: record all nearby points by octree search
        let mut pts_maybe: Vec<usize> = vec![];
        let mut nodes_to_visit = vec![self.root()];
        loop {
            let mut todo = vec![];
            for &parent in nodes_to_visit.iter() {
                let octant = &self[parent];
                match query.relation(&octant) {
                    QORelation::Overlaps | QORelation::Within => {
                        // debug!("overlaps");
                        if octant.children.is_empty() {
                            // is a leaf node: save points
                            pts_maybe.extend(octant.ipoints.iter());
                        } else {
                            // not a leaf node: go down to follow children
                            todo.extend(octant.children.iter());
                        }
                    }

                    QORelation::Contains => {
                        pts_maybe.extend(octant.ipoints.iter());
                    }

                    QORelation::Disjoint => {}
                };
            }

            if todo.is_empty() {
                break;
            }

            nodes_to_visit.clear();
            nodes_to_visit.extend(todo.iter());
        }

        // step 2: linear search
        let (qx, qy, qz) = (query.center[0], query.center[1], query.center[2]);
        let radius = query.radius as f64;
        let rsqr = radius * radius;
        pts_maybe.into_iter().filter_map(move |i| {
            let (px, py, pz) = (self.points[i][0], self.points[i][1], self.points[i][2]);
            let dsqr = (px - qx) * (px - qx) + (py - qy) * (py - qy) + (pz - qz) * (pz - qz);
            if dsqr < rsqr {
                Some((i, dsqr.sqrt()))
            } else {
                None
            }
        })
    }
}
