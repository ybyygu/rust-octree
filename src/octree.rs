use std::collections::HashMap;
use std::ops::{Index, IndexMut};

use super::get_positions_from_xyz_stream;
use super::get_positions_from_xyzfile;

use crate::types::*;

#[derive(Clone, Debug)]
pub struct Octree<'a> {
    /// adjustable parameter for min number points of octant while building octree
    pub bucket_size: usize,
    /// adjustable paramter for min octant extent while building octree
    pub min_extent: f64,

    /// reference points in 3D space
    pub points: &'a Points,
    /// private data storing all octants in octree
    pub octants: Vec<Octant>,
    /// root octant index to Octree.octans
    pub root: OctantId,

    /// for quick access octant containing certain point
    mapping_octants: HashMap<usize, usize>,
}

impl<'a> Octree<'a> {
    /// Construct octree from points in 3D space
    pub fn new(points: &'a Points) -> Self {
        let octant = Octant::from_points(&points);
        let octants = vec![octant];
        let root = OctantId(0);

        Octree {
            points: points,
            octants: octants,
            root: root,

            bucket_size: 8,
            min_extent: 2.0,
            mapping_octants: HashMap::new(),
        }
    }

    pub fn root(&self) -> OctantId {
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

impl<'a> Index<OctantId> for Octree<'a> {
    type Output = Octant;

    fn index(&self, node: OctantId) -> &Octant {
        &self.octants[node.0]
    }
}

impl<'a> IndexMut<OctantId> for Octree<'a> {
    fn index_mut(&mut self, node: OctantId) -> &mut Octant {
        &mut self.octants[node.0]
    }
}

impl<'a> Octree<'a> {
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
fn octree_create_child_octants(octant: &Octant, points: &Points) -> Vec<Octant> {
    let extent = octant.extent as f64 / 2f64;

    let mut octants = vec![];

    // initialize 8 child octants
    // 1. update center
    for i in 0..8 {
        let mut o = Octant::new(extent);
        let factors = get_octant_cell_factor(i);
        // j = 0, 1, 2 => x, y, z
        for j in 0..3 {
            o.center[j] += extent * factors[j] + octant.center[j]
        }
        octants.push(o);
    }

    // 2. update point indices
    if octant.ipoints.len() > 1 {
        let (x0, y0, z0) = (octant.center[0], octant.center[1], octant.center[2]);
        // 1. scan xyz
        for &i in octant.ipoints.iter() {
            let p = points[i];
            let (x, y, z) = (p[0] - x0, p[1] - y0, p[2] - z0);
            let index = get_octant_cell_index(x, y, z);
            octants[index].ipoints.push(i);
        }
    }

    octants
}

// zyx: +++ => 0
// zyx: ++- => 1
// zyx: --- => 7
// morton encode
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
fn get_octant_cell_factor(index: usize) -> Point {
    debug_assert!(index < 8 && index >= 0);
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

impl<'a> Octree<'a> {
    /// build octree by recursively creating all octants
    pub fn build(&mut self) {
        // calculate max allowed depth according min octant extent
        let max_extent = self.octants[0].extent as f64;
        let min_extent = self.min_extent as f64;
        let max_depth = ((max_extent / min_extent).ln() / 2f64.ln()).floor() as usize;

        let root = self.root();
        let npoints = self.points.len();

        if npoints > self.bucket_size {
            let mut depth = 0;
            let mut need_split = vec![root];
            loop {
                // 1. split into child octants
                let mut remained = vec![];
                for &parent_node in need_split.iter() {
                    self.split_octant(parent_node);
                    for &child_node in &self[parent_node].children {
                        let octant = &self[child_node];
                        let n = octant.ipoints.len();
                        if n > self.bucket_size {
                            remained.push(child_node);
                        }
                    }
                }

                // 2. drill down to process child octants
                need_split.clear();
                need_split.extend(remained);

                // 3. loop control
                if need_split.is_empty() {
                    println!("octree built after {:?} cycles.", depth);
                    break;
                }
                depth += 1;
                if depth >= max_depth {
                    eprintln!("octree build: max allowed depth {} reached.", depth);
                    break;
                }
            }
        }

        // cache octants
        // create mapping of point => octant
        for (i, ref octant) in self.octants.iter().enumerate() {
            for &j in octant.ipoints.iter() {
                self.mapping_octants.insert(j, i);
            }
        }
    }
}

#[test]
fn test_octree_struct() {
    let points = get_positions_from_xyz_stream(&XYZ_TXT).unwrap();
    let mut octree = Octree::new(&points);

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
    let points = get_positions_from_xyz_stream(&XYZ_TXT).unwrap();
    let mut octree = Octree::new(&points);
    let root = octree.root();
    octree.split_octant(root);

    // reborrow as immutable
    let octree = octree;
    // root octant
    let octant = &octree[root];
    println!("{:?}", octree);

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
    println!("{:?}", octant);
    println!("{:?}", child7);

    assert!(child7.ipoints.contains(&2));
    assert_eq!(child7.parent, Some(root));
}

impl<'a> Octree<'a> {
    /// Search nearby points within radius of center.
    /// Return
    /// ------
    /// indices of nearby points and distances
    pub fn search(&self, p: Point, radius: f64) -> Vec<(usize, f64)> {
        let mut query = Query::new(radius);
        query.center = p;

        let mut pts_maybe: Vec<usize> = vec![];

        // step 1: record all nearby points by octree search
        let mut nodes_to_visit = vec![self.root()];
        loop {
            let mut todo = vec![];
            for &parent in nodes_to_visit.iter() {
                let octant = &self[parent];
                match query.relation(&octant) {
                    QORelation::Overlaps | QORelation::Within => {
                        // println!("overlaps");
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

        let mut neighbors = vec![];
        for &i in pts_maybe.iter() {
            let (px, py, pz) = (self.points[i][0], self.points[i][1], self.points[i][2]);
            let dsqr = (px - qx) * (px - qx) + (py - qy) * (py - qy) + (pz - qz) * (pz - qz);
            if dsqr < rsqr {
                neighbors.push((i, dsqr.sqrt()));
            }
        }

        neighbors
    }
}

impl<'a> Octree<'a> {
    /// Find neighboring points
    ///
    /// Parameters
    /// ----------
    /// radius: the query ball radius
    ///
    /// Return
    /// ------
    /// indices of neighboring points in pairs
    ///
    pub fn neighbors(&self, radius: f64) -> Vec<(usize, usize, f64)> {
        let mut pairs = vec![];

        for (i, &p) in self.points.iter().enumerate() {
            let neighbors = self.search(p, radius);
            for &(j, d) in neighbors.iter() {
                if j != i {
                    pairs.push((i, j, d));
                }
            }
        }

        pairs
    }
}
