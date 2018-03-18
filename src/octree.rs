// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::7711fb40-175f-4198-bff1-71c5fe1d7bd3][7711fb40-175f-4198-bff1-71c5fe1d7bd3]]
use std::collections::HashMap;

use get_positions_from_xyz_stream;
use get_positions_from_xyzfile;
use Point;
use Points;
// 7711fb40-175f-4198-bff1-71c5fe1d7bd3 ends here

// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::d602663f-9f66-4e18-a538-e60b12985df3][d602663f-9f66-4e18-a538-e60b12985df3]]
#[derive(PartialEq, Eq, Copy, Clone, Debug, Default, Hash)]
pub struct OctantId (usize);

#[derive(Clone, Debug, Default)]
/// A node within a particular octree
pub struct Octant {
    // tree attributes
    parent: Option<OctantId>,
    children: Vec<OctantId>,

    /// The actual data which will be stored within the tree
    center: Point,
    /// The extent of cube (radius)
    extent: f64,
    /// indices of the points in a public array
    ipoints: Vec<usize>,
    /// the ranking in sibling octant
    ranking: usize,
}

impl Octant {
    fn new(extent: f64) -> Self {
        assert!(extent > 0.0, "extent has to be positive: {}", extent);
        Octant {
            extent: extent,
            ..Default::default()
        }
    }

    /// initialize octant struct from point cloud
    fn from_points(points: &Points) -> Self {
        let mut p_min = points[0];
        let mut p_max = points[0];

        for p in points {
            if p[0] > p_max[0] {
                p_max[0] = p[0];
            } else if p[0] < p_min[0] {
                p_min[0] = p[0];
            }

            if p[1] > p_max[1] {
                p_max[1] = p[1];
            } else if p[1] < p_min[1] {
                p_min[1] = p[1];
            }

            if p[2] > p_max[2] {
                p_max[2] = p[2];
            } else if p[2] < p_min[2] {
                p_min[2] = p[2];
            }
        }

        let mut pe = [p_max[0] - p_min[0], p_max[1] - p_min[1], p_max[2] - p_min[2]];
        let mut extent = 0.;
        for &v in pe.iter() {
            if v > extent {
                extent = v;
            }
        }

        let mut octant = Octant::new(extent/2.);
        octant.center = [(p_max[0] + p_min[0])/2., (p_max[1] + p_min[1])/2., (p_max[2] + p_min[2])/2.,];

        let n = points.len();
        octant.ipoints = (0..n).collect();

        octant
    }
}
// d602663f-9f66-4e18-a538-e60b12985df3 ends here

// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::68bdbfaf-0d07-40c4-a77c-5c6b43ab440e][68bdbfaf-0d07-40c4-a77c-5c6b43ab440e]]
#[derive(Debug)]
pub struct Query {
    center : Point,
    radius : f64,
}

/// Four possible relations of a query ball with an octant
#[derive(Debug, PartialEq)]
enum QORelation {
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
        assert!(r > 0.0, "radius has to be positive: {}", r);
        Query {
            center : [0.0; 3],
            radius : r,
        }
    }

    /// get the relation of query ball with octant
    fn relation(&self, octant: &Octant) -> QORelation {
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
                    let r_sqr = radius*radius;
                    let e = extent;
                    let d_sqr = (x+e)*(x+e) + (y+e)*(y+e) + (z+e)*(z+e);
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
        let r_sqr = radius*radius;
        let e = extent;
        let d_sqr = (x-e)*(x-e) + (y-e)*(y-e) + (z-e)*(z-e);
        if d_sqr > r_sqr {
            return QORelation::Disjoint;
        }

        QORelation::Overlaps
    }
}
// 68bdbfaf-0d07-40c4-a77c-5c6b43ab440e ends here

// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::15e377a2-f1f4-483a-a91b-5ddf7f335cb0][15e377a2-f1f4-483a-a91b-5ddf7f335cb0]]
use std::ops::{Index, IndexMut};

#[derive(Clone, Debug)]
pub struct Octree<'a> {
    /// adjustable parameter for min number points of octant while building octree
    pub bucket_size : usize,
    /// adjustable paramter for min octant extent while building octree
    pub min_extent  : f64,

    /// reference points in 3D space
    pub points  : &'a Points,
    /// private data storing all octants in octree
    pub octants : Vec<Octant>,
    /// root octant index to Octree.octans
    pub root    : OctantId,

    /// for quick access octant containing certain point
    mapping_octants : HashMap<usize, usize>,
}

impl<'a> Octree<'a> {
    /// initialize octree from points in 3D space
    pub fn new(points: &'a Points) -> Self {
        let octant = Octant::from_points(&points);
        let octants = vec![octant];
        let root = OctantId(0);

        Octree {
            points    : points,
            octants   : octants,
            root      : root,

            bucket_size : 8,
            min_extent  : 2.0,
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

    // /// Return sibling octants other than current node
    // pub fn siblings(&self, node: OctantId) -> Vec<&Octant> {
    //     let mut sibling_octants = vec![];
    //     let octant = &self[node];
    //     if let Some(parent) = octant.parent {
    //         let parent_octant = &self[parent];
    //         for &n in parent_octant.children.iter() {
    //             let o = &self[n];
    //             if o.ranking != octant.ranking {
    //                 sibling_octants.push(o);
    //             }
    //         }
    //     }

    //     sibling_octants
    // }
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
// 15e377a2-f1f4-483a-a91b-5ddf7f335cb0 ends here

// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::81167b8a-bac9-4a8e-a6c9-56e48dcd6e79][81167b8a-bac9-4a8e-a6c9-56e48dcd6e79]]
#[test]
fn test_octree_query_relations() {
    let octant = Octant::new(2.5);
    let mut query = Query::new(1.4);
    let r = query.relation(&octant);
    assert_eq!(r, QORelation::Within);

    query.radius = 4.4;         // 2.5*sqrt(3)
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
        radius: 3.0
    };

    let mut octant = Octant::new(22.525501499999997);
    octant.center = [67.32353549999999, 40.977877, 144.2379345];
    let r = query.relation(&octant);
    assert_eq!(r, QORelation::Disjoint);
}
// 81167b8a-bac9-4a8e-a6c9-56e48dcd6e79 ends here

// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::85a1bbdb-53b6-4dff-89e2-1ceba40b3c02][85a1bbdb-53b6-4dff-89e2-1ceba40b3c02]]
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

#[test]
fn test_octree_init() {
    let points = get_positions_from_xyz_stream(&XYZ_TXT).unwrap();
    let octant = Octant::from_points(&points);
    assert_relative_eq!(octant.center[0], 0.103126, epsilon=1e-4);
    assert_relative_eq!(octant.center[1], -5.717873145, epsilon=1e-4);
    assert_relative_eq!(octant.center[2], -2.75229595, epsilon=1e-4);
    assert_relative_eq!(octant.extent, 2.145741195, epsilon=1e-4);
}
// 85a1bbdb-53b6-4dff-89e2-1ceba40b3c02 ends here

// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::90433ce9-a63e-4f8e-b497-6cdd3bb88ca8][90433ce9-a63e-4f8e-b497-6cdd3bb88ca8]]
impl<'a> Octree<'a> {
    /// Add octant as orphon node in tree, return OctantId for further operation
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
            o.center[j] += extent*factors[j] + octant.center[j]
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
// 90433ce9-a63e-4f8e-b497-6cdd3bb88ca8 ends here

// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::89da6ad4-0055-4246-84c7-9d19194c5405][89da6ad4-0055-4246-84c7-9d19194c5405]]
// zyx: +++ => 0
// zyx: ++- => 1
// zyx: --- => 7
// morton encode
fn get_octant_cell_index(x: f64, y: f64, z: f64) -> usize {
    // create lookup table, which could be faster
    match (z.is_sign_positive(), y.is_sign_positive(), x.is_sign_positive()) {
        (true, true, true)    => 0,
        (true, true, false)   => 1,
        (true, false, true)   => 2,
        (true, false, false)  => 3,
        (false, true, true)   => 4,
        (false, true, false)  => 5,
        (false, false, true)  => 6,
        (false, false, false) => 7,
    }

    // another way: bit shift
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
        }
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
// 89da6ad4-0055-4246-84c7-9d19194c5405 ends here

// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::9db18239-7b01-48a3-aedc-7bcc082e7949][9db18239-7b01-48a3-aedc-7bcc082e7949]]
impl<'a> Octree<'a> {
    /// build octree by recursively creating all octants
    pub fn build(&mut self) {
        // calculate max allowed depth according min octant extent
        let max_extent = self.octants[0].extent as f64;
        let min_extent = self.min_extent as f64;
        let max_depth = ((max_extent/min_extent).ln()/2f64.ln()).floor() as usize;

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
// 9db18239-7b01-48a3-aedc-7bcc082e7949 ends here

// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::9317478e-996f-4323-9310-e1ca841b8832][9317478e-996f-4323-9310-e1ca841b8832]]
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
// 9317478e-996f-4323-9310-e1ca841b8832 ends here

// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::ea2c2276-5aaa-406e-9d5f-11a258f38cc0][ea2c2276-5aaa-406e-9d5f-11a258f38cc0]]
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
    assert_relative_eq!(x, y, epsilon=1e-4);
    assert_relative_eq!(x, z, epsilon=1e-4);

    assert_eq!(octant.extent, child.extent * 2.0);

    let child = &octree[children[1]];
    let x = child.center[0] - octant.center[0];
    let y = child.center[1] - octant.center[1];
    let z = child.center[2] - octant.center[2];
    assert_relative_eq!(x, child.extent * -1., epsilon=1e-4);
    assert_relative_eq!(y, child.extent * 1., epsilon=1e-4);
    assert_relative_eq!(z, child.extent * 1., epsilon=1e-4);

    let child7 = &octree[children[7]];
    println!("{:?}", octant);
    println!("{:?}", child7);

    assert!(child7.ipoints.contains(&2));
    assert_eq!(child7.parent, Some(root));
}
// ea2c2276-5aaa-406e-9d5f-11a258f38cc0 ends here

// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::bbcfff81-6ec6-4e9e-a787-8641691e6435][bbcfff81-6ec6-4e9e-a787-8641691e6435]]
impl<'a> Octree<'a> {
    /// Search nearby points within radius of center.
    /// Return
    /// ------
    /// indices of nearby points
    pub fn search(&self, p: Point, radius: f64) -> Vec<usize> {

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
                    },

                    QORelation::Contains => {
                        pts_maybe.extend(octant.ipoints.iter());
                    },

                    QORelation::Disjoint => {
                        ;
                    },
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
        let rsqr = radius*radius;

        let mut neighbors = vec![];
        for &i in pts_maybe.iter() {
            let (px, py, pz) = (self.points[i][0], self.points[i][1], self.points[i][2]);
            let dsqr = (px-qx)*(px-qx) + (py-qy)*(py-qy) + (pz-qz)*(pz-qz);
            if dsqr < rsqr {
                neighbors.push(i);
            }
        }

        neighbors
    }
}
// bbcfff81-6ec6-4e9e-a787-8641691e6435 ends here

// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::7e3b12c9-d3f8-4bfc-8ed0-46e2644660d3][7e3b12c9-d3f8-4bfc-8ed0-46e2644660d3]]
impl<'a> Octree<'a> {
    /// Find neighboring points
    ///
    /// Parameters
    /// ----------
    /// radius: the searching radius
    ///
    /// Return
    /// ------
    /// indices of neighboring points in pairs
    ///
    pub fn neighbors(&self, radius: f64) -> Vec<(usize, usize)>{
        let radius2 = radius*radius;

        let mut pairs = vec![];
        let mut count = 0;
        for octant in self.octants.iter() {
            if octant.children.is_empty() {
                count += 1;
                // println!("octant ranking = {:?}", octant.ranking);
                for (i, &pi) in octant.ipoints.iter().enumerate() {
                    for (j, &pj) in octant.ipoints.iter().enumerate().skip(i+1) {
                        let pi = self.points[pi];
                        let pj = self.points[pj];
                        let px = pj[0] - pi[0];
                        let py = pj[1] - pi[1];
                        let pz = pj[2] - pi[2];
                        let d2 = px*px + py*py + pz*pz;
                        if d2 < radius2 {
                            pairs.push((i, j));
                        }
                    }
                }
            }
        }
        println!("non-empty octants = {:?}", count);

        pairs
    }
}
// 7e3b12c9-d3f8-4bfc-8ed0-46e2644660d3 ends here
