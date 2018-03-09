// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::b9a6ab02-bfca-4340-8aec-16fe3c042a9d][b9a6ab02-bfca-4340-8aec-16fe3c042a9d]]
#[macro_use]
extern crate timeit;
extern crate octree;

fn main() {
    // external xyz file
    use octree::*;

    let stream = include_str!("data/pdb4rhv.xyz");
    let points = get_positions_from_xyz_stream(stream).unwrap();

    let q = points[0];
    let mut tree = Octree::new(points);
    tree.bucket_size = 8*8;
    tree.build();

    let x = tree.neighbors(q, 3.0);
    assert!(x.contains(&0));
    assert!(x.contains(&1241));

    println!("neighbors: {:?}", x);

    timeit!({
        for &q in tree.points.iter() {
            tree.neighbors(q, 3.0);
        }
    });
}
// b9a6ab02-bfca-4340-8aec-16fe3c042a9d ends here
