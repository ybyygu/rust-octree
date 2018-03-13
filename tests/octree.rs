// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::2f321e55-2849-4b73-aaf9-cd3271843de0][2f321e55-2849-4b73-aaf9-cd3271843de0]]
#[macro_use]
extern crate timeit;
extern crate octree;
// 2f321e55-2849-4b73-aaf9-cd3271843de0 ends here

// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::feb2e7b9-8cca-4210-a89d-a7f1d2a40d9e][feb2e7b9-8cca-4210-a89d-a7f1d2a40d9e]]
#[test]
fn test_octree() {
    use octree::*;

    let stream = include_str!("data/test.xyz");
    let points = get_positions_from_xyz_stream(stream).unwrap();
    let q = points[0];
    let mut tree = Octree::new(&points);
    tree.bucket_size = 1;
    tree.build();
    let x = tree.search(q, 2.2);
    println!("neighbors: {:?}", x);
}
// feb2e7b9-8cca-4210-a89d-a7f1d2a40d9e ends here
