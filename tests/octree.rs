// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::2f321e55-2849-4b73-aaf9-cd3271843de0][2f321e55-2849-4b73-aaf9-cd3271843de0]]
#[macro_use]
extern crate octree;
// 2f321e55-2849-4b73-aaf9-cd3271843de0 ends here

// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::feb2e7b9-8cca-4210-a89d-a7f1d2a40d9e][feb2e7b9-8cca-4210-a89d-a7f1d2a40d9e]]
#[test]
fn test_octree() {
    use octree::*;

    let stream = include_str!("data/pdb4rhv.xyz");
    let points = get_positions_from_xyz_stream(stream).unwrap();

    let mut tree = Octree::new(&points);
    tree.bucket_size = 8*8;
    tree.build();

    let stream = include_str!("data/result.txt");
    for (line, &p) in stream.lines().zip(points.iter()) {
        let mut expected: Vec<_> = line.split_whitespace().map(|x| x.parse().unwrap()).collect();
        let mut x = tree.search(p, 3.0);
        x.sort();
        assert_eq!(x, expected);
    }
}
// feb2e7b9-8cca-4210-a89d-a7f1d2a40d9e ends here
