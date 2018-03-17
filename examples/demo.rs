// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::b9a6ab02-bfca-4340-8aec-16fe3c042a9d][b9a6ab02-bfca-4340-8aec-16fe3c042a9d]]
#[macro_use]
extern crate timeit;
extern crate octree;

fn main() {
    // external xyz file
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

    timeit!({
        for &q in tree.points.iter() {
            tree.search(q, 3.0);
        }
    });
}
// b9a6ab02-bfca-4340-8aec-16fe3c042a9d ends here
