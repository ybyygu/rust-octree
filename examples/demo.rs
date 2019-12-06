// demo.rs
// :PROPERTIES:
// :header-args: :tangle examples/demo.rs
// :END:

// [[file:~/Workspace/Programming/rust-libs/rust-octree/rust-octree.note::*demo.rs][demo.rs:1]]
#[macro_use]
extern crate timeit;
extern crate octree;

fn main() {
    // external xyz file
    use octree::*;

    let stream = include_str!("data/pdb4rhv.xyz");
    let points = get_positions_from_xyz_stream(stream).unwrap();

    let mut tree = Octree::new(&points);
    tree.bucket_size = 8 * 8;
    tree.build();

    let stream = include_str!("data/result.txt");
    for (line, &p) in stream.lines().zip(points.iter()) {
        let mut expected: Vec<_> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let mut x = tree.search(p, 3.0);
        let mut y: Vec<_> = x.iter().map(|v| v.0).collect();
        y.sort();
        assert_eq!(y, expected);
    }

    timeit!({
        for &q in tree.points.iter() {
            tree.search(q, 3.0);
        }
    });

    timeit!({
        tree.neighbors(3.0);
    });
}
// demo.rs:1 ends here
