// demo.rs
// :PROPERTIES:
// :header-args: :tangle examples/demo.rs
// :END:

// [[file:~/Workspace/Programming/gchemol-rs/octree/octree.note::*demo.rs][demo.rs:1]]
#[macro_use]
extern crate timeit;

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

fn main() {
    use octree::*;

    // external xyz file
    let stream = include_str!("data/pdb4rhv.xyz");
    let points = read_points(stream);

    let mut tree = Octree::new(points.clone());
    let bucket_size = 8 * 8;
    tree.build(bucket_size);

    let stream = include_str!("data/result.txt");
    for (line, &p) in stream.lines().zip(points.iter()) {
        let expected: Vec<_> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let x = tree.search(p, 3.0);
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
