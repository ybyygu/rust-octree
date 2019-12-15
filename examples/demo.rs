#[macro_use]
extern crate timeit;

use rayon::prelude::*;

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

fn main() {
    use octree::*;

    // external xyz file containing 51053 points
    let stream = include_str!("data/3wu2.xyz");
    let points = read_points(stream);

    let mut tree = Octree::new(points.clone());
    let bucket_size = 64;
    tree.build(bucket_size);

    timeit!({
        points.par_iter().for_each(|&q| {
            let neighbors: Vec<_> = tree.search(q, 5.0).collect();
        })
    });
}
