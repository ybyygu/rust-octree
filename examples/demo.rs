#[macro_use]
extern crate timeit;

use rayon::prelude::*;

fn read_points_xyz(txt: &str) -> Vec<[f64; 3]> {
    let mut positions = Vec::new();
    for line in txt.lines().skip(2) {
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
    let points = read_points_xyz(stream);

    let x = timeit_loops!(10, {
        let mut tree = Octree::new(points.clone());
        let bucket_size = 64;
        tree.build(bucket_size);

        points.par_iter().for_each(|&q| {
            let neighbors: Vec<_> = tree.search(q, 8.0).collect();
        })
    });
    dbg!(x);
}
