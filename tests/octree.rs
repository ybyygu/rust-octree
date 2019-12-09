// test_octree.rs
// :PROPERTIES:
// :header-args: :tangle tests/octree.rs
// :ID:       38431531-4955-4c81-9570-86c776602192
// :END:

// [[file:~/Workspace/Programming/gchemol-rs/octree/octree.note::*test_octree.rs][test_octree.rs:1]]
#[test]
fn test_octree() {
    use octree::*;

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

    let stream = include_str!("data/pdb4rhv.xyz");
    let points = read_points(stream);

    let mut tree = Octree::new(&points);
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
}
// test_octree.rs:1 ends here
