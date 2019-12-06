// test_octree.rs
// :PROPERTIES:
// :header-args: :tangle tests/octree.rs
// :ID:       38431531-4955-4c81-9570-86c776602192
// :END:

// [[file:~/Workspace/Programming/octree/octree.note::*test_octree.rs][test_octree.rs:1]]
#[test]
fn test_octree() {
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
}
// test_octree.rs:1 ends here
