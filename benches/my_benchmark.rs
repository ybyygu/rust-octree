use criterion::{criterion_group, criterion_main, Criterion};

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

fn octree_search(points: Vec<[f64; 3]>) {
    use octree::*;

    let mut tree = Octree::new(points);
    let bucket_size = 8;
    tree.build(bucket_size);

    let cutoff = 3.0;
    for &q in tree.points.iter() {
        tree.search(q, cutoff);
    }
}

fn octree_build(points: Vec<[f64; 3]>) {
    use octree::Octree;

    let mut tree = Octree::new(points);
    let bucket_size = 1;
    tree.build(bucket_size);
}

fn criterion_benchmark(c: &mut Criterion) {
    // external xyz file
    let stream = include_str!("data/1crn.xyz");
    let points = read_points(stream);

    c.bench_function("octree build", |b| {
        b.iter(|| octree_build(points.clone()))
    });
    c.bench_function("octree search", |b| {
        b.iter(|| octree_search(points.clone()))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
