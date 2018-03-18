// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::7896c853-796e-4f97-b1d4-45546d4491db][7896c853-796e-4f97-b1d4-45546d4491db]]
#[macro_use]
extern crate timeit;

#[macro_use]
extern crate approx;
// 7896c853-796e-4f97-b1d4-45546d4491db ends here

// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::ca2234bb-c5da-477d-8a3b-d85edc46ddf1][ca2234bb-c5da-477d-8a3b-d85edc46ddf1]]
use std::error;
use std::io::{self, BufReader};
use std::io::prelude::*;
use std::fs::File;

pub mod octree;

pub use octree::Octree;
// ca2234bb-c5da-477d-8a3b-d85edc46ddf1 ends here

// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::9a8ceae5-46dd-4db5-a20d-2d78766568a3][9a8ceae5-46dd-4db5-a20d-2d78766568a3]]
// set up aliases for convenience
pub type Point = [f64; 3];
pub type Points = Vec<Point>;

pub type Result<T> = std::result::Result<T, Box<error::Error>>;
// 9a8ceae5-46dd-4db5-a20d-2d78766568a3 ends here

// [[file:~/Workspace/Programming/rust-octree/rust-octree.note::ddd95be2-3bf5-4478-a21b-ea0c8742f5cb][ddd95be2-3bf5-4478-a21b-ea0c8742f5cb]]
pub fn get_positions_from_xyz_stream(txt: &str) -> Result<Points> {
    let mut positions = Vec::new();

    for line in txt.lines() {
        let attrs: Vec<_> = line.split_whitespace().collect();
        let (symbol, position) = attrs.split_first().ok_or("encountering empty line")?;
        if position.len() != 3 {
            let msg = format!("informal xyz records: {}", line);
            Err(msg)?;
        }

        let p: Vec<f64> = position.iter().map(|x| x.parse().unwrap()).collect();
        positions.push([p[0], p[1], p[2]]);
    }

    Ok(positions)
}

// in a simple and dirty way
pub fn get_positions_from_xyzfile(filename: &str) -> Result<Points> {
    let mut buffer = String::new();
    let f = File::open(filename)?;
    let mut f = BufReader::new(f);

    f.read_line(&mut buffer);
    f.read_line(&mut buffer);
    buffer.clear();

    f.read_to_string(&mut buffer)?;
    get_positions_from_xyz_stream(&buffer)
}
// ddd95be2-3bf5-4478-a21b-ea0c8742f5cb ends here
