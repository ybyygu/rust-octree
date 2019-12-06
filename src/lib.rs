#[cfg(test)]
#[macro_use]
extern crate timeit;
#[cfg(test)]
#[macro_use]
extern crate approx;

use std::error;
use std::fs::File;
use std::io::prelude::*;
use std::io::{self, BufReader};

mod octree;
mod types;
pub use self::octree::Octree;

use crate::types::{Point, Points};

pub type Result<T> = std::result::Result<T, Box<error::Error>>;

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
