# TODO [2019-12-15 Sun] 测试benchmark数据                          :ATTACH:
# :PROPERTIES:
# :ID:       26fc4bcb-66e3-43aa-93f2-3e154212bd95
# :END:

# [[file:~/Workspace/Programming/gchemol-rs/octree/octree.note::*[2019-12-15 Sun] 测试benchmark数据][[2019-12-15 Sun] 测试benchmark数据:1]]
import ase.io
import numpy as np

def load_pts(filename):
    return ase.io.read(filename).positions

test_file_path="examples/data/3wu2.xyz"
points = load_pts(test_file_path)

def get_ckdtree(pts):
    import scipy.spatial as spatial
    tree = spatial.cKDTree(pts, leafsize=64)
    return tree

def test_ckdtree(pts, tree, cutoff):
    x = tree.query_ball_point(pts, cutoff)
    return x

def run_test():
    tree = get_ckdtree(points)
    test_ckdtree(points, tree, cutoff=5)
# [2019-12-15 Sun] 测试benchmark数据:1 ends here
