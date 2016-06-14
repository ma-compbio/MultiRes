#!/usr/bin/env

import sys
import DataStructures.markers as markers

if __name__ == "__main__":
    blocks = markers.read_hom_families_file(sys.argv[1])
    doubled_blocks,identifier = markers.double_hom_families(blocks,1) 
    markers.write_hom_families_file(sys.argv[2],doubled_blocks) 
