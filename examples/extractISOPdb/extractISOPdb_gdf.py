#!/usr/bin/env python 

import argparse
import numpy as np
from gridData import Grid, OpenDX

def parse_args():
    parser = argparse.ArgumentParser(description='Save grid to PDB a given isovalue')
    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('iso', type=float)
    return parser.parse_args()    

def extractISOPdb(input, output, iso):
    g = Grid(input)
    
    # JD: this is clunky but i'm not sure what's a better way
    data = np.array(OpenDX.array(3,g.grid).array.flat)
    data = data.reshape((len(data)/3, 3))

    x, y, z = 0, 0, 0
    counter = 1
    lines = []
    for point in data:
        for pos in point:
            if (iso<0 and pos < iso) or (iso > 0 and pos > iso) :
                line = 'ATOM  %5d  C   PTH     1    %8.3f%8.3f%8.3f%6.2f%6.2f\n'%(counter,g.origin[0]+float(x)*g.delta[0,0],
                                                                                          g.origin[1]+float(y)*g.delta[1,1],
                                                                                          g.origin[2]+float(z)*g.delta[2,2], 0.0,0.0)
                lines.append(line)
                counter += 1
            z+=1
            if z >= g.grid.shape[2]:
                z=0
                y+=1
                if y >=g.grid.shape[1]:
                    y=0
                    x+=1
    
    with open(output, "w") as f:
        f.writelines(lines)            
    
    return g
if __name__ == "__main__":
    args = parse_args()
    g = extractISOPdb(args.input, args.output, args.iso)