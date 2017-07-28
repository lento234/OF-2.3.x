"""
python script to make blockMeshDict
"""

import numpy as np
import os

header =\
r"""
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.4.1                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"""

convertToMeters = 1


H = 4

# Vertices
#xvertices = np.array([-H/2-10*H, -H/2, H/2,   H/2+10*H])
#yvertices = np.array([-H/2-10*H, -H/2, H/2,   H/2+10*H])
#zvertices = np.array([0,          H, H+H, H+H+10*H])

#xvertices = np.array([-1.5*H, -0.5*H, 0.5*H, 1.5*H])
#yvertices = np.array([-1.5*H, -0.5*H, 0.5*H, 1.5*H])
#zvertices = np.array([0,1.0*H, 2.0*H, 3.0*H])

xvertices = np.array([-25, -5, 5, 50])
yvertices = np.array([-25, -5, 5, 25])
zvertices = np.array([0, 5, 10, 30])
# Blocks
nx = np.array([20,20,20])
ny = np.array([20,20,20])
nz = np.array([10,10,20])

# Resolutions
xgrad = np.array([1,1,1])
ygrad = np.array([1,1,1])
zgrad = np.array([1,1,1])

# Nblocks
nvertices = (xvertices.shape[0],yvertices.shape[0],zvertices.shape[0])
nblock    = (nx.shape[0], ny.shape[0], nz.shape[0])

outputfile = 'constant/polyMesh/blockMeshDict'


with open(outputfile,'w') as f:
    # Header text
    f.write(header)
    f.write('\n\n') # space

    # convert to meters
    f.write('convertToMeters %g;' % convertToMeters)
    f.write('\n\n') # space

    # vertices
    f.write('vertices\n(\n')
    for k,kvertices in enumerate(zvertices):
        for j,jvertices in enumerate(yvertices):
            for i,ivertices in enumerate(xvertices):
                f.write('\t\t(%g %g %g)' % (ivertices,jvertices,kvertices))
                f.write(' // vertex #%d\n' % ((nvertices[0]*nvertices[1])*k + nvertices[0]*j + i))

    f.write(');')
    f.write('\n\n') # space

    # blocks
    f.write('blocks\n(\n')

    iblock = 0
    for k, (knz,kzgrad) in enumerate(zip(nz,zgrad)):
        for j, (jny,jygrad) in enumerate(zip(ny,ygrad)):
            for i, (inx,ixgrad) in enumerate(zip(nx,xgrad)):

                # Points
                p1 = (nvertices[0]*nvertices[1])*k + nvertices[0]*j + i
                p2 = (nvertices[0]*nvertices[1])*k + nvertices[0]*j + i+1
                p3 = (nvertices[0]*nvertices[1])*k + nvertices[0]*(j+1) + i+1
                p4 = (nvertices[0]*nvertices[1])*k + nvertices[0]*(j+1) + i
                p5 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*j + i
                p6 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*j + i+1
                p7 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*(j+1) + i+1
                p8 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*(j+1) + i

                # Hex box
                f.write('\t\thex (%d %d %d %d %d %d %d %d)' % (p1, p2, p3, p4, p5, p6, p7, p8))

                # Fluid description
                f.write('\t(%d %d %d) simpleGrading (%g %g %g)' % (inx, jny, knz, ixgrad, jygrad, kzgrad))
                f.write(' // block #%d\n' % ((nblock[0]*nblock[1])*k + nblock[0]*j + i) )
                iblock += 1

    f.write(');')
    f.write('\n\n')

    # edges
    f.write('edges\n(\n);\n\n')

    # boundary
    f.write('boundary\n(\n')

    ################# GROUND
    f.write(2*'\t'+'ground\n'+2*'\t'+'{\n')
    f.write(4*'\t'+'type wall;\n')
    f.write(4*'\t'+'faces\n'+4*'\t'+'(\n')
    k = 0
    for j in range(nblock[1]):
        for i in range(nblock[0]):
            p1 = (nvertices[0]*nvertices[1])*k + nvertices[0]*j + i
            p2 = (nvertices[0]*nvertices[1])*k + nvertices[0]*j + i+1
            p3 = (nvertices[0]*nvertices[1])*k + nvertices[0]*(j+1) + i+1
            p4 = (nvertices[0]*nvertices[1])*k + nvertices[0]*(j+1) + i
            f.write(6*'\t'+'(%d %d %d %d)\n' % (p1, p2, p3, p4))
    f.write(4*'\t'+');\n')
    f.write(2*'\t'+'}\n')

    ################# OUTLET
    f.write(2*'\t'+'outlet\n'+2*'\t'+'{\n')
    f.write(4*'\t'+'type patch;\n')
    f.write(4*'\t'+'faces\n'+4*'\t'+'(\n')
    i = nblock[0]
    for k in range(nblock[2]):
        for j in range(nblock[1]):
            p1 = (nvertices[0]*nvertices[1])*k + nvertices[0]*j + i
            p2 = (nvertices[0]*nvertices[1])*k + nvertices[0]*(j+1) + i
            p3 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*(j+1) + i
            p4 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*j + i
            f.write(6*'\t'+'(%d %d %d %d)\n' % (p1, p2, p3, p4))
    f.write(4*'\t'+');\n')
    f.write(2*'\t'+'}\n')

    ################# SKY
    f.write(2*'\t'+'sky\n'+2*'\t'+'{\n')
    f.write(4*'\t'+'type wall;\n')
    f.write(4*'\t'+'faces\n'+4*'\t'+'(\n')
    k = nblock[2]
    for j in range(nblock[1]):
        for i in range(nblock[0]):
            p1 = (nvertices[0]*nvertices[1])*k + nvertices[0]*j + i
            p2 = (nvertices[0]*nvertices[1])*k + nvertices[0]*j + i+1
            p3 = (nvertices[0]*nvertices[1])*k + nvertices[0]*(j+1) + i+1
            p4 = (nvertices[0]*nvertices[1])*k + nvertices[0]*(j+1) + i
            f.write(6*'\t'+'(%d %d %d %d)\n' % (p1, p2, p3, p4))
    f.write(4*'\t'+');\n')
    f.write(2*'\t'+'}\n')

    ################# INLET
    f.write(2*'\t'+'inlet\n'+2*'\t'+'{\n')
    f.write(4*'\t'+'type patch;\n')
    f.write(4*'\t'+'faces\n'+4*'\t'+'(\n')
    i = 0
    for k in range(nblock[2]):
        for j in range(nblock[1]):
            p1 = (nvertices[0]*nvertices[1])*k + nvertices[0]*j + i
            p2 = (nvertices[0]*nvertices[1])*k + nvertices[0]*(j+1) + i
            p3 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*(j+1) + i
            p4 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*j + i
            f.write(6*'\t'+'(%d %d %d %d)\n' % (p1, p2, p3, p4))
    f.write(4*'\t'+');\n')
    f.write(2*'\t'+'}\n')

    ################# FRONT
    f.write(2*'\t'+'front\n'+2*'\t'+'{\n')
    f.write(4*'\t'+'type empty;\n')
    f.write(4*'\t'+'faces\n'+4*'\t'+'(\n')
    j = 0
    for k in range(nblock[2]):
        for i in range(nblock[0]):
            p1 = (nvertices[0]*nvertices[1])*k + nvertices[0]*j + i
            p2 = (nvertices[0]*nvertices[1])*k + nvertices[0]*j + i+1
            p3 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*j + i+1
            p4 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*j + i
            f.write(6*'\t'+'(%d %d %d %d)\n' % (p1, p2, p3, p4))
    f.write(4*'\t'+');\n')
    f.write(2*'\t'+'}\n')


    ################# BACK
    f.write(2*'\t'+'back\n'+2*'\t'+'{\n')
    f.write(4*'\t'+'type empty;\n')
    f.write(4*'\t'+'faces\n'+4*'\t'+'(\n')
    j = nblock[1]
    for k in range(nblock[2]):
        for i in range(nblock[0]):
            p1 = (nvertices[0]*nvertices[1])*k + nvertices[0]*j + i
            p2 = (nvertices[0]*nvertices[1])*k + nvertices[0]*j + i+1
            p3 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*j + i+1
            p4 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*j + i
            f.write(6*'\t'+'(%d %d %d %d)\n' % (p1, p2, p3, p4))
    f.write(4*'\t'+');\n')
    f.write(2*'\t'+'}\n')

    f.write(');')
    f.write('\n\n')

    # mergePatchPairs
    f.write('mergePatchPairs\n(\n);\n')

    # write
    f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n')
