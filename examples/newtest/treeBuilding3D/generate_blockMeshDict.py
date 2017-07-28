"""
python script to make blockMeshDict
"""

import numpy as np
import os

################################################################################

convertToMeters = 1

xvertices = np.array([-25, -5, 5, 15, 50])
yvertices = np.array([-25, -5, 5, 25])
zvertices = np.array([0, 5, 10, 30])

# Blocks
nx = np.array([20, 20, 30, 20])
ny = np.array([20, 20, 20])
nz = np.array([10, 20, 20])

# Resolutions
xgrad = np.array([1,-3,1,1])
ygrad = np.array([1,-3,1])
zgrad = np.array([1,1,1])

outputfile = 'constant/polyMesh/blockMeshDict'

left = 'patch', 'inlet'
right = 'patch', 'outlet'
front = 'symmetryPlane', 'front'
back  = 'symmetryPlane', 'back'
top = 'patch', 'sky'
bottom = 'wall', 'ground'

holes = [(1,1,0), (1,1,1)]
#holes = [(0, 0, 1), (0, 0, 2)]

custompatchIDs = [['wall', 'cube']]

custompatchVertices = [ [ ( (1,1,0), (1,2,0), (1,2,1), (1,1,1) ), # left
                          ( (1,1,1), (1,2,1), (1,2,2), (1,1,2) ), # left
                          ( (1,1,0), (2,1,0), (2,1,1), (1,1,1) ), # front
                          ( (1,1,1), (2,1,1), (2,1,2), (1,1,2) ), # front
                          ( (2,1,0), (2,2,0), (2,2,1), (2,1,1) ), # right
                          ( (2,1,1), (2,2,1), (2,2,2), (2,1,2) ), # right
                          ( (1,2,0), (2,2,0), (2,2,1), (1,2,1) ), # back
                          ( (1,2,1), (2,2,1), (2,2,2), (1,2,2) ), # back
                          ( (1,1,2), (2,1,2), (2,2,2), (1,2,2) ), # top
                           ]  ]


################################################################################

# Nblocks
nvertices = (xvertices.shape[0],yvertices.shape[0],zvertices.shape[0])
nblock    = (nx.shape[0], ny.shape[0], nz.shape[0])

iholes = [(nblock[0]*nblock[1])*k + nblock[0]*j + i for (i,j,k) in holes]

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

################################################################################

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
                if iblock in iholes:
                    f.write('\t\t//')
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
    f.write('patches\n(\n')

    ################# left
    f.write(2*'\t'+'%s %s\n' % (left[0], left[1]))
    f.write(2*'\t'+'(\n')
    i = 0
    for k in range(nblock[2]):
        for j in range(nblock[1]):
            if ((nblock[0]*nblock[1])*k + nblock[0]*j + i) in iholes:
                f.write(4*'\t'+'//')
            p1 = (nvertices[0]*nvertices[1])*k + nvertices[0]*j + i
            p2 = (nvertices[0]*nvertices[1])*k + nvertices[0]*(j+1) + i
            p3 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*(j+1) + i
            p4 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*j + i
            f.write(4*'\t'+'(%d %d %d %d)\n' % (p1, p2, p3, p4))
    f.write(2*'\t'+')\n')

    ################# right
    f.write(2*'\t'+'%s %s\n' % (right[0], right[1]))
    f.write(2*'\t'+'(\n')
    i = nblock[0]-1
    for k in range(nblock[2]):
        for j in range(nblock[1]):
            if ((nblock[0]*nblock[1])*k + nblock[0]*j + i) in iholes:
                f.write(4*'\t'+'//')
            p1 = (nvertices[0]*nvertices[1])*k + nvertices[0]*j + i+1
            p2 = (nvertices[0]*nvertices[1])*k + nvertices[0]*(j+1) + i+1
            p3 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*(j+1) + i+1
            p4 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*j + i+1
            f.write(4*'\t'+'(%d %d %d %d)\n' % (p1, p2, p3, p4))
    f.write(2*'\t'+')\n')

    ################# top
    f.write(2*'\t'+'%s %s\n' % (top[0], top[1]))
    f.write(2*'\t'+'(\n')
    k = nblock[2]-1
    for j in range(nblock[1]):
        for i in range(nblock[0]):
            if ((nblock[0]*nblock[1])*k + nblock[0]*j + i) in iholes:
                f.write(4*'\t'+'//')
            p1 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*j + i
            p2 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*j + i+1
            p3 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*(j+1) + i+1
            p4 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*(j+1) + i
            f.write(4*'\t'+'(%d %d %d %d)\n' % (p1, p2, p3, p4))
    f.write(2*'\t'+')\n')

    ################# bottom
    f.write(2*'\t'+'%s %s\n' % (bottom[0], bottom[1]))
    f.write(2*'\t'+'(\n')
    k = 0
    for j in range(nblock[1]):
        for i in range(nblock[0]):
            if ((nblock[0]*nblock[1])*k + nblock[0]*j + i) in iholes:
                f.write(4*'\t'+'//')
            p1 = (nvertices[0]*nvertices[1])*k + nvertices[0]*j + i
            p2 = (nvertices[0]*nvertices[1])*k + nvertices[0]*j + i+1
            p3 = (nvertices[0]*nvertices[1])*k + nvertices[0]*(j+1) + i+1
            p4 = (nvertices[0]*nvertices[1])*k + nvertices[0]*(j+1) + i
            f.write(4*'\t'+'(%d %d %d %d)\n' % (p1, p2, p3, p4))
    f.write(2*'\t'+')\n')

    ################# front
    f.write(2*'\t'+'%s %s\n' % (front[0], front[1]))
    f.write(2*'\t'+'(\n')
    j = 0
    for k in range(nblock[2]):
        for i in range(nblock[0]):
            if ((nblock[0]*nblock[1])*k + nblock[0]*j + i) in iholes:
                f.write(4*'\t'+'//')
            p1 = (nvertices[0]*nvertices[1])*k + nvertices[0]*j + i
            p2 = (nvertices[0]*nvertices[1])*k + nvertices[0]*j + i+1
            p3 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*j + i+1
            p4 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*j + i
            f.write(4*'\t'+'(%d %d %d %d)\n' % (p1, p2, p3, p4))
    f.write(2*'\t'+')\n')


    ################# back
    f.write(2*'\t'+'%s %s\n' % (back[0], back[1]))
    f.write(2*'\t'+'(\n')
    j = nblock[1]-1
    for k in range(nblock[2]):
        for i in range(nblock[0]):
            if ((nblock[0]*nblock[1])*k + nblock[0]*j + i) in iholes:
                f.write(4*'\t'+'//')
            p1 = (nvertices[0]*nvertices[1])*k + nvertices[0]*(j+1) + i
            p2 = (nvertices[0]*nvertices[1])*k + nvertices[0]*(j+1) + i+1
            p3 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*(j+1) + i+1
            p4 = (nvertices[0]*nvertices[1])*(k+1) + nvertices[0]*(j+1) + i
            f.write(4*'\t'+'(%d %d %d %d)\n' % (p1, p2, p3, p4))
    f.write(2*'\t'+')\n')

    ################# custom patches
    for patchID,vertices in zip(custompatchIDs, custompatchVertices):
        f.write(2*'\t'+'%s %s\n' % (patchID[0], patchID[1]))
        f.write(2*'\t'+'(\n')
        for ivertex in vertices:
            p = []
            for i,j,k in ivertex:
                p.append((nvertices[0]*nvertices[1])*k + nvertices[0]*j + i)
            f.write(4*'\t'+'(%d %d %d %d)\n' % (p[0], p[1], p[2], p[3]))
        f.write(2*'\t'+')\n')


    f.write(');')
    f.write('\n\n')

    # mergePatchPairs
    f.write('mergePatchPairs\n(\n);\n')

    # write
    f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n')
