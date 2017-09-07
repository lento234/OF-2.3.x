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

"""
xvertices = np.array([0,1,2,3])
yvertices = np.array([0,1,2])
zvertices = np.array([0,1,2,3])

nx = np.array([10,10,10,10,10])
ny = np.array([10,10])
nz = np.array([10,10,10])

xgrad = np.array([1,1,1])
ygrad = np.array([1,1])
zgrad = np.array([1,1,1])
"""
"""
xvertices = np.array([-9.2,-0.5,0.5,25.5])
yvertices = np.array([-0.5, 0.5])
zvertices = np.array([0,0.5,1.5,11.5])

nx = np.array([81, 82, 99])
ny = np.array([1])
nz = np.array([35, 82, 36])

xgrad = np.array([1, 1, 1])
ygrad = np.array([1])
zgrad = np.array([1, 1, 1])
"""

# Vertices
xvertices = np.array([0,1,2,3])
yvertices = np.array([0,1,2])
zvertices = np.array([0,1,2,3])

# Blocks
nx = np.array([10,10,10])
ny = np.array([10,10])
nz = np.array([10,10,10])

# Resolutions
xgrad = np.array([1,1,1])
ygrad = np.array([1,1])
zgrad = np.array([1,1,1])

"""
xvertices = np.array([0,1])
yvertices = np.array([0,1])
zvertices = np.array([0,1])

nx = np.array([10])
ny = np.array([10])
nz = np.array([10])

xgrad = np.array([1])
ygrad = np.array([1])
zgrad = np.array([1])
"""

# Faces
#faces = ('ground',  'outlet',   'sky',  'inlet', 'front', 'back')
#types = ('wall',    'patch',    'wall', 'patch', 'empty', 'empty')


# Nblocks
nvertices = (xvertices.shape[0],yvertices.shape[0],zvertices.shape[0])
nblock    = (nx.shape[0], ny.shape[0], nz.shape[0])

with open('blockMeshDict','w') as f:
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
                f.write('\t\t(%g %g %g)' % (i,j,k))
                f.write(' // vertices #%d\n' % ((nvertices[0]*nvertices[1])*k + nvertices[0]*j + i))

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
    """
    for iface, itype in zip(faces,types):
        f.write('\t\t%s\n\t\t{\n' % iface)
        f.write('\t\t\t\ttype %s;\n' % itype)
        f.write('\t\t\t\tfaces\n\t\t\t\t(\n')
        f.write('\t\t\t\t)\n')
        f.write('\t\t}\n')
    """
    # ground
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

    f.write(');')
    f.write('\n\n')

    # mergePatchPairs
    f.write('mergePatchPairs\n(\n);\n')

    # write
    f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n')
