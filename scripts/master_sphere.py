# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                                         #
# Copyright (c) 2018, Marc De Graef Research Group/Carnegie Mellon University             #
# All rights reserved.                                                                    #
#                                                                                         #
# Author William C. Lenthe                                                                #
#                                                                                         #
# Redistribution and use in source and binary forms, with or without modification, are    #
# permitted provided that the following conditions are met:                               #
#                                                                                         #
#     - Redistributions of source code must retain the above copyright notice, this list  #
#        of conditions and the following disclaimer.                                      #
#     - Redistributions in binary form must reproduce the above copyright notice, this    #
#        list of conditions and the following disclaimer in the documentation and/or      #
#        other materials provided with the distribution.                                  #
#     - Neither the names of Marc De Graef, Carnegie Mellon University nor the names      #
#        of its contributors may be used to endorse or promote products derived from      #
#        this software without specific prior written permission.                         #
#                                                                                         #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"             #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE               #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE              #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE               #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL       #
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR              #
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER              #
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,           #
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE               #
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                #
#                                                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#@brief           : construct a spherical surface mesh of a master pattern for paraview visualization
#@param inputFile : EMsoft master pattern file (.h5)
#@param outputFile: prefix for output file names (outputFile.h5 and outputFile.xdmf)
#@note            : drag outputFile.xdmf into paraview and select 'Xdmf Reader' (not 'Xdmf3 Reader')

import h5py
import numpy as np

inputFile = 'Fo-master.h5' # EMsoft output file
outputFile = 'output' # what should we name the created h5/xdmf file

# inverse square lambert projection (square (-1 <= x,y <= 1) to unit sphere)
def equalAreaSquareToSphere(X, Y):
	# restrict to +/-1 box and determine largest absolute value coordinate
	X = np.clip(X, -1, 1)
	Y = np.clip(Y, -1, 1)
	aX = np.abs(X)
	aY = np.abs(Y)
	vMax = np.maximum(aX, aY)

	# compute x and y
	X0 = np.where(X == 0, 1, X) # don't divide by 0
	Y0 = np.where(Y == 0, 1, Y) # don't divide by 0
	q1  = Y * np.sqrt(2 - Y * Y) # aX <= aY
	q2  = X * np.sqrt(2 - X * X) # aX >  aY
	qq1 = (X * np.pi) / (Y0 * 4) # aX <= aY
	qq2 = (Y * np.pi) / (X0 * 4) # aX >  aY
	xyz = np.zeros(X.shape + (3,))
	xyz[...,0] = np.where(aX <= aY, q1 * np.sin(qq1), q2 * np.cos(qq2))
	xyz[...,1] = np.where(aX <= aY, q1 * np.cos(qq1), q2 * np.sin(qq2))

	# compute z and normalize
	xyz[...,2] = 1 - vMax * vMax
	mag = np.linalg.norm(xyz, axis = -1)
	xyz[...,0] /= mag
	xyz[...,1] /= mag
	xyz[...,2] /= mag
	return xyz;

# read mc and crystal data
f = h5py.File(inputFile)
oc = f['/CrystalData/AtomData'][3] # read occupancies
en = np.sum(f['/EMData/MCOpenCL/accum_e'], axis=(0,1)).astype('float32')

# read master patterns for each atom and accumulate occupancy weighted value
nh = f['/EMData/EBSDmaster/mLPNH'][0] * oc[0]
sh = f['/EMData/EBSDmaster/mLPSH'][0] * oc[0]
for i in range(1, len(oc)):
	nh += f['/EMData/EBSDmaster/mLPNH'][i] * oc[i]
	sh += f['/EMData/EBSDmaster/mLPSH'][i] * oc[i]

# convert energy bins to percentages and compute weighted average of patterns
en = en / np.sum(en)
nh = np.average(nh, axis = 0, weights = en)
sh = np.average(sh, axis = 0, weights = en)

# compute back projected coordinates of each point and vectorize
X = np.linspace(-1, 1, nh.shape[1]).astype('float32')
Y = np.linspace(-1, 1, nh.shape[0]).astype('float32')
[X, Y] = np.meshgrid(X, Y);
xyz = equalAreaSquareToSphere(X, Y)
xyz = np.reshape(xyz, (nh.shape[0]*nh.shape[1], 3))#(x,y,3) -> (x*y,3)

# build quad mesh (squares in ccw order from top left as v0, v1, v2, v3)
v0 = np.tile(np.arange(0, nh.shape[1]-1, 1, dtype = 'uint32'), nh.shape[0]-1) # 0, 1, 2, ..., x-2, x-1, 0, 1, 2, ... y-1 times
v0 = v0 + np.repeat(np.arange(0, nh.shape[0]-1, 1), nh.shape[1]-1) * nh.shape[1] # indicies of all points except for last row/col in row major order
v1 = v0 + nh.shape[1]
v2 = v1 + 1
v3 = v0 + 1

# compute shortest diagonal of each quad and split into trangles
d02 = np.linalg.norm(xyz[v0] - xyz[v2], axis = -1)
d13 = np.linalg.norm(xyz[v1] - xyz[v3], axis = -1)
use02 = d02 <= d13
tris = np.empty((v0.shape[0], 6), dtype = 'uint32')
tris[:,0] = v0
tris[:,1] = v1
tris[:,2] = np.where(use02, v2, v3)
tris[:,3] = v2
tris[:,4] = v3
tris[:,5] = np.where(use02, v0, v1)
tris = tris.reshape((v0.shape[0] * 2, 3))

# add southern hemisphere
trisSh = np.empty_like(tris)
trisSh[:,0], trisSh[:,1], trisSh[:,2] = tris[:,0], tris[:,2], tris[:,1] # flip winding
trisSh = trisSh + xyz.shape[0]
xyzSh = np.copy(xyz)
xyzSh[:,2] = -xyzSh[:,2] # move to southern hemisphere
tris = np.concatenate((tris, trisSh))
xyz = np.concatenate((xyz, xyzSh))

# write mesh to hdf5
f = h5py.File(outputFile + '.h5', 'w')
f['/verts'] = xyz
f['/tris'] = tris
f['/scalar'] = np.concatenate((nh.reshape((nh.shape[0] * nh.shape[1])), sh.reshape((sh.shape[0] * sh.shape[1]))))

# write wrapper for paraview
with open(outputFile + '.xdmf', 'w') as file:
	file.write('<?xml version="1.0"?>\n')
	file.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd"[]>\n')
	file.write('<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">\n')
	file.write(' <Domain>\n')
	file.write('  <Grid Name="MasterPatternMesh" GridType="Uniform">\n')

	# triangles
	file.write('   <Topology TopologyType="Triangle" NumberOfElements="%i">\n' % tris.shape[0])
	file.write('    <DataItem Format="HDF" NumberType="Int" Dimensions="%i 3">\n' % tris.shape[0])
	file.write('     ' + outputFile + '.h5:/tris\n')
	file.write('    </DataItem>\n')
	file.write('   </Topology>\n')

	# points
	file.write('   <Geometry Type="XYZ">\n')
	file.write('    <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="%i 3">\n' % xyz.shape[0])
	file.write('     ' + outputFile + '.h5:/verts\n')
	file.write('    </DataItem>\n')
	file.write('   </Geometry>\n')

	# pattern
	file.write('   <Attribute Name="MasterPattern" AttributeType="Scalar" Center="Node">\n') # node/cell for vertex/triangle attributes
	file.write('    <DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="%i 1" >\n' % xyz.shape[0])
	file.write('     ' + outputFile + '.h5:/scalar\n')
	file.write('    </DataItem>\n')
	file.write('   </Attribute>\n')

	file.write('  </Grid>\n')
	file.write(' </Domain>\n')
	file.write('</Xdmf>\n')
