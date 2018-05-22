####################################################################################
# This is free software; you can redistribute it and/or modify it under the		     #
# terms of the GNU General Public License as published by the Free Software        #
# Foundation; either version 2 of the License, or any later version.               #
#                                                                                  #
# This code is distributed in the hope that it will be useful, but WITHOUT ANY     #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  #
# PARTICULAR PURPOSE. See the GNU General Public License for more details.         #
#                                                                                  #
# You should have received a copy of the GNU General Public License along with     #
# this code; if not, write to the Free Software Foundation, Inc., 51 Franklin      #
# Street, Fifth Floor, Boston, MA 02110-1301, USA                                  #
####################################################################################

This is the reference implementation of the semi-convex hull algorithm proposed in
'Displets: Resolving Stereo Ambiguities using Object Knowledge' as described in 

@INPROCEEDINGS{Guney2015CVPR,
  author = {Fatma Güney and Andreas Geiger},
  title = {Displets: Resolving Stereo Ambiguities using Object Knowledge},
  booktitle = {Conference on Computer Vision and Pattern Recognition (CVPR)},
  year = {2015}
}

It can be used to simplify 3D CAD models to a coarse low-face representation which
renders blazingly fast on a GPU and can be used for inverse graphics or visualization.

Compilation
===========

Run make -j12 in the remesher directory from the Linux command line.

Running the Demo
================

run the command: run_semiConvexHull from within Matlab (tested on Matlab 2014b)

This will show a demo, reading the mesh stored in data and writing the semi-convex
hull at the end. The 3D mesh is visualized during optimization.

Input Data
==========

You can use the code provided with librender to convert your own 3D Warehouse models
or other CAD models in object format to simple vertex/face meshes which can be read:

http://www.cvlibs.net/software/librender/

Dependencies
============

To make the zip file self-contained, it includes copies of

a) QSlim: http://www.cs.cmu.edu/afs/cs/Web/People/garland/quadrics/qslim.html
    
b) Remesher: http://www.gris.informatik.tu-darmstadt.de/~sfuhrman/remesher.html
   
Citation
========

If you find this code useful, we would be happy if you cite us:

@INPROCEEDINGS{Guney2015CVPR,
  author = {Fatma Güney and Andreas Geiger},
  title = {Displets: Resolving Stereo Ambiguities using Object Knowledge},
  booktitle = {Conference on Computer Vision and Pattern Recognition (CVPR)},
  year = {2015}
}

