```@raw html
---
layout: doc
title: Getting Started
---
```

# Getting Started

Determine whether 3D test points lie inside or outside a watertight mesh defined by vertices and faces

## Installation

```julia
using Pkg
Pkg.add("InTriGeom")
```

## Basic Usage

```julia
using InTriGeom.jl
```


intriangulation: Test points in 3d wether inside or outside a (closed) triangulation
usage: in = intriangulation(vertices,faces,testp,heavytest)

arguments: (input)
vertices   - points in 3d as matrix with three columns

faces      - description of triangles as matrix with three columns.
               Each row contains three indices into the matrix of vertices
               which gives the three cornerpoints of the triangle.

  testp      - points in 3d as matrix with three columns

  heavytest  - int n >= 0. Perform n additional randomized rotation tests.

 IMPORTANT: the set of vertices and faces has to form a watertight surface!

arguments: (output)
  in - a vector of length size(testp,1), containing 0 and 1.
       in(nr) =  0: testp(nr,:) is outside the triangulation
       in(nr) =  1: testp(nr,:) is inside the triangulation
       in(nr) = -1: unable to decide for testp(nr,:) 

Thanks to Adam A for providing the FEX submission voxelise. The
algorithms of voxelise form the algorithmic kernel of intriangulation.

Thanks to Sven to discussions about speed and avoiding problems in
special cases.

Author: Johannes Korsawe, heavily based on voxelise from Adam A.
E-mail: johannes.korsawe@volkswagen.de
Release: 1.3
Release date: 25/09/2013

# heavytest = 0  → 1 test,  fast,  rare chance of -1 results
# heavytest = 1  → 2 tests, slower, very rare -1 results  
# heavytest = 3  → 4 tests, slowest, almost no -1 results


[intriangulation](https://es.mathworks.com/matlabcentral/fileexchange/43381-intriangulation-vertices-faces-testp-heavytest?s_tid=prof_contriblnk)
[Mesh voxelisation](https://es.mathworks.com/matlabcentral/fileexchange/27390-mesh-voxelisation)