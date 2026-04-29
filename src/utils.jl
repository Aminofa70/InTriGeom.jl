"""
   rotmatrix(v, deg)
This function creates a matrix that rotates a 3D vector around a given axis `v` by a specified angle `deg`.

* The angle is converted from degrees to radians
* The axis `v` is normalized
* A rotation matrix is built using that axis and angle
* If the angle is 0, it returns the identity matrix (no rotation)

The output matrix `D` can be used to rotate 3D points.
"""
function rotmatrix(v, deg)
    deg = deg / 180 * π
    if deg != 0
        v = v / norm(v)
        v1, v2, v3 = v[1], v[2], v[3]
        ca, sa = cos(deg), sin(deg)
        D = [ca+v1*v1*(1-ca)      v1*v2*(1-ca)-v3*sa    v1*v3*(1-ca)+v2*sa;
             v2*v1*(1-ca)+v3*sa   ca+v2*v2*(1-ca)        v2*v3*(1-ca)-v1*sa;
             v3*v1*(1-ca)-v2*sa   v3*v2*(1-ca)+v1*sa     ca+v3*v3*(1-ca)]
    else
        D = Matrix{Float64}(I, 3, 3)
    end
    return D
end

"""

    VOXELISEinternal(testx, testy, testz, meshXYZ)
This function checks whether a set of 3D points lies inside a closed surface mesh.

For each test point `(testx, testy, testz)`, it casts a vertical ray (along the z-direction) and counts how many times this ray intersects the mesh triangles.

Steps:

* First, it finds nearby candidate triangles using bounding boxes (fast filtering)
* Then, it checks if the ray actually intersects each triangle
* For valid intersections, it computes the exact z-coordinates where the ray crosses the mesh
* These intersection points are sorted and paired

If the test point’s z-value lies between a pair of intersections, it means the point is inside the mesh.

The function returns:

* `OUTPUT`: a boolean array (true = inside, false = outside)
* `correctionLIST`: points with ambiguous intersections (numerical edge cases)
"""
function VOXELISEinternal(testx, testy, testz, meshXYZ)
    # meshXYZ is a (nfaces, 3, 3) array: dim1=face, dim2=coord(x,y,z), dim3=vertex(1,2,3)

    n = length(testx)
    OUTPUT = falses(n)
    correctionLIST = Int[]

    # Min/max z of entire mesh
    meshZmin = minimum(meshXYZ[:, 3, :])
    meshZmax = maximum(meshXYZ[:, 3, :])

    # Per-facet bounding boxes
    nm = size(meshXYZ, 1)
    meshXYZmin = minimum(meshXYZ, dims=3)[:, :, 1]  # (nm, 3)
    meshXYZmax = maximum(meshXYZ, dims=3)[:, :, 1]  # (nm, 3)

    facetCROSSLIST = zeros(Int, 1000)

    for loop in 1:n
        nf = 0

        # Step 1: Find candidate facets whose bounding box contains (testx, testy)
        possibleCROSSLISTy = findall(
            (testy[loop] .- meshXYZmin[:, 2]) .* (meshXYZmax[:, 2] .- testy[loop]) .> 0
        )
        possibleCROSSLISTx = (testx[loop] .- meshXYZmin[possibleCROSSLISTy, 1]) .*
                              (meshXYZmax[possibleCROSSLISTy, 1] .- testx[loop]) .> 0
        possibleCROSSLIST = possibleCROSSLISTy[possibleCROSSLISTx]

        isempty(possibleCROSSLIST) && continue

        # Step 2: Check if ray actually passes through each candidate facet
        for loopCHECKFACET in possibleCROSSLIST

            # Edge 2-3, opposing vertex 1
            denom = meshXYZ[loopCHECKFACET, 1, 2] - meshXYZ[loopCHECKFACET, 1, 3]
            Y1predicted = meshXYZ[loopCHECKFACET, 2, 2] -
                          (meshXYZ[loopCHECKFACET, 2, 2] - meshXYZ[loopCHECKFACET, 2, 3]) *
                          (meshXYZ[loopCHECKFACET, 1, 2] - meshXYZ[loopCHECKFACET, 1, 1]) / denom
            YRpredicted = meshXYZ[loopCHECKFACET, 2, 2] -
                          (meshXYZ[loopCHECKFACET, 2, 2] - meshXYZ[loopCHECKFACET, 2, 3]) *
                          (meshXYZ[loopCHECKFACET, 1, 2] - testx[loop]) / denom

            cond1 = (Y1predicted > meshXYZ[loopCHECKFACET, 2, 1] && YRpredicted > testy[loop]) ||
                    (Y1predicted < meshXYZ[loopCHECKFACET, 2, 1] && YRpredicted < testy[loop]) ||
                    (meshXYZ[loopCHECKFACET, 2, 2] - meshXYZ[loopCHECKFACET, 2, 3]) *
                    (meshXYZ[loopCHECKFACET, 1, 2] - testx[loop]) == 0
            cond1 || continue

            # Edge 3-1, opposing vertex 2
            denom = meshXYZ[loopCHECKFACET, 1, 3] - meshXYZ[loopCHECKFACET, 1, 1]
            Y2predicted = meshXYZ[loopCHECKFACET, 2, 3] -
                          (meshXYZ[loopCHECKFACET, 2, 3] - meshXYZ[loopCHECKFACET, 2, 1]) *
                          (meshXYZ[loopCHECKFACET, 1, 3] - meshXYZ[loopCHECKFACET, 1, 2]) / denom
            YRpredicted = meshXYZ[loopCHECKFACET, 2, 3] -
                          (meshXYZ[loopCHECKFACET, 2, 3] - meshXYZ[loopCHECKFACET, 2, 1]) *
                          (meshXYZ[loopCHECKFACET, 1, 3] - testx[loop]) / denom

            cond2 = (Y2predicted > meshXYZ[loopCHECKFACET, 2, 2] && YRpredicted > testy[loop]) ||
                    (Y2predicted < meshXYZ[loopCHECKFACET, 2, 2] && YRpredicted < testy[loop]) ||
                    (meshXYZ[loopCHECKFACET, 2, 3] - meshXYZ[loopCHECKFACET, 2, 1]) *
                    (meshXYZ[loopCHECKFACET, 1, 3] - testx[loop]) == 0
            cond2 || continue

            # Edge 1-2, opposing vertex 3
            denom = meshXYZ[loopCHECKFACET, 1, 1] - meshXYZ[loopCHECKFACET, 1, 2]
            Y3predicted = meshXYZ[loopCHECKFACET, 2, 1] -
                          (meshXYZ[loopCHECKFACET, 2, 1] - meshXYZ[loopCHECKFACET, 2, 2]) *
                          (meshXYZ[loopCHECKFACET, 1, 1] - meshXYZ[loopCHECKFACET, 1, 3]) / denom
            YRpredicted = meshXYZ[loopCHECKFACET, 2, 1] -
                          (meshXYZ[loopCHECKFACET, 2, 1] - meshXYZ[loopCHECKFACET, 2, 2]) *
                          (meshXYZ[loopCHECKFACET, 1, 1] - testx[loop]) / denom

            cond3 = (Y3predicted > meshXYZ[loopCHECKFACET, 2, 3] && YRpredicted > testy[loop]) ||
                    (Y3predicted < meshXYZ[loopCHECKFACET, 2, 3] && YRpredicted < testy[loop]) ||
                    (meshXYZ[loopCHECKFACET, 2, 1] - meshXYZ[loopCHECKFACET, 2, 2]) *
                    (meshXYZ[loopCHECKFACET, 1, 1] - testx[loop]) == 0
            cond3 || continue

            nf += 1
            if nf > length(facetCROSSLIST)
                append!(facetCROSSLIST, zeros(Int, 1000))
            end
            facetCROSSLIST[nf] = loopCHECKFACET
        end

        nf == 0 && continue
        activeFacets = facetCROSSLIST[1:nf]

        # Step 3: Find z coordinate where ray crosses each facet (plane equation)
        gridCOzCROSS = zeros(Float64, nf)
        for (idx, loopFINDZ) in enumerate(activeFacets)
            planecoA = meshXYZ[loopFINDZ,2,1]*(meshXYZ[loopFINDZ,3,2]-meshXYZ[loopFINDZ,3,3]) +
                       meshXYZ[loopFINDZ,2,2]*(meshXYZ[loopFINDZ,3,3]-meshXYZ[loopFINDZ,3,1]) +
                       meshXYZ[loopFINDZ,2,3]*(meshXYZ[loopFINDZ,3,1]-meshXYZ[loopFINDZ,3,2])
            planecoB = meshXYZ[loopFINDZ,3,1]*(meshXYZ[loopFINDZ,1,2]-meshXYZ[loopFINDZ,1,3]) +
                       meshXYZ[loopFINDZ,3,2]*(meshXYZ[loopFINDZ,1,3]-meshXYZ[loopFINDZ,1,1]) +
                       meshXYZ[loopFINDZ,3,3]*(meshXYZ[loopFINDZ,1,1]-meshXYZ[loopFINDZ,1,2])
            planecoC = meshXYZ[loopFINDZ,1,1]*(meshXYZ[loopFINDZ,2,2]-meshXYZ[loopFINDZ,2,3]) +
                       meshXYZ[loopFINDZ,1,2]*(meshXYZ[loopFINDZ,2,3]-meshXYZ[loopFINDZ,2,1]) +
                       meshXYZ[loopFINDZ,1,3]*(meshXYZ[loopFINDZ,2,1]-meshXYZ[loopFINDZ,2,2])
            planecoD = -meshXYZ[loopFINDZ,1,1]*(meshXYZ[loopFINDZ,2,2]*meshXYZ[loopFINDZ,3,3]-meshXYZ[loopFINDZ,2,3]*meshXYZ[loopFINDZ,3,2]) -
                        meshXYZ[loopFINDZ,1,2]*(meshXYZ[loopFINDZ,2,3]*meshXYZ[loopFINDZ,3,1]-meshXYZ[loopFINDZ,2,1]*meshXYZ[loopFINDZ,3,3]) -
                        meshXYZ[loopFINDZ,1,3]*(meshXYZ[loopFINDZ,2,1]*meshXYZ[loopFINDZ,3,2]-meshXYZ[loopFINDZ,2,2]*meshXYZ[loopFINDZ,3,1])

            abs(planecoC) < 1e-14 && (planecoC = 0.0)

            gridCOzCROSS[idx] = (-planecoD - planecoA * testx[loop] - planecoB * testy[loop]) / planecoC
        end

        # Remove z crossings outside mesh z bounds
        gridCOzCROSS = gridCOzCROSS[gridCOzCROSS .>= meshZmin - 1e-12 .&& gridCOzCROSS .<= meshZmax + 1e-12]
        isempty(gridCOzCROSS) && continue

        # Round and unique
        gridCOzCROSS = round.(gridCOzCROSS .* 1e10) ./ 1e10
        gridCOzCROSS = unique(sort(gridCOzCROSS))

        # Step 4: Check if testz is between a crossing pair (inside)
        n_cross = length(gridCOzCROSS)
        if iseven(n_cross)
            for loopASSIGN in 1:(n_cross ÷ 2)
                if testz[loop] > gridCOzCROSS[2*loopASSIGN-1] && testz[loop] < gridCOzCROSS[2*loopASSIGN]
                    OUTPUT[loop] = true
                    break
                end
            end
        elseif n_cross != 0
            push!(correctionLIST, loop)
        end
    end

    return OUTPUT, correctionLIST
end

"""
    intriangulation(vertices, faces, testp, heavytest)
This function checks whether test points are inside or outside a triangular surface mesh.

It takes:

* `vertices`: mesh vertex coordinates
* `faces`: triangle connectivity
* `testp`: points to test
* `heavytest`: optional extra randomized tests for difficult cases

The function builds triangle coordinate data, then uses `VOXELISEinternal` to test point inclusion by ray casting.

It first tests in the z-direction. If some points are ambiguous, it repeats the test in the x- and y-directions.

If `heavytest` is greater than 0, the mesh and points are randomly rotated and tested again to reduce numerical problems.

The output is:

* `1` = inside
* `0` = outside
* `-1` = unresolved / ambiguous
""" 
function intriangulation(vertices, faces, testp, heavytest::Int=0)
    # vertices: (nv, 3)
    # faces:    (nf, 3) - 1-indexed
    # testp:    (np, 3)
    # heavytest: number of additional randomized rotation tests

    np = size(testp, 1)
    inreturn = zeros(Int, np)
    VER   = copy(vertices)
    TESTP = copy(testp)

    for n in 1:(heavytest + 1)

        # Randomize rotation for heavytest iterations
        if n > 1
            v = rand(3)
            v = v / norm(v)
            D = rotmatrix(v, rand() * 180 / π)
            vertices = VER * D
            testp    = TESTP * D
        else
            vertices = VER
            testp    = TESTP
        end

        # Build meshXYZ: (nfaces, 3, 3) -> dim1=face, dim2=xyz, dim3=vertex
        nf = size(faces, 1)
        meshXYZ = zeros(Float64, nf, 3, 3)
        for v in 1:3
            for fi in 1:nf
                vi = faces[fi, v]
                meshXYZ[fi, :, v] = vertices[vi, :]
            end
        end

        # z-direction test
        in1, cl = VOXELISEinternal(testp[:, 1], testp[:, 2], testp[:, 3], meshXYZ)
        in_result = Int.(in1)  # false->0, true->1

        # x-direction test (only for unresolved points in cl)
        if !isempty(cl)
            meshXYZ_x = meshXYZ[:, [2, 3, 1], :]
            in2, cl2 = VOXELISEinternal(testp[cl, 2], testp[cl, 3], testp[cl, 1], meshXYZ_x)
            in_result[cl[in2]] .= 1
            cl = cl[cl2]
        end

        # y-direction test (only for still unresolved points)
        if !isempty(cl)
            meshXYZ_y = meshXYZ[:, [3, 1, 2], :]
            in3, cl3 = VOXELISEinternal(testp[cl, 3], testp[cl, 1], testp[cl, 2], meshXYZ_y)
            in_result[cl[in3]] .= 1
            cl = cl[cl3]
        end

        # Mark unresolvable points as -1
        in_result[cl] .= -1

        if n == 1
            inreturn = in_result
        else
            # If AT LEAST ONCE inside, mark as inside
            for i in 1:np
                if inreturn[i] != in_result[i] && in_result[i] == 1
                    inreturn[i] = 1
                end
            end
        end

    end

    return inreturn
end