using CellListMap
using StaticArrays
using LinearAlgebra


function _count_pairs_r!(i, j, d2, weights1, weights2, rbins, counts)
    d = sqrt.(d2)
    ibin = searchsortedfirst(rbins, d) - 1
    if ibin > 0
        counts[ibin] += weights1[i] * weights2[j]
    end
    return counts
end 


function count_pairs_r(
    positions1, positions2, weights1, weights2, boxsize, rbins
)
    positions1 = convert(Array{Float64}, positions1)
    positions2 = convert(Array{Float64}, positions2)
    weights1 = convert(Array{Float64}, weights1)
    weights2 = convert(Array{Float64}, weights2)
    boxsize = convert(Array{Float64}, boxsize)
    rbins = convert(Array{Float64}, rbins)
    D1D2 = zeros(Int, length(rbins) - 1)
    
    rmax = maximum(rbins)
    Lbox = [boxsize[1], boxsize[2], boxsize[3]]
    box = Box(Lbox, rmax)

    cl = CellList(positions1, positions2, box)

    map_pairwise!(
        (x, y, i, j, d2, D1D2) ->
        _count_pairs_r!(i, j, d2, weights1, weights2, rbins, D1D2),
        D1D2, box, cl,
        parallel=true
    )

    return D1D2
end


function count_pairs_2d_r(
    positions1, positions2, boxsize, rbins
)
    positions1 = convert(Array{Float64}, positions1)
    positions2 = convert(Array{Float64}, positions2)
    boxsize = convert(Array{Float64}, boxsize)
    rbins = convert(Array{Float64}, rbins)
    D1D2 = zeros(Int, length(rbins) - 1)
    
    rmax = maximum(rbins)
    Lbox = [boxsize[1], boxsize[2]]
    box = Box(Lbox, rmax)

    cl = CellList(positions1, positions2, box)

    map_pairwise!(
        (x, y, i, j, d2, D1D2) ->
        _count_pairs_r!(d2, rbins, D1D2),
        D1D2, box, cl,
        parallel=true
    )
    return D1D2
end


function _count_pairs_rmu!(x, y, i, j, d2, rbins, mubins, com, counts)
    r = y - x
    d = sqrt.(d2)
    mu = LinearAlgebra.dot(r, com) / (d*sqrt(LinearAlgebra.dot(com, com)))
    irbin = searchsortedfirst(rbins, d) - 1
    imubin = searchsortedfirst(mubins, mu) - 1
    if irbin > 0
        counts[irbin, imubin] += 1
    end
    return counts
end 


function count_pairs_rmu(
    positions1, positions2, boxsize, rbins, mubins
)

    positions1 = convert(Array{Float64}, positions1)
    positions2 = convert(Array{Float64}, positions2)
    boxsize = convert(Array{Float64}, boxsize)
    rbins = convert(Array{Float64}, rbins)
    mubins = convert(Array{Float64}, mubins)
    com = [0, 0, 1]
    D1D2 = zeros(Int, length(rbins) - 1, length(mubins) - 1)
    
    rmax = maximum(rbins)
    Lbox = [boxsize[1], boxsize[2], boxsize[3]]

    box = Box(Lbox, rmax)

    cl = CellList(positions1, positions2, box)

    D1D2 = map_pairwise!(
        (x, y, i, j, d2, D1D2) ->
        _count_pairs_rmu!(x, y, i, j, d2, rbins, mubins, com, D1D2),
        D1D2, box, cl,
        parallel=true,
        reduce=reduce_matrix
    )
    return D1D2
end


function reduce_matrix(output, output_threaded)
    output = output_threaded[1]
    for i in 2:length(output_threaded)
        output .+= output_threaded[i]
    end
    return output
end
