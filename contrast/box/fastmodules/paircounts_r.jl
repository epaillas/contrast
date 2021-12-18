using CellListMap
using StaticArrays
using LinearAlgebra
using PyCall
using DelimitedFiles


function _count_pairs_r!(i, j, d2, rbins, counts)
    r = sqrt.(d2)
    ibin = searchsortedfirst(rbins, r) - 1
    if ibin > 0
        counts[ibin] += 1
    end
    return counts
end 

function count_pairs_r(
    positions1, positions2, box_size, rbins
)

    positions1 = convert(Array{Float64}, positions1)
    positions2 = convert(Array{Float64}, positions2)
    npos1 = size(positions1)[1]
    npos2 = size(positions2)[1]

    rbins = convert(Array{Float64}, rbins)
    D1D2 = zeros(Int, length(rbins) - 1)
    R1R2 = zeros(Int, length(rbins) - 1)
    
    rmax = maximum(rbins)
    Lbox = [box_size, box_size, box_size]
    box = Box(Lbox, rmax)

    cl = CellList(positions1, positions2, box)

    D1D2 = map_pairwise!(
        (x, y, i, j, d2, output) ->
        _count_pairs_r!(i, j, d2, rbins, D1D2),
        D1D2, box, cl,
        parallel=true
    )

    return D1D2

end

function count_pairs_2d_r(
    positions1, positions2, box_size, rbins
)

    positions1 = convert(Array{Float64}, positions1)
    positions2 = convert(Array{Float64}, positions2)
    npos1 = size(positions1)[1]
    npos2 = size(positions2)[1]

    rbins = convert(Array{Float64}, rbins)
    D1D2 = zeros(Int, length(rbins) - 1)
    R1R2 = zeros(Int, length(rbins) - 1)
    
    rmax = maximum(rbins)
    Lbox = [box_size, box_size]
    box = Box(Lbox, rmax)

    cl = CellList(positions1, positions2, box)

    D1D2 = map_pairwise!(
        (x, y, i, j, d2, output) ->
        _count_pairs_r!(i, j, d2, rbins, D1D2),
        D1D2, box, cl,
        parallel=true
    )

    return D1D2

end