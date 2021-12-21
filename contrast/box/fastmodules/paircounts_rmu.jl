using CellListMap
using StaticArrays
using LinearAlgebra
using PyCall
using DelimitedFiles


function _count_pairs_rmu!(x, y, i, j, d2, rbins, mubins, com, counts)
    r = x - y
    d = sqrt.(d2)
    mu = (r[1]*com[1] + r[2]*com[2] + r[3]*com[3]) / (d*sqrt(com[1]*com[1] + com[2]*com[2] + com[3]*com[3]))
    irbin = searchsortedfirst(rbins, d) - 1
    imubin = searchsortedfirst(mubins, mu) - 1
    if irbin > 0
        counts[irbin, imubin] += 1
    end
    return counts
end 

function count_pairs_rmu(
    positions1, positions2, box_size, rbins, mubins
)

    positions1 = convert(Array{Float64}, positions1)
    positions2 = convert(Array{Float64}, positions2)

    rbins = convert(Array{Float64}, rbins)
    mubins = convert(Array{Float64}, mubins)
    com = [0, 0, 1]
    D1D2 = zeros(Int, length(rbins) - 1, length(mubins) - 1)
    
    rmax = maximum(rbins)
    Lbox = [box_size, box_size, box_size]
    box = Box(Lbox, rmax)

    cl = CellList(positions1, positions2, box)

    D1D2 = map_pairwise!(
        (x, y, i, j, d2, output) ->
        _count_pairs_rmu!(x, y, i, j, d2, rbins, mubins, com, D1D2),
        D1D2, box, cl,
        parallel=true
    )

    return D1D2

end

