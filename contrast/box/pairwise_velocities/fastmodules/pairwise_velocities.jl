using CellListMap
using StaticArrays
using LinearAlgebra


function _pdf_radial_velocity_r!(x, y, i, j, d2, velocities1, velocities2, rbins, counts)
    d = sqrt.(d2)
    ibin = searchsortedfirst(rbins, d) - 1
    if ibin > 0
        r = y - x
        dv = velocities2[:, j] - velocities1[:, i] 
        v_r = LinearAlgebra.dot(dv, r) / d

        counts[1][i, ibin] += 1
        counts[2][i, ibin] += v_r
    end
    return counts
end 


function _mean_radial_velocity_r!(x, y, i, j, d2, velocities1, velocities2, rbins, counts)
    d = sqrt.(d2)
    ibin = searchsortedfirst(rbins, d) - 1
    if ibin > 0
        r = y - x
        dv = velocities2[:, j] - velocities1[:, i] 
        v_r = LinearAlgebra.dot(dv, r) / d

        counts[1][ibin] += 1
        counts[2][ibin] += v_r
    end
    return counts
end 


function mean_radial_velocity_r(
    positions1, positions2, velocities1, velocities2, box_size, rbins
)
    positions1 = convert(Array{Float64}, positions1)
    positions2 = convert(Array{Float64}, positions2)
    velocities1 = convert(Array{Float64}, velocities1)
    velocities2 = convert(Array{Float64}, velocities2)

    rbins = convert(Array{Float64}, rbins)
    counts = (
        zeros(Int, length(rbins) - 1),
        zeros(Float64, length(rbins) - 1),
    )

    rmax = maximum(rbins)
    Lbox = [box_size, box_size, box_size]
    box = Box(Lbox, rmax)

    cl = CellList(positions1, positions2, box)

    counts = map_pairwise!(
        (x, y, i, j, d2, counts) ->
        _mean_radial_velocity_r!(x, y, i, j, d2, velocities1, velocities2, rbins, counts),
        counts, box, cl,
        parallel=true,
        reduce=reduce_hist
    )
    D1D2 = counts[1]
    V1V2 = counts[2]

    return D1D2, V1V2
end


function pdf_radial_velocity_r(
    positions1, positions2, velocities1, velocities2, box_size, rbins
)
    positions1 = convert(Array{Float64}, positions1)
    positions2 = convert(Array{Float64}, positions2)
    velocities1 = convert(Array{Float64}, velocities1)
    velocities2 = convert(Array{Float64}, velocities2)
    npos1 = size(positions1)[2]

    rbins = convert(Array{Float64}, rbins)
    counts = (
        zeros(Int, npos1, length(rbins) - 1),
        zeros(Float64, npos1, length(rbins) - 1),
    )

    rmax = maximum(rbins)
    Lbox = [box_size, box_size, box_size]
    box = Box(Lbox, rmax)

    cl = CellList(positions1, positions2, box)

    counts = map_pairwise!(
        (x, y, i, j, d2, counts) ->
        _pdf_radial_velocity_r!(x, y, i, j, d2, velocities1, velocities2, rbins, counts),
        counts, box, cl,
        parallel=true,
        reduce=reduce_hist
    )
    D1D2 = counts[1]
    V1V2 = counts[2]

    return D1D2, V1V2
end


function reduce_hist(hist, hist_threaded)
    hist = hist_threaded[1]
    for i in 2:length(hist_threaded)
        hist[1] .+= hist_threaded[i][1]
        hist[2] .+= hist_threaded[i][2]
    end
    return hist
end
