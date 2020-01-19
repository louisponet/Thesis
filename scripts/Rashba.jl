#%%
using DrWatson
quickactivate(@__DIR__)
using Revise
using DFWannier
const DFW = DFWannier
using LinearAlgebra
cd(datadir("Rashba/GeTe"))
#%%
job = DFJob("NSOC")

wgrid = DFW.read_points_from_xsf(joinpath(job, "wan_00001.xsf"))
wfuncs = [DFW.WannierFunction(joinpath(job, "wan_0000$i.xsf"), wgrid) for i = 1:8]


L_up = DFW.WannierFunction(wgrid, wfuncs[6]+1.0im*wfuncs[8])


#all shifting is done in the negative direction
function find_shift_indices(wfunction, cell)
    origin = wfunction.points[1,1,1]
    id1 = Tuple(findmin(map(x -> norm(x - (origin .+ ustrip.(cell[:,1]))), wfunction.points))[2]).-1
    id2 = Tuple(findmin(map(x -> norm(x - (origin .+ ustrip.(cell[:,2]))), wfunction.points))[2]).-1
    id3 = Tuple(findmin(map(x -> norm(x - (origin .+ ustrip.(cell[:,3]))), wfunction.points))[2]).-1
    return (id1, id2, id3)
end


function create_shifted_neg(wfunction, axis_id, cell)
    origin  = wfunction.points[1,1,1]
    id      = find_shift_indices(wfunction, cell)[axis_id]
    shift   = ustrip.(cell[:,axis_id])
    shifted = zeros(wfunction)

    shift_norm = norm(shift)
    for i in CartesianIndices(wfunction.points)
        if i[axis_id] <= id[axis_id]
            shifted.values[(Tuple(i) .+ 2 .* id)...] = wfunction.values[i]
        else
            shifted.values[(Tuple(i) .- id)...] = wfunction.values[i]
        end
    end
    return shifted
end

function create_shifted_pos(wfunction, axis_id, cell)
    origin  = wfunction.points[end,end,end]
    id      = find_shift_indices(wfunction, cell)[axis_id]
    shift   = ustrip.(cell[:,axis_id])
    shifted = similar(wfunction)

    shift_norm = norm(shift)
    for i in CartesianIndices(wfunction.points)
        if i[axis_id] > 2*id[axis_id]
            shifted.values[(Tuple(i) .- 2 .* id)...] = wfunction.values[i]
        else
            shifted.values[(Tuple(i) .+ id)...] = wfunction.values[i]
        end
    end
    return shifted
end


shifted_neg = create_shifted_neg(L_up, 2, cell(job))
sum(shifted_neg.values .- )
sum(L_up.values)
dot(L_up, L_up)
dot(shifted,L_up)
shifted_pos = create_shifted_pos(L_up, 2, cell(job))
sum(shifted_pos)
DFW.write_xsf_file(joinpath(job, "testwan2.xsf"), shifted_pos, job.structure)
DFW.calc_dip(shifted, L_up)

norm(DFW.calc_dip(shifted_neg, shifted_neg) .- DFW.calc_dip(L_up, L_up))
norm(cell(job)[:,2])
DFW.calc_dip(shifted_pos, L_up) + DFW.calc_dip(shifted_neg, L_up)
