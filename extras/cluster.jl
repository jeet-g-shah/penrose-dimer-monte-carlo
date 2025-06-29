# Comment
using DelimitedFiles
using Plots
using Colors
using Random
using Statistics
using LinearAlgebra
#using LsqFit
using Graphs
using JLD
# using CurveFit
using GLM
using DataFrames
using DataStructures
include("common_functions.jl")
include("common_functions2.jl")
default( framestyle = :box, aspect_ratio=:auto)


NUM_SUBDIVISIONS = 7
ngen = NUM_SUBDIVISIONS

# Read the tiling
const ROT = exp(1im * 2 * pi/5)
file_name = "./data/center_"
gph = open(readdlm,file_name*string(ngen)*"_graph.txt")
graph = Array{Int}(undef,size(gph)[1], size(gph)[2]) 
for i in 1:size(gph)[1]
    for j in 1:size(gph)[2]
        graph[i,j] = Int(gph[i,j])
    end
end

pts_x = vec(open(readdlm,file_name*string(ngen)*"_pts_re.txt"))
pts_y = vec(open(readdlm,file_name*string(ngen)*"_pts_im.txt"))
pts = [complex(pts_x[i], pts_y[i]) for i in 1:length(pts_x)];

ind_A = round.(Int, vec(open(readdlm,file_name*string(ngen)*"_ind_A.txt")))
ind_B = round.(Int, vec(open(readdlm,file_name*string(ngen)*"_ind_B.txt")))

AtoG = round.(Int, vec(open(readdlm,file_name*string(ngen)*"_AtoG.txt")))
BtoG = round.(Int, vec(open(readdlm,file_name*string(ngen)*"_BtoG.txt")))

max_match = round.(Int, open(readdlm,file_name*string(ngen)*"_max_match.txt"))
plaq = round.(Int, open(readdlm,file_name*string(ngen)*"_plaq.txt"))
plaq_spins = round.(Int, open(readdlm,file_name*string(ngen)*"_plaq_spins.txt"))
# plaq_spins has the format [s1,s2,s3,s4] where 
# s1, s2 share a vertex
# s1, s3 share a vertex
# s4, s2 share a vertex
# s4, s3 share a vertex
# going cyclically, the order or spins is - s1 -- s2 -- s4 -- s3 -

ABedge_no = round.(Int, open(readdlm,file_name*string(ngen)*"_ABedge_no.txt"))
spins_to_edges = round.(Int, open(readdlm,file_name*string(ngen)*"_spins_to_edges.txt"))

monomer_A = vec(round.(Int, open(readdlm,file_name*string(ngen)*"_monomer_A.txt")))
monomer_B = vec(round.(Int, open(readdlm,file_name*string(ngen)*"_monomer_B.txt")))

const nspins = size(spins_to_edges)[1]
const npts = size(graph)[1];
const nplaqs = length(plaq_spins[:,1])
const ymax = maximum(pts_y)


nbrG = []
for i in 1:npts
    nnn = []
    for j in 1:npts
        if graph[i,j] == 1
            push!(nnn, j)
        end
    end
    push!(nbrG, nnn)
end
nbrG = convert_to_int(nbrG)

Gedge_no = fill(0,npts, npts);
for i in 1:npts
    for j in 1:npts
        if graph[i,j] == 1
            if ind_A[i] < npts+5
                Gedge_no[i,j] = ABedge_no[ind_A[i],ind_B[j]]
            elseif ind_B[i] < npts+5
                Gedge_no[i,j] = ABedge_no[ind_A[j],ind_B[i]]
            else
                println("Error!")
            end
        end
    end
end

spins_pos = generate_spins_pos()
config = reset_config(max_match);


# Generate Graph objects
penrose = Graph(pts, graph, npts, plaq, nplaqs)
spins = Graph(spins_pos, generate_spins_adj(plaq_spins, nspins), nspins, plaq_spins, nplaqs)
dcov = Graph(pts, generate_max_match_adj(max_match, penrose), npts, fill(0,1,1), 0);

mc_steps = 5*1000*100
Nmix = 100
spin_origin = 300 # 4500 # 651
config = reset_config(max_match)
@time sz, szsz = mc_immobile_no_avg(mc_steps, Nmix, config ,  plaq_spins, spin_origin);

mc_steps = 2*10*1000
Nmix = 100
spin_origin = 300
config = reset_config(max_match)
@time sz_disg, szsz_disg = mc_mobile_full_no_avg(mc_steps, Nmix, config ,  plaq_spins, spin_origin, Gedge_no, convert_to_int(nbrG));

monomer_membrane = get_monomer_membranes_numerics(sz);

disg = generate_disconnected_graph(sz_disg)

g = SimpleGraph(disg.adj)
ccg = connected_components(g)
;

# # p = plot_generate(1000)
# # size = 1000;
# legend = false;
# ylim = (-ymax-0.1,ymax+0.1);
# p = plot(size = (1000,1000), aspect_ratio = :equal , legend = legend, ylim = ylim, axis=([], false))
# p = plot_edges(p, penrose)
# # p = plot_sites(p, spins, :blue, 0; annotate = false, markersize=  2)
# p = plot_membrane(p,get_monomer_membranes_numerics(sz_disg) )
# fs = 35
# p = annotate!(7, 2, text("\$ \\mathcal{R}_1 \$",  fs))
# p = annotate!(10,-4, text("\$ \\mathcal{R}_2 \$",  fs))
# p = annotate!(28,0, text("\$ \\mathcal{R}_3 \$",  fs))

# # p = plot!([(7,2)];texts=["\$ \\mathcal{R}_1 \$"])
# # p = plot!([ (10,-4)];texts=["\$ \\mathcal{R}_2 \$"])
# # p = plot!([(28,0)];texts=["\$ \\mathcal{R}_3 \$"])
# # display(p)
# savefig(p,"./images/temp/monomer-membranes.pdf")

# ;

disg = generate_disconnected_graph(sz_disg)

p = plot_generate(1000)
p = plot_edges(p, disg)
# p = plot_sites(p, spins, :blue, 0; annotate = false, markersize=  2)
for i in 1:length(ccg)
    jj = -1
    for k in 1:disg.n
        if disg.adj[ccg[i][1],k] == 1
            jj = k
            break
        end
    end
    spin_number = Gedge_no[ccg[i][1],jj]
    # spin_number = Gedge_no[ccg[i][1], nbrG[ccg[i][1]][1]]
    point = spins_pos[spin_number]
    scatter!(p, [real(point)], [imag(point)], mc=:pink, markerstrokewidth=0, series_annotations = string(spin_number), label = "",  markersize=6)
end
# display(p)
;

# SETTING SPIN ORIGIN HERE -------------------------------------------------
# ##################################################################################################################################################
# Different regions for ngen=8: 921, 707, 1
# spin_origin = 921
spin_origin = 1
# ##################################################################################################################################################

p = plot_generate(500)
plot_edges(p,penrose)
# plot_sites(p,spins, "orange",0, annotate=false)
#config = reset_config(max_match)
#plot_config(p,config)
plot_single_point(p,spins_pos[spin_origin])
# plot_membrane(p,get_monomer_membranes_numerics(sz))
savefig("./images/temp/monomer_membrane.pdf")
;


sub_lat = 0
a, b = spins_to_edges[spin_origin,:]
p1 = AtoG[a]
p2 = BtoG[b]
for i in 1:length(ccg)
    if p1 in ccg[i] && p2 in ccg[i]
        global sub_lat = ccg[i]
        break
    end
end

sub_npts = length(sub_lat)
sub_pts = [pts[sub_lat[i]] for i in 1:length(sub_lat)]

lat_to_sublat = fill(0, npts)
for i in 1:sub_npts
    lat_to_sublat[sub_lat[i]] = i
end

sub_spins = Array{Int64}(undef, 0)
sub_graph = fill(0, (sub_npts, sub_npts)) 

sub_spins_to_full_spins = fill(0, nspins)
sub_spins_lat_mat = fill(0,sub_npts, sub_npts)

for i in 1:sub_npts-1
    for j in i+1:sub_npts
        if graph[sub_lat[i], sub_lat[j]] == 1
            sub_graph[i,j] = 1
            sub_graph[j,i] = 1
            push!(sub_spins, Gedge_no[sub_lat[i], sub_lat[j]])
            sub_spins_to_full_spins[Gedge_no[sub_lat[i], sub_lat[j]]] = length(sub_spins)
            sub_spins_lat_mat[i,j] = length(sub_spins)
            sub_spins_lat_mat[j,i] = length(sub_spins)
        end
    end
end

sub_spins_to_edges = fill(0,2, length(sub_spins))
count = 1
for i in 1:sub_npts-1
    for j in i+1:sub_npts
        if graph[sub_lat[i], sub_lat[j]] == 1
            sub_spins_to_edges[1, count] = i
            sub_spins_to_edges[2, count] = j
            global count += 1
        end
    end
end

sub_spins_pos = fill(0.0+0.0im, length(sub_spins))
for i in 1:length(sub_spins)
    sub_spins_pos[i] = spins_pos[sub_spins[i]]
end

sub_plaq = [] # These are spin plaquettes
for i in 1:size(plaq_spins)[1]
    s1, s2, s3, s4 = plaq_spins[i,:]
    if s1 in sub_spins && s2 in sub_spins && s3 in sub_spins && s4 in sub_spins
        push!(sub_plaq, [sub_spins_to_full_spins[s1], sub_spins_to_full_spins[s2], sub_spins_to_full_spins[s3], sub_spins_to_full_spins[s4]])
    end
end
sub_plaq_temp = fill(0,length(sub_plaq),4)
for i in 1:length(sub_plaq)
    for j in 1:4
        sub_plaq_temp[i,j] = sub_plaq[i][j]
    end
end
sub_plaq = sub_plaq_temp
sub_nspins = length(sub_spins)
sub_nbr = []
for i in 1:sub_npts
    temp = []
    for j in 1:sub_npts
        if sub_graph[i, j] == 1
            push!(temp, j)
        end
    end
    push!(sub_nbr, temp)
end

sub_nplaq = size(sub_plaq)[1]
sub_spins_graph2 = fill(0,sub_nspins, sub_nspins);
for i in 1:sub_nplaq
    s1, s2, s3, s4 = sub_plaq[i,:]
    sub_spins_graph2[s1, s3] = 1
    sub_spins_graph2[s1, s2] = 1
    sub_spins_graph2[s3, s4] = 1
    sub_spins_graph2[s4, s2] = 1
    
    sub_spins_graph2[s3, s1] = 1
    sub_spins_graph2[s2, s1] = 1
    sub_spins_graph2[s4, s3] = 1
    sub_spins_graph2[s2, s4] = 1
end

sub_spins_graph = fill(0, sub_nspins, sub_nspins)
for i in 1:nspins
    if sub_spins_to_full_spins[i] > 0
        for j in 1:nspins
            if spins.adj[i,j] == 1 && sub_spins_to_full_spins[j] > 0
                sub_spins_graph[sub_spins_to_full_spins[i] , sub_spins_to_full_spins[j]] = 1
                sub_spins_graph[sub_spins_to_full_spins[j] , sub_spins_to_full_spins[i]] = 1
            end
        end
    end
end

sub_penrose = Graph(sub_pts, sub_graph, sub_npts, fill(0,1,1), length(sub_plaq[:,1]))
sub_penrose_spin = Graph(sub_spins_pos, sub_spins_graph,length(sub_spins_pos), sub_plaq , length(sub_plaq[:,1]))

sub_lat_nbr_spin = []
sub_penrose_nbr = get_sub_neighbors(sub_penrose)
for i in 1:sub_penrose.n
    temp = []
    for neighbor in sub_penrose_nbr[i]
        append!(temp, sub_spins_lat_mat[i,neighbor])
    end
    push!(sub_lat_nbr_spin, temp)
end
sub_lat_nbr_spin = convert_to_int(sub_lat_nbr_spin)
;

p = plot_generate(1000)
plot_edges(p,sub_penrose)
plot_sites(p,sub_penrose, "orange",0, annotate=false)
display(p)
;

spin_to_plaq = []
 
for i in 1:sub_penrose_spin.n
    temp = []
    for j in 1:size(sub_penrose_spin.plaq)[1]
        flag =  false
        if i in sub_penrose_spin.plaq[j,:]
            push!(temp, j[1])
        end
    end
    push!(spin_to_plaq, temp)
end

spin_to_plaq = convert_to_int(spin_to_plaq);
cplaq = length.(spin_to_plaq);

sub_spins_graph_nbr = []
for i in 1:sub_nspins
    temp = []
    for j in 1:sub_nspins
        if sub_spins_graph[i,j] == 1
            push!(temp, j)
        end
    end
    push!(sub_spins_graph_nbr, temp)
end

vison_spins_graph_nbr = []
for i in 1:sub_nspins
    temp = []
    pq = spin_to_plaq[i]
    for pl in pq
        temp = vcat(temp, sub_plaq[pl,:])
    end
    temp = unique(temp)
    filter!(e->e≠i,temp)
    push!(vison_spins_graph_nbr, temp)
end
for i in 1:length(sub_spins_graph_nbr)
    for xx in sub_spins_graph_nbr[i]
        if !( xx ∈ vison_spins_graph_nbr[i])
            push!(vison_spins_graph_nbr[i], xx)
        end
    end
end
vison_spins_graph_nbr = convert_to_int(convert_to_int(vison_spins_graph_nbr))
;

spin_cno = sum(sub_penrose_spin.adj, dims=1)
# spin_nbr = get_sub_neighbors(sub_penrose_spin)

spins_to_plaqs = [[] for i in 1:sub_nspins]
for i in 1:sub_nplaq    
    for j in sub_penrose_spin.plaq[i,:]
        push!(spins_to_plaqs[j], i)
    end
end

spins_to_plaqs = [[] for i in 1:sub_nspins]
for i in 1:sub_nplaq    
    for j in sub_penrose_spin.plaq[i,:]
        push!(spins_to_plaqs[j], i)
    end
end
spins_to_plaqs = convert_to_int(spins_to_plaqs)

## Finding Vison paths

npaths = 50 # 200
lpath =  40
paths = Vector{Int64}[]
for i in 1:npaths
    pp = []
    while true
        pp = get_good_path_maxl(rand(1:sub_penrose_spin.n), lpath,spins_to_plaqs, sub_nplaq,vison_spins_graph_nbr,3*lpath, 1000)# 100,500)
        if length(pp) > 2
            break
        end
    end
    push!(paths, pp)
end
# paths = paths
# paths = convert_to_int(convert_to_int(paths))

l_paths = [length(paths[i]) for i in 1:length(paths)]
# println(l_paths)
mega_path::Vector{Matrix{Int64}} = [generate_symm_path(pth, sub_penrose_spin.pts) for pth in paths]
;

# #plot paths
# for j in 1:5
#     # p = plot_generate(1000)
#     siz = 6000
#     p = plot(size = (siz,siz) , aspect_ratio = :equal, legend = false, ylim = (-50,50), xlim = (-50,50))
#     plot_edges(p,sub_penrose)
#     # plot_edges(p,sub_penrose_spin, "purple", 2)
#     plot_sites(p,sub_penrose_spin,:blue,6; annotate=false, markersize=0 )
#     # i = rand(1:npaths)
#     # sub_plot_path(p,mega_path[i][:,2])
#     # sub_plot_path(p,custom_path)
#     sub_plot_path(p,mega_path[j][:,3])
#     # sub_plot_path(p,mega_path[1][:,5])
#     savefig(p,"./images/del/foo_vison"*string(j)*".png")
#     # display(p)
# end



# p = plot_generate(1000)
# p = plot_edges(p,sub_penrose, "pink",2)

# p = plot_sites(p, sub_penrose,:red,0;annotate=true)
# display(p)

cno = sum(sub_penrose.adj;dims=1)[1,:];

lat_path_to_spin_path(ppp) = [sub_spins_lat_mat[ppp[i],ppp[i%length(ppp) + 1]] for i in 1:length(ppp)]

sub_spin_nbr = get_spin_nbr(sub_penrose, sub_spins_lat_mat);




my_path = generate_monomer_path(100,get_nbrG(sub_penrose), 200 )
# p = plot_generate(700)
siz = 2000
p = plot(size = (siz,siz), aspect_ratio = :equal , legend = false) #, xlim=(-25,-15), ylim = (-20, -15))
# plot_path(p, bs[4][1], sub_penrose )
plot_path(p, my_path, sub_penrose )
plot_edges(p, sub_penrose)
# plot_single_point(p,sub_penrose.pts[1014])
# plot_single_point(p,sub_penrose.pts[1009])
# plot_single_point(p,sub_penrose.pts[1013])
# savefig(p,"./images/del/foo1.png")
;
 display(p)

println("Vison paths generated")


abasis = generate_mbasis()
length.(abasis)

bulkbasis = []
for i in 1:length(abasis)
    temp = []
    for j in 1:length(abasis[i])
        if boundary_path_check(abasis[i][j], sub_penrose_spin) == false
            push!(bulkbasis, abasis[i][j])
        end
    end
    # push!(bulkbasis, temp)
    # push!(bulkbasis, temp)
    # push!(bulkbasis, convert_to_int(temp))
end

 # save("./data/paper/mfm_R2_path20.jld", "path20", convert_to_int(bulkbasis[1]))
 # save("./data/monomer_paths_$ngen $spin_origin.jld", "bulkbasis", bulkbasis);
 # bulkbasis = load("./data/monomer_paths_$ngen $spin_origin.jld")["bulkbasis"]
# save("./data/R3_bulkbasis.jld", "bulkbasis", bulkbasis)

mfmnum_basis = convert_to_int(bulkbasis)
# temp_mfmnum_basis = copy(mfmnum_basis)
# mfmnum_basis = []
# for i in 1:length(temp_mfmnum_basis)
#     for j in 0:2:10
#         push!(mfmnum_basis, shift(temp_mfmnum_basis[i],j))
#     end
# end
# mfmnum_basis = [
#     shift(mfmnum_basis[end],i) for i in 1:10
# ]
mfmnum_basis1 = [split_loop(lp)[1] for lp in mfmnum_basis]
mfmnum_basis2 = [split_loop(lp)[2] for lp in mfmnum_basis]
for i in 1:length(mfmnum_basis1)
    if isodd( length(mfmnum_basis1[i]) +length(mfmnum_basis2[i]) )
        println("Loop length not even")
    end
end
mfmnum1::Vector{Matrix{Int64}} = basis_path_to_all_path(mfmnum_basis1, sub_penrose_spin.pts)
mfmnum2::Vector{Matrix{Int64}} = basis_path_to_all_path(mfmnum_basis2, sub_penrose_spin.pts)
npath = length(mfmnum1)
lpath1 = [size(mfmnum1[i])[1] for i in 1:npath]
lpath2 = [size(mfmnum2[i])[1] for i in 1:npath]
# println([lpath1 lpath2])
mfmpath = [mfmnum1, mfmnum2]

cap = fill(0,2, 10, npath)
for i in 1:npath
    for j in 1:10
        s1 = mfmnum1[i][1,j]
        s2 = mfmnum2[i][end,j]
        temp = intersect( sub_spins_to_edges[:,s1], sub_spins_to_edges[:,s2])
        if length(temp) == 1
            cap[1,j,i] = temp[1]
        else
            println("ERROR")
        end

        s1 = mfmnum1[i][end,j]
        s2 = mfmnum2[i][1,j]
        temp = intersect( sub_spins_to_edges[:,s1], sub_spins_to_edges[:,s2])
        if length(temp) == 1
            cap[2,j,i] = temp[1]
        else
            println("ERROR")
        end
    end
end

clean_mfmnum1 = []
clean_mfmnum2 = []
clean_cap = []
for i in 1:npath
    c1 = cap[1,1,i]
    c2 = cap[2,1,i]
    if (cno[c1] == 3 || cno[c1] == 5) && (cno[c2] == 3 || cno[c1] == 5)
        push!(clean_mfmnum1, mfmnum1[i])
        push!(clean_mfmnum2, mfmnum2[i])
        push!(clean_cap, cap[:,:,i])
    end
end
npath = length(clean_mfmnum1)
cap = fill(0,2, 10, npath)
for i in 1:npath
    for j in 1:10
        cap[2,j,i] = clean_cap[i][2,j]
        cap[1,j,i] = clean_cap[i][1,j]
    end
end
mfmnum1 = convert_to_int(clean_mfmnum1)
mfmnum2 = convert_to_int(clean_mfmnum2)
mfmpath = [mfmnum1, mfmnum2]
lpath1 = [size(mfmnum1[i])[1] for i in 1:npath]
lpath2 = [size(mfmnum2[i])[1] for i in 1:npath]  
; 

@show counter(lpath1)
println("MFM paths generated")

p = plot_generate(1000)
p = plot_edges(p,sub_penrose, "pink",2)

p = plot_sites(p,sub_penrose;markersize=3)
i = rand(1:8)
p = plot_path(p,mfmpath[1][i][:,1], sub_penrose_spin)
p = plot_path(p,mfmpath[2][i][:,1], sub_penrose_spin)

#vfmnum1::Vector{Matrix{Int64}} = basis_path_to_all_path(vfmnum_basis1, sub_penrose_spin.pts)
#vfmnum2::Vector{Matrix{Int64}} = basis_path_to_all_path(vfmnum_basis2, sub_penrose_spin.pts)
#npath = length(vfmnum1)
#lpath1 = [size(vfmnum1[i])[1] for i in 1:npath]
#lpath2 = [size(vfmnum2[i])[1] for i in 1:npath]
#vfmpath = [vfmnum1, vfmnum2]


# @time vbasis = generate_vbasis(vison_spins_graph_nbr, convert_to_int(spin_to_plaq), length.(spin_to_plaq))
#  save("./data/vfm_paths_$ngen $spin_origin.jld","vbasis",vbasis)

# vfmnum_basis1, vfmnum_basis2, vfmpath, lpath1, lpath2, npath = generate_vfnumbasis(vbasis, sub_penrose_spin)
# vfmnum1, vfmnum2 = vfmpath;

#  save("./data/paper/vfm_paths_$ngen $spin_origin.jld","vfmpath",vfmpath)
# vfmpath = load("./data/paper/vfm_paths_$ngen $spin_origin.jld")["vfmpath"]
# vfmnum1, vfmnum2 = vfmpath;


siz = 3000
p = plot(size = (siz,siz), aspect_ratio = :equal , legend = false)#, ylim = (-3, 20), xlim=(5,25))
# p = plot_generate(700)
plot_edges(p,sub_penrose)
plot_sites(p, sub_penrose_spin, :blue, 0; annotate = false, markersize=  1)
# plot_edges(p,sub_penrose_spin, "purple", 2)
j = 4
for i in length(vfmnum1)-4:length(vfmnum1)-3
    # sub_plot_path(p,vfmnum_basis1[i])
    sub_plot_path(p,vfmnum1[i][:,j])
    sub_plot_path(p,vfmnum2[i][:,j])
    # sub_plot_path(p,vfmnum_basis2[i])
end
# plot_sub_config(p,sub_config)
# display(p)
savefig(p,"./images/del/foovfm1.png")


function mc_mobile_sub_all(mc_steps::Int64, Nmix::Int64, config::Vector{Int64} ,  plaq_spins::Matrix{Int64}, npts::Int64, sub_spins_lat_mat::Matrix{Int64}, nbr::Vector{Vector{Int64}} , paths::Vector{Vector{Int64}}, mega_path::Vector{Matrix{Int64}}, vfmpath::Vector{Vector{Matrix{Int64}}}, mfmpath::Vector{Vector{Matrix{Int64}}}, cap::Array{Int64, 3}, sub_spin_nbr::Vector{Vector{Int64}})
    nplaq::Int64 = size(plaq_spins)[1]
    nspins::Int64 = length(config)
    sz_corr_array::Matrix{Int64} = fill(0,nspins, nspins)
    sz_avg::Vector{Int64} = fill(0,nspins)
    npaths::Int64 = length(paths)
    lpath::Vector{Int64} = length.(paths)
    vis_corr_array::Vector{Matrix{Int64}} = [fill(0,lpath[i], lpath[i]) for i in 1:npaths]
    vis_corr_array_2::Vector{Matrix{Int64}} = [ fill(0,size(mega_path[i])[1],10) for i in 1:length(mega_path) ]
    nmegapaths::Int64 = length(mega_path)
    lmegapath::Vector{Int64} = [size(mega_path[i])[1] for i in 1:nmegapaths]
    
    vfmnum1 = vfmpath[1]
    vfmnum2 = vfmpath[2]
    vfmnpath = length(vfmnum1)
    vfmlpath1 = [size(vfmnum1[i])[1] for i in 1:vfmnpath]
    vfmlpath2 = [size(vfmnum2[i])[1] for i in 1:vfmnpath]
    num1 = fill(0,10,vfmnpath)
    num2 = fill(0,10,vfmnpath)
    den = fill(0,10,vfmnpath)

    mfmnum1 = mfmpath[1]
    mfmnum2 = mfmpath[2]
    mnpath = length(mfmnum1)
    m_num1 = fill(0,10,mnpath)
    m_num2 = fill(0,10,mnpath)
    m_den = fill(0,10, mnpath)
    m_lpath1 = [size(mfmnum1[i])[1] for i in 1:mnpath]
    m_lpath2 = [size(mfmnum2[i])[1] for i in 1:mnpath]

    s1::Int64 = 0
    s2::Int64 = 0
    s3::Int64 = 0
    s4::Int64 = 0

    random_plaq::Int64 = 0

    for i in 1:mc_steps
        for j in 1:Nmix
            while !monomer_move(config, npts, sub_spins_lat_mat,nbr)
            end
            while !plaq_move(config, nplaq, plaq_spins )
            end
        end
        # sz_avg .= sz_avg .+ config
        # sz_correlator_mat(config, sz_corr_array, nspins)
        # vis_partial(config, vis_corr_array_2, nmegapaths,lmegapath,  mega_path)
        # vis_fm(config, num1, num2, den, vfmnum1, vfmnum2,vfmnpath, vfmlpath1, vfmlpath2 )
        mon_fm(config, m_num1, m_num2, m_den, mfmnum1, mfmnum2, mnpath, m_lpath1, m_lpath2, cap, sub_spin_nbr)
    end

    sz_corr_array_float= sz_corr_array/mc_steps
    sz_avg_float = sz_avg/mc_steps
    vis_corr_array_float = vis_corr_array_2/mc_steps
    vfm_arr_float = [num1, num2, den] ./ mc_steps
    mfm_arr_float = [m_num1, m_num2, m_den] ./ mc_steps
    return sz_avg_float, (sz_corr_array_float+sz_corr_array_float')+I, vis_corr_array_float, (vfm_arr_float), mfm_arr_float, config
end
;
mc_steps = 1000*1
Nmix = 500
sub_config = get_sub_config( reset_config(max_match) )
println("The main run")
@time sz1, szsz1, vison1, vfm_arr, mfm_arr, sub_config1 = mc_mobile_sub_all(mc_steps, Nmix, sub_config ,  sub_penrose_spin.plaq , sub_npts, sub_spins_lat_mat, get_sub_neighbors(sub_penrose), paths, mega_path, vfmpath, mfmpath, cap, convert_to_int(sub_spin_nbr));
# symmetrise_vison(vison1, length.(paths), length(paths))

data_to_save = [[sz1, szsz1, vison1, vfm_arr, mfm_arr], [mc_steps, Nmix, sub_config ,  sub_penrose_spin.plaq , sub_npts, sub_spins_lat_mat, get_sub_neighbors(sub_penrose), paths, mega_path, vfmpath, mfmpath, cap, sub_spin_nbr]]
# save("./data/paper/paper_mc_run_all_corr_R3.jld","data_to_save",data_to_save)

# loaded_data = load("./data/mc_run_all_corr.jld")["data_to_save"]
# sz1, szsz1, vison1, vfm_arr, mfm_arr = loaded_data[1]
# mc_steps, Nmix, sub_config ,  temp1 , sub_npts, sub_spins_lat_mat, temp2, paths, mega_path, vfmpath, mfmpath, cap, sub_spin_nbr = loaded_data[2];
;
