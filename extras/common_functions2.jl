function sz_correlator(config::Vector{Int64} , running_sum::Vector{Int64}, spin_origin::Int64)
    if config[spin_origin] == 1
        for i in eachindex(running_sum)
            running_sum[i] = running_sum[i] + config[i]
        end
    else
        for i in eachindex(running_sum)
            running_sum[i] = running_sum[i] - config[i]
        end
    end
end

function add_vec!(vec1::Vector{Int64}, vec2::Vector{Int64})
    for i in eachindex(vec1)
        vec1[i] = vec1[i] + vec2[i]
    end
end

function mc_immobile_no_avg(mc_steps::Int64, Nmix::Int64, config::Vector{Int64} ,  plaq_spins::Matrix{Int64}, spin_origin::Int64 )

    nplaq::Int64 = size(plaq_spins)[1]
    nspins::Int64 = length(config)
    szsz::Vector{Int64} = fill(0,nspins)
    sz_avg::Vector{Int64} = fill(0,nspins)

    s1::Int64 = 0
    s2::Int64 = 0
    s3::Int64 = 0
    s4::Int64 = 0

    random_plaq::Int64 = 0

    for i in 1:mc_steps
        for j in 1:Nmix
            while true
                random_plaq = rand(1:nplaq)
                s1 = plaq_spins[random_plaq,1]
                s2 = plaq_spins[random_plaq,2]
                s3 = plaq_spins[random_plaq,3]
                s4 = plaq_spins[random_plaq,4]

                if config[s1] == 1 && config[s4] == 1
                    config[s1] = -1
                    config[s4] = -1
                    config[s2] = 1
                    config[s3] = 1
                    break
                elseif config[s3] == 1 && config[s2] == 1
                    config[s1] = 1
                    config[s4] = 1
                    config[s2] = -1
                    config[s3] = -1
                    break
                end
            end
        end
        sz_avg .= sz_avg .+ config
        #add_vec!(sz_avg, config)
        sz_correlator(config, szsz, spin_origin)

    end

    sz_corr_array_float= szsz/mc_steps
    sz_avg_float = sz_avg/mc_steps

    return sz_avg_float, sz_corr_array_float
end

function monomer_move(config::Vector{Int64}, npts::Int64, Gedge_no::Matrix{Int64}, nbrG::Vector{Vector{Int64}})
    initial_site::Int64 = rand(1:npts)
    nbr_initial::Vector{Int64} = nbrG[initial_site]

    for i in nbr_initial
        if config[Gedge_no[initial_site, i]] == 1
            # there is a dimer
            return false
        end
    end

    # at this point, we know initial_site is a monomer
    special_nbr::Int64 = rand(nbr_initial)
    nnbr::Vector{Int64} = nbrG[special_nbr]
    for i in nnbr
        if config[Gedge_no[i,special_nbr]] == 1
            config[Gedge_no[i,special_nbr]] = -1
            config[Gedge_no[initial_site, special_nbr]] = 1
            return true
        end
    end
    return false
end

function plaq_move(config::Vector{Int64}, nplaqs::Int64, plaq_spins::Matrix{Int64})
    random_plaq = rand(1:nplaqs)
    s1 = plaq_spins[random_plaq,1]
    s2 = plaq_spins[random_plaq,2]
    s3 = plaq_spins[random_plaq,3]
    s4 = plaq_spins[random_plaq,4]

    if config[s1] == 1 && config[s4] == 1
        config[s1] = -1
        config[s4] = -1
        config[s2] = 1
        config[s3] = 1
        return true
    elseif config[s3] == 1 && config[s2] == 1
        config[s1] = 1
        config[s4] = 1
        config[s2] = -1
        config[s3] = -1
        return true
    end
    return false
end
    
function mc_mobile_full_no_avg(mc_steps::Int64, Nmix::Int64, config::Vector{Int64} ,  plaq_spins::Matrix{Int64}, spin_origin::Int64, Gedge_no::Matrix{Int64} , nbrG::Vector{Vector{Int64}})
    
    nplaqs::Int64 = size(plaq_spins)[1]
    nspins::Int64 = length(config)
    sz_corr_array::Vector{Int64} = fill(0,nspins)
    sz_avg::Vector{Int64} = fill(0,nspins)
    actual_steps = 0


    for i in 1:mc_steps
        for j in 1:Nmix
            while !monomer_move(config, npts, Gedge_no, nbrG)
            end
            while !plaq_move(config, nplaqs, plaq_spins)
            end
        end
        add_vec!(sz_avg, config)
        sz_correlator(config, sz_corr_array, spin_origin)
    end

    sz_corr_array_float= sz_corr_array/mc_steps
    sz_avg_float = sz_avg/mc_steps

    return sz_avg_float, sz_corr_array_float
end

function get_monomer_membranes_numerics(sz)
    membrane_edges = []
    for i in 1:nspins
        if abs(sz[i] +1.0)<1e-3
            a, b = spins_to_edges[i,:]
            push!(membrane_edges, [AtoG[a], BtoG[b]])
        end
    end
    return convert_to_int(membrane_edges)
end

function generate_disconnected_graph(sz)
    disg =  Graph(copy(pts), copy(graph), copy(npts), copy(plaq), copy(nplaqs))
    monomer_membrane = get_monomer_membranes_numerics(sz)
    for i in monomer_membrane
        p1, p2 = i
        disg.adj[p1,p2] = 0
        disg.adj[p2,p1] = 0
    end
    return disg
end

function sz_correlator_mat(config::Vector{Int64}, running_sum::Matrix{Int64}, nspins::Int64)
    for origin in 1:nspins-1
        if config[origin] == 1
            for pt2 in origin+1:nspins
                running_sum[pt2,origin] = running_sum[pt2,origin] + config[pt2]
            end
        else
            for pt2 in origin+1:nspins
                running_sum[pt2,origin] = running_sum[pt2,origin] - config[pt2]
            end 
        end
    end
end

function get_good_path(src::Int64, len::Int64, spin_cno::Matrix{Int64},vison_spins_graph_nbr:: Vector{Vector{Int64}}, spins_to_plaqs::Vector{Vector{Int64}}, sub_nplaq::Int64, Ntries::Int64 = 10)
    path = [src]
    plaq_visited = fill(0, sub_nplaq)
    for i in spins_to_plaqs[src]
        plaq_visited[i] += 1
    end
    for i in 1:Ntries
        if len<0
            return path
        end
        next_pt = 0
        for j in shuffle(vison_spins_graph_nbr[path[end]])
            # println("j = ",j)
            if j ∈ path
                # backtracking
                continue
            end
            plaqs = spins_to_plaqs[j]
            if length(plaqs) == 1 && plaq_visited[plaqs[1]] <= 1
                push!(path, j)
                len -= 1
                plaq_visited[plaqs[1]] += 1
                break
            elseif length(plaqs) == 2 && (plaq_visited[plaqs[1]] == 0 && plaq_visited[plaqs[2]] == 1) || (plaq_visited[plaqs[1]] == 1 && plaq_visited[plaqs[2]] == 0)
                push!(path, j)
                len -= 1
                plaq_visited[plaqs[1]] += 1
                plaq_visited[plaqs[2]] += 1
                break
            end
        end
    end
    return path
end

function get_good_path_maxl(src::Int64, len::Int64, spins_to_plaqs::Vector{Vector{Int64}}, sub_nplaq::Int64,vison_spins_graph_nbr:: Vector{Vector{Int64}}, Ntries, Nsamp = 1000)
    spin_cno = sum(sub_penrose_spin.adj, dims=1)
    # spin_nbr = convert_to_int(get_sub_neighbors(sub_penrose_spin))
    list_paths = []
    for i in 1:Nsamp
        new_path = get_good_path(src, len, spin_cno, vison_spins_graph_nbr,spins_to_plaqs, sub_nplaq, Ntries)
        push!(list_paths,new_path )
        if length(new_path) == len
            break
        end
    end
    mxval, mxindex = findmax(length.(list_paths))
    return list_paths[mxindex]
end

function generate_symm_path(path::Vector{Int64}, pts::Vector{ComplexF64})
    symm_path = []
    for i in 0:4
        curr_path::Vector{Int64} = []
        curr_path_r::Vector{Int64} = []
        for j in path
            symm_pt = ROT^i * pts[j]
            symm_pt_r = reflectY( ROT^i * pts[j] )
            push!(curr_path, pos_to_spin_index(symm_pt, pts))
            push!(curr_path_r, pos_to_spin_index(symm_pt_r, pts))
        end
        push!(symm_path, curr_path)
        push!(symm_path, curr_path_r)
    end
    return transp( convert_to_matrix_int(symm_path))
end

function vis(sub_config::Vector{Int64}, running_sum::Vector{Matrix{Int64}}, npaths::Int64,lpath::Vector{Int64},  paths::Vector{Vector{Int64}})
    for i in 1:npaths
        for j in 1:lpath[i]-1
            prod = sub_config[paths[i][j]]
            for k in j+1:lpath[i]
                prod = prod * sub_config[paths[i][k]]
                running_sum[i][k,j] = running_sum[i][k,j] + prod
            end
        end
    end
end
function vis_partial(config::Vector{Int64}, running_sum::Vector{Matrix{Int64}}, npaths::Int64,lpath::Vector{Int64},  mega_path::Vector{Matrix{Int64}})
    for i in 1:npaths
        for j in 1:10
            prod = config[mega_path[i][1,j]]
            for k in 2:lpath[i]
                prod = prod * config[mega_path[i][k,j]]
                running_sum[i][k,j] = running_sum[i][k,j] + prod
            end
        end
    end
end
function symmetrise_vison(vison::Vector{Matrix{Float64}}, lpath::Vector{Int64}, npaths::Int64)
    for i in 1:npaths
        for j in 1:lpath[i]-1
            for k in j+1:lpath[i]
                vison[i][j,k] = vison[i][k,j]
            end
            vison[i][j,j] = 1
        end
        vison[i][lpath[i], lpath[i]] = 1
    end
end

function post_process_vison_lite(vison, paths)
    distance = []
    vis_vector = []
    for i in 1:length(paths)
        current_path = paths[i]
        for j in 1:length(current_path)-1
            for k in j+1:length(current_path)
                push!(distance, abs(sub_spins_pos[current_path[j]] - sub_spins_pos[current_path[k]]))
                push!(vis_vector, abs(vison[i][j,k]))
            end
        end
    end
    return distance, vis_vector
end

function post_process_vison(vison, paths;  averaging_radius = 0.4)    
    distance = []
    vis_vector = []
    for i in 1:length(paths)
        current_path = paths[i]
        for j in 1:length(current_path)-1
            for k in j+1:length(current_path)
                push!(distance, abs(sub_spins_pos[current_path[j]] - sub_spins_pos[current_path[k]]))
                push!(vis_vector, abs(vison[i][j,k]))
            end
        end
    end
    #display(scatter(distance, (vis_vector)))

    clean_dist = []
    clean_vis = []
    visited = fill(0,length(distance))
    for i in 1:length(distance)
        if visited[i] == 1 
            continue
        end
        push!(clean_dist, distance[i])
        temp = []
        for j in 1:length(distance)
            if abs(distance[i] - distance[j]) < averaging_radius && visited[j] == 0
                push!(temp, vis_vector[j])
                visited[j] = 1
            end
        end
        # visited[i] = 1
        push!(clean_vis,  mean(temp))
        if i== 1
            println(length(temp))
        end
    end
    # display( scatter( clean_dist, clean_vis ) )
    display( scatter( (clean_dist), log.(clean_vis) ) )
    return clean_dist, clean_vis
end

function generate_monomer_path(l::Int64,sub_nbrG::Vector{Vector{Int64}}, Ntries::Int64=100)
    npts = sub_penrose.n
    counter = l
    start = rand(1:npts)
    path = [start, rand(sub_nbrG[start])]
    done = false
    for i in 1:Ntries
        current_pt = path[end]
        previous_pt = path[end-1]
        rnbr = shuffle(sub_nbrG[current_pt])
        closing = false
        for next_pt in rnbr
            if next_pt != previous_pt
                if next_pt in path
                    push!(path, next_pt)
                    # println("Found loop")
                    temp = [pop!(path)]
                    while length(path) > 0
                        new_pt = pop!(path)
                        if new_pt == next_pt
                            # if (cno[temp[1]] == 3 || cno[temp[1]] == 5 ) && (cno[temp[end]] == 3 || cno[temp[end]] == 5)
                            #     return temp
                            # else
                            #     return []
                            # end 
                            return temp
                        else
                            push!(temp, new_pt)
                        end
                    end
                    return path
                else
                    push!(path, next_pt)
                    break
                end
            end
        end
        if length(path) >= l
            return []
        end
        
    end
    return  []
end

function split_loop(loop)
    len = length(loop)
    # if isodd(len)
    #     println("Not a closed loop")
    #     return []
    # end
    return [loop[1:len ÷ 2], loop[len÷2 + 1 : end] ]
end
function boundary_path_check(path, spin_graph)
    for sp in path
        count = 0
        for j in 1:spin_graph.nplaq
            if sp in sub_penrose_spin.plaq[j,:]
                count = count + 1
            end
        end
        if count == 1
            return true
        end
    end
    return false
end

# OLD FUNCTION
# function generate_mbasis()
#     lval = [4, 6, 8, 10, 12, 14, 16, 18, 20]
#     basis = [[] for i in 1:length(lval)]
    
#     np = 2000
#     counter = 0
#     while counter < 1000
#         curr = generate_monomer_path(100,get_nbrG(sub_penrose), 200)
#         if minimum(lval) <= length(curr) <= maximum(lval)
#             ind = (length(curr) - lval[1]) ÷ 2 +1
#             # ind = length(curr) ÷ lval[1] # -2
#             if length(basis[ind]) < np
#                 push!(basis[ind], curr)
#             end
#         end
#         if norm( length.(basis) .- np ) < 1e-5
#             spin_basis = [[lat_path_to_spin_path(i) for i in basis[j]] for j in 1:length(basis)]
#             return [b for b in spin_basis]
#         end
#         counter += 1
#     end
#     spin_basis = [[lat_path_to_spin_path(i) for i in basis[j]] for j in 1:length(basis)]
#     return [b for b in spin_basis]
#     # return nothing
# end

function generate_mbasis(sub_penrose::Graph)
    lval = [4, 6, 8, 10, 12, 14, 16, 18, 20]
    basis = [[] for i in 1:length(lval)]
    
    np = 2000
    counter = 0
    nbrG_sub_penrose = get_nbrG(sub_penrose)
    while counter < 400 * 10
        curr = generate_monomer_path(100,nbrG_sub_penrose, 200)
        if minimum(lval) <= length(curr) <= maximum(lval)
            ind = (length(curr) - lval[1]) ÷ 2 +1
            # ind = length(curr) ÷ lval[1] # -2
            if length(basis[ind]) < np
                push!(basis[ind], curr)
            end
        end
        if norm( length.(basis) .- np ) < 1e-5
            spin_basis = [[lat_path_to_spin_path(i) for i in basis[j]] for j in 1:length(basis)]
            return [b for b in spin_basis]
        end
        counter += 1
    end
    spin_basis = [[lat_path_to_spin_path(i) for i in basis[j]] for j in 1:length(basis)]
    return [b for b in spin_basis]
    # return nothing
end


function shift(bbb, i)
    aaa = copy(bbb)
    for j in 1:i
        x = pop!(aaa)
        pushfirst!(aaa, x)
    end
    return aaa
end
function get_nbrG(g)
    n = g.n
    graph = g.adj
    nbrG = []
    for i in 1:n
        temp = []
        for j in 1:n
            if graph[i,j] == 1
                push!(temp, j)
            end
        end
        push!(nbrG, temp)
    end
    return convert_to_int(convert_to_int(nbrG))
end
function get_spin_nbr(g, edges_to_spins)
    n = g.n
    graph = g.adj
    nbr = []
    for i in 1:n
        temp = []
        for j in 1:n
            if graph[i,j] == 1
                push!(temp, edges_to_spins[j,i])
            end
        end
        push!(nbr, temp)
    end
    return convert_to_int(nbr)
end

function ismonomer(config::Vector{Int64}, site::Int64, sub_spin_nbr::Vector{Vector{Int64}})
    for k in sub_spin_nbr[site]
        if config[k] > 0
            return false
        end
    end
    return true
end
function foo1(config::Vector{Int64},mfmnum::Vector{Matrix{Int64}}, i::Int64, j::Int64 )
    siz = 700
    p = plot(size = (siz,siz), aspect_ratio = :equal , legend = false, ylim = (-3, 20), xlim=(5,25))
    plot_edges(p,sub_penrose)
    # plot_sites(p, sub_penrose_spin, :blue, 0; annotate = false, markersize=  1)
    sub_plot_path(p,mfmnum[i][:,j])
    plot_sub_config(p,config)
    display(p)
end
function foo2(config::Vector{Int64},mfmnum1::Vector{Matrix{Int64}}, mfmnum2::Vector{Matrix{Int64}}, i::Int64, j::Int64 )
    siz = 700
    p = plot(size = (siz,siz), aspect_ratio = :equal , legend = false, ylim = (-3, 20), xlim=(5,25))
    plot_edges(p,sub_penrose)
    # plot_sites(p, sub_penrose_spin, :blue, 0; annotate = false, markersize=  1)
    sub_plot_path(p,mfmnum1[i][:,j])
    sub_plot_path(p,mfmnum2[i][:,j])
    plot_sub_config(p,config)
    display(p)
end
function mon_fm(config::Vector{Int64}, num1::Matrix{Int64}, num2::Matrix{Int64}, den::Matrix{Int64}, mfmnum1::Vector{Matrix{Int64}}, mfmnum2::Vector{Matrix{Int64}}, npath::Int64, lpath1::Vector{Int64}, lpath2::Vector{Int64}, cap::Array{Int64, 3}, sub_spin_nbr::Vector{Vector{Int64}})
    for i in 1:npath
        for j in 1:10
            # First thing: check if the path is alternating
            alternating_path1 = true
            s1 = config[mfmnum1[i][1,j]]
            for k in 2:lpath1[i]
                s2 = config[mfmnum1[i][k,j]]
                if s1 == s2
                    alternating_path1 = false
                    break
                end
                s1 = s2
            end
            if alternating_path1
                monomer1 = ismonomer(config, cap[1,j,i], sub_spin_nbr)
                if monomer1
                    num1[j,i] = num1[j,i] + 1
                    # println("here")
                else
                    monomer2 = ismonomer(config, cap[2,j,i], sub_spin_nbr)
                    if monomer2
                        num1[j,i] = num1[j,i] + 1
                        # println("-------------")
                        # for k in 1:lpath1[i]
                        #     println(config[mfmnum1[i][k,j]])
                        # end
                        # println("-------------")
                        # foo1(config, mfmnum1, i, j)
                        # println("Here")
                    end
                end
            end

            alternating_path2 = true
            s1 = config[mfmnum2[i][1,j]]
            for k in 2:lpath2[i]
                s2 = config[mfmnum2[i][k,j]]
                if s1 == s2
                    alternating_path2 = false
                    break
                end
                s1 = s2
            end
            if alternating_path2
                monomer1 = ismonomer(config, cap[1,j,i], sub_spin_nbr) # repeated calculation, potential to improve
                if monomer1
                    num2[j,i] = num2[j,i] + 1
                else
                    monomer2 = ismonomer(config, cap[2,j,i], sub_spin_nbr) # repeated calculation, potential to improve
                    if monomer2
                        num2[j,i] = num2[j,i] + 1
                    # foo1(config, mfmnum2, i, j)
                    end
                end
            end
            # Now take care of the denominator
            if alternating_path1 && alternating_path2
                if config[mfmnum1[i][1,j]] != config[mfmnum2[i][end,j]]
                    den[j,i] = den[j,i] + 1
                    # foo2(config,mfmnum1, mfmnum2, i, j)
                end
            end
        end
    end
end

function post_process_mfm(mfm_arr, mfmpath)
    mfmnum1, mfmnum2 = mfmpath
    mat = (mfm_arr[1] .* mfm_arr[2] ./ mfm_arr[3])
    num1 = []
    num2 = []
    den = []
    corr = []
    dist = []
    npath = size(mfm_arr[1])[2]
    for i in 1:npath
        n1 = mean(mfm_arr[1][:,i])
        n2 = mean(mfm_arr[2][:,i])
        de = mean(mfm_arr[3][:,i])
        push!(num1, n1)
        push!(num2, n2)
        push!(den, de)
        push!(dist, size(mfmpath[1][i])[1] + size(mfmpath[2][i])[1])
        push!(corr, n1*n2/de)
    end
    return dist, num1, num2, den, corr
end

function vis_fm(config::Vector{Int64}, num1::Matrix{Int64}, num2::Matrix{Int64}, 
    den::Matrix{Int64}, vfmnum1::Vector{Matrix{Int64}}, vfmnum2::Vector{Matrix{Int64}}, 
    npath::Int64, lpath1::Vector{Int64}, lpath2::Vector{Int64})
    for i in 1:npath
        for j in 1:10
            c1 = 1
            for k in 1:lpath1[i]
                if config[vfmnum1[i][k,j]] < 0
                    c1 = -c1
                end
            end
            num1[j,i] = num1[j,i] + c1
            c2 = 1
            for k in 1:lpath2[i]
                if config[vfmnum2[i][k,j]] < 0
                    c2 = -c2
                end
            end
            num2[j,i] = num2[j,i] + c2
            if c1 == c2 
                den[j,i] = den[j,i] + 1
            else
                den[j,i] = den[j,i] - 1
            end
        end
    end
end

function check_good_vison_path(path, plaq)
    nplaq = size(plaq)[1]
    visit = fill(0,nplaq)
    for pt in path
        for j in 1:nplaq
            if pt in plaq[j,:]
                visit[j] = visit[j] + 1
            end
        end
    end
    c = counter(visit)
    if c[1] != 2
        println("Error: More than 2 ones")
        return false
    elseif c[3] != 0
        println("Error: Three exists")
        return false
    elseif c[4] != 0
        println("Error: Four exists")
        return false
    end
    return true
end

function check_good_vison_loop(visit)
    c = counter(visit)
    if c[1] != 0
        println("Error: One exists")
        return false
    elseif c[3] != 0
        println("Error: Three exists")
        return false
    elseif c[4] != 0
        println("Error: Four exists")
        return false
    end
    return true
end

function generate_vfm_loop(l::Int64,vison_spins_graph_nbr::Vector{Vector{Int64}},spin_to_plaq::Vector{Vector{Int64}}, cplaq::Vector{Int64}, nspins::Int64, nplaq::Int64, visit::Vector{Int64}, Ntries::Int64=100)
    start = rand(1:nspins)
    counter = l
    path = [start, rand(vison_spins_graph_nbr[start])]
    if cplaq[path[1]] <2 || cplaq[path[2]] < 2
        return []
    end
    for i in 1:nplaq
        visit[i] = 0
    end
    for p in spin_to_plaq[path[1]]
        visit[p] += 1
    end
    for p in spin_to_plaq[path[2]]
        visit[p] += 1
    end
    for i in 1:Ntries
        current_pt = path[end]
        previous_pt = path[end-1]
        rnbr = shuffle(vison_spins_graph_nbr[current_pt])
        for next_pt in rnbr
            if next_pt != previous_pt
                if next_pt in path
                    push!(path, next_pt)
                    counter = counter - 1
                    temp = [pop!(path)]
                    while length(path) > 0
                        new_pt = pop!(path)
                        if new_pt == next_pt
                            return temp
                            # return temp, visit
                        else
                            push!(temp, new_pt)
                        end
                    end
                    return path
                    # return path, visit
                else
                    p = spin_to_plaq[next_pt]
                    if length(p)<2
                        break
                    end
                    if visit[p[1]] > 1 || visit[p[2]] > 1
                        break
                    end
                    push!(path, next_pt)
                    counter = counter - 1
                    visit[p[1]] = visit[p[1]] + 1
                    visit[p[2]] = visit[p[2]] + 1
                    break
                end
            end
        end
        if counter < 0 
            return []
        end
        
    end
    return  []
end

function generate_vbasis(vison_spins_graph_nbr::Vector{Vector{Int64}}, spin_to_plaq::Vector{Vector{Int64}}, len::Vector{Int64})
    Ncalls = 1000*1000*50
    result_list = []
    nspins = length(vison_spins_graph_nbr)
    nplaq = size(plaq)[1]
    visit = fill(0,nplaq)

    for i in 1:Ncalls
        result = generate_vfm_loop(100, vison_spins_graph_nbr, spin_to_plaq,  len, nspins, nplaq, visit)
        if result != []
            push!(result_list, result)
        end
    end
    return [res for res in result_list]
end

function generate_vfnumbasis(vbasis, sub_penrose_spin)
    Npth = 110
    vfmnum_basis1 = []
    vfmnum_basis2 = []
    cntr = counter(length.(vbasis))
    good_length = []
    for i in 4:maximum(length.(vbasis))
        if cntr[i] > Npth
            push!(good_length, i)
        end
    end
    println(good_length)
    acc = fill(0,maximum(length.(vbasis)))
    for i in 1:length(vbasis)
        lcurrpath = length(vbasis[i])
        if lcurrpath in good_length
            if acc[lcurrpath] < Npth
                pth1, pth2 = split_loop(vbasis[i])
                push!(vfmnum_basis1, pth1)
                push!(vfmnum_basis2, pth2)
                acc[lcurrpath] += 1
                if length(pth1) + length(pth2) != lcurrpath
                    println("cry")
                end
            end
        end
    end
    vfmnum_basis1 = convert_to_int(vfmnum_basis1)
    vfmnum_basis2 = convert_to_int(vfmnum_basis2)
    
    vfmnum1::Vector{Matrix{Int64}} = basis_path_to_all_path(vfmnum_basis1, sub_penrose_spin.pts)
    vfmnum2::Vector{Matrix{Int64}} = basis_path_to_all_path(vfmnum_basis2, sub_penrose_spin.pts)
    npath = length(vfmnum1)
    lpath1 = [size(vfmnum1[i])[1] for i in 1:npath]
    lpath2 = [size(vfmnum2[i])[1] for i in 1:npath]
    vfmpath = [vfmnum1, vfmnum2]
    return vfmnum_basis1, vfmnum_basis2, vfmpath, lpath1, lpath2, npath
end

function post_process_vfm(vfm_arr, vfmpath)
    cbin = []
    fmbin = []
    dbin = []
    ebin = []

    num1, num2, den = vfm_arr
    num_path1, num_path2 = vfmpath
    fm = sqrt.( abs.(num1 .* num2 ./ den) )
    println(size(fm))
    mean1, mean2, mean_den, mean_fm = mean(abs.(num1), dims=1), mean(abs.(num2), dims=1), mean(abs.(den), dims=1), mean(fm, dims=1)
    # std_fm = std(fm, dims=1)
    disnum = []
    cnum = [] 
    # cnum1 = []
    # cnum2 = []
    disden = []
    cden = []

    for i in 1:length(num_path1)
        push!(disnum, size(num_path1[i])[1])
        # push!(cnum1, mean1[1,i])
        push!(cnum, mean1[1,i])
        push!(disnum, size(num_path2[i])[1])
        # push!(cnum2, mean2[1,i])
        push!(cnum, mean2[1,i])

        push!(disden, size(num_path1[i])[1]+ size(num_path2[i])[1])
        push!(cden, mean_den[1,i])
    end
    return disnum, cnum,  disden, cden
    # return disnum, cnum1, cnum2,  disden, cden
end

function post_process_vfm_asymmetric(vfm_arr, vfmpath)
    cbin = []
    fmbin = []
    dbin = []
    ebin = []

    num1, num2, den = vfm_arr
    num_path1, num_path2 = vfmpath
    fm =  sqrt.( abs.(num1 .* num2 ./ den) )
    println(size(fm))
    # mean1, mean2, mean_den, mean_fm = mean(abs.(num1), dims=1), mean(abs.(num2), dims=1), mean(abs.(den), dims=1), mean(fm, dims=1)
    # std_fm = std(fm, dims=1)
    disnum = []
    cnum = [] 
    cnum1 = []
    cnum2 = []
    disden = []
    cden = []
    dis = []
    c = []

    for i in 1:length(num_path1)
        for j in 1:10
            # push!(disnum, size(num_path1[i])[1])
            # # push!(cnum1, mean1[1,i])
            # push!(cnum, num1[j,i])
            # push!(disnum, size(num_path2[i])[1])
            # # push!(cnum2, mean2[1,i])
            # push!(cnum, num2[j,i])

            # push!(disnum, size(num_path1[i])[1])
            push!(cnum1, num1[j,i] )
            # push!(disnum, size(num_path2[i])[1])
            push!(cnum2, num2[j,i] )

            push!(disden, size(num_path1[i])[1]+ size(num_path2[i])[1])
            push!(cden, den[j,i])

            push!(dis, size(num_path1[i])[1]+ size(num_path2[i])[1])
            push!(c, sqrt( abs( num1[j,i] * num2[j,i] / den[j,i])) )
            # push!(c, fm[j,i])
        end
    end
    # return dis, abs.(c), disnum, abs.(cnum),  disden, abs.(cden)
    return dis, abs.(c), disden, abs.(cnum1), abs.(cnum2), abs.(cden)
end

function post_process_vfm_euclidean(vfm_arr, vfmpath, sub_penrose_spin)
    cbin = []
    fmbin = []
    dbin = []
    ebin = []

    num1, num2, den = vfm_arr
    num_path1, num_path2 = vfmpath
    fm = sqrt.( abs.(num1 .* num2 ./ den) )
    println(size(fm))
    mean1, mean2, mean_den, mean_fm = mean(abs.(num1), dims=1), mean(abs.(num2), dims=1), mean(abs.(den), dims=1), mean(fm, dims=1)
    # std_fm = std(fm, dims=1)
    disnum = []
    cnum = [] 
    # cnum1 = []
    # cnum2 = []
    disden = []
    cden = []
    pts = sub_penrose_spin.pts

    for i in 1:length(num_path1)
        p_start1 = num_path1[i][1,1]
        p_end1 = num_path1[i][end,1]
        d1 = abs(pts[p_start1] - pts[p_end1])
        push!(disnum, d1)
        push!(cnum, mean1[1,i])

        p_start2 = num_path2[i][1,1]
        p_end2 = num_path2[i][end,1]
        d2 = abs(pts[p_start2] - pts[p_end2])
        push!(disnum, d2)
        push!(cnum, mean2[1,i])

        push!(disden, (d1+d2)/2)
        push!(cden, mean_den[1,i])
    end
    return disnum, cnum,  disden, cden
    # return disnum, cnum1, cnum2,  disden, cden
end


function post_process_vfm_ratio(vfm_arr, vfmpath, sub_penrose_spin)
    num1, num2, den = vfm_arr
    num_path1, num_path2 = vfmpath
    # println(size(fm))
    mean1, mean2, mean_den = mean(abs.(num1), dims=1), mean(abs.(num2), dims=1), mean(abs.(den), dims=1)
    println(size(mean1), "\t sadfasd")
    fm =  abs.(mean1 .* mean2 ./ mean_den) 
    dis  = []
    c = []
    disnum = []
    cnum = [] 
    disden = []
    cden = []

    pts = sub_penrose_spin.pts
    println(size(fm), "\t size fm")

    for i in 1:length(num_path1)
        push!(disnum, size(num_path1[i])[1])
        # push!(cnum1, mean1[1,i])
        push!(cnum, mean1[1,i])
        push!(disnum, size(num_path2[i])[1])
        # push!(cnum2, mean2[1,i])
        push!(cnum, mean2[1,i])

        push!(disden, size(num_path1[i])[1]+ size(num_path2[i])[1])
        push!(cden, mean_den[1,i])

        push!(c, fm[1,i])
        push!(dis, size(num_path1[i])[1]+ size(num_path2[i])[1])
    end

    return dis,  c, disnum, cnum, disden, cden
    # return disnum, cnum1, cnum2,  disden, cden
end



