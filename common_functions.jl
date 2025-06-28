using Statistics

logdev(corr, dev) = abs.( log.(corr) .- log.(corr .+ dev)  )

function swapcol!(x,i,j)
    for k in axes(x, 1)  # <- give dimension as input to axes function
        x[k, i], x[k, j] = x[k, j], x[k, i]
    end
end

function swaprow!(x,i,j)
    for k in axes(x, 2)  # <- give dimension as input to axes function
        x[ i, k], x[ j, k] = x[ j, k], x[ i, k]
    end
end

function linear_fit1(x, y)
    bhat = [x ones(length(x))]\ y
    return bhat
end

function linear_fit(x, y)
    Lx = length(x)
    Ly = length(y)
    L = Lx
    if Lx != Ly
        println("X and Y have different dimensions")
    end
    x_bar = sum(x)/L
    y_bar = sum(y)/L
    b = sum( (x .- x_bar) .* (y .- y_bar) ) / sum( (x .- x_bar) .^ 2 )
    a = y_bar - b * x_bar
    return [a, b]
end 

function print_mat(mat)
    for i in 1:size(mat)[1]
        for j in 1:size(mat)[2]
            print(mat[i,j])
            print(", ")
        end
        print("\n")
    end
end

struct Graph
    pts::Vector{ComplexF64}
    adj::Matrix{Int64}
    n::Int64
    plaq::Matrix{Int64}
    nplaq::Int64
end
Base.copy(g::Graph) = Graph(g.pts, g.adj, g.n, g.plaq, g.nplaq)

function convert_to_int(nbr::Vector{Vector{Any}})
    nbrG = Vector{Int64}[]
    for i in 1:length(nbr)
        temp = Int64[]
        for j in 1:length(nbr[i])
            push!(temp, nbr[i][j])
        end
        push!(nbrG, temp)
    end
    return nbrG
end
convert_to_int(vec::Vector{Any}) = [v for v in vec]

function convert_to_matrix_int(A)
    # assuming A is a list of lists
    A_mat = fill(0, length(A), length(A[1]))
    for i in 1:length(A)
        for j in 1:length(A[i])
            A_mat[i,j] = A[i][j]
        end
    end
    return A_mat
end
function transp(M::Matrix{Int64})
    a,b = size(M)
    A = fill(0,b, a)
    for i in 1:a
        for j in 1:b
            A[j,i] = M[i,j]
        end
    end
    return A
end

function comma_separate(name::String)
    mega_list::Vector{String} = []
    list::Vector{Char} = []
    for i in 1:length(name)
        if name[i] == ','
            push!(mega_list, String(list))
            list = []
        elseif name[i] == ' '
            nothing
        else
            push!(list, name[i])
        end
    end
    push!(mega_list, String(list))
    return mega_list
end

function pos_to_spin_index(position, spins_pos)
    for i in 1:length(spins_pos)
        if abs(position - spins_pos[i]) < 1e-5
            return i
        end
    end
    return nothing
end
reflectY(z::ComplexF64) = -conj(z)

# --------------------------------- PLOTTING FUNCTIONS ---------------------------- 

function plot_generate(size = 500; legend = false, ylim = false)
    # function plot_generate(size = 500; legend = false, ylim = (-ymax-0.1,ymax+0.1))
    if ylim == false
        return plot(size = (size,size), aspect_ratio = :equal , legend = legend)
    else
        return plot(size = (size,size), aspect_ratio = :equal , legend = legend, ylim = ylim)
    end
end

function plot_edges(p, g, linecolor = "pink", linewidth = 2)
    x = real(g.pts)
    y = imag(g.pts)
    for i in 1:g.n-1
        for j in i+1:g.n
            if g.adj[i,j] == 1
                plot!(p, [x[i], x[j]],[y[i] , y[j]], linecolor=linecolor, linewidth=linewidth, label="")
            end
        end
    end
    return p
end

function plot_sites(p, g, mc=:blue, markerstrokewidth=0; annotate = false, markersize=  6)
    if annotate
        scatter!(p, real(g.pts), imag(g.pts), mc=mc, markerstrokewidth=markerstrokewidth, series_annotations = text.(1:g.n), label = "",  markersize=markersize)
    else
        scatter!(p, real(g.pts), imag(g.pts), mc=mc, markerstrokewidth=markerstrokewidth, label = "",  markersize=markersize)
    end
    return p
end

function plot_single_point(p, point, mc=:magenta, markersize=6)
    scatter!(p, [real(point)], [imag(point)], mc=mc, markerstrokewidth=0, label="", markersize=markersize)
    return p
end

function plot_graph(p,g, edge_color="pink", linewidth = 2, point_color=:blue)
        p = plot_edges(p,g, edge_color, linewidth)
        p = plot_sites(p, g, point_color)
    return p
end

function plot_membrane(p, membrane_edges) 
    #membrane_edges = get_membrane_edges_theory()
    for l in membrane_edges
        i, j = l
        plot!(p, [pts_x[i], pts_x[j]],[pts_y[i] , pts_y[j]], linecolor="brown", linewidth=5, label="")
    end
    return p
end

function plot_config(p, config, size = 1000, annotate = false)
    ms = 6    
    for i in 1:nspins
        if config[i] == 1
            a, b = spins_to_edges[i,:]
            ptA = AtoG[a]
            ptB = BtoG[b]
            plot!(p, [pts_x[ptA], pts_x[ptB]],[pts_y[ptA] , pts_y[ptB]], linecolor="black", linewidth=4)
            scatter!(p, [pts_x[ptA]], [pts_y[ptA]], mc=:blue, markerstrokewidth=0, markersize = ms)
            scatter!(p, [pts_x[ptB]], [pts_y[ptB]], mc=:green, markerstrokewidth=0, markersize = ms)
        end
    end
    m_list = monomer_list(config)
    scatter!(p, real(pts[m_list]), imag(pts[m_list]), mc=:red, markerstrokewidth=0, markersize = ms)#, series_annotations = text.(1:npts))
    if annotate
        scatter!(p, real(pts), imag(pts), mc=:blue, markerstrokewidth=0, series_annotations = text.(1:npts))
    end
    #display(p)
    return p
end



function plot_sub_lat_edges(p)
    for i in 1:length(sub_lat)-1
        for j in i+1:length(sub_lat)
            if sub_graph[i,j] == 1
                plot!(p, [pts_x[sub_lat[i]], pts_x[sub_lat[j]]],[pts_y[sub_lat[i]] , pts_y[sub_lat[j]]], linecolor="pink", linewidth=2)
            end
        end
    end
    return p
end

function plot_sub_lat(p, annotate=false)
    
    if annotate
        scatter!(p,real(sub_pts), imag(sub_pts), mc=:red, markerstrokewidth=0, series_annotations = text.(1:sub_npts))
    else
        scatter!(p,real(sub_pts), imag(sub_pts), mc=:red, markerstrokewidth=0)
    end
    return p
end

function plot_sub_config(p, sub_config, size = 1000, annotate = false)
    ms = 6

    for i in 1:sub_nspins
        if sub_config[i] == 1
            a, b = spins_to_edges[sub_spins[i],:]
            ptA = AtoG[a]
            ptB = BtoG[b]
            plot!(p, [pts_x[ptA], pts_x[ptB]],[pts_y[ptA] , pts_y[ptB]], linecolor="black", linewidth=4)
            scatter!(p, [pts_x[ptA]], [pts_y[ptA]], mc=:blue, markerstrokewidth=0, markersize = ms)
            scatter!(p, [pts_x[ptB]], [pts_y[ptB]], mc=:green, markerstrokewidth=0, markersize = ms)
        end
    end
    return p
end

function plot_sub_spins(p)
    scatter!(p,real(sub_spins_pos), imag(sub_spins_pos), mc=:red, markerstrokewidth=0)
    return p
end

function sub_plot_path(p,path)
    for i in 1:length(path)-1
        pos1 = sub_spins_pos[path[i]]
        pos2 = sub_spins_pos[path[i+1]]
        plot!(p, [real(pos1), real(pos2)],[imag(pos1) ,imag(pos2)], linecolor="orange", linewidth=7)
    end
    return p
end

function plot_path(p,path, g)
    for i in 1:length(path)-1
        pos1 = g.pts[path[i]]
        pos2 = g.pts[path[i+1]]
        plot!(p, [real(pos1), real(pos2)],[imag(pos1) ,imag(pos2)], linecolor="orange", linewidth=7)
    end
    return p
end
# ----------------------------------------------------------------------
# Plotting functions

function plot_generate1(size = 500; legend = false)
    return plot(size = (size,size), aspect_ratio = :equal , legend = legend)
end


function plot_edges(p, g, linecolor = "pink", linewidth = 2)
    x = real(g.pts)
    y = imag(g.pts)
    for i in 1:g.n-1
        for j in i+1:g.n
            if g.adj[i,j] == 1
                plot!(p, [x[i], x[j]],[y[i] , y[j]], linecolor=linecolor, linewidth=linewidth, label="")
            end
        end
    end
    return p
end

function plot_sites(p, g, mc=:blue, markerstrokewidth=0; annotate = false, markersize=  6)
    if annotate
        scatter!(p, real(g.pts), imag(g.pts), mc=mc, markerstrokewidth=markerstrokewidth, series_annotations = text.(1:g.n), label = "",  markersize=markersize)
    else
        scatter!(p, real(g.pts), imag(g.pts), mc=mc, markerstrokewidth=markerstrokewidth, label = "",  markersize=markersize)
    end
    return p
end

function plot_single_point(p, point, mc=:magenta, markersize=6)
    scatter!(p, [real(point)], [imag(point)], mc=mc, markerstrokewidth=0, label="", markersize=markersize)
    return p
end


function plot_graph(p,g, edge_color="pink", linewidth = 2, point_color=:blue)
        p = plot_edges(p,g, edge_color, linewidth)
        p = plot_sites(p, g, point_color)
    return p
end


# --------------------------

function distance_trim(d, c, dev, dmin, dmax)
    d_trim = []
    c_trim = []
    dev_trim = []
    for i in 1:length(d)
        if dmin <= d[i] <= dmax && !(isnan(c[i])) && (c[i] != Inf) && (c[i] != -Inf)
            push!(d_trim, d[i])
            push!(c_trim, c[i])
            push!(dev_trim, dev[i])
        end
    end
    return d_trim, c_trim, dev_trim
end 

function distance_trim(d, c, dmin, dmax)
    d_trim = []
    c_trim = []
    for i in 1:length(d)
        if dmin <= d[i] <= dmax && !(isnan(c[i])) && (c[i] != Inf) && (c[i] != -Inf)
            push!(d_trim, d[i])
            push!(c_trim, c[i])
        end
    end
    return d_trim, c_trim
end 

function distance_average_int(distance, correlator)
    # p = sortperm(distance)
    # dist = distance[p]
    # corr = correlator[p]
    d = []
    c = []
    dev = []
    count = 0
    for i in 0:maximum(distance)
        temp = []
        for j in 1:length(distance)
            if abs(i-distance[j]) < 1/2 && correlator[j] != 0.0 && correlator[j] != Inf && !(isnan(correlator[j]))
                push!(temp, correlator[j])
            end
        end
        # println(temp)
        if length(temp) > 0
            push!(d, i)
            push!(c, mean(temp))
            push!(dev, std(temp))
        end
    end
    return d, c, dev
end

function distance_average_int_env(distance, correlator, fraction = 1.0)
    # p = sortperm(distance)
    # dist = distance[p]
    # corr = correlator[p]
    d = []
    c = []
    dev = []
    log_dev = []
    env_up = []
    env_dn = []
    bins = []
    count = 0
    for i in 0:maximum(distance)
        temp = []
        for j in 1:length(distance)
            if abs(i-distance[j]) < 1/2 && correlator[j] != 0.0 && correlator[j] != Inf && !(isnan(correlator[j]))
                push!(temp, correlator[j])
            end
        end
        # println(temp)
        if length(temp) > 0
            push!(d, i)
            push!(c, mean(temp))
            push!(dev, std(temp))
            push!(log_dev, std(log10.(abs.(temp))))
            nsamp = length(temp)
            push!(bins, temp)
            dat = sort(temp)
            ndat = length(dat)
            bot = Int64(round((1-fraction)/2 * ndat ; digits=0))
            top = Int64(round( (1 - (1-fraction)/2) * ndat ; digits = 0))
            interval = dat[bot+1:top]
            
            push!(env_up, maximum(interval))
            push!(env_dn, minimum(interval))
        end
    end
    return d, c, dev, env_up, env_dn, log_dev, bins
end

function distance_average(distance, corr, averaging_radius=0.4)
    clean_dist = []
    clean_corr = []
    dev = []
    visited = fill(0,length(distance))
    for i in 1:length(distance)
        if visited[i] == 1 
            continue
        end
        temp = []
        for j in 1:length(distance)
            if abs(distance[i] - distance[j]) < averaging_radius && visited[j] == 0 && (corr[j] != Inf) && (!(isnan(corr[j])))
                push!(temp, corr[j])
                visited[j] = 1
            end
        end
        visited[i] = 1
        
        if length(temp) > 0
            # println("In")
            push!(clean_dist, distance[i])
            push!(clean_corr, mean(temp))
            push!(dev, std(temp))
        end
        if i == 1
            println(length(temp))
        end
    end
    return clean_dist, clean_corr, dev
end

function sf(p,path)
    savefig(p,path*".png")
    savefig(p,path*".pdf")
end

function make_flat(arr)
    temp = []
    for i in arr
        temp = vcat(temp, i)
    end
    return temp
end

function post_process_sz(szsz::Matrix{Float64} , sz::Vector{Float64}, sub_spins_pos::Vector{ComplexF64}, averaging_radius=0.4, fraction=1.0)
    distance = []
    szsz_conn = []
    sub_nspins = length(sub_spins_pos)
    for i in 1:sub_nspins-1
        # if abs(sz[i]) == 1.0
        #     continue
        # end

        for j in i+1:sub_nspins
            # if abs(sz[j]) == 1.0
            #     continue
            # end
            append!(distance, abs(sub_spins_pos[i] - sub_spins_pos[j]) )
            append!(szsz_conn, abs(szsz[i,j] - sz[i]*sz[j]))
        end
    end
    # display(scatter(distance, (abs.(szsz_conn))))
    clean_dist = []
    clean_szsz_conn = []
    logcorr = []
    logerr = []
    dev = []
    env_up = []
    env_dn = []
    visited = fill(0,length(distance))
    for i in 1:length(distance)
        if visited[i] == 1 
            continue
        end
        
        push!(clean_dist, distance[i])
        temp = []
        for j in 1:length(distance)
            if abs(distance[i] - distance[j]) < averaging_radius && visited[j] == 0
                push!(temp, szsz_conn[j])
                visited[j] = 1
            end
        end
        visited[i] = 1

        dat = sort(temp)
        ndat = length(dat)
        bot = Int64(round((1-fraction)/2 * ndat ; digits=0))
        top = Int64(round( (1 - (1-fraction)/2) * ndat ; digits = 0))
        interval = dat[bot+1:top]
        
        push!(env_up, maximum(interval))
        push!(env_dn, minimum(interval))

        push!(clean_szsz_conn,  mean(temp))
        # push!(env_up, maximum(temp))
        # push!(env_dn, minimum(temp))
        push!(dev, std(temp))
        push!(logcorr, log.(mean(temp)))
        push!(logerr, std(log.(temp)))
    end
    return clean_dist, clean_szsz_conn, dev, logcorr, logerr, env_up, env_dn
    # return distance, szsz_conn
end

function post_process_basic_vison(vison::Vector{Matrix{Float64}}, mega_path::Vector{Matrix{Int64}}, pts; rad = 0.4)
    v = []
    d = []
    for i in 1:length(mega_path)
        for j in 1:size(mega_path[i])[1]
            for k in 1:size(mega_path[i])[2]
                push!(v, abs(vison[i][j,k]))
                p1 = mega_path[i][1,k]
                p2 = mega_path[i][j,k]
                push!(d, abs(pts[p1] - pts[p2]))
            end
        end
    end
    # display(scatter(d,v))
    return distance_average(d, v, rad)
end    

function post_process_basic_vison_path_length(vison::Vector{Matrix{Float64}}, mega_path::Vector{Matrix{Int64}}, pts; rad = 0.4)
    v = []
    d = []
    for i in 1:length(mega_path)
        for j in 1:size(mega_path[i])[1]
            for k in 1:size(mega_path[i])[2]
                push!(v, abs(vison[i][j,k]))
                push!(d, abs(j-k))
            end
        end
    end
    # display(scatter(d,v))
    return d, v # distance_average(d, v, rad)
end   

function generate_spins_pos()
    spins_pos = fill(0.0 + 1im*0., nspins);
    for i in 1:nspins
        A = spins_to_edges[i,1]
        B = spins_to_edges[i,2]
        spins_pos[i] = (pts[AtoG[A]] + pts[BtoG[B]])/2
    end
    return spins_pos
end

function reset_config(max_match)
    config = fill(-1, nspins)
    for i in 1:size(max_match)[1]
        config[ABedge_no[max_match[i,1],max_match[i,2]]] = 1
    end
    return config
end
# Generate the spin adjacancy matrix

function generate_spins_adj(plaq_spins, nspins)
    spins_adj = fill(0,nspins, nspins)
    for i in 1:nplaqs
        s1, s2, s3, s4 = plaq_spins[i,:]
        spins_adj[s1, s2] = 1
        spins_adj[s2, s4] = 1
        spins_adj[s4, s3] = 1
        spins_adj[s3, s1] = 1

        spins_adj[s2, s1] = 1
        spins_adj[s4, s2] = 1
        spins_adj[s3, s4] = 1
        spins_adj[s1, s3] = 1
    end
    return spins_adj
end

function generate_max_match_adj(max_match, penrose)
    max_match_adj = fill(0, penrose.n, penrose.n)
    A, B = 0,0
    for i in 1:size(max_match)[1]
        ptA = AtoG[max_match[i,1]]
        ptB = BtoG[max_match[i,2]]
        if ptA > npts+1 || ptB > npts+1
            println(i)
        end

        for i in 1:penrose.n
            if abs(pts[ptA] - pts[i]) < 1e-5
                A = i
            end
            if abs(pts[ptB] - pts[i]) < 1e-5
                B = i
            end
        end
        max_match_adj[A,B] = 1
        max_match_adj[B, A] = 1
    end
    return max_match_adj
end

function get_membrane_edges_theory()
    coordination_no = sum(graph; dims =1)
    even_sites = [ i for i in 1:npts if coordination_no[i]%2 == 0];
    even_visited = fill(0,length(even_sites));
    to_visit = copy(even_sites)

    membrane_edges = []

    for l in even_sites
        nbr_start = nbrG[l]
        for j in nbr_start
            if coordination_no[j] == 5
                push!(membrane_edges, [l, j])
            end
        end
    end

    # Remove the (4-)^5 membranes
    membrane_edges_copy = copy(membrane_edges)
    elist = [ (membrane_edges[i][1], membrane_edges[i][2]) for i in 1:length(membrane_edges)];
    g = SimpleGraph(Graphs.SimpleEdge.(elist));
    loops = cycle_basis(g)
    for i in 1:length(loops)
        loop = loops[i]
        if length(loop) == 10
            to_remove = []
            l_loop = length(loop)
            for k in 1:l_loop
                if coordination_no[loop[k]] == 4
                    push!(to_remove, [loop[k], loop[(k)%l_loop+1]])
                elseif coordination_no[loop[k]] == 5
                    push!(to_remove, [loop[(k)%l_loop+1], loop[k]])
                else
                    println("Error!!!")
                end
            end
            filter!(e->e ∉ to_remove,membrane_edges_copy)
        end
    end
    
    return membrane_edges_copy
end

function monomer_list(config)
    monomer_list = Int64[]
    for i in 1:npts
        is_monomer = true
        for j in nbrG[i]
            if config[Gedge_no[i,j]] == 1
                # there is a dimer so not monomer
                is_monomer = false
                break
            end
        end
        if is_monomer
            push!(monomer_list, i)
        end
    end
    return monomer_list
end

function get_connected_corr(sz, szsz, spin_origin)
    return abs.(szsz .- sz[spin_origin] .* sz)
end

# ---------------------- AUTO CORRELATION STUFF ---------------------
# nchanges = 0
# for i in 1:MC_steps-1
#     if s0_array[i]*s0_array[i+1] < 0
#         nchanges += 1
#         #println(i)
#     end
# end
# println("1 change every ", MC_steps/nchanges)

# function acf(dt , s0_array)
#     ar1 = copy(s0_array[1:end-dt])
#     ar2 = copy(s0_array[dt+1:end])
    
#     return (dot(ar1, ar2)/length(ar1)-mean(ar1)*mean(ar2))/(mean(s0_array.^2) - mean(s0_array)^2)
# end

# ti = 1:1000
# y = []
# for dt in ti
#     push!(y,acf(dt,s0_array))
# end

# model(t, p) = p[1] * exp.(-p[2] * t)
# p0 = [1,0.002]
# fit = curve_fit(model, ti, y, p0)

# p = plot(ti,y)
# plot!(p, ti, model(ti,fit.param), label="exponential fit")
# display(p)
# println("corrrelation time from the fit: ", 1/fit.param[2])
# -----------------------------------------------------------------------

function get_sub_neighbors(g)
    adj = g.adj
    nbr = []
    for i in 1:g.n
        temp = []
        for j in 1:g.n
            if g.adj[i,j] == 1
                push!(temp, j)
            end
        end
        push!(nbr, temp)
    end
    return convert_to_int(convert_to_int(nbr))
end


function get_sub_config(config)
    sub_config = fill(0, length(sub_spins))
    config = reset_config(max_match);
    for i in 1:length(sub_spins)
        sub_config[i] = config[sub_spins[i]]
    end
    return sub_config
end

basis_path_to_all_path(basis_path::Vector{Vector{Int64}}, pts::Vector{ComplexF64}) = [generate_symm_path(pth, pts) for pth in basis_path]



    