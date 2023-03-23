
function pre_defined_params()
    Δt = 0.0001
    Δx = 1
    T = 10
    X = 20
    aspect = 1
    boundary = "no_flux"
    simulation_iterations = 4
    parameter_iterations = 30
    λ₀,λ₁ = 0.017,0.045
    h = 1
    dif_min = 0.0049
    dif_max = 0.1
    λmax = 100
    hmax = 10
    amax = 4
    data_root = "data/"
    input_values = T6ssDataGenerateInput(
            Δt,Δx,T,X,aspect,boundary,
            simulation_iterations,parameter_iterations,
            λ₀,λ₁,h,dif_min,dif_max,λmax,hmax,amax,data_root)
    return input_values
end

function pre_defined_params(parameter_iterations)
    Δt = 0.0001
    Δx = 1
    T = 10
    X = 20
    aspect = 1
    boundary = "no_flux"
    simulation_iterations = 4
    parameter_iterations = parameter_iterations
    λ₀,λ₁ = 0.017,0.045
    h = 1
    dif_min = 0.0049
    dif_max = 0.1
    λmax = 100
    hmax = 10
    amax = 4
    data_root = "data/"
    input_values = T6ssDataGenerateInput(
            Δt,Δx,T,X,aspect,boundary,
            simulation_iterations,parameter_iterations,
            λ₀,λ₁,h,dif_min,dif_max,λmax,hmax,amax,data_root)
    return input_values
end

function data_path_filenames(input_values,data_numeric_label)
    @unpack Δt, Δx, T, X, aspect,boundary,simulation_iterations,parameter_iterations,λ₀,λ₁,h,dif_min,dif_max,λmax,hmax,amax,data_root = input_values
    iter = Iterations(simulation_iterations,parameter_iterations)
    parameter_span = "Span"
    iterations = lpad.(1:iter.parameter_iterations,2,"0")
    data_path_json_vec = [
        data_root*
        join_figure_number_letter(x,parameter_span)*
        "ParamIter"*
        string(z) 
            for x in data_numeric_label 
            for z in iterations] .* 
            ".json"
    return data_path_json_vec
end


function get_grouped_file_paths(input_values)

    grouped_file_paths = @chain input_values begin
            _.data_root
            find(_,1,2,"json") 
            sort(_)
        end
    return grouped_file_paths
end

function generate_iters_get_values() 
    iters_get_values = [
    (group = 1, data_type = "parameters", data_value = "λ₀", data_label = L"\lambda_0", data_measure = "",xaxis = L"d_T",yaxis = L"P(d_T)", dist_axis = L"\bar{d}_T"),
    (group = 2, data_type = "parameters", data_value = "h", data_label = L"h_0", data_measure = "",xaxis = L"d_T",yaxis = L"P(d_T)", dist_axis = L"\bar{d}_T"),
    (group = 3, data_type = "variables", data_value = "aspect", data_label = L"aspect", data_measure = "",xaxis = L"d_T",yaxis = L"P(d_T)", dist_axis = L"\bar{d}_T"),
    (group = 4, data_type = "parameters", data_value = "h", data_label = L"D_0 \> (\mu m^2/s)", data_measure = L"\mu m^2/s",xaxis = L"d_{T_{OC}}\>\mu m",yaxis = L"P(d_{T_{OC}})", dist_axis = L"\bar{d}_{T_{OC}}\>\mu m")]
return iters_get_values
end

function generate_telegraph_times_one_cycle(var,Δp)
    st_one = get_state_time_series(T6ssBiologicalInspired(),var,Δp)
    lt = length(st_one)
    time_vec = range(0,lt*var.Δt,lt)
    lines(time_vec,Int.(st_one))
    gi = get_indices(st_one,time_vec)

    gi |> typeof
    if gi.i₁ > gi.i₀
    push!(gi.i₁,length(time_vec))       
    elseif gi.i₀ > gi.i₁
        push!(gi.i₀,length(time_vec))       
    else
        #huh
    end
    gi = telegraph_points(gi.i₀,gi.i₁)
    tele_times = get_times(gi,time_vec,st_one[1])

    t₀ = filter(x->x!=0.0,tele_times.t₀)
    t₁ = filter(x->x!=0.0,tele_times.t₁)
    return telegraph_time_dist(t₀,t₁)
end

function generate_figure_3_data()

    input_values = pre_defined_params()
    @unpack Δt, Δx, T, X, aspect,boundary,simulation_iterations,parameter_iterations,λ₀,λ₁,h,dif_min,dif_max,λmax,hmax,amax,data_root = 
                input_values
    T = 200         
    Δt = 1E-2   

    aspect = 1.6
    var = Variables(Δt, Δx, T, X, aspect,boundary)
    par =  Parameters(λ₀,λ₁,h)
    Δp = Δ(par,var) 
    Y_length = Int(var.aspect*var.X)
    dom = Domain(var.X,Y_length)
    time_vec = t⃗(var)

    """
        Generate telegram and random walk
    """
    telegraph_walker = get_all_walker_position(var,Δp,dom)


    """
        Compute data for telegraph plot
    """
    st = Int.(telegraph_walker.states)
    indices_per_state = get_indices(st,time_vec)
    t₀ = time_vec[indices_per_state.i₀]
    t₁ = time_vec[indices_per_state.i₁]
    s₀ = st[indices_per_state.i₀]
    s₁ = st[indices_per_state.i₁]

    """
        Compute random walk time series
    """

    Δx = 1E-1
    pos = RandomWalker.position.(telegraph_walker.walk)
    x_pos = [p[1]*Δx for p in pos]
    y_pos = [p[2]*Δx for p in pos]
    rwx = x_pos[indices_per_state.i₀]
    rwy = y_pos[indices_per_state.i₁]



    tele_times = [generate_telegraph_times_one_cycle(var,Δp) for i in 1:10^simulation_iterations]
    not_firing_time = reduce(vcat,[t.t₀ for t in tele_times])
    firing_time = reduce(vcat,[t.t₁ for t in tele_times])





    data = (
        time_vec = time_vec,
        firing = indices_per_state.i₁,
        not_firing = indices_per_state.i₀,
        telegraph = (
            states = st,
            indices = 
                (t₀ = t₀,
                t₁ = t₁,
                s₀ = s₀,
                s₁=s₁),
            waiting_times = 
                (firing_time = firing_time,
                not_firing_time = not_firing_time),
            parameters = 
                (λ₀ = λ₀,
                λ₁ = λ₁)),
        walker = (
            x = x_pos,
            y = y_pos,
            indices = (
                rwx = rwx,
                rwy = rwy))
    )

    return data
end

function generate_figure_3_a(data)
    time_vec = data[:time_vec]
    st = data[:telegraph][:states]
    t₀ = data[:telegraph][:indices][:t₀]
    t₁ = data[:telegraph][:indices][:t₁]
    s₀ = data[:telegraph][:indices][:s₀]
    s₁ = data[:telegraph][:indices][:s₁]
    fire = colorant"#d4af37"
    not_fire = colorant"#2b50c8"
    line =colorant"#222a5f"

    figsa = Figure(resolution=(800,800))
        fontsize_theme = Theme(fontsize = 35)
        set_theme!(fontsize_theme)
        ax = Axis(figsa[1,1],
        width = 512, height = 512,aspect=1,
        xlabel = L"t",ylabel = L"s(t)",xlabelsize=35,ylabelsize=35)
        lines!(ax,time_vec,st,color=line)
        scatter!(ax,t₀,s₀,label = L"s=0",color=not_fire,markersize=20)
        scatter!(ax,t₁,s₁,label = L"s=1",color=fire,markersize=20)
        axislegend()

    return figsa
end

function generate_figure_3_b(data)
    firing_time = data[:telegraph][:waiting_times][:firing_time]
    not_firing_time = data[:telegraph][:waiting_times][:not_firing_time]
    λ₀ = data[:telegraph][:parameters][:λ₀]
    λ₁ = data[:telegraph][:parameters][:λ₁]
    fire = colorant"#d4af37"
    not_fire = colorant"#2b50c8"
    bins = 100
    largest_time=maximum([maximum(firing_time),maximum(not_firing_time)])
    exp_time_vec = range(0,largest_time,length=bins)
    λ₀_exp = generate_exponential.(λ₀,exp_time_vec)
    λ₁_exp = generate_exponential.(λ₁,exp_time_vec)
    figsb = Figure(resolution=(800,800))
        fontsize_theme = Theme(fontsize = 35)
        set_theme!(fontsize_theme)
        ax = Axis(figsb[1,1],
        width = 512, height = 512,aspect=1,
        xlabel = L"\tau",ylabel = L"P(\tau)",xlabelsize=35,ylabelsize=35)
        hist!(ax,not_firing_time,label=L"\lambda_0",bins=bins,normalization=:pdf,color=not_fire)
        hist!(ax,firing_time,label=L"\lambda_1",bins=bins,normalization=:pdf,color=fire)
        lines!(ax,exp_time_vec,λ₀_exp,color=:red)
        lines!(ax,exp_time_vec,λ₁_exp,color=:red)
        axislegend()

    return figsb
end

function generate_figure_3_c(data)
    time_vec = data[:time_vec]
    x_pos = data[:walker][:x]
    y_pos = data[:walker][:y]
    t₀ = data[:telegraph][:indices][:t₀]
    t₁ = data[:telegraph][:indices][:t₁]
    rwx = data[:walker][:indices][:rwx]
    rwy = data[:walker][:indices][:rwy]
    fire = colorant"#d4af37"
    not_fire = colorant"#2b50c8"
    figsc = Figure(resolution=(800,800))
        fontsize_theme = Theme(fontsize = 35)
        set_theme!(fontsize_theme)
        ax = Axis(figsc[1,1],
        width = 512, height = 512,aspect=1,
        xlabel = L"t",ylabel = "",xlabelsize=35,ylabelsize=35)
        lines!(ax,time_vec,x_pos,label = L"x(t)",color=not_fire)    
        lines!(ax,time_vec,y_pos,label = L"y(t)",color=fire)
        scatter!(ax,t₀,rwx,color=not_fire,markersize=20)
        scatter!(ax,t₁,rwy,color=fire,markersize=20)
        axislegend()
    return figsc
end

function generate_figure_3_d(data)
    x_pos = data[:walker][:x]
    y_pos = data[:walker][:y]
    firing = data[:firing]
    not_firing = data[:not_firing]
    fire = colorant"#d4af37"
    not_fire = colorant"#2b50c8"
    line =colorant"#222a5f"
    figsd = Figure(resolution=(800,800))
        fontsize_theme = Theme(fontsize = 35)
        set_theme!(fontsize_theme)
        ax = Axis(figsd[1,1],
        width = 512, height = 512,aspect=1,
        xlabel = L"x",ylabel = L"y",xlabelsize=35,ylabelsize=35)
        lines!(ax,x_pos,y_pos,label = L"RW",color=line)    
        scatter!(ax,x_pos[firing],y_pos[firing],color=fire,markersize=20)    
        scatter!(ax,x_pos[not_firing],y_pos[not_firing],color=not_fire,markersize=20)    
        axislegend()
    return figsd
end

function generate_figure_4_data(input_values)

    @unpack Δt, Δx, T, X, aspect,boundary,simulation_iterations,parameter_iterations,λ₀,λ₁,h,dif_min,dif_max,λmax,hmax,amax,data_root = input_values
    var = Variables(Δt, Δx, T, X, aspect,boundary)
    par =  Parameters(λ₀,λ₁,h)
    Δp = Δ(par,var)    
    dom = dimensions(var)
    iter = Iterations(simulation_iterations,parameter_iterations)
    λ₀ₘᵢₙ = 0.01
    λ₀ₘₐₓ = 10
    figure_number = 4
    λ₀_vec = range(λ₀ₘᵢₙ,λ₀ₘₐₓ,iter.parameter_iterations)
    λ₁_vec = ones(iter.parameter_iterations) .* λ₁
    parλ₀ = [Parameters(λ₀_vec[i],λ₁_vec[i],h) for i in 1:iter.parameter_iterations]
    Δparλ₀ = [Δ(p,var) for p in parλ₀]       



    data_path_json_vec = data_path_filenames(input_values,figure_number)


    function sol(var, Δp, dom, iter,data_path_json_vec)
        res = get_experimental_sample_dist_vecs(var, Δp, dom, iter)
        solu = SolutionVarParDom(var,Δp,dom,iter,res.experimental,res.sample)
        return write_solution_to_file(solu,data_path_json_vec)
    end

    for i in 1:iter.parameter_iterations
        sol(var, Δparλ₀[i], dom, iter,data_path_json_vec[i])
    end
end

function generate_figure_4_a(input_values)
    @unpack Δt, Δx, T, X, aspect,boundary,simulation_iterations,parameter_iterations,λ₀,λ₁,h,dif_min,dif_max,λmax,hmax,amax,data_root = 
    input_values
    vars = Variables(Δt, Δx, T, X, aspect,boundary)
    par =  Parameters(λ₀,λ₁,h)
    Δp = Δ(par,vars)
    iter = Iterations(simulation_iterations,parameter_iterations)
    λ₀_vec = range(0,.05,iter.parameter_iterations)
    parameter_colours = [colorant"#d471d1", colorant"#60dce5"]
    x = t⃗(vars)

    # Show two generic time series
    par =  map(x -> Parameters(x,λ₁,h),[λ₀_vec[10],λ₀_vec[end-1]])
    Δp = map(x -> Δ(x,vars),par)
    y = map(x->get_state_time_series(vars,x),Δp)
    size_x = size(x,1)
    xy12 = [(x[i],y[j][i]) for i in 1:size_x, j in 1:2]

    fig4a = Figure(resolution=(800,800))
    fontsize_theme = Theme(fontsize = 35)
    set_theme!(fontsize_theme)
    ga = fig4a[2, 1] = GridLayout()
    axtop = Axis(ga[1, 1],ylabel = L"s(t)",width = 512, height = 245,aspect=2)
    axbottom = Axis(ga[2, 1],
    width = 512, height = 256,aspect=2,
    xlabel = L"t",ylabel = L"s(t)",xlabelsize=35,ylabelsize=35)
    linkxaxes!(axbottom, axtop)
    lines!(axtop,xy12[:,1],color=parameter_colours[1])
    lines!(axbottom,xy12[:,2],color=parameter_colours[2])

    return fig4a
end

function generate_figure_4_b(iters_get_values)
    Δx = 0.1
    first_plot = "02"
    figure_4 = 4

    parameter_colours = [colorant"#d471d1", colorant"#60dce5"]

    # For everything except the last figure
    data_figsb = @chain figure_4 begin
        data_path_filenames(input_values,_)
        filter(x -> occursin("Iter$(first_plot).",x),_)[1]
        read_solution_from_memory(_,SolutionVarParDom)
    end


    data_figsb = @set data_figsb.experimental = data_figsb.experimental ./ (1/Δx)
    data_figsb = @set data_figsb.sample = data_figsb.sample ./ (1/Δx)
    data_figsb = @set data_figsb.domain = data_figsb.domain * Δx
    data_figsb = @set data_figsb.variables.Δx = data_figsb.variables.Δx ./ (1/Δx)


    figsb = view_distance_and_mean(data_figsb,iters_get_values[1],parameter_colours[1])
    return figsb
end

function generate_figure_4_c(input_values,iters_get_values)
    Δx = 0.1
    second_plot = lpad(input_values.parameter_iterations - 1,2,"0")
    figure_4 = 4

    parameter_colours = [colorant"#d471d1", colorant"#60dce5"]


    # For everything except the last figure
    data_figsc = @chain figure_4 begin
        data_path_filenames(input_values,_)
        filter(x -> occursin("Iter$(second_plot).",x),_)[1]
        read_solution_from_memory(_,SolutionVarParDom)
    end



    data_figsc = @set data_figsc.experimental = data_figsc.experimental ./ (1/Δx)
    data_figsc = @set data_figsc.sample = data_figsc.sample ./ (1/Δx)
    data_figsc = @set data_figsc.domain = data_figsc.domain * Δx
    data_figsc = @set data_figsc.variables.Δx = data_figsc.variables.Δx ./ (1/Δx)


    figsc = view_distance_and_mean(data_figsc,iters_get_values[1],parameter_colours[2])
    return figsc
end

function generate_figure_4_d(input_values,iters_get_values)

    Δx = 0.1
    λ₀ₘᵢₙ = 0.01
    λ₀ₘₐₓ = 10
    figure_number = 4
    λ₀_vec = range(λ₀ₘᵢₙ,λ₀ₘₐₓ,input_values.parameter_iterations)
    figsd = @chain figure_number begin
        data_path_filenames(input_values,_) 
        read_solution_from_memory.(_,SolutionVarParDom)
    end

    theoretical = get_t6ss_walk_theoretical_dist.(figsd) 
    [i.parameters.λ₀ for i in figsd] |> sort
    param = λ₀_vec 
    exper = [mean(i.experimental) for i in figsd] ./ (1/Δx)
    sampl = [mean(i.sample) for i in figsd] ./ (1/Δx)
    theor = theoretical ./ (Δx)
    distance = (parameter = param, experimental = exper, sample = sampl,theoretical = theor)        


    control_colour = colorant"#d5b670"  # Medium yellow 
    simulation_color = colorant"#443b28" # Sangria
    theoretical_color = colorant"#d92121" # Maximum red
    small_parameter_color = colorant"#d471d1" # Magenta
    large_parameter_color = colorant"#60dce5" # Light cyan


    x = distance[:parameter]
    l = length(x)-1
    x2 = [x[2]]
    x3 = [x[l]]
    yexp = distance[:experimental]
    yexp2 = [yexp[2]]
    yexp3 = [yexp[l]]
    ysam = distance[:sample]
    theoretical = distance[:theoretical]
    figsdscat = Figure(resolution=(800,800))
    fontsize_theme = Theme(fontsize = 35)
    set_theme!(fontsize_theme)
    ax = Axis(figsdscat[2,1],
        width = 512, height = 512,aspect=1,
        xlabel = iters_get_values[1][:data_label],ylabel = iters_get_values[1][:dist_axis], xlabelsize=35, ylabelsize=35)
    scatter!(ax,x,ysam,color=control_colour)
    scatter!(ax,x,yexp,color=simulation_color)
    lines!(ax,x,theoretical,color=theoretical_color)
    scatter!(ax,x2,yexp2,markersize=30,color=(small_parameter_color,0.5))
    scatter!(ax,x3,yexp3,markersize=30,alpha=.1,color=(large_parameter_color,0.5))

    return figsdscat
end

function generate_figure_5_data(input_values)
    
    @unpack Δt, Δx, T, X, aspect,boundary,simulation_iterations,parameter_iterations,λ₀,λ₁,h,dif_min,dif_max,λmax,hmax,amax,data_root = 
    input_values
    figure_number = 5
    
    
    
    Δt = 0.0001
    Δx = 1
    T = 10
    X = 20
    aspect = 1.6
    boundary = "no_flux"
    
    iter = Iterations(simulation_iterations,parameter_iterations) 
    parameter_iterations_save = iter.parameter_iterations+1
    data_path_json_vec = 
        data_path_filenames(
            pre_defined_params(parameter_iterations_save),
            figure_number) # sort different itrerations
    vars = Variables(Δt, Δx, T, X, aspect,boundary) # set variables to Variable struct
    dom = dimensions(vars) # get dimensions of domain stored to Domain struct
    
    Literature_scale = 2.7 # Ratio 1:2.7 of TssB:TssL
    
    λ₀ = 1/60
    λ₁ = λ₀*Literature_scale
    dif_min = 1E-3#0.0049
    dif_max = 0.1
    Δx = 0.1  # need to add in later
    D⃗ = @chain range(dif_min,dif_max,iter.parameter_iterations) begin
        vcat(_,0.0049)
        sort(_)    
    end # Set range of the diffusion
    
    h⃗ = round.(LitRate.(D⃗,Δx),digits=4) # Set range for h based on diffusion (D⃗) from literature
    par = map(h -> Parameters(λ₀,λ₁,h),h⃗)
    Δp = map(p -> Δ(p,vars),par) 
    
    
    
    # Update Δt according to Db and Δx
    function Δt_func(φ,Δx,Db,var::Variables)
        ϵ = 1e-6
        var = @set var.Δt = round((φ*Δx^2)/(Db) - ϵ,digits = 5)
        return var
    end


    φ = 1e-2
    vars = map(Db -> Δt_func(φ,Δx,Db,vars),D⃗)
    results = map(
        (v,p) -> 
            get_experimental_sample_dist_vecs(T6ssBiologicalInspired(),v, p, dom, iter), 
            vars, Δp) .* 
        Δx 
    vars = [@set i.Δx = Δx for i in vars] # Re-set Δx
    vars = [@set i.X = Int(X*Δx) for i in vars] # Re-set X
    dom = [dimensions(i) for i in vars] # Upate dimensions
    par = map(D -> Parameters(λ₀,λ₁,D),D⃗) # Update parameters
    
    sols = [
        SolutionVarParDom(
            vars[i],
            par[i],
            dom[i],
            iter,
            results[i].experimental,
            results[i].sample) 
            for i in 1:parameter_iterations_save] # Convert to struct # account for special D
    

    write_solution_to_file.(sols,data_path_json_vec)
end

function generate_figure_5_a(iters_get_values)

    first_plot = "02"
    figure_number = 5

    parameter_colours = [colorant"#d471d1", colorant"#60dce5"]

    figsa = @chain figure_number begin
        data_path_filenames(input_values,_)
        filter(x -> occursin("Iter$(first_plot).",x),_)[1]
        read_solution_from_memory(_,SolutionVarParDom)
        view_distance_and_mean(_,iters_get_values[4],parameter_colours[1])
    end

    return figsa

end

function generate_figure_5_b(input_values,iters_get_values)
    second_plot = lpad(input_values.parameter_iterations - 1,2,"0")
    figure_number = 5

    parameter_colours = [colorant"#d471d1", colorant"#60dce5"]

    figsb = @chain figure_number begin
        data_path_filenames(input_values,_)
        filter(x -> occursin("Iter$(second_plot).",x),_)[1]
        read_solution_from_memory(_,SolutionVarParDom)
        view_distance_and_mean(_,iters_get_values[4],parameter_colours[2])
    end

    return figsb

end

function generate_figure_5_c(iters_get_values)
    figure_number = 5
    first_D_value = 0.0049

    figsd_data = @chain figure_number begin
        data_path_filenames(input_values,_)
        read_solution_from_memory.(_,SolutionVarParDom)
    end


    D_param = @chain figsd_data begin
                map(x-> x.parameters.h,_)
        end

    l = length(D_param)
    first_point = findall(d -> d==first_D_value,D_param)
    second_point = l-3

    theoretical = get_t6ss_walk_theoretical_dist.(Ref(T6ssBiologicalInspired()),figsd_data)

    sols_plot = @chain figsd_data begin
            map(x -> (
                    get_value(
                    x,
                    iters_get_values[4][:data_type], 
                    iters_get_values[4][:data_value]),
                    mean(x.experimental),
                    mean(x.sample)),
                    _)
                    get_solutions_vector(_,theoretical)
    end


    control_colour = colorant"#d5b670"  # Medium yellow 
    simulation_color = colorant"#443b28" # Sangria
    theoretical_color = colorant"#d92121" # Maximum red
    small_parameter_color = colorant"#d471d1" # Magenta
    large_parameter_color = colorant"#60dce5" # Light cyan



    x2 = D_param[first_point]
    x3 = D_param[second_point]
    yexp = sols_plot[:experimental]
    yexp2 = yexp[first_point]
    yexp3 = yexp[second_point]

    ysam = sols_plot[:sample]
    theoretical = sols_plot[:theoretical]
    figsd = Figure(resolution=(800,800))
            fontsize_theme = Theme(fontsize = 35)
            set_theme!(fontsize_theme)
            ax = Axis(figsd[2,1],
            width = 512, height = 512,aspect=1,
            xlabel = iters_get_values[4][:data_label],ylabel = iters_get_values[4][:dist_axis],xlabelsize=35,ylabelsize=35)
            scatter!(ax,D_param,ysam,color=control_colour)
            scatter!(ax,D_param,yexp,color=simulation_color)
            lines!(ax,D_param,theoretical,color=theoretical_color)
            scatter!(ax,x2,yexp2,markersize=20,color=(small_parameter_color,0.5))
            scatter!(ax,x3,yexp3,markersize=20,alpha=.1,color=(large_parameter_color,0.5))

    return figsd
end



"""
Figure 3:
1. Time series of a single simulation
2. Waiting time histogram
3. Time series random walk in x and y
4. Random walk 
"""
function generate_figure_3()

    data = generate_figure_3_data()
    figsa = generate_figure_3_a(data)
    figsb = generate_figure_3_b(data)
    figsc = generate_figure_3_c(data)
    figsd = generate_figure_3_d(data)

    return (a = figsa,b = figsb,c=figsc,d=figsd)
end

"""
Figure 4:
1. Generic time series
2. Small value of λ₀
3. Large value of λ₀
4. Range of λ₀ and extpected distance travelled
"""
function generate_figure_4()
    input_values = pre_defined_params()
    iters_get_values = generate_iters_get_values()
    generate_figure_4_data(input_values)
    figsa = generate_figure_4_a(input_values)
    figsb = generate_figure_4_b(iters_get_values)
    figsc = generate_figure_4_c(input_values,iters_get_values)
    figsd = generate_figure_4_d(input_values,iters_get_values)


return (a = figsa,b = figsb,c=figsc,d=figsd)
end

"""
Figure 5:
1. Small value of D
2. Large value of D
3. Range of D and extpected distance travelled
"""
function generate_figure_5()
    parameter_iterations = 49
    input_values = pre_defined_params(parameter_iterations)
    iters_get_values = generate_iters_get_values()
    generate_figure_5_data(input_values)
    figsa = generate_figure_5_a(iters_get_values)
    figsb = generate_figure_5_b(input_values,iters_get_values)
    figsc = generate_figure_5_c(iters_get_values)

    return (a = figsa,b = figsb,c=figsc)
end