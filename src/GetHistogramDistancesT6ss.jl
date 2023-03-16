"""
Take type `SolutionVarParDom` and returns a `Figure`. The figure is a bar plot of sample and experimental simulations as well as the means.
""" 
# add color option for distinguishing which plot I care about
function view_distance_and_mean(data::SolutionVarParDom,label::NamedTuple,simulation_color)
    normal(y) = y ./ sum(y)
    control_colour = colorant"#d5b670"  # Medium yellow 
    control_mean_colour = colorant"#005432" # Green Cyan
    simulation_mean_colour = colorant"#92000a" # Sangria
    
    x = set_range_iter(data.domain,data.variables)
    yexp = normal(get_count(data.experimental,x))
    ysam = normal(get_count(data.sample,x))
    fig = Figure(resolution=(800,800))
    fontsize_theme = Theme(fontsize = 35)
    set_theme!(fontsize_theme)
    ax = Axis(fig[2,1],
        width = 512, height = 512,aspect=1,
        xlabel = label[:xaxis],ylabel = label[:yaxis],xlabelsize=35,ylabelsize=35)
    barplot!(ax,x,ysam,color=control_colour)
    vlines!(ax,mean(data.sample),color = control_colour)
    barplot!(ax,x,yexp,color=simulation_color)
    vlines!(ax,mean(data.experimental),color=simulation_color)
    #Legend(fig[1,1],[[control_barplot,control_mean_vline],[simulation_barplot,simulation_mean_barplot]],["Control","Simulation"], label[:data_label],orientation=:horizontal)
    return fig
end


"""
Take type `SolutionVarParDom` and returns a `Figure`. The figure is a bar plot of sample and experimental simulations as well as the means.
""" 
# add color option for distinguishing which plot I care about
function view_mean_heatmap(data::SolutionVarParDom,label::NamedTuple,simulation_color)
    normal(y) = y ./ sum(y)
    control_colour = colorant"#d5b670"  # Medium yellow 
    control_mean_colour = colorant"#005432" # Green Cyan
    simulation_mean_colour = colorant"#92000a" # Sangria
    
    x = set_range_iter(data.domain,data.variables)
    yexp = normal(get_count(data.experimental,x))
    ysam = normal(get_count(data.sample,x))
    fig = Figure(resolution=(800,800))
    fontsize_theme = Theme(fontsize = 35)
    set_theme!(fontsize_theme)
    ax = Axis(fig[2,1],
        width = 512, height = 512,aspect=1,
        xlabel = L"d",ylabel = L"P(d)",xlabelsize=35,ylabelsize=35)
    control_barplot = barplot!(ax,x,ysam,color=control_colour)
    control_mean_vline = vlines!(ax,mean(data.sample),color = control_mean_colour)
    simulation_barplot = barplot!(ax,x,yexp,color=simulation_color)
    simulation_mean_barplot = vlines!(ax,mean(data.experimental),color = simulation_mean_colour)
    Legend(fig[1,1],[[control_barplot,control_mean_vline],[simulation_barplot,simulation_mean_barplot]],["Control","Simulation"], label[:data_label],orientation=:horizontal)
    return fig
end





"""
Take a `NamedTuple` of distances to plot scatter.
"""
function view_scatter_mean(distance::NamedTuple,label::NamedTuple)
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
    theoretical = [i[2] for i in distance[:theoretical]]
    fig = Figure(resolution=(800,800))
    fontsize_theme = Theme(fontsize = 35)
    set_theme!(fontsize_theme)
    ax = Axis(fig[2,1],
        width = 512, height = 512,aspect=1,
        xlabel = label[:data_label],ylabel = label[:dist_axis], xlabelsize=35, ylabelsize=35)
    scatter!(ax,x,ysam,color=control_colour)
    scatter!(ax,x,yexp,color=simulation_color)
    lines!(ax,x,theoretical,color=theoretical_color)
    scatter!(ax,x2,yexp2,markersize=30,color=(small_parameter_color,0.5))
    scatter!(ax,x3,yexp3,markersize=30,alpha=.1,color=(large_parameter_color,0.5))
    #=
    Legend(fig[1,1],\
            [control_scatter,simulation_scatter,theoretical_scatter],
            ["Control","Simulation","Theoretical"],
            label[:data_label], 
            orientation=:horizontal)
            =#
    return fig    
end


"""
Takes type `Figure` and `String` and saves a figre to the path, which is the string.
"""
function save_figure(figure::Figure,path::String)
    save(path,figure,pt_per_unit = 1)
end





