@unpack Δt, Δx, T, X, aspect,boundary,simulation_iterations,parameter_iterations,λ₀,λ₁,h,dif_min,dif_max,λmax,hmax,amax,data_root = 
            input_values
vars = Variables(Δt, Δx, T, X, aspect,boundary)
par =  Parameters(λ₀,λ₁,h)
Δp = Δ(par,vars)
dom = dimensions(vars)
iter = Iterations(simulation_iterations,parameter_iterations)
λ₀_vec = range(0,.05,iter.parameter_iterations)
h_vec = range(0,hmax,iter.parameter_iterations)
a_vec = range(0,amax,iter.parameter_iterations)
D_vec = range(dif_min,dif_max,iter.parameter_iterations) 
parameter_colours = [colorant"#d471d1", colorant"#60dce5"]

# Figure 3
# Show two time series 
x = t⃗(vars)
par =  map(x -> Parameters(x,λ₁,h),[λ₀_vec[10],λ₀_vec[end-1]])
Δp = map(x -> Δ(x,vars),par)
y = map(x->get_state_time_series(vars,x),Δp)
xy12 = [(x[i],y[j][i]) for i in 1:size(x,1), j in 1:2]

fig3a = Figure(resolution=(800,800))
fontsize_theme = Theme(fontsize = 35)
set_theme!(fontsize_theme)
ga = fig3a[2, 1] = GridLayout()
axtop = Axis(ga[1, 1],ylabel = L"s(t)",width = 512, height = 245,aspect=2)
axbottom = Axis(ga[2, 1],
    width = 512, height = 256,aspect=2,
    xlabel = L"t",ylabel = L"s(t)",xlabelsize=35,ylabelsize=35)
linkxaxes!(axbottom, axtop)
λ₀₁ = lines!(axtop,xy12[:,1],color=parameter_colours[1])
λ₀₂ = lines!(axbottom,xy12[:,2],color=parameter_colours[2])

fig3a


save_figure(fig3a,figure_root*"/Figure3/Figure3a.pdf")
#=
Legend(fig3a[1,1],
        [λ₀₁,λ₀₂],
        [L"\mathrm{Small} \> \lambda_0", L"\mathrm{Large} \> \lambda_0"], 
        iters_get_values[1][:data_label],
        orientation=:horizontal)=#


 




# Figure 4
# Two random walkers
par =  map(x -> Parameters(λ₀,λ₁,x),[h_vec[2],h_vec[end-1]])
Δp = map(x -> Δ(x,vars),par)
rw = @chain Δp begin
    map(x->get_all_walker_position(vars,x,dom),_)
    map(x->RandomWalker.position.(x.walk),_)
end
dim = (Int(-(1/2)*dom.width),Int((1/2)*dom.width),Int(-(1/2)*dom.length),Int((1/2)*dom.length))

fig4a = Figure(resolution=(800,800))
fontsize_theme = Theme(fontsize = 35)
set_theme!(fontsize_theme)
ga = fig4a[2, 1] = GridLayout()
ax = Axis(ga[1, 1],
    width = 512, height = 512,aspect=1,
    xlabel = L"x",ylabel = L"y",xlabelsize=35,ylabelsize=35,
    limits = dim)
rw₁ = lines!(ax,rw[1],color=parameter_colours[1],linewidth = 5)
rw₂ = lines!(ax,rw[2],color=parameter_colours[2])
#=Legend(fig4a[1,1],
        [rw₁,rw₂],
        [L"\mathrm{Small} \> h", L"\mathrm{Large} \> h"], 
        iters_get_values[2][:data_label],
        orientation=:horizontal)=#



    save_figure(fig4a,figure_root*"/Figure4/Figure4a.pdf")


        


# Figure 5
vars = map(x->Variables(Δt, Δx, T, X, x,boundary),[a_vec[2],a_vec[end-1]])
function dims(dom::Domain) 
    x₁ = Int(-(1/2)*dom.width)
    x₂ = Int((1/2)*dom.width)
    y₁ = Int(-(1/2)*dom.length)
    y₂ = Int((1/2)*dom.length)
    return Point2f[(x₁,y₁),(x₂,y₁),(x₂,y₂),(x₁,y₂)]
end
dom = map(x->dimensions(x),vars)
rect_points = map(dims,dom)

fig5a = Figure(resolution=(800,800))
fontsize_theme = Theme(fontsize = 35)
set_theme!(fontsize_theme)
ga = fig5a[2, 1] = GridLayout()
axleft = Axis(ga[1, 1],
    width = 256, height = 512,aspect=.5,
    xlabel = L"x",ylabel = L"y",xlabelsize=35,ylabelsize=35)
axright = Axis(ga[1, 2],
            width = 256, height = 512,aspect=.5,xlabel = L"x")
            #xlabel = L"x",ylabel = L"y",xlabelsize=30,ylabelsize=30)
linkyaxes!(axleft, axright)
p₁ = poly!(axleft,rect_points[1], color = parameter_colours[1], strokecolor = :black, strokewidth = 1)
p₂ = poly!(axright,rect_points[2], color = parameter_colours[2], strokecolor = :black, strokewidth = 1)
#=Legend(fig5a[1,1],
        [p₁,p₂],
        ["Small aspect", "Large aspect"], 
        iters_get_values[3][:data_label],
        orientation=:horizontal)=#


    save_figure(fig5a,figure_root*"/Figure5/Figure5a.pdf")