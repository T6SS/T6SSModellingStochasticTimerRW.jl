# Parameters

xwidth = 20
aspect = 1
ylength = dimensions(xwidth,aspect = aspect)[:ylength] 
Δx = 0.1
Δt = 0.001
duration = 60
t⃗ = range(0,duration,step = Δt)




λ₀ = 1/20
λ₁ = 1/20
h = 10
λ₀Δt = XΔt(λ₀,Δt)
λ₁Δt = XΔt(λ₁,Δt)
hΔt = XΔt(h,Δt)
boundary = "no_flux"
N = 3
x₀ = 0
y₀ = 0
walker = Walker2D(x₀,y₀)



fig_num = 5
fig_sub_num = ("a","b","c")
order = ifelse(fig_sub_num == ("a","b","c"),true,false)



fig_folder = "/home/pushinglimits/Projects/Dundee_PhD/documentation/ModellingStochasticTimerRandomMotionT6SS/images/FiguresPaper/"
data_root = "/home/pushinglimits/Projects/Dundee_PhD/data/02_telegraph_data/"







