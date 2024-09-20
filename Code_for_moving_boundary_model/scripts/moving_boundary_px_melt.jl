using OrdinaryDiffEq
using DiffEqCallbacks
using Symbolics
using LinearSolve
using Plots
using Parameters
using Plots:mm
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

cd(@__DIR__)

"""
    semi_discretisation!(du, u, p, t)

Semi-discretisation function for the diffusion equation using finite differences.

# Arguments
- `du::Vector{Float64}`: The derivative of `u` with respect to time, to be updated in-place.
- `u::Vector{Float64}`: The current state vector.
- `p::Tuple`: A tuple containing the parameters:
    - `D_all::Dict{Symbol, Float64}`: A dictionary of diffusion coefficients for different elements.
    - `R::Float64`: The reaction rate.
    - `Δx::Float64`: The spatial step size.
    - `C_pyro::Float64`: The pyroxene concentration.
    - `C_factor::Float64`: A scaling factor for the pyroxene concentration.
    - `el_name::Symbol`: The name of the element for which the diffusion coefficient is to be used.
- `t::Float64`: The current time (not used in this function but typically required for time-stepping functions).

# Description
This function computes the semi-discretisation of the diffusion equation using finite differences. It updates the `du` vector in-place based on the current state `u` and the parameters `p`.

The function handles both the inner points and the boundary points separately. For the inner points, it uses a finite difference scheme to compute the second spatial derivative and the first spatial derivative. For the boundary points, it applies specific boundary conditions.
"""
function semi_discretisation!(du, u, p, t)

    # unpack the parameters
    D_all, R, Δx, C_pyro, C_factor, el_name = p

    # get the diffusion coefficient
    D = D_all[el_name]

    # compute the inner points
    @inbounds for i in 2:length(u)-1
        du[i] = D * (u[i-1] - 2*u[i] + u[i+1]) / (Δx^2) + R * (u[i+1] - u[i]) / Δx
    end

    # compute the boundary points
    C_left = - R * 2 * Δx / D * (C_pyro* C_factor - u[1]) + u[2]
    du[1] = D * (C_left - 2*u[1] + u[2]) / (Δx^2) + R * (u[2] - u[1]) / Δx
    du[end] = 0.0
end

function main(el_name, tmax, time_unit, length_corona, C_factor_melt, C_pyro)

    if time_unit == "days"
        tmax = tmax
    elseif time_unit == "years"
        tmax = tmax * 365.25
    end

    tdis = tmax  # time of dissolution, in days

    R = - length_corona / tdis  # Dissolution rate in μm/day, depends on total time of dissolution
    C_melt = 0.0  # concentration in the melt (ug/g)

    # create named tuple with the elements and their diffusion coefficients
    D_all = (
        Lu = 3246.9,  #  diffusion coefficient in μm²/day (value from Jiahui)
        Yb = 5021.0,  #  diffusion coefficient in μm²/day (value from Jiahui)
        Er = 3595.7,  #  diffusion coefficient in μm²/day (value from Jiahui)
        Ho = 5248.6,  #  diffusion coefficient in μm²/day (value from Jiahui)
        Dy = 5580.6,  #  diffusion coefficient in μm²/day (value from Jiahui)
        Y = 6011.4,  #  diffusion coefficient in μm²/day (value from Jiahui)
    )

    # print all the parameters
    println("Start new model...")

    # print tmax
    println("tmax = $tmax $(time_unit)")

    # print R, C_pyro, C_melt, D
    println("R = $R μm/$(time_unit)")
    println("Initial concentration in the pyroxene = $C_pyro")
    println("D($(el_name)) = $(D_all[el_name]) μm²/$(time_unit)")

    u0 = zeros(nx) .+ C_melt  # initial concentration in the melt

    # define the parameters of the model
    p = (D_all, R, Δx, C_pyro, C_factor_melt, el_name)

    # define the time span
    tspan = (0.0, tmax)

    du0 = copy(u0)
    # calculate jacobian sparsity pattern with symbolics
    jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> semi_discretisation!(du, u, p, 0.0), du0, u0)

    f = ODEFunction(semi_discretisation!; jac_prototype = float.(jac_sparsity))

    # define the ODE problem
    prob = ODEProblem(f, u0, tspan, p)

    # solve the ODE problem
    @time sol = solve(prob, TRBDF2(linsolve=KLUFactorization()), abstol=1e-10, reltol=1e-8, progress=true); #, callback=update_R_factor);

    println("...Model finished!")

    return sol, R
end

# define the model
const lx = 10000.0  # in μm (garnet + melt domain in the thin section)
const nx = 3000  # number of nodes
const Δx = lx / nx  # in μm
const x = range(Δx/2, stop=lx-Δx/2, length=nx)  # in μm, with cell centered positions

# define the initial condition
# size of the pyroxene in μm
# it is half of the real size as it is being dissolved from both sides

length_corona = 400.0  # in μm
length_garnet = 1700.0  # in μm

C_factor_melt = 3.09  # factor to multiply the concentration in the melt, based on partitional calculations
C_pyro = 0.30  # estimated concentration in the igneous pyroxene (μg/g)
tmax = 1
time_unit = "days"
sol_1_day, R_1_day = main(:Lu, tmax, time_unit, length_corona, C_factor_melt, C_pyro);
tmax = 10
time_unit = "days"
sol_10_day, R_10_day = main(:Lu, tmax, time_unit, length_corona, C_factor_melt, C_pyro);
tmax = 100
time_unit = "days"
sol_100_day, R_100_day = main(:Lu, tmax, time_unit, length_corona, C_factor_melt, C_pyro);
tmax = 1
time_unit = "years"
sol_1_yrs, R_1_yrs = main(:Lu, tmax, time_unit, length_corona, C_factor_melt, C_pyro);
tmax = 10
time_unit = "years"
sol_10_yrs, R_10_yrs = main(:Lu, tmax, time_unit, length_corona, C_factor_melt, C_pyro);
tmax = 100
time_unit = "years"
sol_100_yrs, R_100_yrs = main(:Lu, tmax, time_unit, length_corona, C_factor_melt, C_pyro);

# find index corresponding to length=3000 μm
idx = findfirst(x -> x >= 3000, x)

# start plotting
layout = @layout [a b]

plt1 = plot(legend=:topright, xlabel="Distance from the pyroxene-melt interface [μm]", ylabel="Lu mass fraction [μg/g]", ylim=(-0.005, 2.29), dpi=200, left_margin = 5mm)

# calculate distance from interface
grt_dist = 0
length_pyro = 400.0  # in μm

# plot initial conditions
plt1 = plot!(plt1, x, sol_1_day(0)[:], label="Melt", linewidth=1.5)

# add pyroxene position
plt1 = plot!(plt1, [(grt_dist-length_pyro+Δx/2), Δx/2], [C_pyro, C_pyro], label="Pyroxene ", color=:green, linetype=:steppost, xlim=((grt_dist-length_pyro+Δx/2), x[idx]), ylim=(-0.005, 1), linewidth=2)
# add vertical line at the initial pyroxene-melt interface
plt1 = vline!(plt1, [grt_dist], label="Initial pyroxene-melt\ninterface", color=:red, linestyle=:dot, linewidth=1.5)

colors_plot = colormap("Blues",8)

plt2 = plot(legend=:topright, xlabel="Distance from the pyroxene-melt interface [μm]", ylabel="", ylim=(-0.005, 1), dpi=200, label="Melt", xlim=(0, x[idx]), left_margin = 0mm,right_margin=5mm)
plt2 = plot!(x, sol_1_day[end], label="Duration = 1 day,\nR = $(R_1_day) μm/day", linewidth=2, color = colors_plot[8])
plt2 = plot!(x, sol_10_day[end], label="Duration = 10 days,\nR = $(R_10_day) μm/day", linewidth=2, color = colors_plot[7])
plt2 = plot!(x, sol_100_day[end], label="Duration = 100 days,\nR = $(R_100_day) μm/day", linewidth=2, color = colors_plot[6])
plt2 = plot!(x, sol_1_yrs[end], label="Duration = 1 year,\nR = $(round(R_1_yrs;digits=2)) μm/year", linewidth=2, color = colors_plot[5])
plt2 = plot!(x, sol_10_yrs[end], label="Duration = 10 years,\nR = $(round(R_10_yrs;digits=2)) μm/year", linewidth=2, color = colors_plot[4])
plt2 = plot!(x, sol_100_yrs[end], label="Duration = 100 years,\nR = $(round(R_100_yrs;digits=2)) μm/year", linewidth=2, color = colors_plot[3])

grt_dist = -sol_1_day.t[end]*R_1_day

plt2 = vline!(plt2, [grt_dist], label="", color=:red, linestyle=:dot, linewidth=1.5)

plt = plot(plt1, plt2, layout=layout,title = ["A" "B"], titlelocation = :left, guidefont=7, plot_title=" ", legendfontsize=6, ylims=(0,1))

display(plt)

# save fig
savefig(plt, "Figure_9_paper.png")
savefig(plt, "Figure_9_paper.pdf")

