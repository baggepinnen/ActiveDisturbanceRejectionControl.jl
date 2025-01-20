
cd(@__DIR__)
using Pkg
Pkg.activate(".")

using Symbolics
# NOTE: extreme hacks to workaround https://github.com/JuliaSymbolics/Symbolics.jl/issues/1404
Base.zero(::Type{Any}) = Num(0)
Base.one(::Type{Any}) = Num(1)
Base.oneunit(::Type{Any}) = Num(1)
Base.inv(A::Matrix{Any}) = inv(identity.(A))
Base.convert(::Type{T}, x::AbstractArray{Num}) where T <: Array{Num} = T(map(Num, x)) # https://github.com/JuliaSymbolics/Symbolics.jl/issues/1405

using ControlSystemsBase, Plots, RobustAndOptimalControl, Test, LinearAlgebra
default(margin=4Plots.mm, l=3, titlefontsize=12)
t = 0:0.002:15

# The plant model used in experiments
K = 1
D = 1
T = 1
P = tf([K],[T^2, 2D*T, 1])

"""
    adrc(Tsettle, ogain)

Construct an ADRC controller with settling time `Tsettle` and observer poles that are `ogain` times faster than the closed-loop poles. The returned controller has two inputs, `r` and `y`, and one output, `u`.
"""
function adrc(Tsettle, ogain, b0 = 1)
    T = promote_type(typeof(Tsettle), typeof(ogain))
    Pdes1 = ss(tf([one(T)],[1, 0, 0]))*b0
    Pdes2 = add_low_frequency_disturbance(Pdes1, 2)
    sCL = -6/Tsettle

    Kp = sCL^2
    Kd = -2sCL
    sESO = ogain*sCL
    k1 = -3sESO
    k2 = 3sESO^2
    k3 = -sESO^3
    K = T[k1; k2; k3;;]

    obs = observer_filter(Pdes2, K, output_state=true)

    obsn = named_ss(ss(diagm([one(T), Kd, one(T)]))*obs, u=[:u, :y], y=[:yh, :yh2, :fh])
    Kpn = named_ss(ss(Kp), u=:e, y=:up)
    b0n = named_ss(ss(1/b0), u=:ub0, y=:u)
    s1 = sumblock("e = r-yh")
    s2 = sumblock("u0 = up - yh2")
    sd = sumblock("ub0 = u0 - fh")

    systems = [
        obsn,
        Kpn,
        b0n,
        s1,
        s2,
        sd,
    ]
    connections = [
        :e => :e
        :up => :up
        :u0 => :u0
        :ub0 => :ub0
        :u => :u
        :fh => :fh
        :yh => :yh
        :yh2 => :yh2
    ]
    
    external_inputs = [:r, :y]
    external_outputs = [:u]
    
    connect(systems, connections; external_inputs, external_outputs, unique=true) #* diagm([1, -1])
end

function pid_2dof_2filt(kp, ki, kd, Tf, d, b, c)
    # r through filter
    tempA = [0 0 1; 1 0 0; -1 / (Tf^2) 0 (-2d) / Tf]
    tempB = [0 0; 1 0; c / (Tf^2) -1 / (Tf^2)]
    tempC = [kp ki kd]
    tempD = [b*kp 0]

    named_ss(ss(tempA, tempB, tempC, tempD), "PID", u=[:r, :y], y=:u)
end
filt(Tf, d) = tf(1, [Tf^2, 2d*Tf, 1])

##

"""
    equivalent_pid(Tsettle, ogain; simplified_r = true)

Construct a PID controller that is equivalent to the ADRC controller with settling time `Tsettle` and observer poles that are `ogain` times faster than the closed-loop poles. The returned controller has two inputs, `r` and `y`, and one output, `u`.

If `simplified_r` is true, the controller is a PI controller with set-point weighting on the proportional term and a first-order lowpass filter on the measurement. If `simplified_r` is false, the controller exactly matches the ADRC contorller, which is a filtered PID controller from the reference signal.
"""
function equivalent_pid(Tsettle, ogain, b0=1; simplified_r = true)
    # The expressions for the PID parameters are obtained symbolically in the bottom of the script
    kpy = (72.0*ogain^3 + 108.0*ogain^2)/(3.0*Tsettle^2*ogain^2 + 6.0*Tsettle^2*ogain + Tsettle^2)
    kiy = 216.0*ogain^3/(3.0*Tsettle^3*ogain^2 + 6.0*Tsettle^3*ogain + Tsettle^3)    
    kdy = (6.0*ogain^3 + 36.0*ogain^2 + 18.0*ogain)/(3.0*Tsettle*ogain^2 + 6.0*Tsettle*ogain + Tsettle)    
    Ty = -0.166666666666667*Tsettle*sqrt(1/(3.0*ogain^2 + 6.0*ogain + 1.0))    
    dy = -0.5*(3.0*ogain + 2.0)*sqrt(1/(3.0*ogain^2 + 6.0*ogain + 1.0))  
  
    b = (6/Tsettle)^2 / kpy # Found by matching asymptotes
    C = if simplified_r
        pid_2dof_2filt(kpy, kiy, kdy, Ty, dy, b, 0)
    else
        error("Not implemented")
        # Cpidr = (b*kpy + kiy/s)
        # Cpidy = (kpy + kiy/s + kdy*s)/(s^2*Ty^2 + 2*dy*s*Ty + 1)
        # [Cpidr Cpidy]
    end
    named_ss(1/b0 * ss(C) * diagm([1, 1]), y=:u, u=[:r, :y])
end

Tsettle = 5 # Parameters suggested in the paper
ogain = 10
b0 = 1
Ca = adrc(Tsettle, ogain, b0)
Cr = Ca[:u,:r]
Cy = Ca[:u,:y]

s = tf('s')
Ki = 0.6
TZ1 = TZ2 = T
T1 = 0.2
C_suggested_pid = Ki * (1+TZ1*s)*(1 + TZ2*s) / (s*(1 + T1*s))
C_equivalent_pid = equivalent_pid(Tsettle, ogain, b0; simplified_r=true)
label = ["ADRC" "Suggested PID" "Equivalent PIDF"]

w = exp10.(LinRange(-2, 3, 200))
gangoffourplot(P, [Ca[:u,:y], C_suggested_pid, C_equivalent_pid[:u,:y]]; label)

##
feedback2d(P, Ca) = feedback(P, Ca[:u, :y], pos_feedback=true)*(Ca[:u, :r])


## Analysis
# The analysis below compares the ADRC controller to the in the paper suggested PID controller, as well as the approximately equivalent PI controller with set-point weighting

# ADRC controller has overall much higher gain
# It looks like a PI controller from r, and a filtered PID controller from y
# Overall, it has a much higher gain but rolloff from measurements
bodeplot([Ca, [C_suggested_pid -C_suggested_pid], C_equivalent_pid], w; label=repeat(label, inner=(1,4)), background_color_legend=nothing, foreground_color_legend=nothing, title=["\$C_{ur}\$" "\$C_{uy}\$" "" ""], legend=[true false false false], linestyle=repeat([:solid :solid :dash], inner=(1,4)))
savefig("paper/figures/second_order_bode_C.pdf")

# Step responses from r are almost identical
plot(step.([feedback2d(P, Ca), feedback(P*C_suggested_pid), feedback2d(P, C_equivalent_pid)], Ref(t)); label)

# So are closed-loop tf from r -> y 
bodeplot([feedback2d(P, Ca), feedback(P*C_suggested_pid), feedback2d(P, C_equivalent_pid)], w; label=repeat(label, inner=(1,2)), title="\$G_{yr}\$", linestyle=repeat([:solid :solid :dash], inner=(1,2)))
savefig("paper/figures/second_order_bode_ry.pdf")

# Bode plots from y -> u look different, ADRC is tuned much more aggressively but uses rolloff. The equivalent PID is identical to ADRC
bodeplot([G_CS(P, -Ca[:u, :y]), G_CS(P, C_suggested_pid), G_CS(P, -C_equivalent_pid[:u, :y])]; label=repeat(label, inner=(1,2)), title="\$G_{yu}\$", legend=[:topleft false], background_color_legend=nothing, foreground_color_legend=nothing, linestyle=repeat([:solid :solid :dash], inner=(1,2)), l=3)
savefig("paper/figures/second_order_bode_uy.pdf")

## Reproduce response-plots from 
using MonteCarloMeasurements
K = 1
T = 1
Ku = Particles([0.1, 0.2, 0.5, 1, 2, 5])
Tu = Particles([0.1, 0.2, 0.5, 1, 2, 5])
Pu = tf([Ku],[1, T])
plot(step.([
    feedback2d(Pu, Ca),
    feedback(Pu*C_suggested_pid),
    feedback2d(Pu, C_equivalent_pid)
], Ref(t)); label, ri=false, layout=(1,3), sp=(1:3)', size=(800,400), ylabel="y", c=[1 2 3])
savefig("paper/figures/second_order_K.pdf")

Pu = tf([K],[1, Tu])
plot(step.([
    feedback2d(Pu, Ca),
    feedback(Pu*C_suggested_pid),
    feedback2d(Pu, C_equivalent_pid)
], Ref(t)); label, ri=false, layout=(1,3), sp=(1:3)', size=(800,400), ylabel="y", c=[1 2 3])
savefig("paper/figures/second_order_T.pdf")


## To figure this out, we propagate symbolic variables through the adrc constructor
using Symbolics, SymbolicControlSystems
@variables Tsettle ogain b0
K = adrc(Tsettle, ogain)
K = ss(identity.(K.A), identity.(K.B), identity.(K.C), identity.(K.D))
ex = [to_num(K[1,1]) to_num(K[1,2])]
tf.(ex)
# We then simply match the coefficients for each order of s, obtaining a system of equations that we can solve for the PID parameters
# The PID controllers on symbolic form are constructed below


using SymbolicControlSystems
using SymbolicControlSystems: @syms
Base.active_repl.options.hint_tab_completes = false # This messes with sympy https://discourse.julialang.org/t/sympy-makes-repl-to-stuck/124814/6
@syms kp ki kd Tf d Tsettle ogain b c
# Cy

tf_des = tf((kp + ki/s + kd*s)/(s^2*Tf^2 + 2*d*s*Tf + 1))
tf_act = tf(symbolics_to_sympy(Symbolics.simplify(ex[2], expand=true)))

numvec1 = numvec(tf_des)[]
numvec2 = numvec(tf_act)[]
denvec1 = denvec(tf_des)[][1:end-1] 
denvec2 = denvec(tf_act)[][1:end-1] 


eqs = [
    numvec1./denvec1[end] - numvec2./denvec2[end],
    denvec1[1:end-1]./denvec1[end]  - denvec2[1:end-1]./denvec2[end] # Make the last element 1
]
vars = [kp, ki, kd, Tf, d]

sols = sp.solve(eqs, vars)

soli = 2
@show kpy = sols[soli][1]
@show kiy = sols[soli][2]
@show kdy = sols[soli][3]
@show Ty =  sols[soli][4]
@show dy =  sols[soli][5]

## Latex representation of the state-space realization of the equivalent PID controller
@syms k_p k_i k_d T_f d b
C_equiv_sym = pid_2dof_2filt(k_p, k_i, k_d, T_f, d, b, 0)
function show_latex_ss(sys::AbstractStateSpace)
    A,B,C,D = to_latex.(ssdata(sys))
    println("\\begin{align}")
    println("\\dot{x} &= $(A)x + $(B)u \\\\")
    println("y &= $(C)x + $(D)u")
    println("\\end{align}")
end

show_latex_ss(C_equiv_sym)

## Figure out a state-space realization of the equivalent PID controller
# The result of these symbolic computations are palced inside the function pid_2dof_2filt above
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using RobustAndOptimalControl


function assemblesys(;name)
    t = ModelingToolkit.t_nounits
    D = ModelingToolkit.D_nounits
    pars = @parameters kp ki kd T d b c
    vars = @variables y(t) r(t) u(t) e(t) fy(t) fu(t) ei(t)
    eqs = [
        u ~ kp*(b*r + fy) + ki*ei + kd*D(fy) # The filtered measurement is already negated
        D(ei) ~ r + fy # Integral term
        fu ~ c*r - y 
        T^2*D(D(fy)) + 2*d*T*D(fy) + fy ~ fu # Lowpass filter
    ]
    ODESystem(eqs, t, vars, pars; name)
end

@named model = assemblesys()
cm = complete(model)
mats, ssys = ModelingToolkit.linearize_symbolic(model, [cm.r, cm.y], [cm.u])
RobustAndOptimalControl.show_construction(ss(mats.A, mats.B, mats.C, mats.D))

## Get nice expressions for the state-space realization
sp.latex(sp.simplify.(symbolics_to_sympy.(mats.A)))
sp.latex(sp.simplify.(symbolics_to_sympy.(mats.B)))
sp.latex(sp.simplify.(symbolics_to_sympy.(mats.C)))
sp.latex(sp.simplify.(symbolics_to_sympy.(mats.D)))

## Get nice expressions for the state-space realization in original parameters
@syms T_s g
sp.latex(sp.simplify.(equivalent_pid(T_s, g, simplified_r=true).A))
sp.latex(sp.simplify.(equivalent_pid(T_s, g, simplified_r=true).B))
sp.latex(sp.simplify.(equivalent_pid(T_s, g, simplified_r=true).C))
sp.latex(sp.simplify.(equivalent_pid(T_s, g, simplified_r=true).D))

##

function gangofsevenplot(P, C, F, args...; c, name="", kwargs...)
    S,D,CS,T = gangoffour(P,C)
    RY = T*F
    RU = CS*F
    RE = S*F
    bodeplot!(S, args...; show=false, title="\$S = 1/(1+PC)\$", lab="$name: \$S\$", c, sp=1, plotphase=false, legend=:bottomright, kwargs...)
    bodeplot!(D, args...; show=false, title="\$PS = P/(1+PC)\$", lab="$name: \$PS\$", c, sp=2, plotphase=false, legend=:bottom, kwargs...)
    bodeplot!(CS, args...; show=false, title="\$CS = C/(1+PC)\$", lab="$name: \$CS\$", c, sp=3, plotphase=false, legend=:topleft, kwargs...)
    bodeplot!(T, args...; show=false, title="\$T = PC/(1+PC)\$", lab="$name: \$T\$", c, sp=4, plotphase=false, legend=:bottomleft, kwargs...)
    Plots.hline!([1], l=(:black, :dash, 1), primary=false, sp=4)
    bodeplot!(RE, args...; show=false, title="\$S = 1/(1+PC)\$", lab="$name: \$SF\$", l=(:dash,), c, sp=1, plotphase=false, kwargs...)
    bodeplot!(RY, args...; show=false, title="\$T = PC/(1+PC)\$", lab="$name: \$TF = r\\to y\$", l=(:dash,), c, sp=4, plotphase=false, kwargs...)
    bodeplot!(RU, args...; show=false, title="\$CS = C/(1+PC)\$", lab="$name: \$CSF\$", l=(:dash,), c, sp=3, plotphase=false, kwargs...)
end
default(titlefontsize=14, legendfontsize=8)
w = exp10.(LinRange(-2, 4, 200))
F = tf(Cr) / tf(-Cy) # This computes the equivalent reference prefilter appearing before the error calculation
plot(; layout=4, ticks=:default, xscale=:log10, size=(1200,700))#, link=:both)
gangofsevenplot(P, -tf(Cy), F, w; name="ADRC", c=1, background_color_legend=nothing, foreground_color_legend=nothing)
gangofsevenplot(P, C_suggested_pid, tf(1), w; name="Suggested PID", label, c=2, background_color_legend=nothing, foreground_color_legend=nothing)

F_equivalent_pid = tf(C_equivalent_pid[:u,:r]) / tf(-C_equivalent_pid[:u,:y])
gangofsevenplot(P, -tf(C_equivalent_pid[:u,:y]), F_equivalent_pid, w; name="Equivalent PID", c=3, background_color_legend=nothing, foreground_color_legend=nothing, linestyle=:dot)
savefig("paper/figures/second_order_7.pdf")