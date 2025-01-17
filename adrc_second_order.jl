
cd(@__DIR__)
using Pkg
Pkg.activate(".")

using Symbolics
# NOTE: extreme hack to workaround https://github.com/JuliaSymbolics/Symbolics.jl/issues/1404
Base.zero(::Type{Any}) = Num(0)
Base.one(::Type{Any}) = Num(1)
Base.convert(::Type{T}, x::AbstractArray{Num}) where T <: Array{Num} = T(map(Num, x)) # https://github.com/JuliaSymbolics/Symbolics.jl/issues/1405

using ControlSystemsBase, Plots, RobustAndOptimalControl, Test, LinearAlgebra
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
    # L = [Kp/b0 1/b0]

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

"""
Convert from ``1 / ((sT_f)^2/(4d^2) + sT_f + 1)`` to ``ω^2 / (s^2 + 2ζω s + ω^2)``
"""
function Tfd2wd(Tf, d)
    0 ≤ d ≤ 2 || error("The damping ratio must be between 0 and 2, but got $d")
    ω = 2d/Tf
    ζ = d
    ω, ζ
end

function wd2Tfd(w, d)
    Tf = 2d/w
    Tf, d
end

function pid_2dof_2filt(kp, ki, kd, T, d, b, c)
    # r through filter
    tempA = [0 0 1; 1 0 0; -1 / (T^2) 0 (-2d) / T]
    tempB = [0 0; 1 0; c / (T^2) -1 / (T^2)]
    tempC = [kp ki kd]
    tempD = [b*kp 0]

    named_ss(ss(tempA, tempB, tempC, tempD), "PID", u=[:r, :y], y=:u)
end
filt(Ty, dy) = tf(1, [Ty^2, 2dy*Ty, 1])

##

"""
    equivalent_pid(Tsettle, ogain; simplified_r = true)

Construct a PID controller that is equivalent to the ADRC controller with settling time `Tsettle` and observer poles that are `ogain` times faster than the closed-loop poles. The returned controller has two inputs, `r` and `y`, and one output, `u`.

If `simplified_r` is true, the controller is a PI controller with set-point weighting on the proportional term and a first-order lowpass filter on the measurement. If `simplified_r` is false, the controller exactly matches the ADRC contorller, which is a filtered PID controller from the reference signal.
"""
function equivalent_pid(Tsettle, ogain, b0=1; simplified_r = true)
    # The expressions for the PID parameters are obtained symbolically in the bottom of the script
    kpy = (-72.0*ogain^3 - 108.0*ogain^2)/(3.0*Tsettle^2*ogain^2 + 6.0*Tsettle^2*ogain + Tsettle^2)
    kiy = -216.0*ogain^3/(3.0*Tsettle^3*ogain^2 + 6.0*Tsettle^3*ogain + Tsettle^3)    
    kdy = (-6.0*ogain^3 - 36.0*ogain^2 - 18.0*ogain)/(3.0*Tsettle*ogain^2 + 6.0*Tsettle*ogain + Tsettle)    
    Ty = -0.166666666666667*Tsettle*sqrt(1/(3.0*ogain^2 + 6.0*ogain + 1.0))    
    dy = -0.5*(3.0*ogain + 2.0)*sqrt(1/(3.0*ogain^2 + 6.0*ogain + 1.0))  
  
    b = -(6/Tsettle)^2 / kpy # Found by matching asymptotes
    C = if simplified_r
        -pid_2dof_2filt(kpy, kiy, kdy, Ty, dy, b, 0)
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
Ca = adrc(Tsettle, ogain)
Cr = Ca[:u,:r]
Cy = Ca[:u,:y]

s = tf('s')
Ki = 0.6
TZ1 = TZ2 = T
T1 = 0.2
C_suggested_pid = Ki * (1+TZ1*s)*(1 + TZ2*s) / (s*(1 + T1*s))
C_equivalent_pid = equivalent_pid(Tsettle, ogain; simplified_r=true)
label = ["ADRC" "Suggested PID" "Equivalent PID (simp.)"]

gangoffourplot(P, [Ca[:u,:y], C_suggested_pid, C_equivalent_pid[:u,:y]]; label)

##
feedback2d(P, Ca) = feedback(P, Ca[:u, :y], pos_feedback=true)*(Ca[:u, :r])


## Analysis
# The analysis below compares the ADRC controller to the in the paper suggested PID controller, as well as the approximately equivalent PI controller with set-point weighting

# ADRC controller has overall much higher gain
# It looks like a PI controller from r, and a filtered PID controller from y
# Overall, it has a much higher gain but rolloff from measurements
bodeplot([Ca, [C_suggested_pid -C_suggested_pid], C_equivalent_pid], w; label=repeat(label, inner=(1,4)), background_color_legend=nothing, foreground_color_legend=nothing)

# Step responses from r are almost identical
plot(step.([feedback2d(P, Ca), feedback(P*C_suggested_pid), feedback2d(P, C_equivalent_pid)], Ref(t)); label)

# So are closed-loop tf from r -> y 
bodeplot([feedback2d(P, Ca), feedback(P*C_suggested_pid), feedback2d(P, C_equivalent_pid)], w; label=repeat(label, inner=(1,2)), title="Gry", background_color_legend=nothing, foreground_color_legend=nothing)

# Bode plots from y -> u look different, ADRC is tuned much more aggressively but uses rolloff. The equivalent PID is identical to ADRC
bodeplot([G_CS(P, -Ca[:u, :y]), G_CS(P, C_suggested_pid), G_CS(P, -C_equivalent_pid[:u, :y])]; label=repeat(label, inner=(1,2)), title="Gyu", legend=:topleft, background_color_legend=nothing, foreground_color_legend=nothing)

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
], Ref(t)); label, ri=false, layout=(1,3), sp=(1:3)', size=(800,400), ylabel="y")

Pu = tf([K],[1, Tu])
plot(step.([
    feedback2d(Pu, Ca),
    feedback(Pu*C_suggested_pid),
    feedback2d(Pu, C_equivalent_pid)
], Ref(t)); label, ri=false, layout=(1,3), sp=(1:3)', size=(800,400), ylabel="y")


## To figure this out, we propagate symbolic variables through the adrc constructor
using Symbolics, SymbolicControlSystems
@variables Tsettle ogain b0
K = adrc(Tsettle, ogain)#, b0)
K = ss(identity.(K.A), identity.(K.B), identity.(K.C), identity.(K.D))
ex = Num(K)
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