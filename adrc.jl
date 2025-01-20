cd(@__DIR__)
using Pkg
Pkg.activate(".")

using Symbolics
# NOTE: extreme hacks to workaround https://github.com/JuliaSymbolics/Symbolics.jl/issues/1404
Base.zero(::Type{Any}) = Num(0)
Base.one(::Type{Any}) = Num(1)
Base.oneunit(::Type{Any}) = Num(1)
Base.inv(A::Matrix{Any}) = inv(identity.(A))
isinteractive() && (Base.active_repl.options.hint_tab_completes = false) # This messes with sympy https://discourse.julialang.org/t/sympy-makes-repl-to-stuck/124814/6
using ControlSystemsBase, Plots, RobustAndOptimalControl, Test, LinearAlgebra
default(margin=4Plots.mm, l=3, titlefontsize=12, legendfontsize=10, background_color_legend=nothing, foreground_color_legend=nothing)
t = 0:0.001:2

# The plant model used in experiments
K = 1
T = 1
P = tf([K],[T, 1])
s = tf('s')
w = exp10.(LinRange(-3, 3, 200))


"""
    adrc(Tsettle, ogain)

Construct an ADRC controller with settling time `Tsettle` and observer poles that are `ogain` times faster than the closed-loop poles. The returned controller has two inputs, `r` and `y`, and one output, `u`.
"""
function adrc(Tsettle, ogain)
    # See second-order analysis for how b0 is incorporated here
    Pdes1 = ss(tf([1],[1, 0]))
    Pdes2 = add_low_frequency_disturbance(Pdes1, 1)

    Kp = 4/Tsettle
    sCL = -Kp
    sESO = ogain*sCL
    k1 = -2sESO
    k2 = sESO^2
    K = [k1; k2;;]

    obs = observer_filter(Pdes2, K, output_state=true)

    obsn = named_ss(obs, u=[:u, :y], y=[:yh, :fh])
    Kpn = named_ss(ss(Kp), u=:e, y=:u0)
    s1 = sumblock("e = r-yh")
    s2 = sumblock("u = u0 - fh")

    systems = [
        obsn,
        Kpn,
        s1,
        s2,
    ]
    connections = [
        :e => :e
        :u0 => :u0
        :u => :u
        :fh => :fh
        :yh => :yh
    ]
    
    external_inputs = [:r, :y]
    external_outputs = [:u]
    
    connect(systems, connections; external_inputs, external_outputs, unique=true)
end

"""
    equivalent_pid(Tsettle, ogain; simplified_r = true)

Construct a PID controller that is equivalent to the ADRC controller with settling time `Tsettle` and observer poles that are `ogain` times faster than the closed-loop poles. The returned controller has two inputs, `r` and `y`, and one output, `u`.

If `simplified_r` is true, the controller is a PI controller with set-point weighting on the proportional term and a first-order lowpass filter on the measurement. If `simplified_r` is false, the controller exactly matches the ADRC contorller, which is a filtered PID controller from the reference signal.
"""
function equivalent_pid(Tsettle, ogain; simplified_r = true)
    den = Tsettle^2*(4 + 8ogain)
    ki = 64ogain^2 / den
    kpy = (32Tsettle*ogain + 16Tsettle*ogain^2) / den
    T = Tsettle^3 / den
    
    C = if simplified_r
        b = 4/Tsettle / kpy # Found by matching asymptotes
        pid_2dof(kpy, ki; b, form=:parallel) * [tf(1) 0; 0 tf(1, [T, 1])]
    else
        Cpidy = pid(kpy, ki, form=:parallel)*tf(1, [T, 1])
        kpr = 32Tsettle*ogain / den
        kdr = 4Tsettle^2 / den
        b = kpr / kpy
        Cpidr = pid(kpr, ki, kdr, form=:parallel)*tf(1, [T, 1])
        [Cpidr -Cpidy]
    end
    named_ss(C, y=:u, u=[:r, :y])
end

Tsettle = 1 # Parameters suggested in the paper
ogain = 10
Ca = adrc(Tsettle, ogain)
Cr = Ca[:u,:r]
Cy = Ca[:u,:y]

C_suggested_pid = pid(3.85, 3.85, form=:parallel)
C_equivalent_pid = equivalent_pid(Tsettle, ogain, simplified_r=true)
label = ["ADRC" "Suggested PI" "Equivalent PIF"]

feedback2d(P, Ca) = feedback(P, -Ca[:u, :y])*(Ca[:u, :r])
gangoffourplot(P, [Ca[:u,:y], C_suggested_pid, C_equivalent_pid[:u,:y]]; label, background_color_legend=nothing, foreground_color_legend=nothing, Ms_lines=[])
plot!(ylims=(-Inf, Inf), legend=[true false false false])



## Analysis
# The analysis below compares the ADRC controller to the in the paper suggested PID controller, as well as the approximately equivalent PI controller with set-point weighting

# ADRC controller has overall higher gain
# It looks like a PI controller from r, and a filtered PI controller from y
# Overall, it has a much higher low-frequency gain but rolloff from measurements
bodeplot([Ca, [C_suggested_pid -C_suggested_pid], C_equivalent_pid]; label=repeat(label, inner=(1,4)), background_color_legend=nothing, foreground_color_legend=nothing, title=["\$C_{ur}\$" "\$C_{uy}\$" "" ""], legend=[true false false false], linestyle=repeat([:solid :solid :dash], inner=(1,4)))
savefig("paper/figures/first_order_bode_C.pdf")

# Cr looks _almost_ like a non-filtered PI contorller (as opposed to a filtered PID contorller) We can identify the parameters of this simplified PI controller by matching the asymptotes of the two controllers. The low-frequency asymptote is given by ki, and the high-frequency asymptote is given by kpr. The high-frequency asymptote is given by the limit of Cr when s → ∞, which simlifies to the Kp used in the adrc controller.
# Cpidr2 = pid(4/Tsettle, kir, form=:parallel)


# Step responses from r are identical
plot(step.([feedback2d(P, Ca), feedback(P*C_suggested_pid), feedback2d(P, C_equivalent_pid)], Ref(t)); label)

# So are closed-loop tf from r -> y 
bodeplot([feedback2d(P, Ca), feedback(P*C_suggested_pid), feedback2d(P, C_equivalent_pid)]; label=repeat(label, inner=(1,2)), title="\$G_{yr}\$", linestyle=repeat([:solid :solid :dash], inner=(1,2)))
plot!(legend=[:bottomleft false])
savefig("paper/figures/first_order_bode_ry.pdf")

# Bode plots from y -> u look different, ADRC is tuned much more aggressively but uses rolloff
bodeplot([G_CS(P, -Ca[:u, :y]), G_CS(P, C_suggested_pid), G_CS(P, -C_equivalent_pid[:u, :y])]; label=repeat(label, inner=(1,2)), title="\$G_{yu}\$", legend=[:topleft false], background_color_legend=nothing, foreground_color_legend=nothing, linestyle=repeat([:solid :solid :dash], inner=(1,2)), l=3)
savefig("paper/figures/first_order_bode_uy.pdf")

## Reproduce response-plots from 
using MonteCarloMeasurements
K = 1
T = 1
Ku = Particles([0.1, 0.2, 0.5, 1, 2, 5, 10])
Tu = Particles([0.1, 0.2, 0.5, 1, 2, 5, 10])
Pu = tf([Ku],[T, 1])
plot(step.([
    feedback2d(Pu, Ca),
    feedback(Pu*C_suggested_pid),
    feedback2d(Pu, C_equivalent_pid)
], Ref(t)); label, ri=false, layout=(1,3), sp=(1:3)', size=(800,400), ylabel="y", c=[1 2 3], legend=:bottomright)
savefig("paper/figures/first_order_K.pdf")

Pu = tf([K],[Tu, 1])
plot(step.([
    feedback2d(Pu, Ca),
    feedback(Pu*C_suggested_pid),
    feedback2d(Pu, C_equivalent_pid)
], Ref(0:0.001:2)); label, ri=false, layout=(1,3), sp=(1:3)', size=(800,400), ylabel="y", c=[1 2 3], legend=:bottomright)
savefig("paper/figures/first_order_T.pdf")

# ## 




## To figure this out, we propagate symbolic variables through the adrc constructor
s = tf('s')
using Symbolics, SymbolicControlSystems
@variables Tsettle ogain vb0
K = adrc(Tsettle, ogain)#, vb0)
K = ss(identity.(K.A), identity.(K.B), identity.(K.C), identity.(K.D))
ex = [to_num(K[1,1]) to_num(K[1,2])]
tf.(ex)
# We then simply match the coefficients for each order of s, obtaining a system of equations that we can solve for the PID parameters
# The PID controllers on symbolic form are constructed below


Base.active_repl.options.hint_tab_completes = false # This messes with sympy https://discourse.julialang.org/t/sympy-makes-repl-to-stuck/124814/6
using SymbolicControlSystems
using SymbolicControlSystems: @syms
@syms kp ki kd Tf d Tsettle ogain b c
# Cy

tf_des = tf((kp + ki/s)/(s*Tf + 1))
tf_act = tf(symbolics_to_sympy(Symbolics.simplify(ex[2], expand=true)))

numvec1 = numvec(tf_des)[]
numvec2 = numvec(tf_act)[]
denvec1 = denvec(tf_des)[][1:end-1] 
denvec2 = denvec(tf_act)[][1:end-1] 


eqs = [
    numvec1./denvec1[end] - numvec2./denvec2[end],
    denvec1[1:end-1]./denvec1[end]  - denvec2[1:end-1]./denvec2[end] # Make the last element 1
]
vars = [kp, ki, Tf]

sol = sp.solve(eqs, vars)

@show kpy = sol[kp]
@show kiy = sol[ki]
@show Ty =  sol[Tf]

## Get nice expressions for the state-space realization
# Requires a temporary change to an if-statement in pid_ss_2dof
@variables T_s g
to_latex(x) = sp.latex(sp.simplify.(symbolics_to_sympy.(x)))
to_latex(equivalent_pid(T_s, g, simplified_r=true).A) |> println 
to_latex(equivalent_pid(T_s, g, simplified_r=true).B) |> println 
to_latex(equivalent_pid(T_s, g, simplified_r=true).C) |> println 
to_latex(equivalent_pid(T_s, g, simplified_r=true).D) |> println 

##
@syms k_p k_i T_f b
C_equiv_sym = pid_2dof(k_p, k_i; b, form=:parallel) * [tf(1) 0; 0 tf(1, [T_f, 1])]
function show_latex_ss(sys::AbstractStateSpace)
    A,B,C,D = to_latex.(ssdata(sys))
    println("\\begin{align}")
    println("\\dot{x} &= $(A)x + $(B)u \\\\")
    println("y &= $(C)x + $(D)u")
    println("\\end{align}")
end

show_latex_ss(C_equiv_sym)
## Gang of seven

function gangofsevenplot(P, C, F, args...; c, name="", kwargs...)
    S,D,CS,T = gangoffour(P,C)
    RY = T*F
    RU = CS*F
    RE = S*F
    bodeplot!(S, args...; show=false, title="\$S = 1/(1+PC)\$", lab="$name: \$S\$", c, sp=1, plotphase=false, legend=:bottomright, kwargs...)
    bodeplot!(D, args...; show=false, title="\$PS = P/(1+PC)\$", lab="$name: \$PS\$", c, sp=2, plotphase=false, legend=:bottom, kwargs...)
    bodeplot!(CS, args...; show=false, title="\$CS = C/(1+PC)\$", lab="$name: \$CS\$", c, sp=3, plotphase=false, legend=:bottom, kwargs...)
    bodeplot!(T, args...; show=false, title="\$T = PC/(1+PC)\$", lab="$name: \$T\$", c, sp=4, plotphase=false, legend=:bottomleft, kwargs...)
    Plots.hline!([1], l=(:black, :dash, 1), primary=false, sp=4)
    bodeplot!(RE, args...; show=false, title="\$S = 1/(1+PC)\$", lab="$name: \$SF = r\\to e\$", l=(:dash,), c, sp=1, plotphase=false, kwargs...)
    bodeplot!(RY, args...; show=false, title="\$T = PC/(1+PC)\$", lab="$name: \$TF = r\\to y\$", l=(:dash,), c, sp=4, plotphase=false, kwargs...)
    bodeplot!(RU, args...; show=false, title="\$CS = C/(1+PC)\$", lab="$name: \$CSF = r\\to u\$", l=(:dash,), c, sp=3, plotphase=false, kwargs...)
end
default(titlefontsize=14, legendfontsize=9)
w = exp10.(LinRange(-2, 4, 200))
F = tf(Cr) / tf(-Cy) # This computes the equivalent reference prefilter appearing before the error calculation
plot(; layout=4, ticks=:default, xscale=:log10, size=(1200,700))#, link=:both)
gangofsevenplot(P, -tf(Cy), F, w; name="ADRC", c=1, background_color_legend=nothing, foreground_color_legend=nothing)
gangofsevenplot(P, C_suggested_pid, tf(1), w; name="Suggested PID", label, c=2, background_color_legend=nothing, foreground_color_legend=nothing)

F_equivalent_pid = tf(C_equivalent_pid[:u,:r]) / tf(-C_equivalent_pid[:u,:y])
gangofsevenplot(P, -tf(C_equivalent_pid[:u,:y]), F_equivalent_pid, w; name="Equivalent PID", c=3, background_color_legend=nothing, foreground_color_legend=nothing, linestyle=:dot)
savefig("paper/figures/first_order_7.pdf")


## Complete simulation in MTK
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Der

@mtkmodel MTKCL begin
    @parameters begin
        K = 1.0
        T = 1.0
        Tsettle = 1.0
        ogain = 10.0
        b0 = 1
        Kp = 4/Tsettle
    end
    begin
        sCL = -Kp
        sESO = ogain*sCL
        k1 = -2sESO
        k2 = sESO^2
    end
    @variables begin
        x1(t) = 0
        x2(t) = 0
        Px(t) = 0
        u(t)
        y(t)
        r(t)
    end
    @equations begin
        r ~ (t >= 0.01)
        y ~ Px
        T*Der(Px) + Px ~ K*u
        Der(x1) ~ x2 + b0*u + k1*(y-x1)
        Der(x2) ~ k2*(y-x1)
        u ~ (Kp*(r-x1) - x2) / b0
    end
end

@named model = MTKCL()
model = complete(model)

ssys = structural_simplify(model)
prob = ODEProblem(ssys, [model.T => 0.1], (0.0, 2.0))

using OrdinaryDiffEqTsit5, Plots
sol = solve(prob, Tsit5())
plot(sol)

mats, ssys = ModelingToolkit.linearize_symbolic(model, [model.r], [model.y])
RobustAndOptimalControl.show_construction(ss(mats.A, mats.B, mats.C, mats.D))

T,K,ogain,Tsettle,Kp,b0 = 1,1,10,1,4,1
temp = let
    tempA = [-1 / T (-K*Kp) / (T*b0) (-K) / (T*b0); 2Kp*ogain -Kp - 2Kp*ogain 0; (Kp^2)*(ogain^2) -(Kp^2)*(ogain^2) 0]
    tempB = [(K*Kp) / (T*b0); Kp; 0;;]
    tempC = [1 0 0]
    tempD = [0;;]
    ss(tempA, tempB, tempC, tempD)
end

