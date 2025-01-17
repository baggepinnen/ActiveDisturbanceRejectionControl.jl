
cd(@__DIR__)
using Pkg
Pkg.activate(".")

using Symbolics
# NOTE: extreme hacks to workaround https://github.com/JuliaSymbolics/Symbolics.jl/issues/1404
Base.zero(::Type{Any}) = Num(0)
Base.one(::Type{Any}) = Num(1)
Base.oneunit(::Type{Any}) = Num(1)
Base.inv(A::Matrix{Any}) = inv(identity.(A))

using ControlSystemsBase, Plots, RobustAndOptimalControl, Test, LinearAlgebra
t = 0:0.001:2

# The plant model used in experiments
K = 1
T = 1
P = tf([K],[1, T])
s = tf('s')


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
label = ["ADRC" "Suggested PID" "Equivalent PID (simp.)"]

feedback2d(P, Ca) = feedback(P, -Ca[:u, :y])*(Ca[:u, :r])
gangoffourplot(P, [Ca[:u,:y], C_suggested_pid, C_equivalent_pid[:u,:y]]; label)



## Analysis
# The analysis below compares the ADRC controller to the in the paper suggested PID controller, as well as the approximately equivalent PI controller with set-point weighting

# ADRC controller has overall higher gain
# It looks like a PI controller from r, and a filtered PI controller from y
# Overall, it has a much higher low-frequency gain but rolloff from measurements
bodeplot([Ca, [C_suggested_pid -C_suggested_pid], C_equivalent_pid]; label=repeat(label, inner=(1,2)))

# Cr looks _almost_ like a non-filtered PI contorller (as opposed to a filtered PID contorller) We can identify the parameters of this simplified PI controller by matching the asymptotes of the two controllers. The low-frequency asymptote is given by ki, and the high-frequency asymptote is given by kpr. The high-frequency asymptote is given by the limit of Cr when s → ∞, which simlifies to the Kp used in the adrc controller.
# Cpidr2 = pid(4/Tsettle, kir, form=:parallel)


# Step responses from r are identical
plot(step.([feedback2d(P, Ca), feedback(P*C_suggested_pid), feedback2d(P, C_equivalent_pid)], Ref(t)); label)

# So are closed-loop tf from r -> y 
bodeplot([feedback2d(P, Ca), feedback(P*C_suggested_pid), feedback2d(P, C_equivalent_pid)]; label=repeat(label, inner=(1,2)), title="Gry")

# Bode plots from y -> u look different, ADRC is tuned much more aggressively but uses rolloff
bodeplot([G_CS(P, -Ca[:u, :y]), G_CS(P, C_suggested_pid), G_CS(P, -C_equivalent_pid[:u, :y])]; label=repeat(label, inner=(1,2)), title="Gyu", legend=:topleft, background_color_legend=nothing, foreground_color_legend=nothing)

## Reproduce response-plots from 
using MonteCarloMeasurements
K = 1
T = 1
Ku = Particles([0.1, 0.2, 0.5, 1, 2, 5, 10])
Tu = Particles([0.1, 0.2, 0.5, 1, 2, 5, 10])
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


# ## 




## To figure this out, we propagate symbolic variables through the adrc constructor
s = tf('s')
using Symbolics, SymbolicControlSystems
@variables vTsettle vogain vb0
K = adrc(vTsettle, vogain)#, vb0)
K = ss(identity.(K.A), identity.(K.B), identity.(K.C), identity.(K.D))
ex = [Num(K[1,1]) Num(K[1,2])]
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