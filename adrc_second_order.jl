
cd(@__DIR__)
using Pkg
Pkg.activate(".")

# NOTE: extreme hack to workaround https://github.com/JuliaSymbolics/Symbolics.jl/issues/1404
Base.zero(::Type{Any}) = Num(0)
Base.one(::Type{Any}) = Num(1)
Base.convert(::Type{T}, x::AbstractArray{Num}) where T <: Array{Num} = T(map(Num, x)) # https://github.com/JuliaSymbolics/Symbolics.jl/issues/1405

using ControlSystemsBase, Plots, MonteCarloMeasurements, RobustAndOptimalControl, Test, LinearAlgebra
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
    s1 = sumblock("e = r-yh")
    s2 = sumblock("u0 = up - yh2")
    sd = sumblock("u = u0 - fh")

    systems = [
        obsn,
        Kpn,
        s1,
        s2,
        sd,
    ]
    connections = [
        :e => :e
        :up => :up
        :u0 => :u0
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
    

"""
    equivalent_pid(Tsettle, ogain; simplified_r = true)

Construct a PID controller that is equivalent to the ADRC controller with settling time `Tsettle` and observer poles that are `ogain` times faster than the closed-loop poles. The returned controller has two inputs, `r` and `y`, and one output, `u`.

If `simplified_r` is true, the controller is a PI controller with set-point weighting on the proportional term and a first-order lowpass filter on the measurement. If `simplified_r` is false, the controller exactly matches the ADRC contorller, which is a filtered PID controller from the reference signal.
"""
function equivalent_pid(Tsettle, ogain; simplified_r = true)
    kpy, kiy, kdy, Ty, dy = 
    #((-72*ogain^3 - 108*ogain^2)/(3*Tsettle^2*ogain^2 + 6*Tsettle^2*ogain + Tsettle^2), -216*ogain^3/(3*Tsettle^3*ogain^2 + 6*Tsettle^3*ogain + Tsettle^3), (-6*ogain^3 - 36*ogain^2 - 18*ogain)/(3*Tsettle*ogain^2 + 6*Tsettle*ogain + Tsettle), Tsettle*sqrt(1/(3*ogain^2 + 6*ogain + 1))/6, (3*ogain + 2)*sqrt(1/(3*ogain^2 + 6*ogain + 1))/2)
    
    ((-72*ogain^3 - 108*ogain^2)/(3*Tsettle^2*ogain^2 + 6*Tsettle^2*ogain + Tsettle^2), -216*ogain^3/(3*Tsettle^3*ogain^2 + 6*Tsettle^3*ogain + Tsettle^3), (-6*ogain^3 - 36*ogain^2 - 18*ogain)/(3*Tsettle*ogain^2 + 6*Tsettle*ogain + Tsettle), -Tsettle*sqrt(1/(3*ogain^2 + 6*ogain + 1))/6, -(3*ogain + 2)*sqrt(1/(3*ogain^2 + 6*ogain + 1))/2)

    # Cpidy = (kpy + kiy/s + kdy*s)/(s^2*Ty^2 + 2*dy*s*Ty + 1)
    
    C = if simplified_r
        b = -(6/Tsettle)^2 / kpy # Found by matching asymptotes
        # Cpidr = pid(4/Tsettle, kir, form=:parallel) # Equivalent if concatenated with Cpidy below
        # Tf, _ = wd2Tfd(1/Ty, dy)
        # We create the PID below with a tiny filter constant Tf=1e-9 since the constructor can only construct state space systems and those need a filter on the derivative term
        -pid_2dof(kpy, kiy, kdy; b, c = 0, form=:parallel, state_space=false, Tf=1e-6) * [tf(1) 0; 0 tf(1, [Ty^2, 2dy*Ty, 1])]
    else
        Cpidy = pid(kpy, ki, form=:parallel)*tf(1, [T, 1])
        kpr = 32*Tsettle*ogain / den
        kdr = 4Tsettle^2 / den
        b = kpr / kpy
        Cpidr = pid(kpr, ki, kdr, form=:parallel)*tf(1, [T, 1])
        [Cpidr -Cpidy]
    end
    named_ss(C, y=:u, u=[:r, :y])
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
C_equivalent_pid = equivalent_pid(Tsettle, ogain)
label = ["ADRC" "Suggested PID" "Equivalent PID (simp.)"]

gangoffourplot(P, [Ca[:u,:y], C_suggested_pid, C_equivalent_pid[:u,:y]]; label)

feedback2d(P, Ca) = feedback(P, -Ca[:u, :y])*(Ca[:u, :r])


## Analysis
# The analysis below compares the ADRC controller to the in the paper suggested PID controller, as well as the approximately equivalent PI controller with set-point weighting

# ADRC controller has overall much higher gain
# It looks like a PI controller from r, and a filtered PID controller from y
# Overall, it has a much higher low-frequency gain but rolloff from measurements
bodeplot([Ca, [C_suggested_pid -C_suggested_pid], C_equivalent_pid]; label=repeat(label, inner=(1,2)))


# Step responses from r are identical
plot(step.([feedback2d(P, Ca), feedback(P*C_suggested_pid), feedback2d(P, C_equivalent_pid)], Ref(t)); label)

# So are closed-loop tf from r -> y 
bodeplot([feedback2d(P, Ca), feedback(P*C_suggested_pid), feedback2d(P, C_equivalent_pid)]; label=repeat(label, inner=(1,2)), title="Gry")

# Bode plots from y -> u look different, ADRC is tuned much more aggressively but uses rolloff
bodeplot([G_CS(P, -Ca[:u, :y]), G_CS(P, C_suggested_pid), G_CS(P, -C_equivalent_pid[:u, :y])]; label=repeat(label, inner=(1,2)), title="Gyu", legend=:topleft, background_color_legend=nothing, foreground_color_legend=nothing)

## Reproduce response-plots from 
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
@variables Tsettle ogain
K = adrc(Tsettle, ogain)
K = ss(identity.(K.A), identity.(K.B), identity.(K.C), identity.(K.D))
ex = Num(K)
tf.(ex)
# We then simply match the coefficients for each order of s, obtaining a system of equations that we can solve for the PID parameters
# The PID controllers on symbolic form are constructed below

using SymbolicControlSystems: @syms
@syms kp ki kd T d Tsettle ogain

# Cy

numvec1 = [
    -216*Tsettle^2*ogain^3 - 1296*Tsettle^2*ogain^2 - 648*Tsettle^2*ogain
    -2592*Tsettle*ogain^3 - 3888*Tsettle*ogain^2
                                     -7776*ogain^3
]

numvec2 = [
    kd
    kp
    ki
]

denvec1 = [
    1*Tsettle^5 # obtained from tf(ex[2])
    18*Tsettle^4*ogain + 12*Tsettle^4
    108*Tsettle^3*ogain^2 + 216*Tsettle^3*ogain + 36*Tsettle^3
]
denvec2 = [
    T^2
    2T*d
    1
]

eqs = [
    numvec1./denvec1[end] .- numvec2./denvec2[end]
    denvec1[1:end-1]./denvec1[end] .- denvec2[1:end-1]./denvec2[end]
]
vars = [kp, ki, kd, T, d]

sols = sp.solve(eqs, vars)


# obtained from sol
kpy, kiy, kdy, Ty, dy = ((-72*ogain^3 - 108*ogain^2)/(3*Tsettle^2*ogain^2 + 6*Tsettle^2*ogain + Tsettle^2), -216*ogain^3/(3*Tsettle^3*ogain^2 + 6*Tsettle^3*ogain + Tsettle^3), (-6*ogain^3 - 36*ogain^2 - 18*ogain)/(3*Tsettle*ogain^2 + 6*Tsettle*ogain + Tsettle), Tsettle*sqrt(1/(3*ogain^2 + 6*ogain + 1))/6, (3*ogain + 2)*sqrt(1/(3*ogain^2 + 6*ogain + 1))/2)
# ((-72*ogain^3 - 108*ogain^2)/(3*Tsettle^2*ogain^2 + 6*Tsettle^2*ogain + Tsettle^2), -216*ogain^3/(3*Tsettle^3*ogain^2 + 6*Tsettle^3*ogain + Tsettle^3), (-6*ogain^3 - 36*ogain^2 - 18*ogain)/(3*Tsettle*ogain^2 + 6*Tsettle*ogain + Tsettle), -Tsettle*sqrt(1/(3*ogain^2 + 6*ogain + 1))/6, -(3*ogain + 2)*sqrt(1/(3*ogain^2 + 6*ogain + 1))/2)

Cpidy = (kpy + kiy/s + kdy*s)/(s^2*Ty^2 + 2*dy*s*Ty + 1)

