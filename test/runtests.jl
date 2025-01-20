using Test, ActiveDisturbanceRejectionControl, ControlSystemsBase, RobustAndOptimalControl

if isinteractive()
    using Plots
end

## order = 1
order = 1
t = 0:0.001:2

# The plant model used in experiments
K = 1
T = 0.1
P = tf([K],[1, T])
s = tf('s')
w = exp10.(LinRange(-3, 3, 200))

Tsettle = 1.0 # Parameters suggested in the paper
ogain = 5.0
b0 = 1.0
Ca = adrc(Tsettle, ogain, b0; order)
C_equivalent_pid = equivalent_pid(Tsettle, ogain, b0; simplified_r=true, order)
label = ["ADRC" "Equivalent PIF"]

feedback2d(P, Ca) = feedback(P, Ca[:u, :y], pos_feedback=true)*(Ca[:u, :r])


# gangoffourplot(P, [Ca[:u,:y], C_equivalent_pid[:u,:y]]; label, background_color_legend=nothing, foreground_color_legend=nothing, Ms_lines=[])
# plot!(ylims=(-Inf, Inf), legend=[true false false false])
# bodeplot([Ca, C_equivalent_pid]; label=repeat(label, inner=(1,4)), background_color_legend=nothing, foreground_color_legend=nothing, title=["\$C_{ur}\$" "\$C_{uy}\$" "" ""], legend=[false true false false], linestyle=repeat([:solid :dash], inner=(1,4)))

@test hinfnorm(minreal(Ca[:u, :y] - C_equivalent_pid[:u, :y], 1e-6))[1] < 1e-8
@test hinfnorm(minreal(Ca[:u, :r] - C_equivalent_pid[:u, :r], 1e-6))[1] < 1.1

Gcl = feedback(P, Ca[:u, :y], pos_feedback=true)
@test isstable(Gcl)


# Test with different parameters
Tsettle = 2.0 
ogain = 8.0
b0 = 10.0
Ca = adrc(Tsettle, ogain, b0; order)
C_equivalent_pid = equivalent_pid(Tsettle, ogain, b0; simplified_r=true, order)
@test hinfnorm(minreal(Ca[:u, :y] - C_equivalent_pid[:u, :y], 1e-6))[1] < 1e-8
@test hinfnorm(minreal(Ca[:u, :r] - C_equivalent_pid[:u, :r], 1e-6))[1] < 0.1


# Without approximation in reference response
C_equivalent_pid = equivalent_pid(Tsettle, ogain, b0; simplified_r=false, order)
@test hinfnorm(minreal(Ca[:u, :y] - C_equivalent_pid[:u, :y], 1e-6))[1] < 1e-8
@test hinfnorm(minreal(Ca[:u, :r] - C_equivalent_pid[:u, :r], 1e-6))[1] == 0


## Order 2

order = 2
# The plant model used in experiments
K = 1
D = 1
T = 1
P = tf([K],[T^2, 2D*T, 1])

Tsettle = 5.0 # Parameters suggested in the paper
ogain = 10.0
b0 = 1.0
Ca = adrc(Tsettle, ogain, b0; order)
C_equivalent_pid = equivalent_pid(Tsettle, ogain, b0; simplified_r=true, order)
@test hinfnorm(minreal(Ca[:u, :y] - C_equivalent_pid[:u, :y], 1e-6))[1] < 1e-8
@test hinfnorm(minreal(Ca[:u, :r] - C_equivalent_pid[:u, :r], 1e-6))[1] < 0.6


# Test with different parameters
Tsettle = 2.0 
ogain = 8.0
b0 = 10.0
Ca = adrc(Tsettle, ogain, b0; order)
C_equivalent_pid = equivalent_pid(Tsettle, ogain, b0; simplified_r=true, order)
@test hinfnorm(minreal(Ca[:u, :y] - C_equivalent_pid[:u, :y], 1e-6))[1] < 1e-8
@test hinfnorm(minreal(Ca[:u, :r] - C_equivalent_pid[:u, :r], 1e-6))[1] < 0.4
