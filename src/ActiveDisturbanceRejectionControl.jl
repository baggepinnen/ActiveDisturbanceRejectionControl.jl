module ActiveDisturbanceRejectionControl
using LinearAlgebra
using ControlSystemsBase, RobustAndOptimalControl

export adrc, equivalent_pid

"""
    adrc(Tsettle, ogain)

Construct an ADRC controller with settling time `Tsettle` and observer poles that are `ogain` times faster than the closed-loop poles. The returned controller has two inputs, `r` and `y`, and one output, `u`. The controller is to be used with positive feedback, the closed-loop transfer function from `r` to `y` may be obtained with
```
feedback2d(P, C) = feedback(P, C[:u, :y], pos_feedback=true)*C[:u, :r]
```
"""
function adrc(Tsettle, ogain, b0 = 1; order = 1)
    T = promote_type(typeof(Tsettle), typeof(ogain))
    Pdes1 = ss(tf([one(T)],[1; zeros(order)]))*b0
    Pdes2 = add_low_frequency_disturbance(Pdes1, order)
    b0n = named_ss(ss(1/b0), u=:ub0, y=:u)
    s1 = sumblock("e = r-yh")

    if order == 1
        Kp = 4/Tsettle
        sCL = -Kp
        sESO = ogain*sCL
        k1 = -2sESO
        k2 = sESO^2
        K = [k1; k2;;]
        obs = observer_filter(Pdes2, K, output_state=true)
        obsn = named_ss(obs, u=[:u, :y], y=[:yh, :fh])
        Kpn = named_ss(ss(Kp), u=:e, y=:u0)
        s2 = sumblock("ub0 = u0 - fh")
    
        systems = [
            obsn,
            Kpn,
            b0n,
            s1,
            s2,
        ]
        connections = [
            :e => :e
            :u0 => :u0
            :ub0 => :ub0
            :u => :u
            :fh => :fh
            :yh => :yh
        ]
    elseif order == 2
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
    else
        error("Only order 1 and 2 are supported")
    end

    external_inputs = [:r, :y]
    external_outputs = [:u]
    
    connect(systems, connections; external_inputs, external_outputs, unique=true)
end


"""
    equivalent_pid(Tsettle, ogain; simplified_r = true)

Construct a PID controller that is equivalent to the ADRC controller with settling time `Tsettle` and observer poles that are `ogain` times faster than the closed-loop poles. The returned controller has two inputs, `r` and `y`, and one output, `u`.

If `simplified_r` is true, the controller is a PI controller with set-point weighting on the proportional term and a first-order lowpass filter on the measurement. If `simplified_r` is false, the controller exactly matches the ADRC contorller, which is a filtered PID controller from the reference signal.
"""
function equivalent_pid(Tsettle, ogain, b0; simplified_r = true, order=1)
    if order == 1
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
        return named_ss(1/b0 * C, y=:u, u=[:r, :y])

    elseif order == 2
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
        end
        return named_ss(1/b0 * ss(C) * diagm([1, 1]), y=:u, u=[:r, :y])
    else
        error("Only order 1 and 2 are supported")
    end
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

end