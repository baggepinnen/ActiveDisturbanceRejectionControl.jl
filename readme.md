# Active Disturbance Rejection Control (ADRC)
ADRC is a recent trend (fad?) in the control-systems community. A recent paper, reference below, demonstrated that linear ADRC is a rebranding of a well established concept, state feedback with state provided by a disturbance observer, where the disturbance model is $1/s$. This kind approach is a common way of providing integral action in a state-feedback controller.

This repo demonstrates that linear first-order ADRC is in fact equivalent also to a filtered 2DOF PID controller, and very close to equivalent to a PI controller with set-point weighting and first-order low-pass filtering of the measurement.

The proof makes use of the formulation of ADRC provided in 

> "A Simulative Study on Active Disturbance Rejection Control (ADRC) as a Control Tool for Practitioners", Gernot Herbst

that is, with the "bandwidth parametrization" of the disturbance observer. The paper suggests tuning the PI controller to meet a certain reference step response, which is oftentimes a practice recommended against: 

> It may be strongly misleading to only show properties of a few systems for example the response of the output to command signals. A common omission in many papers and books. "Feedback Fundamentals", Karl Johan Åström

> “The user should not test the loop using set-point changes if the set point is to remain constant most of the time. To tune for fast recovery from load changes, a load disturbance should be simulated by stepping the controller output in manual, and then transferring to auto. For lag-dominant processes, the two responses are markedly different.”, - Shinskey 1993

This is also why the paper makes it look like the PID controller is much worse than the ADRC controller. If the PI(D) controller is instead tuned to perform well for disturbance rejection, one can later tune the reference response by adjusting the set-point weight (or performing more elaborate reference prefiltering). This little analysis demonstrates that the linear first-order ADRC controller is in fact _completely equivalent_ to a 2DOF PID controller, and with only very minor approximation error, equivalent also to a PI controller with set-point weighting and first-order low-pass filtering of the measurement.



For second-order linear ADRC, the conclusions are vey similar. The controller transfer function from the measurement to the control signal is (at least very close to) a second-order filtered PID controller, and the response from reference to control signal is approximately a PI controller. The expressions for the parameters are much more complex, but the conclusion that a 2DOF PID controller can do the job equally well remains.