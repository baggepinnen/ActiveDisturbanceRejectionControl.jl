# Active Disturbance Rejection Control (ADRC)
ADRC is a recent trend (fad?) in the control-systems community. A recent paper, reference below, demonstrated that linear ADRC is a rebranding of a well established concept, state feedback with state provided by a disturbance observer, where the disturbance model is $1/s$. This kind approach is a common way of providing integral action in a state-feedback controller.

This repo demonstrates that linear first-order ADRC is in fact equivalent also to a filtered 2DOF PID controller, and very close to equivalent to a PI controller with set-point weighting and first-order low-pass filtering of the measurement.

The proof makes use of the formulation of ADRC provided in 

> "A Simulative Study on Active Disturbance Rejection Control (ADRC) as a Control Tool for Practitioners", Gernot Herbst

that is, with the "bandwidth parametrization" of the disturbance observer. The paper suggests tuning the PI controller to meet a certain reference step response, which is oftentimes a practice recommended against. This is also why the paper makes it look like the PID controller is much worse than the ADRC controller. If the PI(D) controller is instead tuned to perform well for disturbance rejection, one can later tune the reference response by adjusting the set-point weight (or performing more elaborate reference prefiltering). This little analysis demonstrates that the linear first-order ADRC controller is in fact _completely equivalent_ to a 2DOF PID controller, and with only very minor approximation error, equivalent also to a PI controller with set-point weighting and first-order low-pass filtering of the measurement.