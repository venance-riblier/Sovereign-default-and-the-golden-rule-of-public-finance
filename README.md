# Sovereign-default-and-the-golden-rule-of-public-finance
 Replication files for my master thesis: Sovereign default and the golden rule of public finance.

The baseline model is solved by running `Solve_baseline.m` and the corresponding policy functions are stored in `policies.mat`. 
The results of the simulations are not stored here but they can be replicated by running `Simulations_baseline.m`. The analysis of the baseline
model then can be replicated by running `Analysis_baseline.m`.

The models with enforceable debt contracts, the golden rule, and the debt ceiling rule experiments are solved by running `Solve_experiments.m` in the folder Experiments.
The corresponding policies are stored in the subfolder Experiments/Policies. The simulations can be replicated by running `Simulations_experiments.m`,
the results will store in the folder Experiments/Simulations. The analysis of the experiments then can be replicated by running `Analysis_experiments.m`.

Please note that the code requires two routines that are provided in the corresponding folder. Please consider including them in your own MATLAB routine path, or adding them on the top of the MATLAB search path using ``addpath()``. The simulations also require the ``hpfilter()`` function from the [Econometrics Toolbox](https://www.mathworks.com/help/econ/index.html?s_tid=CRUX_lftnav).
