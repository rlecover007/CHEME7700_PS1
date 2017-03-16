## Gene regulatory network (GRN) model ##
This repository contains the model code for an effective GRN model implemented in the [Julia](http://julialang.org) programming language.

### Installation and Requirements
You can download this repository as a zip file, or `clone`/`pull` it by using the command (from the command-line):

	$ git pull https://github.com/varnerlab/(your repository name here).git

or

	$ git clone https://github.com/varnerlab/(your repository name here).git

The model code was machine generated using the [Gene Regulatory Network in Julia (JuGRN)](https://github.com/varnerlab/JuGRN-Generator) code generation system.
The model code uses several [Julia](http://julialang.org) packages:

Package | Description | Command
--- | --- | ---
ODE | Contains the ``ode23s`` subroutine to solve the model equations | Pkg.add("ODE")
PyPlot | Used to make figures (assume you have Python installed) | Pkg.add("PyPlot")

### How do I solve the model equations?
The model equations can be solved using the ``Driver.jl`` script, or by calling the ``SolveBalances`` function directly, encoded in ``SolveBalances.jl``:

	(T,X) = SolveBalances(time_start,time_stop,time_step_size,data_dictionary)

where:

Argument | Type | Description
--- | --- | ---
``time_start`` | scalar (float/double) | Initial simulation time (default time unit is min)
``time_stop`` | scalar (float/double) | Final simulation time (default time unit is min)
``time_step`` | scalar (float/double) | Simulation time step (depending upon the solver, the internal time step used during the solution of the GRN equations may be different)
``data_dictionary`` | [Julia dictionary](http://docs.julialang.org/en/stable/stdlib/collections/?highlight=dict#Base.Dict) | Instance of the data dictionary structure created by the ``DataDictionary`` function
``T`` | [array](http://docs.julialang.org/en/stable/stdlib/arrays/?highlight=array) (``number_of_steps x 1``) | Output time array (``time_state:time_step:time_stop``)
``X`` |[array](http://docs.julialang.org/en/stable/stdlib/arrays/?highlight=array) (``number_of_steps x number_of_states``) | Output state array (solution of the model equations)

### How do I solve the adjoint balances?
The adjoint balances can be solved (for a particular parameter index) using the ``AdjDriver.jl`` script, or by calling the ``SolveAdjBalances`` function
encoded in ``SolveBalances.jl``:

	(T,X) = SolveAdjBalances(time_start,time_stop,time_step_size,parameter_index,data_dictionary)

where:

Argument | Type | Description
--- | --- | ---
``time_start`` | scalar (float/double) | Initial simulation time (default time unit is min)
``time_stop`` | scalar (float/double) | Final simulation time (default time unit is min)
``time_step`` | scalar (float/double) | Simulation time step (depending upon the solver, the internal time step used during the solution of the GRN equations may be different)
``parameter_index`` | scalar (int) | index of the parameter that we are calculating the sensitivity for
``data_dictionary`` | [Julia dictionary](http://docs.julialang.org/en/stable/stdlib/collections/?highlight=dict#Base.Dict) | Instance of the data dictionary structure created by the ``DataDictionary`` function
``T`` | [array](http://docs.julialang.org/en/stable/stdlib/arrays/?highlight=array) (``number_of_steps x 1``) | Output time array (``time_state:time_step:time_stop``)
``X`` |[array](http://docs.julialang.org/en/stable/stdlib/arrays/?highlight=array) (``number_of_steps x (2*number_of_states)``) | Output state array (solution of the adjoint model equations). The first number of states is the model solution, while the second state block is the sensitivity of the state with respect to the parameter index
