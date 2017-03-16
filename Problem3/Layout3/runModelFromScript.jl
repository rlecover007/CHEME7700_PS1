function runWashout(time_start, time_stop, time_step_size, data_dictionary)
	# Script to solve the balance equations -
		# First - run the model to steady-state w/no ATRA forcing -
#	println("Estimating SS")
	  XSS = estimate_steady_state(0.001,data_dictionary)

	  # Next, set the IC to the steady-state value -
	  initial_condition_array = XSS;
	  data_dictionary["initial_condition_array"] = initial_condition_array;

	  # Phase 1: Run the model 1/4 of the final time w/o ATRA
	  # Run the model for a section of time w/no ATRA forcing -
	  time_start_phase_1 = 0.0;
	  time_stop_phase_1 = 4.0;

	  # Solve the model equations -
#		println("Solving before induce")
	  (TP1,XP1) = SolveBalances(time_start_phase_1,time_stop_phase_1,time_step_size,data_dictionary);

	  # Phase 2: Express gene_1 and run to end-time (ic = last point phase 1) -
	  # Update the IC again -
	  initial_condition_array = XP1[end,:];
	  data_dictionary["initial_condition_array"] = initial_condition_array;

	  # Grab the control parameters - turn on gene_1 =
	  control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
	  control_parameter_dictionary["W_gene_1_RNAP"] = 1.0;
	  data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary;

	  # Run the model for a section of time w/no ATRA forcing -
	  time_start_phase_2 = time_stop_phase_1+time_step_size
	  time_stop_phase_2 = time_start_phase_2 + 6

	  # Solve the model equations -
#		println("inducing")
	  (TP2,XP2) = SolveBalances(time_start_phase_2,time_stop_phase_2,time_step_size,data_dictionary);

	  # Washout -
	  initial_condition_array = XP2[end,:];
	  data_dictionary["initial_condition_array"] = initial_condition_array;

	  # Grab the control parameters - turn on gene_1 =
	  control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
	  control_parameter_dictionary["W_gene_1_RNAP"] = 0.0;
	  data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary;

	  time_start_phase_3 = time_stop_phase_2+time_step_size
	  time_stop_phase_3 = time_stop
#		println("After inducer removed")
	  (TP3,XP3) = SolveBalances(time_start_phase_3,time_stop_phase_3,time_step_size,data_dictionary);

	  # Package the two phases together -
	  T = [TP1 ; TP2 ; TP3];
	  X = [XP1 ; XP2 ; XP3];
	#makePlot(T,X, "figures/Prob3Layout1")
	return(T,X)
	
end

function runFromCmd()
	include("Include.jl")
	time_start = parse(ARGS[1])
	time_stop = parse(ARGS[2])
	time_step_size = parse(ARGS[3])
	strparams=(ARGS[4:end])
	params = Float64[]
	gene_params = Float64[]
	for j in collect(1:size(strparams,1)-3)
		push!(params, parse(strparams[j]))
	end
	for j in collect(size(strparams,1)-2:size(strparams,1))
		push!(gene_params, parse(strparams[j]))
	end
	data_dictionary =DataDictionary(time_start, time_stop, time_step_size, params, gene_params)
	T,X=runWashout(time_start, time_stop, time_step_size, data_dictionary)
	writedlm("Layout3Output.txt", X)
end

runFromCmd()


