include("Include.jl")
include("Simulations.jl")
using PyPlot

function runModel()
	# Script to solve the balance equations -
	time_start = 0.0
	time_stop = 120.0
	time_step_size = 0.01
	data_dictionary = DataDictionary(time_start,time_stop,time_step_size)
	XSS = estimate_steady_state(0.001,data_dictionary)

	  # Next, set the IC to the steady-state value -
	  initial_condition_array = XSS;
	  data_dictionary["initial_condition_array"] = initial_condition_array;
	  # Grab the control parameters - turn on gene_1 =
	  control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
	  control_parameter_dictionary["W_gene_1_RNAP"] = 1.0;
	  data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary;

	# Solve the model equations -
	(T,X) = SolveBalances(time_start,time_stop,time_step_size,data_dictionary)
	makePlot(T,X)
	return(T,X)
	
end

function runWashout()
	time_start = 0.0
	time_stop = 120.0
	time_step_size = 0.01
	data_dictionary = DataDictionary(time_start,time_stop,time_step_size)
	parameter_index = 1
	T,X=adj_washout_simulation(time_start, time_stop, time_step_size, parameter_index, data_dictionary)
	makePlot(T,X, "figures/WashoutSim_")
	return T,X

end

function runWashout(time_start, time_stop, time_step_size, data_dictionary)
	parameter_index = 1
	T,X=washout_simulation(time_start, time_stop, time_step_size, parameter_index, data_dictionary)
	makePlot(T,X, "figures/WashoutSim_")
	return T,X

end

function makePlot(T,X, savestr)
	close("all")
	fig1=figure(figsize=[15,15])
	fig3=figure(figsize=[15,15])
	fig2=figure(figsize=[15,15])
	for j in collect(1:size(X,2))
		@show j
		if(j<=3)
			PyPlot.figure(1)
			plot(T, X[:,j])
		elseif (j>3 && j<=6)
			PyPlot.figure(2)
			plot(T, X[:,j])
		elseif (j>6 && j<=9)
			PyPlot.figure(3)
			plot(T, X[:,j])
		end
	end
	legarr = ["gene 1", "gene 2", "gene 3", "mRNA 1", "mRNA2", "mRNA3", "protein 1", "protein 2", "protein 3" ]
	PyPlot.figure(1)
	legend(legarr[1:3], loc="best")
	savefig(string(savestr, "genes.pdf"))
	PyPlot.figure(2)
	legend(legarr[4:6], loc="best")
	savefig(string(savestr, "mRNA.pdf"))
	PyPlot.figure(3)
	legend(legarr[7:9], loc="best")
	savefig(string(savestr, "proteins.pdf"))
end
