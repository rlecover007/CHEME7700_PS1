include("Include.jl")
using PyPlot

function calculateAvgSensitivityArr(time_start, time_stop)
	pathToSensitivityFiles = "sensitivity/"
	pattern = "AdjSimulation-P"
	timeSkip = 1.0
#	time_start = 0.0
#	time_stop = 120.0
	time_step_size=0.01
	data_dictionary = DataDictionary(time_start,time_stop,time_step_size)
	parameter_name_mapping_array = data_dictionary["parameter_name_mapping_array"]
#	sensitivity_arr = calculate_sensitivity_array(pathToSensitivityFiles, pattern, timeSkip, data_dictionary)
	avg_sensitivity_arr = calculate_average_scaled_sensitivity_array(pathToSensitivityFiles, pattern, data_dictionary, time_start, time_stop)
	estable_params = estimate_identifiable_parameters(avg_sensitivity_arr, .01)
	measurable_species = hcat(zeros(1,3),ones(1,6)) #we can measures proteins and mRNAs, but not gene levels
#	FIM = calculate_fisher_information_matrix(avg_sensitivity_arr, measurable_species,estable_params)
#	@show FIM
	U,S,V=svd(avg_sensitivity_arr, thin=false)
	#A = U*diagm(S)*V'
	timestr = string("From", time_start, "to", time_stop)
	makeHeatMap(U, timestr)
	makeHeatMap(V,timestr)
	makeHeatMap(avg_sensitivity_arr, timestr)
	close("all")
	@show parameter_name_mapping_array[estable_params]
	return U,S,V, avg_sensitivity_arr
end

function calculateSensitivityArr(timeSkip)
	pathToSensitivityFiles = "sensitivity/"
	pattern = "AdjSimulation-P"
	time_start = 0.0
	time_stop = 120.0
	time_step_size=0.01
	data_dictionary = DataDictionary(time_start,time_stop,time_step_size)
	parameter_name_mapping_array = data_dictionary["parameter_name_mapping_array"]
	time_sample_arr,sensitivity_arr = calculate_sensitivity_array(pathToSensitivityFiles, pattern, timeSkip, data_dictionary)
	estable_params = estimate_identifiable_parameters(sensitivity_arr, .01)
	return size(estable_params,1),parameter_name_mapping_array[estable_params]
end

function calculateSensitivityArr(timeSkip, selectedSpecies)
	#selectedSpecies should be an array containing the indexes of the species we can measure
	pathToSensitivityFiles = "sensitivity/"
	pattern = "AdjSimulation-P"
	time_start = 0.0
	time_stop = 120.0
	time_step_size=0.01
	data_dictionary = DataDictionary(time_start,time_stop,time_step_size)
	parameter_name_mapping_array = data_dictionary["parameter_name_mapping_array"]
	time_sample_arr,sensitivity_arr = calculate_sensitivity_array_selected_species(pathToSensitivityFiles, pattern, timeSkip, data_dictionary, selectedSpecies)
	estable_params = estimate_identifiable_parameters(sensitivity_arr, .001)
	return size(estable_params,1),parameter_name_mapping_array[estable_params]
end

function createW(num_species)
	W = zeros(num_species, num_species)
end

function makeHeatMap(data, timestr)
	figure(figsize=[15,15])
	timeSkip = 1.0
	time_start = 0.0
	time_stop = 120.0
	time_step_size=0.01
	data_dictionary = DataDictionary(time_start,time_stop,time_step_size)
	labels=["gene 1", "gene 2", "gene 3", "mRNA 1", "mRNA2", "mRNA3", "protein 1", "protein 2", "protein 3" ]
	if(size(data,1)==9 && size(data,2)==9)
		x = collect(0:size(labels,1)-1)
		y= collect(0:size(labels,1)-1)
		pcolormesh(x,y,abs(data),cmap="Reds")
		colorbar()
		ax = gca()
		ax[:xaxis][:set_ticks](x-.5)
		ax[:xaxis][:set_ticklabels](labels, rotation = 60, fontsize = 8)
		ax[:yaxis][:set_ticks](y+.5)
		ax[:yaxis][:set_ticklabels](labels, rotation = 0, fontsize = 8)
		savefig(string("figures/ProteinInteractions",timestr, ".pdf"))
	elseif(size(data,1)==24 && size(data,2)==24)
		x = collect(0:size(data,1))
		y= collect(0:size(data,1))
		pcolormesh(x,y,abs(data), cmap="Reds")
		colorbar()
		ax = gca()
		labels = data_dictionary["parameter_name_mapping_array"]
		ax[:xaxis][:set_ticks](x-.5)
		ax[:xaxis][:set_ticklabels](labels, rotation = 60, fontsize = 8)
		ax[:yaxis][:set_ticks](y+.5)
		ax[:yaxis][:set_ticklabels](labels, rotation = 0, fontsize = 8)
		savefig(string("figures/ParameterInteractions",timestr, ".pdf"))
	elseif(size(data,1)==9 && size(data,2)==24)
		x = collect(0:size(data,2))
		y= collect(0:size(data,1))
		pcolormesh(x,y,abs(data), cmap="Reds")
		colorbar()
		ax = gca()
		xlabels = data_dictionary["parameter_name_mapping_array"]
		ax[:xaxis][:set_ticks](x-.5)
		ax[:xaxis][:set_ticklabels](xlabels, rotation = 60, fontsize = 8)
		ax[:yaxis][:set_ticks](y+.5)
		ax[:yaxis][:set_ticklabels](labels, rotation = 0, fontsize = 8)
		savefig(string("figures/ParameterSpeciesInteractions",timestr, ".pdf"))

	end
end

function calculate_params_per_stepsize()
	close("all")
	Skips = [10,25,50,100,200,400,700,1000,1500]
	step_size = .01
	num_params_all = Float64[]
	num_params_selected = Float64[]
	for j in collect(1:size(Skips,1))
		curr_params_all,names_all = calculateSensitivityArr(Skips[j])
		curr_params_sel,names_sel = calculateSensitivityArr(Skips[j], [4,9])
		push!(num_params_all, curr_params_all)
		push!(num_params_selected, curr_params_sel)
		@show Skips[j]*step_size, names_all, names_sel
	end
	plot(Skips*step_size, num_params_all, "kx-")
	plot(Skips*step_size, num_params_selected, "rx-")
	xlabel("Time between measurement in arabitrary time units")
	ylabel("Number of Estimatable Parameters")
	legend(["Measuring All Species", "Only Measuring P3 and mRNA1"], loc="best")
	savefig("figures/NumberOfEstamatableParamsAsAFunctionOfFreq.pdf")
	axis([0,20, 0, 24])
end
