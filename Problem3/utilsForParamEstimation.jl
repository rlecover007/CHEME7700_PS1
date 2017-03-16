using NLopt
using DataStructures
include("Include.jl")
include("runModel.jl")
expdata = readdlm("expdata/Washout-Data-PS1-Q3.dat")
  idx_small = find(expdata.<1E-15) #make anything less than 1E-15 zero
  expdata[idx_small] = 0.0

function calculateMSE(model, measured)
	numPoints = size(model,1)
	sum = 0.0
	for j in collect(1:numPoints)
		sum = sum+(model[j]-measured[j])^2
	end
	MSE = sum/numPoints
	return MSE
	
end

function calculateSSE(model, measured, uncertainty)
	numPoints = size(model,1)
	sum = 0.0
	for j in collect(1:numPoints)
		if(uncertainty[j]>0) #so we don't divide by zero
			sum = sum+(model[j]-measured[j])^2/uncertainty[j]
		end
	end
	SSE = sum
	return SSE
	
end

function objective(params::Vector, grad::Vector)
	@show params
	#need to write something to build dictionary from parameters
	time_start = 0.0
	time_stop = 24.0
	time_step_size = 0.01
	data_dict = DataDictionary(time_start, time_stop, time_step_size, params)
#	@show values(data_dict["binding_parameter_dictionary"])
#	@show values(data_dict["control_parameter_dictionary"])
#	@show values(data_dict["misc_parameter_dictionary"])
	T,X = runWashout(time_start, time_stop, time_step_size, data_dict)
	#not getting any dynamics-why?
	X = X[1:10:end, :] #pull out the points for which we have measured data
	#calculate SSE between P1, P2, P3
	SSEP1 = calculateSSE(X[:,7], expdata[:,2], expdata[:,3])
	SSEP2 = calculateSSE(X[:,8], expdata[:,4], expdata[:,5])
	SSEP3 = calculateSSE(X[:,9], expdata[:,6], expdata[:,7])
	trueSSE = SSEP1+SSEP2+SSEP3
	@show trueSSE
	return trueSSE
end

function runOpt()
	time_start = 0.0
	time_stop = 24.0
	time_step_size = 0.01
	initial_dict = DataDictionary(time_start,time_stop,time_step_size)
	numvars = length(initial_dict["parameter_name_mapping_array"])+1
	@show numvars
	opt = Opt(:LN_NELDERMEAD,numvars)
	lower_bounds!(opt, vec(fill(1E-9,1,numvars)))
	upper_bounds!(opt, vec(fill(1E7,1,numvars)))
	min_objective!(opt, objective)
	binding = DataStructures.SortedDict(initial_dict["binding_parameter_dictionary"])
	control = DataStructures.SortedDict(initial_dict["control_parameter_dictionary"])
	misc = DataStructures.SortedDict(initial_dict["misc_parameter_dictionary"])
	initial_parameter_estimate = Float64[]
	total = DataStructures.SortedDict(merge(binding,control, misc))
	@show total
	j = 1
	for (index, value) in total
		@show j, index
		push!(initial_parameter_estimate, value)
		j = j+1
	end
	
#	@show initial_parameter_estimate
	(minf, minx, ret) = NLopt.optimize(opt, vec(initial_parameter_estimate))
	println("got $minf at $minx after $count iterations (returned $ret)")
#	
end
