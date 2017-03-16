using NLopt
using DataStructures
include("Include.jl")
include("runModel.jl")
expdata = readdlm("../expdata/Washout-Data-PS1-Q3.dat")
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
		else
			sum = sum+(model[j]-measured[j])^2/1.0
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
#	SSEP1 = calculateSSE(X[:,7]/1000, expdata[:,2], expdata[:,3])
#	SSEP2 = calculateSSE(X[:,8]/1000, expdata[:,4], expdata[:,5])
#	SSEP3 = calculateSSE(X[:,9]/1000, expdata[:,6], expdata[:,7])
#	trueSSE = SSEP1+SSEP2+SSEP3
#	@show trueSSE
#	return trueSSE
	MSE1 = calculateMSE(X[:,7]/1000, expdata[:,2])
	MSE2 = calculateMSE(X[:,8]/1000, expdata[:,6])
	MSE3 = calculateMSE(X[:,9]/1000, expdata[:,4])
	totalMSE = MSE1+MSE2+MSE3
	@show totalMSE
	return totalMSE
end

function runOpt()
	time_start = 0.0
	time_stop = 24.0
	time_step_size = 0.01
	initial_dict = DataDictionary(time_start,time_stop,time_step_size)
	numvars = length(initial_dict["parameter_name_mapping_array"])+1
	@show numvars
	opt = Opt(:LN_NELDERMEAD,numvars)
	min_objective!(opt, objective)
	binding = DataStructures.SortedDict(initial_dict["binding_parameter_dictionary"])
	control = DataStructures.SortedDict(initial_dict["control_parameter_dictionary"])
	misc = DataStructures.SortedDict(initial_dict["misc_parameter_dictionary"])
	initial_parameter_estimate = Float64[]
	total = DataStructures.SortedDict(merge(binding,control, misc))
	@show total
	j = 1
	for (index, value) in total
#		@show j, index
		push!(initial_parameter_estimate, value)
		j = j+1
	end
	
	@show size(initial_parameter_estimate)
	lower_bounds!(opt, vec(initial_parameter_estimate)*1E-2)
	upper_bounds!(opt, vec(initial_parameter_estimate)*100)
	@show opt
	(minf, minx, ret) = NLopt.optimize(opt, vec(initial_parameter_estimate))
	#(minf, minx, ret) = NLopt.optimize(opt, vec(ones(numvars,1)))
	println("got $minf at $minx after $count iterations (returned $ret)")
#	
end

function makePostNMplot()
	close("all")
	#params = readdlm("bestparamsAfterNM_MinMSE_NetworkLayout1.txt", ',')
	params = [248.583,40.6087,438.076,1.0e-6,0.0,1.00186e-6,0.194626,800.0,1041.94,424.806,822.441,0.99021,17819.7,1.75725e5,209.662,75.4871,0.724132,3.15911,2.00168e5,14017.6,8878.05,30391.5]
	time_start = 0.0
	time_stop = 24.0
	time_step_size = 0.01
	data_dict = DataDictionary(time_start, time_stop, time_step_size, params)
	numvars = length(data_dict["parameter_name_mapping_array"])+1
	T,X = runWashout(time_start, time_stop, time_step_size, data_dict)
	selectedX = X[1:10:end, :]
	SSEP1 = calculateSSE(selectedX[:,7]/1000, expdata[:,2], expdata[:,3])
	SSEP2 = calculateSSE(selectedX[:,8]/1000, expdata[:,6], expdata[:,7])
	SSEP3 = calculateSSE(selectedX[:,9]/1000, expdata[:,4], expdata[:,5])
	trueSSE = SSEP1+SSEP2+SSEP3
	AIC = 2*numvars-2*log(trueSSE)
	figure()
	plot(T, X[:,7]/1000, "b-")
	plot(T, X[:,8]/1000, "g-")
	plot(T, X[:,9]/1000, "y-")
	plot(expdata[:,1],expdata[:,2], "bx")
	plot(expdata[:,1],expdata[:,6], "gx") #if varner sent wrong data
	plot(expdata[:,1],expdata[:,4], "yx")
	legarr = ["Predicted P1", "Predicted P2", "Predicted P3", "Experimental P1", "Experimental P2", "Experimental P3"]
	legend(legarr, loc="best")
	ylabel("Concentration, microMolar")
	xlabel("Time, hours")
	annotate(string("AIC = ", AIC), xy=[.15, .85], xycoords="figure fraction")
	savefig("../figures/postNMPlotLayout1.pdf")
end
