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
	opt = Opt(:LN_COBYLA,numvars)
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
	
	@show size(initial_parameter_estimate)
	lower_bounds!(opt, vec(initial_parameter_estimate)*1E-2)
	upper_bounds!(opt, vec(initial_parameter_estimate)*100)
	@show opt
	(minf, minx, ret) = NLopt.optimize(opt, vec(initial_parameter_estimate))
	println("got $minf at $minx after $count iterations (returned $ret)")
#	
end

function makePostNMplot()
	close("all")	
	time_start = 0.0
	time_stop = 24.0
	time_step_size = 0.01
	#params =[4155.79,6133.31,532.547,220.169,0.0,2.46427,0.0,2.15024,179.448,0.05,0.015,6.07073,10.2806,0.00424873,216.397,606.734,1.8171,3.98848,2.24524,0.653426,1.34155,38449.3,8137.06,8549.89,92498.4]
	params = [3773.63,5392.53,397.219,188.634,0.0,2.01276,0.0,1.83374,140.389,0.5,0.015,3.35037,6.3188,0.0112672,269.936,365.924,1.29386,4.55393,2.35943,0.703358,1.70803,1.63354e5,6623.33,5054.6,633899.0]	
	data_dict = DataDictionary(time_start, time_stop, time_step_size, params)
	#data_dict = DataDictionary(time_start, time_stop, time_step_size)
	T,X = runWashout(time_start, time_stop, time_step_size, data_dict)
	@show size(X)
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
	savefig("../figures/postNMPlotLayout4.pdf")
end
