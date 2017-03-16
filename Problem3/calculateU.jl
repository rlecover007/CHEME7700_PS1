using NLopt
using DataStructures
include("DataDictionary.jl")

function calculate_u(X1, X2, tstart, tend, tstep)
	t = collect(tstart:tstep:tend)
#	@show size(t)
#	@show any(isnan, X1)
#	@show any(isnan, X2)
	num_species = 3 #we are only going to compare performance on the proteins that we're measuring
	idx_small = find(X1.<1E-9)
  	X1[idx_small] = 0.0
	idx_small = find(X2.<1E-9)
  	X2[idx_small] = 0.0
	P1_1 = X1[:,7]
	P1_2 = X1[:,8]
	P1_3 = X1[:,9]
	P2_1 = X2[:,7]
	P2_2 = X2[:,8]
	P2_3 = X2[:,9]

	#there is a way to generalize this, but I don't presently want to spend a lot of time figuring it out'
#	@show size((P1_1-P2_1).^2./((P1_1+P2_1)/2).^2)#
#	@show size((P1_2-P2_2).^2./((P1_2+P2_2)/2).^2)
#	@show size((P1_3-P2_3).^2./((P1_3+P2_3)/2).^2)
	W = (P1_1-P2_1).^2./((P1_1+P2_1)/2).^2.+(P1_2-P2_2).^2./((P1_2+P2_2)/2).^2.+(P1_3-P2_3).^2./((P1_3+P2_3)/2).^2
	idxNan = find(isnan, W)
	W[idxNan] =0.0
#	@show W
#	@show size(W)
#	@show size(vec(W))
#	@show any(isnan, W), any(isnan, t)
	u=trapz(t, vec(W))
	return u
	
end

function objWrapper(params::Vector, grad::Vector)
	tstart = 0.0
	tstep = .01
	tend = 24.0
	@show params
	cd("Layout1/")
	run(`echo $params`)
	run(`julia runModelFromScript.jl  $tstart $tend $tstep $params`)
	X1 = readdlm("Layout1Output.txt")
	cd("..")
	cd("Layout5/")
	run(`julia runModelFromScript.jl  $tstart $tend $tstep $params`)
	X2 = readdlm("Layout5Output.txt")
	cd("..")
	u = calculate_u(X1, X2, tstart, tend, tstep)
	@show u
	return u
end

function objWrapper(params::Vector)
	tstart = 0.0
	tstep = .01
	tend = 24.0
	cd("Layout1/")
	run(`echo $params`)
	run(`julia runModelFromScript.jl  $tstart $tend $tstep $params`)
	X1 = readdlm("Layout1Output.txt")
	cd("..")
	@show size(X1)
	cd("Layout5/")
	run(`julia runModelFromScript.jl  $tstart $tend $tstep $params`)
	X2 = readdlm("Layout5Output.txt")
	cd("..")
	u = calculate_u(X1, X2, tstart, tend, tstep)
	return u
end

function trapz{Tx<:Number, Ty<:Number}(x::Vector{Tx}, y::Vector{Ty})
    # Trapezoidal integration rule
    local n = length(x)
    if (length(y) != n)
        error("Vectors 'x', 'y' must be of same length")
    end
    r = zero(zero(Tx) + zero(Ty))
    if n == 1; return r; end
    for i in 2:n
        r += (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    #= correction -h^2/12 * (f'(b) - f'(a))
    ha = x[2] - x[1]
    he = x[end] - x[end-1]
    ra = (y[2] - y[1]) / ha
    re = (y[end] - y[end-1]) / he
    r/2 - ha*he/12 * (re - ra)
    =#
    return r/2
end

function runOpt()
	time_start = 0.0
	time_stop = 24.0
	time_step_size = 0.01
	initial_dict = DataDictionary(time_start,time_stop,time_step_size)
	numvars = 6
	@show numvars
	opt = Opt(:LN_NELDERMEAD,numvars)
	max_objective!(opt, objWrapper)
	lower_bounds = [.001, .001, .001, 0, 0, 0]
	upper_bounds = ones(6)
	lower_bounds!(opt, vec(lower_bounds))
	upper_bounds!(opt, vec(upper_bounds))
	@show opt
	(minf, minx, ret) = NLopt.optimize(opt, vec(lower_bounds+upper_bounds)/2)
	#(minf, minx, ret) = NLopt.optimize(opt, vec(ones(numvars,1)))
	println("got $minf at $minx after $count iterations (returned $ret)")
	#got 81.50600675275699 at [0.96765,0.001,0.0940649,0.994224,0.999786,0.892191] after count iterations (returned XTOL_REACHED) for comparing 1 and 3
end


