using ODE
using PyPlot

function BalanceEquations(t,x,data_dict)
	 idx_small = find(x.<0)
 	 x[idx_small] = 0.0
	#amount of gene is constant
	params = data_dict["params"]
	kx = params[1]
	kt =params[2]
	G = params[3]
	kd = params[4]
	kdbar =params[5]
	Lt = params[6]
	Kt = params[7]
	Kx = params[8]
	mubar = params[9]
	Rt =params[10]
	Rx = params[11]
	Lx= params[12]

	specifics = data_dict["specifics"]
	Lt_j = specifics[1]
	k2 = specifics[2]
	k3 = specifics[3]
	ftl = specifics[4]
	k2x = specifics[5]
	k3x = specifics[6]
	ftlX = specifics[7]
	ft = 1-ftl
	#assume 3 bp = 1AA
	Lx_j = Lt_j/3

	m = x[1]
	p =x[2]

	rtbar = kt*Rt*(Lt/Lt_j)*(G/(Kt+G))
	u = k2*ftl/(1+k2*ftl+k3*ft)
	rt = rtbar*u


	w=k2x*ftlX/(1+k2x*ftlX+k3x*(1-ftlX))
	rx = kx*Rx*(Lx/Lx_j)*(m/(Kx+m))*w
	
	dmdt = rt-(kd+mubar)*m
	dpdt = rx-(kdbar+mubar)*p
	return [dmdt;dpdt]
end

function DataDictionary()
	data_dict = Dict()
	params = Float64[];
	# ------------------------------------------------------------------------------------------#
	# constants (from bionumbers)       units
	# ------------------------------------------------------------------------------------------#
	cell_diameter = 1.1                 # mum
	number_of_rnapII = 4600            	# copies/cells
	number_of_ribosome = 50000         	# copies/cells
	mRNA_half_life_TF = 0.083           # hrs
	protein_half_life = 70              # hrs
	doubling_time_cell = 0.33           # hrs
	max_translation_rate = 16.5         # aa/sec
	max_transcription_rate = 60.0       # nt/sec
	average_transcript_length = 1000   	# nt
	average_protein_length = 400       	# aa
	fraction_nucleus = 0.0             	# dimensionless
	av_number = 6.02e23                 # number/mol
	avg_gene_number = 2                 # number of copies of a gene
	polysome_number = 4									# number of ribsomoses per transcript
	# ------------------------------------------------------------------------------------------#
	#
	# ------------------------------------------------------------------------------------------#
	# Calculate constants using bionumber values
	# ------------------------------------------------------------------------------------------#
	# Calculate the volume (convert to L)
	V = ((1-fraction_nucleus)*(1/6)*(3.14159)*(cell_diameter)^3)*(1e-15)

	# Calculate the rnapII_concentration and ribosome_concentration
	rnapII_concentration = number_of_rnapII*(1/av_number)*(1/V)*1e9                   # nM
	ribosome_concentration = number_of_ribosome*(1/av_number)*(1/V)*1e9               # nM

	# degrdation rate constants -
	degradation_constant_mRNA = -(1/mRNA_half_life_TF)*log(0.5)                       # hr^-1
	degradation_constant_protein = -(1/protein_half_life)*log(0.5)                    # hr^-1

	# kcats for transcription and translation -
	kcat_transcription = max_transcription_rate*(3600/average_transcript_length)      # hr^-1
	kcat_translation = polysome_number*max_translation_rate*(3600/average_protein_length)             # hr^-1

	# Maximum specific growth rate -
	maximum_specific_growth_rate = (1/doubling_time_cell)*log(2)                      # hr^-1

	# What is the average gene concentration -
	avg_gene_concentration = avg_gene_number*(1/av_number)*(1/V)*1e9                  # nM

	# How fast do my cells die?
	death_rate_constant = 0.05*maximum_specific_growth_rate                            # hr^-1

	# Saturation constants for translation and trascription -
	saturation_transcription = 4600*(1/av_number)*(1/V)*1e9                           # nM
	saturation_translation = 150000*(1/av_number)*(1/V)*1e9                           # nM
	# -------------------------------------------------------------------------------------------# 
	push!(params, kcat_translation)
	push!(params, kcat_transcription)
	push!(params, avg_gene_concentration)
	push!(params, degradation_constant_mRNA)
	push!(params, degradation_constant_protein*1.1)
	push!(params, average_transcript_length); #Lt
	push!(params, saturation_transcription)
	push!(params, saturation_translation)
	push!(params, maximum_specific_growth_rate)
	push!(params,rnapII_concentration)
	push!(params, ribosome_concentration)
	push!(params, average_protein_length)
	data_dict["params"]=params
	
	specifics = Float64[]
	push!(specifics, 100) #Lt_j
	push!(specifics, .1) #k2
	push!(specifics, .1) #k3
	push!(specifics, .5) #ftl
	push!(specifics, .1) #k2x
	push!(specifics, .1) #k3x
	push!(specifics, .5) #ftlX

	data_dict["specifics"] = specifics
	return data_dict
end

function makePlot(t,x, stylestr)
#	figure()
	plt[:subplot](2,1,1)
	axis([0,10,0,50])
	plot(t, ([a[1] for a in x]), string("r", stylestr))
	title("mRNA")
	plt[:subplot](2,1,2)
	axis([0,10,0,400])
	plot(t, ([a[2] for a in x]), string("b", stylestr))
	title("Protein")
	
end

function main()
	close("all")
	time_start = 0.0
	time_stop = 1.0
	time_step_size = 0.0001
	data_dictionary = DataDictionary()
	fed_equations(t,x) = BalanceEquations(t,x,data_dictionary)
	T1,X1 = ODE.ode23(fed_equations, [1.0;0], collect(time_start:time_step_size:time_stop), abstol = 1E-6, reltol =1E-6)
	makePlot(T1*60,X1, "-")
	data_dictionary["specifics"][3] =data_dictionary["specifics"][3]*10 #make k3 large to correspond to a lot of open DNA 
	data_dictionary["specifics"][6] =data_dictionary["specifics"][6]*10 #make k3x large to correspond to a lot of unbound RNA-p 
	fed_equations(t,x) = BalanceEquations(t,x,data_dictionary)
	T2,X2 = ODE.ode23(fed_equations, [1.0;0], collect(time_start:time_step_size:time_stop), abstol = 1E-6, reltol =1E-6)
	makePlot(T2*60,X2, "-.")
	legend(["Case 1", "Case 2"], loc="best")
	xlabel("Time, Minues")
	ylabel("Concentration, nM")
	savefig("Problem1NotQuiteLinear.pdf")

end
