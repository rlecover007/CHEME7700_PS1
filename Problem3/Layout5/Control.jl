# ----------------------------------------------------------------------------------- #
# Copyright (c) 2017 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #
#
# ----------------------------------------------------------------------------------- #
# Function: Control
# Description: Calculate the transcriptional control array at time t
# Generated on: 2017-03-14T22:23:49.222
#
# Input arguments:
# t::Float64 => Current time value (scalar) 
# x::Array{Float64,1} => State array (number_of_species x 1) 
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model parameters 
#
# Output arguments:
# control_array::Array{Float64,1} => Transcriptional control array (number_of_genes x 1) at time t 
# ----------------------------------------------------------------------------------- #
function Control(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

	# initialize the control - 
	control_array = zeros(3)

	# Alias the species - 
	gene_1 = x[1]
	gene_2 = x[2]
	gene_3 = x[3]
	mRNA_gene_1 = x[4]
	mRNA_gene_2 = x[5]
	mRNA_gene_3 = x[6]
	protein_gene_1 = x[7]
	protein_gene_2 = x[8]
	protein_gene_3 = x[9]

	# Alias the binding parameters - 
	binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
	n_gene_2_gene_1 = binding_parameter_dictionary["n_gene_2_gene_1"]
	K_gene_2_gene_1 = binding_parameter_dictionary["K_gene_2_gene_1"]
	n_gene_3_gene_2 = binding_parameter_dictionary["n_gene_3_gene_2"]
	K_gene_3_gene_2 = binding_parameter_dictionary["K_gene_3_gene_2"]

	# Alias the control function parameters - 
	control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
	W_gene_1_RNAP = control_parameter_dictionary["W_gene_1_RNAP"]
	W_gene_2_RNAP = control_parameter_dictionary["W_gene_2_RNAP"]
	W_gene_2_gene_1 = control_parameter_dictionary["W_gene_2_gene_1"]
	W_gene_3_RNAP = control_parameter_dictionary["W_gene_3_RNAP"]
	W_gene_3_gene_2 = control_parameter_dictionary["W_gene_3_gene_2"]

	# Control function for gene_1 - 
	control_array[1] = (W_gene_1_RNAP)/(1+W_gene_1_RNAP)

	# Transfer function target:gene_2 actor:gene_1
	actor_set_gene_2_gene_1 = [
		protein_gene_1
	]
	actor = prod(actor_set_gene_2_gene_1)
	b_gene_2_gene_1 = (actor^(n_gene_2_gene_1))/(K_gene_2_gene_1^(n_gene_2_gene_1)+actor^(n_gene_2_gene_1))

	# Control function for gene_2 - 
	control_array[2] = (W_gene_2_RNAP+W_gene_2_gene_1*b_gene_2_gene_1)/(1+W_gene_2_RNAP+W_gene_2_gene_1*b_gene_2_gene_1)

	# Transfer function target:gene_3 actor:gene_2
	actor_set_gene_3_gene_2 = [
		protein_gene_2
	]
	actor = prod(actor_set_gene_3_gene_2)
	b_gene_3_gene_2 = (actor^(n_gene_3_gene_2))/(K_gene_3_gene_2^(n_gene_3_gene_2)+actor^(n_gene_3_gene_2))

	# Control function for gene_3 - 
	control_array[3] = (W_gene_3_RNAP)/(1+W_gene_3_RNAP+W_gene_3_gene_2*b_gene_3_gene_2)

	# return - 
	return control_array
end
