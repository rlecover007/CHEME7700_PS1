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
# Function: calculate_transcription_rates
# Description: Calculate the transcriptional rate array at time t
# Generated on: 2017-03-09T13:14:13.798
#
# Input arguments:
# t::Float64 => Current time value (scalar) 
# x::Array{Float64,1} => State array (number_of_species x 1) 
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model parameters 
#
# Output arguments:
# transcription_rate_array::Array{Float64,1} => Transcriptional rate array (number_of_genes x 1) at time t 
# ----------------------------------------------------------------------------------- #
function calculate_transcription_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

	# Alias the species - 
	gene_1 = x[1]
	gene_2 = x[2]
	gene_3 = x[3]

	# Initialize the transcription rate - 
	transcription_rate_array = zeros(3)
	KSAT = data_dictionary["saturation_constant_transcription"]
	kcat_transcription = data_dictionary["kcat_transcription"]
	rnapII_concentration = data_dictionary["rnapII_concentration"]
	average_transcript_length = data_dictionary["average_transcript_length"]
	gene_coding_length_array = data_dictionary["gene_coding_length_array"]

	# Populate the transcription rate array - 
	# Gene: gene_1
	gene_length = gene_coding_length_array[1]
	scale_factor = (average_transcript_length/gene_length)
	transcription_rate_array[1] = scale_factor*kcat_transcription*(rnapII_concentration)*((gene_1)/(KSAT+gene_1))

	# Gene: gene_2
	gene_length = gene_coding_length_array[2]
	scale_factor = (average_transcript_length/gene_length)
	transcription_rate_array[2] = scale_factor*kcat_transcription*(rnapII_concentration)*((gene_2)/(KSAT+gene_2))

	# Gene: gene_3
	gene_length = gene_coding_length_array[3]
	scale_factor = (average_transcript_length/gene_length)
	transcription_rate_array[3] = scale_factor*kcat_transcription*(rnapII_concentration)*((gene_3)/(KSAT+gene_3))


	# return transcription_rate_array - 
	return transcription_rate_array
end

# ----------------------------------------------------------------------------------- #
# Function: calculate_background_transcription_rates
# Description: Calculate the leak transcriptional rate array at time t
# Generated on: 2017-03-09T13:14:13.798
#
# Input arguments:
# t::Float64 => Current time value (scalar) 
# x::Array{Float64,1} => State array (number_of_species x 1) 
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model parameters 
#
# Output arguments:
# background_transcription_rate_array::Array{Float64,1} => Background transcriptional rate array (number_of_genes x 1) at time t 
# ----------------------------------------------------------------------------------- #
function calculate_background_transcription_rates(t::Float64,x::Array{Float64,1},transcription_rate_array::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})
	return zeros(length(x))
end


# ----------------------------------------------------------------------------------- #
# Function: calculate_translation_rates
# Description: Calculate the translation rate array at time t
# Generated on: 2017-03-09T13:14:13.798
#
# Input arguments:
# t::Float64 => Current time value (scalar) 
# x::Array{Float64,1} => State array (number_of_species x 1) 
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model parameters 
#
# Output arguments:
# translation_rate_array::Array{Float64,1} => Translation rate array (number_of_genes x 1) at time t 
# ----------------------------------------------------------------------------------- #
function calculate_translation_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

	# Alias the species - 
	mRNA_gene_1 = x[4]
	mRNA_gene_2 = x[5]
	mRNA_gene_3 = x[6]

	# Initialize the translation rate - 
	translation_rate_array = zeros(3)
	KSAT = data_dictionary["saturation_constant_translation"]
	kcat_translation = data_dictionary["kcat_translation"]
	ribosome_concentration = data_dictionary["ribosome_concentration"]
	average_protein_length = data_dictionary["average_protein_length"]
	protein_coding_length_array = data_dictionary["protein_coding_length_array"]

	# Populate the translation rate array - 
	# Transcript: mRNA_gene_1
	protein_length = protein_coding_length_array[1]
	scale_factor = (average_protein_length/protein_length)
	translation_rate_array[1] = scale_factor*kcat_translation*(ribosome_concentration)*((mRNA_gene_1)/(KSAT+mRNA_gene_1))

	# Transcript: mRNA_gene_2
	protein_length = protein_coding_length_array[2]
	scale_factor = (average_protein_length/protein_length)
	translation_rate_array[2] = scale_factor*kcat_translation*(ribosome_concentration)*((mRNA_gene_2)/(KSAT+mRNA_gene_2))

	# Transcript: mRNA_gene_3
	protein_length = protein_coding_length_array[3]
	scale_factor = (average_protein_length/protein_length)
	translation_rate_array[3] = scale_factor*kcat_translation*(ribosome_concentration)*((mRNA_gene_3)/(KSAT+mRNA_gene_3))


	# return translation array - 
	return translation_rate_array
end

# ----------------------------------------------------------------------------------- #
# Function: calculate_mRNA_degradation_rates
# Description: Calculate the mRNA degradation rate array at time t
# Generated on: 2017-03-09T13:14:13.798
#
# Input arguments:
# t::Float64 => Current time value (scalar) 
# x::Array{Float64,1} => State array (number_of_species x 1) 
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model parameters 
#
# Output arguments:
# mRNA_degradation_rate_array::Array{Float64,1} => mRNA degradation rate array (number_of_genes x 1) at time t 
# ----------------------------------------------------------------------------------- #
function calculate_mRNA_degradation_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

	# Alias the species - 
	mRNA_gene_1 = x[4]
	mRNA_gene_2 = x[5]
	mRNA_gene_3 = x[6]

	# Initialize the degrdation array - 
	degradation_rate_array = zeros(3)
	mRNA_degrdation_constant = data_dictionary["degradation_constant_mRNA"]

	# Calculate the degradation_rate_array - 
	degradation_rate_array[1] = (mRNA_degrdation_constant)*mRNA_gene_1
	degradation_rate_array[2] = (mRNA_degrdation_constant)*mRNA_gene_2
	degradation_rate_array[3] = (mRNA_degrdation_constant)*mRNA_gene_3

	# return the degrdation rate array - 
	return degradation_rate_array
end

# ----------------------------------------------------------------------------------- #
# Function: calculate_protein_degradation_rates
# Description: Calculate the protein degradation rate array at time t
# Generated on: 2017-03-09T13:14:13.798
#
# Input arguments:
# t::Float64 => Current time value (scalar) 
# x::Array{Float64,1} => State array (number_of_species x 1) 
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model parameters 
#
# Output arguments:
# protein_degradation_rate_array::Array{Float64,1} => protein degradation rate array (number_of_proteins x 1) at time t 
# ----------------------------------------------------------------------------------- #
function calculate_protein_degradation_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

	# Alias the species - 
	protein_gene_1 = x[7]
	protein_gene_2 = x[8]
	protein_gene_3 = x[9]

	# Initialize the degrdation array - 
	degradation_rate_array = zeros(3)
	protein_degrdation_constant = data_dictionary["degradation_constant_protein"]

	# Calculate the degradation_rate_array - 
	degradation_rate_array[1] = (protein_degrdation_constant)*protein_gene_1
	degradation_rate_array[2] = (protein_degrdation_constant)*protein_gene_2
	degradation_rate_array[3] = (protein_degrdation_constant)*protein_gene_3

	# return the degrdation rate array - 
	return degradation_rate_array
end


