# FluSubsample_Rshiny
R Shiny App to create stratified subsamples of Genbank-NCBI unaligned Influenza sequence data.
Version: 26 Feb 2018

# Introduction
There is now a wealth of publically available Influenza virus sequence data, 
however to use it for evolutionary analyses it is often advisable to perform stratified subsampling.  
This type of subsampling takes allows a maximum of N sequences per category, and thus avoids having datasets too biased towards hosts/locations etc where there was intensive sampling effort.  
For Influenza a useful subsampling strategy might be to allow a maximum of 3(nper=3) sequences per category of Time,Location,Subtype,Host-type.  
Where Time can be Year, or maybe Year+Month, and Location could be Continent,Geographic Region,Individual Country, or Country or State for large countries (USA,Canada,China,Russian Federation).  
When doing Avian-influenza centric analyses, it is also often useful to subdivide the Host-type of Avian into bird orders, or at least domestic galliformes, domestic anseriformes and wild birds.
This R Shiny App is intended to perform configurable stratified subsampling on sequence data, and will output the subsamples and details of what was done.

# Inputs
This R Shiny App works on publically available Influenza virus sequence data downloaded from the NCBI Influenza Virus Resource 
https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=genomeset

The data must be in FASTA format and fasta header of the data must be:\n
'>{serotype}_{host}_{accession}_{strain}_{country}_{year}/{month}/{day}_{segname}'

i.e. you must click 'Customize FASTA Def line' after composing your query and before downloading.

The data can be either individual segment data, or whole genome data (the App will figure it out from the _{segname} field in the sequence names)

There is a small dataset of H5N5 avian influenza data in the example_data directory.

# Outputs
Sequence data will be processed by the App and Fasta format subsamples are written to the userGenerated directory.

For each subsample replicate, there is one fasta file and one corresponding traits file (tab separated) for each segment, 
and for each subsampling setting there the complete table corresponding to all the sequences is output (*.csv) as well as a table of whole genome only acession numbers and indices into the original sequence file (*.csv).  
The settings for the subsampling are also output as a *_settings.txt file.
