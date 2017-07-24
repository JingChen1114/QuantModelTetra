
# [Download QuantModelTetra](https://github.com/LJLeach/QuantModelTetra/archive/master.zip)

<h2> Contents </h2> 

#### [Introduction](https://github.com/LJLeach/QuantModelTetra#-introduction-)

#### Basic Usage

##### [Usage Instructions](https://github.com/LJLeach/QuantModelTetra#-usage-instructions-)

##### [Input Files](https://github.com/LJLeach/QuantModelTetra#-input-files-)

##### [Output Files](https://github.com/LJLeach/QuantModelTetra#-output-files-)

##### [Example Analysis](https://github.com/LJLeach/QuantModelTetra#-example-analysis-)

##### [Population Simulation](https://github.com/LJLeach/QuantModelTetra#-population-simulation-)

##### [Extras](https://github.com/LJLeach/QuantModelTetra#-extras-)

#### Misc.

##### [References](https://github.com/LJLeach/QuantModelTetra#-references-)


<h2> Introduction </h2>

QuantModelTetra is a series of related programs that illustrate the orthogonal contrast based model for quantitative genetic analysis in autoetraploid species and its practical application in trait segregation analysis to detect segregation of major QTL in an autotetraploid population.

<h2> Usage Instructions </h2>

QuantModelTetra is supplied as a set of 6 Fortran 90 (.for) source files in the folder FortranSource, including:

<b> main.for </b> 
The main fortran program for carrying out trait segregation analysis

<b> gametogenesis.for </b> 
Subroutine for simulating gametogenesis (meiosis) in any given autotetraploid parental genotype

<b> offspring.for </b> 
Subroutines for generating an offspring population from any two given parental autotetraploid genotypes

<b> orthogonalscales.for </b> 
Subroutine for calculating the orthogonal scales in the quantitative genetic model 

<b> ReducedModel.for </b> 
Subroutine for implementing the orthogonal scales model when the number of possible offspring genotypes is less than the maximum (5). 

To use the program an Intel Fortran Compiler and the RogueWave ISML 7 Fortran Numerical Library is required. The programs can be be compiled and built under Linux or Windows operating systems.
An example is given here for how to compile, build and run the programs using the Microsoft Visual Studio IDE in Windows.

In Microsoft Developer Studio, open an empty project, give it a name, for example "MySegregationAnalysis" and save to the desired location. Add all 6 *.for files as source files to the project.
Build the solution from the Build menu and check there are no errors in the output window. 

<h2> Input files </h2>

<b> data_in.txt </b> 
A single input file containing the trait phenotype data for the segregating population. A single header line contains the number of individuals in the population. The trait phenotype score for each individual is given on separate lines. 
Two example input files are provided in the folder "RealPotatoPhenotypeData", including "data_in_FloweringTime" and "data_in_PlantHeight".

<h2> Output files </h2>

6 output files are generated including

<b> data_out.txt </b> 
Contains results of all 600 runs of the trait segregation analysis - across all 12 possible parental genotype configurations and over 50 different levels of double reduction in steps of 0.005, from its minimum value of 0 to the maximum of 0.25.

<b> LOD_scores.txt </b> 
Contains LOD scores for all 600 runs of the trait segregation analysis.

<b> AIC.txt </b> 
Contains the Akaike Information Criterion scores for all 600 runs of the trait segregation analysis.

<b> Likelihood.txt </b> 
Contains likelihood scores for all 600 runs of the trait segregation analysis.

<b> Chi_square.txt </b> 
Contains chi square test values for all 600 runs of the trait segregation analysis.

<b> IOff.txt </b> 
Contains the number of offspring genotypes (i = 2,3,4,5) for all 600 runs of the trait segregation analysis.

<h2> Example Analysis </h2>

<b> Analysis of Flowering Time </b>
The program can be run using the input file "data_in_FloweringTime.txt" which contains a header indicating 297 offspring individuals followed by 297 lines containing the average trait score for 297 individuals.
For example, using Microsoft Developer Studio, from the Debug menu select "Start debugging" to run the code, which should complete within 1 minute on a standard desktop computer.
The output files are included here in the folder "FloweringTimeResults".

<h2> Population Simulation </h2>
Two sets of Fortran 90 source files are provided for simulation of a segregating population under either bivalent oe quadrivalent settings for chromosome pairing in an autotetraploid meiosis.
To simulate offspring populations under either setting, the programs need to be compiled and built as a single fortran project with the main program given as main.for. 
The input file "data_in" contains the simulation parameters. 
The output file "data_out" contains the estimated parameters of the quantitative genetic model.
The output file "info_mapping" contains parental genotypes, recombination frequency (and double reduction in the quadrivalent setting) parameters, and the genotype and phenotype scores for the simulated offspring individuals.

<h2> Extras </h2>

<h6> R_Scripts/FittingModel_FloweringTime.r </h6>A simple script to produce Figure 1e in Chen et al (2017).
<h6> R_Scripts/FittingModel_PlantHeight.r </h6>A simple script to produce Figure 1f in Chen et al (2017).
<h6> R_Scripts/Profile_Of_LOD_Scores.r </h6>A simple script to produce Supplementary Figure S1 in Chen et al (2017).
<h6> Heritability/Heritability.r </h6>A script using REML to calculate the narrow sense heritability of quantitative traits using plant height as an example in RawData_PlantHeight.txt.

<h3> References </h3>

Chen, J, Zhang, F, Wang L, Leach L and Luo, Z. "Orthogonal contrast based models for quantitative genetic analysis ins autotetraploid species." In prep.
