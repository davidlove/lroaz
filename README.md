lroaz
=====

Implementation of a Phi-Divergence Two-Stage Linear Program with Recourse (PhiLP-2) for Tucson water data.  This project originally started life as an implementation fo the Likelihood Robust Two-Stage Linear Program with Recourse, but was later generalized to include other phi-divergences as well.

Classes are as follows:

* *@LPMODEL* contains the structure of a two-stage stochastic linear program with recourse.
* *@LRLP* implements the modified Bender's Decomposition used to solve the PhiLP-2 model, formerly the Likelihood Robust model.
* *@PhiDivergence* is a class to contain the phi-divergences and their important properties.
* *@Solution* class is used to store the candidate and best solution as the modified Bender's Decomposition algorithm runs.

A few functions of interest:

* *SolveLRLP.m* implements the basic Bender's Decomposition loop using the LRLP class
* *GatherLROData.m* loops over a variety of parameter values to collect data on optimal solution of the PhiLP-2 problem.

=====

**DR. LANSEY'S STUDENTS:**

For the code to build an LP based on the excel spreadsheets, you want the *LPModel* class, contained in the folder @LPModel

The file *OneStageAZ.m* shows you the basic use of the LPModel class.

The *doc* folder contains some technical documentation on how LPModel constructs the water model.

I currently have LPModel set to never write to Excel (it doesn't work on Linux).  To enable writing to excel, open @LPModel/LPModel.m, and set the variable writeToExcel = true.  This is line 52 as of the time I'm writing this.

For reference, the files in LPModel:
* *LPModel.m* contains the basic information about LPModel--data and some small methods.
* *BuildA.m* builds the submatrices A, A_st and A_log.  Alicia can tell you all about those.
* *ExpandA.m* creates the full constraint matrix out of A, A_st and A_log.
* *BuildVectors.m* Builds the cost, right hand side, and lower and upper bound vectors.  It should work for both one and two stage problems.
	*NOTE* BuildVectors will only work for Time_Lag = 1 or Time_Lag = 12.  This is based on what I see in sf.xlsx
* *GenerateDiagonal.m* is used by ExpandA.m to build the full constraint matrix.  Don't modify this file.
* *PrintConstraint.m* prints the constraints in a human readable format, which should be helpful for debugging any changes to the code, or creating new spreadsheets.
