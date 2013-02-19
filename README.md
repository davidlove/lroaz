lroaz
=====

Implementation of a Likelihood Robust Two-Stage Linear Program with Recourse for Tucson water data

=====

**DR. LANSEY'S STUDENTS:**

For the code to build an LP based on the excel spreadsheets, you want the *LPModel* class, contained in the folder @LPModel

The file *one_stage_az.m* shows you the basic use of the LPModel class.

I currently have LPModel set to never write to Excel (it doesn't work on Linux).  To enable writing to excel, open @LPModel/LPModel.m, and set the variable writeToExcel = true.  This is line 52 as of the time I'm writing this.

For reference, the files in LPModel:
* LPModel.m contains the basic information about LPModel--data and some small methods.
* BuildA.m builds the submatrices A, A_st and A_log.  Alicia can tell you all about those.
* ExpandA.m creates the full constraint matrix out of A, A_st and A_log.
* BuildVectors.m Builds the cost, right hand side, and lower and upper bound vectors.  It should work for both one and two stage problems.
	*NOTE* BuildVectors will only work for Time_Lag = 1 or Time_Lag = 12.  This is based on what I see in sf.xlsx
* GenerateDiagonal.m is used by ExpandA.m to build the full constraint matrix.  Don't modify this file.
