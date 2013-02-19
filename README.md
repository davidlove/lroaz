lroaz
=====

Implementation of a Likelihood Robust Two-Stage Linear Program with Recourse for Tucson water data

=====

**DR. LANSEY'S STUDENTS:**

For the code to build an LP based on the excel spreadsheets, you want the *LPModel* class, contained in the folder @LPModel

The file *one_stage_az.m* shows you the basic use of the LPModel class.

I currently have LPModel set to never write to Excel (it doesn't work on Linux).  To enable writing to excel, open @LPModel/LPModel.m, and set the variable writeToExcel = true.  This is line 52 as of the time I'm writing this.
