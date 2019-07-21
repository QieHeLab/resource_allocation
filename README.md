Overview

This directory contains code necessary to run all the numerical experiment in our paper: 

	"A fast algorithm for separable convex resource allocation with nested constraints", Zeyang Wu, Qie He, Kameng Nip, In preparation for Submission to INFORMS Journal on Computing

There are three experiments:

1. Numerical experiment on DCA and Gurobi on DRAP-NC with linear and quadratic objectives
	Test file: TestGurobiDCA_Lin.java, TestGurobiDCA_Qua.java, TestSparseGurobiDCA_Lin, TestSparseGurobiDCA_Qua
	Test instances: DRAP-NC with linear and quadratic objectives;
	Benchmark: Gurobi
	Requirements: Gurobi is installed 

2. Numerical experiment on DCA and MDA on DRAP-NC with three benchmark convex objectives function
	Test file: TestFlow1Con
	Test instances: DRAP-NC with three classes of benchmark convex objectives, [F], [CRASH], and [FUEL]
	Benchmark: MDA
	Requirements: None

3. Numerical experiment on DCA and MDA on SOVREX 
	Test file: TestInseparable.java
	Test instances: The RAP-NC projection subproblems in the eight instances of SOVREX
	Benchmark: MDA
	Requirements: all the data files in the folder "data" 


Running the code

Here are the short description on the files


1. Numerical experiment on DCA and Gurobi on DRAP-NC with linear and quadratic objectives
	To be included

2. Numerical experiment on DCA and MDA on DRAP-NC with three benchmark convex objectives function
	To be included

3. Numerical experiment on DCA and MDA on SOVREX 
	Our last numerical experiment is on the Support Vector Ordinal Regression (SOVREX) model. In this algorithm, we apply DCA and MDA to solve RAP-NC (the continuous problem) subproblems in a projected gradient descent method for SOVREX. 

	

	a. Main files:

	TestInseparable.java :  This is the main test class to run the test

	TestSVM_Continuous2.java :  This is a test class to evaluate the numerical performance of DCA, MDA when they are used as a subroutine in the projected gradient descent of SOVREX.

	RAP_Continuous.java : This a class for simple resource allocation with separable convex quadratic objectives. A bisection method is implemented to solve the RAP subproblems in the two algorithms to speed up the performance.

	RAPnc_Continuous.java :  A class for resource allocation problem with lower and upper nested constraints (RAPNC). DCA and MDA for RAP-NC is implemented to solve an RAP-NC instance.

	d. Files for reading and preprocessing the data:
	ReadFile.java
	ReadFileLineByLineUsingScannerAbalone.java
	ReadFileBank32nh.java
	ReadFileBostonHousing.java
	ReadFileCaliforniaHousing.java 
	ReadFileCensus.java 
	ReadFileComputer.java 
	ReadFileMachineCPU.java
	ReadFilePyrimidines.java		
		
	c. ResultType files: 
	ResultTypeRAP_Continuous.java : This is a class for storing the solution to an instance continuous RAP
	ResultTypeMDA_Continuous.java : This is a class for storing the solution to the subproblems of RAP-NC in the MDA method
	ResultTypeRAPNC_Continuous.java : This is a class for storing the solution to an instance continuous RAP-NC

	d. Data files:
	Data folder. The data sets are available at \url{https://www.dcc.fc.up.pt/~ltorgo/Regression/DataSets.html





