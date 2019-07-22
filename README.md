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

Here are the short description on each file

0. Main source files:
	Function.java : This is an abstract class for the function oracles
	
	RAP.java : This a class for discrete simple resource allocation (DRAP) with separable convex objectives. 
	Hochbaum's algorithm (1994) is implemented to solve DRAP with general convex objectives. 
	In addition, for linear objectives, a more efficient O(n) time algorithm is also implemented.
	
	RAPNC.java : This is a class for RAPNC with separable convex objectives. It provides three methods to solve the problem (DCA, MDA, SFA).
	
	RAP_Continuous.java : This is a class for simple resource allocation with separable convex quadratic objectives. A bisection method is implemented to solve the RAP subproblems in the two algorithms to speed up the performance.
	
	RAPNC_Continuous.java : This is class for resource allocation problem with lower and upper nested constraints (RAPNC). DCA and MDA for RAP-NC is implemented to solve an RAP-NC instance.


	Unit test files: 


	ResultType files:
		Discrete problem:
			ResultTypeRAP.java : This is a class for storing the solution to an instance of DRAP
			ResultTypeMDA.java : This is a class for storing the solution to the subproblems of DRAP-NC in the MDA method
			ResultTypeRAPNC.java : This is a class for storing the solution to an instance of DRAP-NC
		Continuous problem: 
			ResultTypeRAP_Continuous.java : This is a class for storing the solution to an instance of continuous RAP 
			ResultTypeMDA_Continuous.java : This is a class for storing the solution to the subproblems of RAP-NC in the MDA method
			ResultTypeRAPNC_Continuous.java : This is a class for storing the solution to an instance of continuous RAP-NC


1. Numerical experiment on DCA and Gurobi on DRAP-NC with linear and quadratic objectives
	a. Main files:
		TestGurobiDCA_Lin.java : This is a test class to evaluate the numerical performance of DCA, and Gurobi for DRAP-NC with linear objectives
		TestGurobiDCA_Qua.java : This is a test class to evaluate the numerical performance of DCA, and Gurobi for DRAP-NC with quadratic objectives
		TestSparseGurobiDCA_Lin.java : This is a test class to evaluate the numerical performance of DCA, and Gurobi for DRAP-NC with linear objectives under a sparse formulation. The running time of Gurobi is sped up by more than 30 times.
		TestSparseGurobiDCA_Qua.java : This is a test class to evaluate the numerical performance of DCA, and Gurobi for DRAP-NC with quadratic objectives under a sparse formulation. The running time of Gurobi is sped up by more than 30 times.


2. Numerical experiment on DCA and MDA on DRAP-NC with three benchmark convex objectives function
	To be included

3. Numerical experiment on DCA and MDA on SOVREX 
	Our last numerical experiment is on the Support Vector Ordinal Regression (SOVREX) model. In this algorithm, we apply DCA and MDA to solve RAP-NC (the continuous problem) subproblems in a projected gradient descent method for SOVREX. 

	a. Main files:

	TestInseparable.java :  This is the main test class to run the test

	TestSVM_Continuous2.java :  This is a test class to evaluate the numerical performance of DCA, MDA when they are used as a subroutine in the projected gradient descent of SOVREX.

	b. Files for reading and preprocessing the data:
	ReadFile.java
	ReadFileLineByLineUsingScannerAbalone.java
	ReadFileBank32nh.java
	ReadFileBostonHousing.java
	ReadFileCaliforniaHousing.java 
	ReadFileCensus.java 
	ReadFileComputer.java 
	ReadFileMachineCPU.java
	ReadFilePyrimidines.java			

	c. Data files:
	Data folder. The data sets are available at \url{https://www.dcc.fc.up.pt/~ltorgo/Regression/DataSets.html





