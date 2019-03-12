/*
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    For the GNU General Public License, see <http://www.gnu.org/licenses/>.
*/
/*
	LPSolver.h:  Contains interface for 4 classes:
	Matrix class is the container for a 2D dynamic array.
	Tableau class is the container for an LP tableau and implements the SIMPLEX solver for canoncal Tableau. It inherits from Matrix.
	TableauSolver class is composed of a Tableau, and can make it canonical (if needed).
	LPSolver class gets an LP and makes the corresponding Tableau by adding necessary slack variables.
*/

#ifndef __LPSOLVER_H__
#define __LPSOLVER_H__
#include <vector>
#include <cfloat>
#include <limits>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <assert.h>

#define MAXVARS 10000	//maximum number of variables
#define MAXCSTS 128		//maximum number of constraints
#define INF DBL_MAX 
#define EPSILON 1.0e-6
#define LPSOLVER_STATUS_SOLVED 0
#define LPSOLVER_STATUS_UNBOUNDED -1
#define LPSOLVER_STATUS_INFEASIBLE -2

#define TABLEAU_STATUS_CANONICAL 0
#define TABLEAU_STATUS_CANONICAL_BUT_NOT_PRICED_OUT 1
#define TABLEAU_STATUS_NOT_CANONICAL -1

typedef double *dblPtr;

class Matrix{
public:
	bool *rowIsSet;
	dblPtr *val;
	int height,width;
	Matrix(int height, int width);
	void setRow(int rowNumber,dblPtr row);	    //Sets a row of Matrix.	
	void rowAdd(int rowIdxTarget,int rowIdxSource, double ratio); //rowIdxTarget += rowIdxSource * ratio;
	void rowMult(int rowIdx,double a); //row *= a;
};

class TableauSolver;
class Tableau: public Matrix{
	friend class TableauSolver;	
private:
	int *locationOfDedicatedVariables;
	/* countNonZeroRowsForVariable
	counts and returns the number of constraint rows i which 
	have non-zero val[i][variableIdx], loc is the index of first such row , or -1.
	*/
	int countNonZeroRowsForVariable(int variableIdx, int *loc); 
	std::vector<int> findConstraintsWithoutDedicatedVariable(); //parse locationOfDedicatedVariables and collect those <0.
	void findLocationOfDedicatedVariables();	//Finds dedicated variables for each constraint row and sets locationOfDedicatedVariables.
public:
	int numberOfConstraints;
	int numberOfVariables;
	bool pricedOut;
	Tableau(int numberOfConstraints, int numberOfVariables);
	~Tableau();
	void setRow(int rowNumber,dblPtr row);	    //Sets a row of Tableau.	
	int status(); //returns 0 if Tableau is canonical and priced out, 1 if Tableau is canonical but not priced out, -1 if it is not canonical, and -2 if some rows are not yet set.
	void priceOut(); //Price out a canonical Tableau. Should be called only when status() > 0.
	/*
	solves a canonical Tableau which is in the form:
	[1 0 c v]
	[0 I A b]
	*/
	int Solve();
	friend std::ostream& operator<<(std::ostream &out, const Tableau &T); 
};	

class TableauSolver{
private:
	Tableau *T;
	Tableau *phaseITableau; //Tableu for solving phase I, if T is not canonical.
	bool needsPhaseI;
	void MakePhaseITableau(); //makes a phase I canonical Tableau from T, if T is not already canonical.

public:
	TableauSolver(Tableau *T);
	~TableauSolver();
	int Solve(); //Solves a Tableau
};

/*
Solve LP of the form:
Maximize w . x
Subject To: 
	A x <= b,		size(A) = [numberOfConstraints, numberOfVariables]
and C x == d,		size(C) = size(A)
and x_i >= 0,		i=0, 1, ..., numberOfVariables-1
*/
class LPSolver{
public:
	dblPtr w;
	Matrix *A, *C;
	dblPtr b,d;
	int numberOfConstraints, numberOfVariables;
	LPSolver(int numberOfVariables, dblPtr w, Matrix *A, dblPtr b, Matrix *C=NULL, dblPtr d=NULL);
	~LPSolver();
	int Solve();
	dblPtr OptimumVariables;
	double optval;
	int status;
};

#endif //__LPSOLVER_H__
