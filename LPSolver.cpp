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

    LPSolver.cpp  Solver of a linear program based on the Simplex Algorithm.
      Implementation of Matrix, Tableau, TableauSolver, and LPSolver classes.
*/
#include "LPSolver.h"

//=== Utility Functions ===
bool all(bool *B,int size){
	for (int i=0;i<size;i++)
		if (!B[i])
			return false;
	return true;
}

bool isZero(double x){
	return abs(x) < EPSILON ? true : false;
}

//=== Matrix Functions ===
Matrix::Matrix(int height,int width){
	this->height = height;
	this->width = width;
	this->val = (dblPtr*) malloc(this->height*sizeof(dblPtr));
	this->rowIsSet = (bool*) malloc(this->height*sizeof(bool));
	for (int i=0; i<this->height; i++){
		this->val[i] = (dblPtr) malloc(this->width*sizeof(double));
		this->rowIsSet[i] = false;
	}
}

void Matrix::setRow(int rowIdx,dblPtr row){
	for (int i=0;i<this->width;i++)
		this->val[rowIdx][i] = row[i];
	this->rowIsSet[rowIdx] = true;
}

void Matrix::rowAdd(int rowIdxTarget ,int rowIdxSource, double ratio){
	for (int j=0;j<this->width;j++)
		this->val[rowIdxTarget][j] = this->val[rowIdxTarget][j] + this->val[rowIdxSource][j] * ratio;
}

void Matrix::rowMult(int rowIdx,double a){
	for (int j=0;j<this->width;j++)
		this->val[rowIdx][j] *= a;
}

//=== Tableau Functions ===
Tableau::Tableau(int numberOfConstraints,int numberOfVariables) : Matrix(numberOfConstraints,numberOfVariables) {
	this->numberOfConstraints = numberOfConstraints;
	this->numberOfVariables = numberOfVariables;
	this->pricedOut = false;
	this->locationOfDedicatedVariables = (int*) malloc(this->numberOfConstraints*sizeof(int));
	for (int i=0; i<this->numberOfConstraints; i++){
		this->locationOfDedicatedVariables[i] = -1;
	}
}

void Tableau::setRow(int rowIdx,dblPtr row) {
	Matrix::setRow(rowIdx,row);
	this->pricedOut = false;
}

int Tableau::countNonZeroRowsForVariable(int variableIdx, int *loc){
	//counts and returns the number of constraint rows i which have non-zero val[i][variableIdx]
	//loc is the index of first such row , or -1.
	int cnt=0;
	*loc=-1;
	for (int i=1;i<this->numberOfConstraints;i++)
		if (!isZero(this->val[i][variableIdx])){
			cnt++;
			if ((*loc)<0)
				*loc = i;
		}
	return cnt;
}

void Tableau::findLocationOfDedicatedVariables(){
	assert(all(this->rowIsSet,this->numberOfConstraints));
	//find dedicated variables for each constraint row.
	for (int i=0;i<this->numberOfVariables;i++){
		int loc;
		int cnt = this->countNonZeroRowsForVariable(i,&loc);
		if (cnt==1) //&& this->locationOfDedicatedVariables[*loc]==-1)
			this->locationOfDedicatedVariables[loc] = i;
		if (i==0){
			assert(cnt==0);
			this->locationOfDedicatedVariables[0] = 0;
		}
	}
}

std::vector<int> Tableau::findConstraintsWithoutDedicatedVariable(){
	findLocationOfDedicatedVariables();
	std::vector<int> constraintsWithoutDedicatedVariable;
	for(int i=0;i<this->numberOfConstraints;i++)
		if (this->locationOfDedicatedVariables[i] < 0)  //need to add a variable for this constraint
			constraintsWithoutDedicatedVariable.push_back(i);
	return constraintsWithoutDedicatedVariable;
}

int Tableau::status(){
	this->findLocationOfDedicatedVariables();
	std::vector<int> constraintsWithoutDedicatedVariable = this->findConstraintsWithoutDedicatedVariable();
	if (constraintsWithoutDedicatedVariable.size() > 0)
		return TABLEAU_STATUS_NOT_CANONICAL;
	else if (this->pricedOut)
		return TABLEAU_STATUS_CANONICAL;
	else
		return TABLEAU_STATUS_CANONICAL_BUT_NOT_PRICED_OUT;
}

void Tableau::priceOut(){
	assert(this->status()==TABLEAU_STATUS_CANONICAL_BUT_NOT_PRICED_OUT);
	for (int i=1;i<this->numberOfConstraints;i++){
		int colIdx = this->locationOfDedicatedVariables[i];
		rowMult(i,1/this->val[i][colIdx]);
		//if ( !isZero(this->val[0][colIdx]))
		rowAdd(0,i,-this->val[0][colIdx]);
	}
	this->pricedOut = true;
}

int Tableau::Solve(){
	assert(this->status()==TABLEAU_STATUS_CANONICAL);
	int round =0;
	while (true){
		round++;
		std::cout << (*this);
		//A. Choose the pivot column
		int pivotIdx = -1;
		for (int i=1;i< this->numberOfVariables;i++){
			if (this->val[0][i] < 0){ //minimize --> set the first negative relative cost as the pivot
				pivotIdx = i;
				break;
			}
		}
		if (pivotIdx <0)
			break; //optimum found
		std::cout << "Round#" << round << "Pivot#" << pivotIdx << "(" << this->val[0][pivotIdx] << ")" << std::endl;

		//B. Choose the pivot row
		double minRatio=INF;
		int minRatioIdx=-1;
		for (int i=1;i< this->numberOfConstraints;i++){
			double ratio;
			if (this->val[i][pivotIdx]>0){ //find the smallest b/d ratio for positive d
				ratio = this->val[i][this->numberOfConstraints-1] / this->val[i][pivotIdx];
				if (ratio<minRatio){
					minRatio = ratio;
					minRatioIdx = i;
				}
			}
		}
		if (minRatioIdx <0) {
			return LPSOLVER_STATUS_UNBOUNDED;
		}
		std::cout << "Row#" << minRatioIdx << "(" << this->val[minRatioIdx][pivotIdx] << ")" << std::endl;

		//C. Reduce the pivot row from other rows, so that there is only one 1 in the pivotIdx column (i.e. the selected row).
		for (int i=0;i< this->numberOfConstraints;i++){
			if (   isZero(this->val[i][pivotIdx])   )
				continue;
			double val_i_piv = this->val[i][pivotIdx];
			double ratio= val_i_piv / this->val[minRatioIdx][pivotIdx];
			if (i==minRatioIdx)
				this->rowMult(i,1/val_i_piv);
			else
				this->rowAdd(i,minRatioIdx,-ratio);
		}
	}
	return LPSOLVER_STATUS_SOLVED;
}


std::ostream& operator<<(std::ostream &outStream, const Tableau &T){
	for (int i=0;i<T.numberOfConstraints;i++){
		outStream << '[';
		for (int j=0;j<T.numberOfVariables;j++)
			outStream << std::setprecision(5) << T.val[i][j] << '\t';
		outStream << ']'<< std::endl ;
	}
	return outStream;
}


//=== TableauSolver Functions ===
TableauSolver::TableauSolver(Tableau *T){
	this->T = T;
	this->MakePhaseITableau();
}

int TableauSolver::Solve(){
	if (needsPhaseI){		
		int status = phaseITableau->Solve();
		int L = phaseITableau->numberOfVariables-1;
		assert(status==LPSOLVER_STATUS_SOLVED && isZero(phaseITableau->val[0][L]));
		//rebuild T based on the solution of Phase I.
		dblPtr row = (dblPtr) malloc(sizeof(double)*T->numberOfVariables);
		for(int i=1;i< T->numberOfConstraints+1;i++){
			for (int j=0;j<T->numberOfVariables-1;j++)
				row[j] = phaseITableau->val[i][j+1];
			row[T->numberOfVariables-1] = phaseITableau->val[i][L];
			T->setRow(i-1,row);
		}
	}
	T->priceOut();
	return T->Solve();
}



void TableauSolver::MakePhaseITableau(){	
	//assert(T->status()==TABLEAU_STATUS_NOT_CANONICAL);
	std::vector<int> constraintsWithoutDedicatedVariable = T->findConstraintsWithoutDedicatedVariable();
	if (constraintsWithoutDedicatedVariable.size() == 0){
		needsPhaseI = false;
		return;
	}

	phaseITableau = new Tableau(T->numberOfConstraints+1, 
											T->numberOfVariables+1+constraintsWithoutDedicatedVariable.size());

	dblPtr row = (dblPtr) malloc(sizeof(double)*phaseITableau->numberOfVariables);
	
	//set the new objective: maximize sum(newly_added_vars) 
	row[0] = 1.0;
	for (int i=1;i<phaseITableau->numberOfVariables-1;i++)
		if (i < T->numberOfVariables)
			row[i] = 0.0;
		else
			row[i] = 1.0;
	row[phaseITableau->numberOfVariables-1] = 0.0;
	phaseITableau->setRow(0,row);

	//set the constraint rows:
	for (int rowIdx=T->numberOfConstraints-1;rowIdx>=0;rowIdx--){
		row[0]=0.0;
		row[phaseITableau->numberOfVariables-1] = T->val[rowIdx][T->numberOfVariables-1];		
		for (int j=1;j<phaseITableau->numberOfVariables-1;j++){
			if (j< T->numberOfVariables)
				row[j] = T->val[rowIdx][j-1];
			else
				row[j]=0;
		}
		if (constraintsWithoutDedicatedVariable.size()>0 && constraintsWithoutDedicatedVariable.back()==rowIdx){
			row[T->numberOfVariables+constraintsWithoutDedicatedVariable.size()-1] =1.0;
			constraintsWithoutDedicatedVariable.pop_back();
		}
		phaseITableau->setRow(rowIdx+1,row);
	}
	phaseITableau->priceOut();
	assert(phaseITableau->status()==TABLEAU_STATUS_CANONICAL);
	needsPhaseI = true;
	return;
}

//=== LPSolver funtions ===
LPSolver::LPSolver(int numberOfVariables, dblPtr w, Matrix *A, dblPtr b, Matrix *C, dblPtr d){
	if (!A){
		assert(C); //both A and C can not be NULL at the same time --> Unbounded.
		A = new Matrix(0,numberOfVariables);
		b = NULL;
	}
	this->numberOfConstraints = (C? (A->height + C->height) : A->height);
	this->numberOfVariables = numberOfVariables; 
	assert(A->width == numberOfVariables);
	if (C) assert(A->width == C->width);
	this->w = w;
	this->A = A;
	this->b = b;
	this->C = C;
	this->d = d;
}

int LPSolver::Solve(){
	Tableau *T = new Tableau(numberOfConstraints+1,numberOfVariables+A->height+2); //for each inequality cst, one new slack variable should be added
	dblPtr row = (dblPtr) malloc(sizeof(double)*T->numberOfVariables);
	//first row:
	row[0] = 1.0;
	row[T->numberOfVariables-1] = 0.0;
	for (int i=1;i<T->numberOfVariables-1;i++)
		if (i<this->numberOfVariables+1)
			row[i] = w[i-1];
		else
			row[i] = 0.0;
	T->setRow(0,row);
	
	//inequality constraint rows:
	for (int rowIdx=0;rowIdx<A->height;rowIdx++){
		row[0] = 0.0;
		row[T->numberOfVariables-1] = b[rowIdx];
		for(int i=1;i<T->numberOfVariables-1;i++){
			if (i<this->numberOfVariables+1)
				row[i]=A->val[rowIdx][i-1];
			else
				row[i] = 0.0;
		}
		row[rowIdx+this->numberOfVariables+1] = 1.0; //slack variable for inequality cst rowIdx
		T->setRow(rowIdx+1,row);
	}

	//equality constraint rows:
	if (C)
		for (int rowIdx=0;rowIdx<C->height;rowIdx++){
			row[0] = 0.0;
			row[T->numberOfVariables-1] = d[rowIdx];
			for(int i=1;i<T->numberOfVariables-1;i++){
				if (i<this->numberOfVariables+1)
					row[i]=C->val[rowIdx][i-1];
				else
					row[i] = 0.0;
			}
			T->setRow(rowIdx+A->height+1,row);
		}
	std::cout << "LPSolver has constructed the following Tableau:" << std::endl << (*T) << std::endl;
	TableauSolver *TS = new TableauSolver(T);
	this->status = TS->Solve();
	this->optval = T->val[0][T->numberOfVariables-1];
	//this->OptimumVariables = ;
	return this->status;
}

