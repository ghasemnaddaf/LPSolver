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

    Main.cpp  Test examples from: http://en.wikipedia.org/wiki/Simplex_algorithm
*/
#include "LPSolver.h"

int main(int argc,char*argv[]){
	std::cout << "--- Example 1 ---" << std::endl;
	Tableau *T1 = new Tableau(3,7);
	double D0[7] = {1,-2,-3,-4,0,0,0};
	double D1[7] = {0,3,2,1,1,0,10};
	double D2[7] = {0,2,5,3,0,1,15};
	T1->setRow(0,D0);
	T1->setRow(1,D1);
	T1->setRow(2,D2);
	TableauSolver *TS1 = new TableauSolver(T1);
	int status = TS1->Solve();
	std::cout << "success: " << status << std::endl;
	std::cout << "optval: " << T1->val[0][6] << std::endl;

	std::cout << "--- Example 2 ---" << std::endl;
	Tableau *T2 = new Tableau(4,8);
	double D0_2[8] = {1,0,-5,-7,-4,0,0,-25};
	double D1_2[8] = {0,1,2,3,4,0,0,0};
	double D2_2[8] = {0,0,3,2,1,1,0,10};
	double D3_2[8] = {0,0,2,5,3,0,1,15};
	T2->setRow(0,D0_2);
	T2->setRow(1,D1_2);
	T2->setRow(2,D2_2);
	T2->setRow(3,D3_2);
	TableauSolver *TS2 = new TableauSolver(T2);
	int status2 = TS2->Solve();
	std::cout << "success: " << status2 << std::endl;
	std::cout << "optval: " << T2->val[0][7] << std::endl;

	std::cout << "--- Example 3: TableauSolver Phase I Needed---" << std::endl;
	Tableau *T3 = new Tableau(3,5);
	double D0_3[5] = {1,-2,-3,-4,0};
	double D1_3[5] = {0,3,2,1,10};
	double D2_3[5] = {0,2,5,3,15};
	T3->setRow(0,D0_3);
	T3->setRow(1,D1_3);
	T3->setRow(2,D2_3);
	TableauSolver *TS3 = new TableauSolver(T3);
	int status3 = TS3->Solve();
	std::cout << "success: " << status3 << std::endl;
	std::cout << "optval: " << T3->val[0][4] << std::endl;

	std::cout << "--- Example 4: LPSolver for Example 1---" << std::endl;
	double w4[3] = {-2,-3,-4};
	double b4[2] = {10,15};
	Matrix *A4 = new Matrix(2,3);
	double A40[3] = {3,2,1};
	double A41[3] = {2,5,3};
	A4->setRow(0,A40);
	A4->setRow(1,A41);
	LPSolver *LS4 = new LPSolver(3,w4,A4,b4);
	LS4->Solve();
	std::cout << "success: " << LS4->status << std::endl;
	std::cout << "optval: "  << LS4->optval << std::endl;

	std::cout << "--- Example 5: LPSolver for Example 3---" << std::endl;
	double w5[3] = {-2,-3,-4};
	double d5[2] = {10,15};
	//dblPtr b5 = NULL;
	//Matrix *A5 = new Matrix(0,3);
	Matrix *C5 = new Matrix(2,3);
	double C50[3] = {3,2,1};
	double C51[3] = {2,5,3};
	C5->setRow(0,C50);
	C5->setRow(1,C51);
	LPSolver *LS5 = new LPSolver(3,w5,NULL,NULL,C5,d5);
	LS5->Solve();
	std::cout << "success: " << LS5->status << std::endl;
	std::cout << "optval: "  << LS5->optval << std::endl;

	return 0;
}
