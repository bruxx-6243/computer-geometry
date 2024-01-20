#include <cstdlib>
#include <iostream>
#include <math.h>
#include <graphics.h>
using namespace std;



class R3_to_Z2 {
private:
	int Wx, Wy;
	double **MV; int **ME; int **MF;
public:
    int phiZ, phiY;
	double Com[4][4], Cmk[4][4], Ckd[4][4], Cmd[4][4], Cod[4][4];
	double Rz[4][4];
	double Ry[4][4];
	bool *Evis, *Fvis;
	R3_to_Z2(int phiz, int phiy, double Xmin, double Xmax, double Ymin, double Ymax, int wx, int wy);
	void SetAngles(int phiz, int phiy);
	void Ob_to_Ekr(double Pobj[4], double Pdisp[4]);
	void DrawSegment(double Pobj1[4], double Pobj2[4], int color = 15);
	void visible();
    void FacetNormal(int f, double N[4]);
	void DrawPolyhedron();
	void DrawAxes();
};


void VecProduct(double a[4], double b[4], double c[4]) {
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

void MultMatr(double A[4][4], double B[4][4], double C[4][4])
{
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			C[i][j] = 0;
			for (int k = 0; k < 4; ++k)
			{
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}

void MultVecMatr(double v[4], double A[4][4], double u[4]) {
	for (int i = 0; i < 4; ++i)
	{
		u[i] = 0;
		for (int k = 0; k < 4; ++k)
		{
			u[i] += v[k] * A[k][i];
		}
	}
}





R3_to_Z2::R3_to_Z2(int phiz, int phiy, double Xmin, double Xmax, double Ymin, double Ymax, int wx, int wy)
{
    this -> phiZ = phiz;
    this -> phiY = phiy;
    this -> Wx = wx;
    this -> Wy = wy;
    

	Rz[2][2] = Rz[3][3] = 1;
	Ry[1][1] = Ry[3][3] = 1;

	double kx = Wx / (Xmax-Xmin), ky = Wy / (Ymax-Ymin), k;
	if( kx < ky) k = kx; else k = ky;
	Cmk[1][0] = Cmk[2][1] = Cmk[3][3] = 1;
	Ckd[0][0] = k; Ckd[1][1] = -k; Ckd[3][3] = 1; Ckd[3][0] = Wx / 2; Ckd[3][1] = Wy / 2;
	MultMatr(Cmk, Ckd, Cmd);
	SetAngles(phiz, phiy);
	DrawAxes();
}

void R3_to_Z2::SetAngles(int phiz, int phiy) {
	phiZ = phiz; phiY = phiy;
	double phiZr = phiZ * M_PI / 180, phiYr = phiY * M_PI / 180;

	Rz[0][0] = Rz[1][1] = cos(phiZr);
	Rz[0][1] = sin(phiZr); Rz[1][0] = -Rz[0][1];

	Ry[0][0] = Ry[2][2] = cos(phiYr);
	Ry[2][0] = sin(phiYr); Ry[0][2] = -Ry[2][0];
	MultMatr(Rz, Ry, Com);
	MultMatr(Com, Cmd, Cod);
}

void R3_to_Z2::Ob_to_Ekr(double Pobj[4], double Pdisp[4]) {
	MultVecMatr(Pobj, Cod, Pdisp);
}


void R3_to_Z2::DrawSegment(double Pobj1[4], double Pobj2[4], int color) {
	double Pd1[4], Pd2[4];
	Ob_to_Ekr(Pobj1, Pd1); Ob_to_Ekr(Pobj2, Pd2);
	setcolor(color);
	line(Pd1[0], Pd1[1], Pd2[0], Pd2[1]);
}

void R3_to_Z2::DrawAxes() {
	double O[4] = {0, 0, 0, 1};
	double x[4] = {2, 0, 0, 1};
	DrawSegment(O, x, RGB(50, 209, 249));
	double y[4] = {0, 2, 0, 1};
	DrawSegment(O, y, RGB(253, 73, 140));
	double z[4] = {0, 0, 2, 1};
	DrawSegment(O, z, RGB(255,255,0));
}


void initMatrix(double **&MV, int **&ME, int **&MF) {
	const int v = 6, e = 10, f = 6; // V, E, F
	
	/*double mv[v][4] = {{2,1,0,1}, {3,2,0,1}, {3,4,0,1}, {2,5,0,1}, 
                        {1,3,0,1}, {2,3,5,1}};
      */              
          
    double mv[v][4] = {{0,0,0,1}, {0,0,1,1}, {1,0,0,1}, {1.5,0.5,0,1},
                       {1,1,0,1}, {0,1,0,1}};
	MV = new double*[v];
	for (int i = 0; i < v; ++i) {
        MV[i] = new double [4];
		for (int j = 0; j < 4; ++j) {
			MV[i][j] = mv[i][j];
		}
	}
	
	int me[e][4] = {{2,1,4,0}, {3,1,0,1}, {4,1,1,2}, {5,1,2,3}, 
                    {0,1,3,4}, {2,0,5,4}, {0,5,5,3}, {5,4,5,2}, 
                    {4,3,5,1}, {3,2,5,0}};
                    
	ME = new int*[e];
	for (int i = 0; i < e; ++i) {
		ME[i] = new int[4];
		for (int j = 0; j < 4; ++j) {
			ME[i][j] = me[i][j];
		}
	}

	int mf[f][3] = {{3,1,2}, {4,1,3}, {5,1,4}, {0,1,5}, {2,1,0}, 
                     {2,0,5}};
	MF = new int*[f];
	for (int i = 0; i < f; ++i) {
		MF[i] = new int[3];
		for (int j = 0; j < 3; ++j) {
			MF[i][j] = mf[i][j];
		}
	}
}

const int V = 6, E = 10, F = 6;



void R3_to_Z2::visible() {
	double N[4], Nm[4];
	for (int i = 0; i < F; ++i)
	{
		FacetNormal(i, N);
		MultVecMatr(N, Com, Nm);
		if(Nm[0] > 0) Fvis[i] = 1;
		else Fvis[i] = 0;
	}
	for (int i = 0; i < E; ++i)
	{
		int Fl = ME[i][2], Fr = ME[i][3];
		Evis[i] = Fvis[Fl] || Fvis[Fr];
	}
}

void R3_to_Z2::FacetNormal(int f, double N[4]) {
	int v0 = int(MF[f][0]), v1 = int(MF[f][1]), v2 = int(MF[f][2]);
	double e1[4] = { MV[v1][0] - MV[v0][0],
				     MV[v1][1] - MV[v0][1],
				     MV[v1][2] - MV[v0][2],
				                        1  };
    double e2[4] = { MV[v2][0] - MV[v0][0],
				     MV[v2][1] - MV[v0][1],
				     MV[v2][2] - MV[v0][2],
				                        1  };
	VecProduct(e1, e2, N);
}

void R3_to_Z2::DrawPolyhedron() {
	int v0, v1;
	initMatrix(MV, ME, MF);
	Evis = new bool[E];
	Fvis = new bool[F];
	visible();
	for (int i = 0; i < E; ++i)
	{
		v0 = ME[i][0];
		v1 = ME[i][1];
		int color = 0;
		if (Evis[i]) color = RGB(192,192,192); else color = RGB(192,192,192);
		v0 = ME[i][0];
		v1 = ME[i][1];
		DrawSegment(MV[v0], MV[v1], color);
	}
}


int main(int argc, char const *argv[])
{

	initwindow(800, 700);
	char key, Esckey, Lkey, Rkey, Upkey, Dkey;
    Esckey = 27;
    Lkey = 75;
    Rkey = 77;
    Upkey = 72;
    Dkey = 80;
    int phiz = -30, phiy = 30;
	R3_to_Z2* RW = new R3_to_Z2(phiz, phiy, -1.5, 1.5, -1.5, 1.5, 800, 700);
	RW->DrawPolyhedron();
    while(true){
        key = getch();
        if(key == Esckey){closegraph();break;}
        else if(key == Lkey) phiz -= 2;
        else if(key == Rkey) phiz += 2; 
        else if(key == Upkey)phiy += 2; 
        else if(key == Dkey) phiy -= 2; 
        RW->SetAngles(phiz, phiy); 
        cleardevice();
        RW->DrawAxes();
        RW->DrawPolyhedron();
    }
    system("pause");
	closegraph();
}
