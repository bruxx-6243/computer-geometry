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
	void DrawPoint(double po[4], int color);
	void DrawAxes();
	void DrawULine(double v0, int color);
	void DrawVLine(double u0, int color);
	void DrawCoverSkeleton(int color);


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
	k = kx < ky ? kx : ky;
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

void R3_to_Z2::DrawPoint(double po[4], int color)
{
    double pe[4];
    Ob_to_Ekr(po, pe);
    putpixel(pe[0], pe[1], color);
}

void R3_to_Z2::DrawSegment(double Pobj1[4], double Pobj2[4], int color) {
	double Pd1[4], Pd2[4];
	Ob_to_Ekr(Pobj1, Pd1); Ob_to_Ekr(Pobj2, Pd2);
	setcolor(color);
	line(Pd1[0], Pd1[1], Pd2[0], Pd2[1]);
}

void R3_to_Z2::DrawAxes() {
	double O[4] = {0, 0, 0, 1};
	double x[4] = {8, 0, 0, 1};
	DrawSegment(O, x, BLUE);
	double y[4] = {0, 8, 0, 1};
	DrawSegment(O, y, RED);
	double z[4] = {0, 0, 8, 1};
	DrawSegment(O, z, GREEN);
}

void initMatrix(double **&MV, int **&ME, int **&MF) {
	const int v = 8, e = 12, f = 12; // V, E, F
	double mv[v][4] = {};
	MV = new double*[v];
	for (int i = 0; i < v; ++i) {
        MV[i] = new double [4];
		for (int j = 0; j < 4; ++j) {
			MV[i][j] = mv[i][j];
		}
	}
	int me[e][4] = {};
	ME = new int*[e];
	for (int i = 0; i < e; ++i) {
		ME[i] = new int[4];
		for (int j = 0; j < 4; ++j) {
			ME[i][j] = me[i][j];
		}
	}

	int mf[f][3] = {};
	MF = new int*[f];
	for (int i = 0; i < f; ++i) {
		MF[i] = new int[3];
		for (int j = 0; j < 3; ++j) {
			MF[i][j] = mf[i][j];
		}
	}
}





const double v0 = 5, u0 = 20;

const double Ua = 0, Ub = M_PI * 4, Va = 0, Vb = v0;

const double a = 1;

const int Nu = 10, Nv = 7, N = 1000;


/*
double sh (double x)
{ return (exp(x) - exp(-x))/2; }

double ch (double x)
{ return (exp(x) + exp(-x))/2; }

*/


double X (double u, double v)
{ return v*cos(u); }


double Y (double u, double v)
{ return v*sin(u); }

double Z (double u, double v)
{ return a*u; }



void R3_to_Z2::DrawULine(double v0, int color)
{
    double h = (Ub-Ua)/N, u, x, y, z;
    for (int k = 0; k <= N; k++)
    {
        u = Ua + h*k;
        x = X(u, v0); y = Y(u,v0); z = Z(u, v0);
        double p[4]={x,y,z,1};
        DrawPoint(p,color);
    }
}

void R3_to_Z2::DrawVLine(double u0, int color)
{
    double h  = (Vb-Va)/N, v, x, y, z;
    for (int k = 0; k <= N; k++)
    {
        v = Va + h*k;
        x = X(u0,v); y = Y(u0,v); z = Z(u0,v);
        double p[4]={x,y,z,1};
        DrawPoint(p,color);
    }
}

void R3_to_Z2::DrawCoverSkeleton(int color)
{
    double hv = (Vb - Va)/Nu, v;
    for (int i = 0; i <= Nu; i++)
    {
        v = Va + hv*i;
        DrawULine(v, color);
    }
    double hu = (Ub - Ua)/Nv, u;
    for(int i = 0; i <= Nv; i++)
    {
        u = Ua + hu*i;
        DrawVLine(u, color);
    }

}








int main(int argc, char const *argv[])
{

	initwindow(800, 600);
	char key, Esckey, Lkey, Rkey, Upkey, Dkey;
    Esckey = 27;
    Lkey = 75;
    Rkey = 77;
    Upkey = 72;
    Dkey = 80;
    int phiz = 20, phiy = 30, D = 10;
	R3_to_Z2* RW = new R3_to_Z2(phiz, phiy, -20, D, -D, D, 800, 600);



    //RW->DrawCoverSkeleton(RGB(253, 73, 140));
    RW->DrawCoverSkeleton(WHITE);

    system("pause");
	closegraph();
}
