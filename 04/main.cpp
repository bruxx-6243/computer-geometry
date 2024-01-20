#include <cstdlib>
#include <iostream>
#include <graphics.h>
#include <math.h>


using namespace std;

class RealWin
{
    private:
        double  Xmin, Xmax, Ymin, Ymax, X0, Y0, K;
        int Wx, Wy, X0e, Y0e;
    public:
        RealWin(double Xmin, double Ymin, double Xmax, double Ymax, int Wx, int Wy);
        ~RealWin();
        void R2_to_Z2(double x, double y, int& xe, int& ye);
        void DrawPoint(double Tx, double Ty, int color);
        void DrawSegment(double x0, double y0, double x1, double y1, int color);
        void DrawAxes();
};


RealWin::RealWin(double Xmin, double Xmax,
                 double Ymin, double Ymax,
                 int Wx, int Wy) 
{
    this->Xmin = Xmin;
    this->Xmax = Xmax;
    this->Ymin = Ymin;
    this->Ymax = Ymax;
    this->Wx = Wx;
    this->Wy = Wy;
    X0 = (Xmin+Xmax)/2; Y0 = (Ymin+Ymax)/2;
    X0e = Wx/2;         Y0e = Wy/2;
    double Kx = Wx/(Xmax-Xmin);
    double Ky = Wy/(Ymax-Ymin);
    if (Kx < Ky) K = Kx; else K = Ky;
    initwindow(Wx,Wy);
    DrawAxes();
} 

RealWin::~RealWin() {
    closegraph();
}

void RealWin::DrawAxes()
{
    DrawSegment(0,Ymin, 0,Ymax, GREEN);
    DrawSegment(Xmin,0, Xmax,0, GREEN);
    int xe, ye;
    R2_to_Z2(0,0, xe,ye);
    setcolor(GREEN);
    outtextxy(xe,ye, "0");
}

void RealWin::R2_to_Z2(double x, double y, int& xe, int& ye) {
    xe = (int)( K*(x-X0) + X0e);
    ye = (int)(-K*(y-Y0) + Y0e);
}

void RealWin::DrawPoint(double Tx, double Ty, int color) {
    int xe, ye;
    R2_to_Z2(Tx, Ty, xe, ye);
    putpixel(xe, ye, color);
}

void RealWin::DrawSegment(double x0, double y0,
                          double x1, double y1, int color)
{
    int xe0, ye0, xe1, ye1;
    R2_to_Z2(x0, y0, xe0, ye0);
    R2_to_Z2(x1, y1, xe1, ye1);
    setcolor(color);
    line(xe0, ye0, xe1, ye1);
}

//=============================================================================

const double a = 2;



double X(double t)
{
    return 4*cos(t);
}

double Y(double t)
{
    return 4*sin(t);
}


void draw_Curve(double a, double b, int n)
{
    double Xmin = X(a), Xmax = X(a), Ymin  = Y(a),  Ymax = Y(a);
    double dt = (b - a)/n;
    for (int i  = 0; i <= n; i++)
    {
        double t = a + dt*i;
        double x = X(t);
        double y = Y(t);
        if(x < Xmin) Xmin = x;
        if(x > Xmax) Xmax = x;
        if(y < Ymin) Ymin = y;
        if(y > Ymax) Ymax = y;
    }
    
    
    int Wx = 800, Wy = 600;
    
    RealWin* RW = new RealWin(-10, 10, -10, 10, Wx, Wy);
    
    for(int i = 0; i <= n; i++)
    {
        double t = a + dt*i;
        double x = X(t);
        double y = Y(t);
        RW -> DrawPoint(x,y, YELLOW);
    }
}

int main()
{
    draw_Curve(0, 2*M_PI, 1000);
    system("PAUSE");
    return 0;
}
