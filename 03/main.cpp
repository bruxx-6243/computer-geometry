#include <cstdlib>
#include <iostream>
#include <graphics.h>


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
        void DrawTrue(int color);
        void DrawFalse(int color);
        void DrawPolygon(int n, double x[], double y[], int color);
        void DrawSegment(double x0, double y0, double x1, double y1, int color);
        void DrawCircle(double x0, double y0, double R, int color);
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


void RealWin::DrawTrue(int color)
{
    setcolor(color);
    outtextxy(50,50, "Convex");
}

void RealWin::DrawFalse(int color)
{
    setcolor(color);
    outtextxy(50,50, "No convex");
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

void RealWin::DrawCircle(double x0, double y0, double R, int color)
{
    int xe, ye, re;;
    R2_to_Z2(x0, y0, xe, ye);
    re = (int)(K*R);
    setcolor(color);
    setfillstyle(1, RED);
    circle(xe, ye, re);
    floodfill(xe, ye, color);
}

void RealWin::DrawPolygon(int n, double x[], double y[], int color) 
{
    int xe, ye;
    for (int i = 0; i < n; i++) 
    {
        DrawSegment(x[i], y[i], x[(i+1)%n], y[(i+1)%n], color);
    }
}


class Polygon_test_2 
{
    private :
    int n;
    double* x; double* y;
public:
    Polygon_test_2 (int n, double x[], double y[]);
   ~Polygon_test_2();
   double Sor(double a1,  double b1, double a2, double b2);
   bool convex ();
   bool PointBelongPolygon();
};

Polygon_test_2::Polygon_test_2(int n,
                         double x[], double y[]) 
                         {
    this -> n = n;
    this -> x = x;
    this -> y = y;
}

Polygon_test_2::~Polygon_test_2() 
{
    delete [] x;
    delete [] y;
}

double Polygon_test_2::Sor(double a1,  double b1, double a2, double b2) 
{
    return (a1 * b2 - b1 * a2);
}


bool Polygon_test_2::convex() 
{
    int pos = 0, neg = 0, zero = 0;
    for(int i = 1; i <= n; i++) 
    {
         double s = Sor(x[i%n] - x[i-1], y[i%n] - y[i-1], x[(i+1)%n] - x[i%n], y[(i+1)%n]-y[i%n]);
         if (s > 0) pos ++;
         else if (s < 0) neg ++;
         else zero++;
         if (pos > 0  &&  neg > 0) break;
    }
    return (pos == 0 || neg == 0);
}



int main()
{
    RealWin* RW = new RealWin(-2, 9, -5, 5, 800, 600);
    const int n = 5;
    double x[n] = {2, 3, 6.5, 7.2, 3.2}; 
    double y[n] = {0.85, -2, -2.3, 1.5, 2};
    RW->DrawPolygon(n, x, y, YELLOW);
    Polygon_test_2* T2 = new Polygon_test_2(n, x, y);
    if(T2->convex())
    RW->DrawTrue(WHITE);    
    else RW->DrawFalse(WHITE);
    
    system ("PAUSE");
    delete[]x;
    delete[]y;
    return 0;
}
