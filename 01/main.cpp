#include <cstdlib>
#include <iostream>
#include <graphics.h>
using namespace std;


class RealWin
{
   public:
        RealWin(double Xmin, double Xmax, double Ymin, double Ymax, int Wx, int Wy);
        ~RealWin();
        void R2_to_Z2(double x, double y, int &xe, int &ye);
        void DrawPoint(double x, double y, int color);
        void DrawSegment(double x0, double y0, double x1, double y1, int color);
        void DrawCircle(double X0, double X0, double R, int color);
        void DrawAxes(double x, double y);
        
   private:
        double Xmin, Xmax, Ymin, Ymax, X0, Y0, k;
        int Wx, Wy, X0e, Y0e;    
};

RealWin::RealWin(double Xmin, double Xmax, double Ymin, double Ymax, int Wx, int Wy)
{
    this -> Xmin = Xmin;
    this -> Xmax = Xmax;
    this -> Ymin = Ymin;
    this -> Ymax = Ymax;
    this -> Wx = Wx;
    this -> Wy = Wy;
    
    X0 = (Xmin + Xmax)/2; Y0 = (Ymin + Ymax)/2;
    X0e = Wx/2; Y0e = Wy/2;
    double kx, ky;
    if(kx > ky) k = kx; else k = ky;
    initindow(Wx, Wy);    
}

RealWin::~RealWin()
{
   closegraph(); 
}

void RealWin::R2_to_Z2(double x, double y, int &xe, int &ye)
{
    xe = int (k*(x - X0) + X0e);
    ye = int (-k*(y - Y0) + Y0e);
}

void RealWin::DrawPoint(double x, double y, int color)
{
   int xe, ye;
   R2_to_Z2(x, y, xe, ye);
   putpixel(xe, ye, color);
}

void RealWin::DrawSegment(double x0, double y0, double x1, double y1, int color)
{
  int xe0, ye0, xe1, ye1;
  R2_to_Z2(x0, y0, xe0, ye0, xe1, ye1);
  line(xe0, ye0, xe1, ye1);  
}

void RealWin::DrawCircle(double X0, double X0, double R, int color)
{
  int xe, ye, re;
  R2_to_z2(x, y, xe, ye);
  re = int (k*r);
  setcolor(color);
  circle(xe, ye, re);  
}

void RealWin::DrawAxes(double x, double y)
{
    
}

int main(int argc, char *argv[])
{
   

    system("PAUSE");
    closegraph();
       
   return EXIT_SUCCESS;
}


