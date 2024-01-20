#include <cstdlib>
#include <iostream>
#include <graphics.h>

using namespace std;

class Ball1 
{
 private:
    int Wx,Wy;
    int R;
    int x,y;
    int Vx, Vy;
 public:
    Ball1(int Wx, int Wy, int R, int x, int y, int Vx, int Vy);
    ~Ball1();
    void BallStep();
    void FrameChange();
    void Go();
};

Ball1::Ball1(int Wx, int Wy, int R, int x, int y, int Vx, int Vy) 
{
    this -> Wx = Wx;
    this -> Wy = Wy;
    this->R = R;
    this -> x = x;
    this -> y = y;
    this -> Vx = Vx;
    this -> Vy = Vy;
    initwindow(Wx,Wy);
}
    
Ball1::~Ball1(){
    system("PAUSE");
    closegraph();
}
    
void Ball1::BallStep()
{
    x = x + Vx;
    y = y + Vy;
    if ((x < R) && (Vx < 0))          Vx = -Vx;
    if ((x > Wx-R) && (Vx > 0))       Vx = -Vx;
    if ((y < R) && Vy < 0)            Vy = -Vy;
    if ((y > Wy-R) && Vy > 0)         Vy = -Vy;
}

void Ball1::FrameChange()
{
    BallStep();
    clearviewport();
    circle(x,y,R);
}

void Ball1::Go()
{
    const int keyEsc =27;
    int key = 0;
    while (key != keyEsc)
    {
        FrameChange();
        if (kbhit()) key = getch();
        delay(20);
    }
}

int main(int argc, char *argv[])
{
   Ball1* b = new Ball1(800, 600, 30, 400, 300, -1, 2);
   b->Go();
   
   return 0;
}

