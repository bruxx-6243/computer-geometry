#include <cstdlib>
#include <iostream>
#include <graphics.h>
#include <math.h>

using namespace std;

class Ball_3
{
    int Wx, Wy;
    int x, y, R;
    double v, phi, Vx, Vy;
    int curKey;
public:
    Ball_3(int Wx, int Wy, int x, int y, int R, double v, double phi);
    ~Ball_3();
    void BallStep();
    void PressKey(int key);
    void FrameChange();
    void GO();
};

Ball_3::Ball_3(int Wx, int Wy, int x, int y, int R, double v, double phi)
{
    this -> Wx = Wx;
    this -> Wy = Wy;
    this -> x = x;
    this -> y = y;
    this -> R = R;
    this -> v = v;
    this -> phi = phi;
    Vx = v * cos(phi);
    Vy = v * sin(phi);
    initwindow(Wx, Wy);
    circle(x, y, R);
}

Ball_3::~Ball_3()
{
    closegraph();
}

void Ball_3::BallStep()
{
    x += Vx;
    y += Vy;
    if ((x < R) && (Vx < 0))         Vx = -Vx;
    if ((x > Wx - R) && (Vx > 0))    Vx = -Vx;
    if ((y < R) && (Vy < 0))         Vy = -Vy;
    if ((y > Wy - R) && (Vy > 0))    Vy = -Vy;
}

void Ball_3::PressKey(int key)
{
    const double dv = 0.5; 
    const double dphi = 5.0 * M_PI/180;
    if (key == KEY_UP)  v += dv;
    if (key == KEY_DOWN)  v -= dv;
    if (key == KEY_LEFT)  phi += dphi;
    if (key == KEY_RIGHT) phi -= dphi;
    if ((key == KEY_UP)||(key == KEY_DOWN)||(key == KEY_LEFT)||(key == KEY_RIGHT))
    {
        Vx = v * cos(phi);
        Vy = v * sin(phi);
    }
}

void Ball_3::FrameChange()
{
    BallStep(); clearviewport();
    circle(x, y, R);
    line(x, y, x + R * cos(phi), y + R * sin(phi));
}

void Ball_3::GO()
{
    curKey = 0;
    while (curKey != 27)
    {
        if (kbhit())
        {
            curKey = getch();
            PressKey(curKey);
        }
        FrameChange();
        delay(20);
    }
}


int main(int argc, char *argv[])
{
    Ball_3 ball(800, 600, 400, 300, 10, 1, 45);
    ball.GO();
    return EXIT_SUCCESS;
}
