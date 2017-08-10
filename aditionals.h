#ifndef ADITIONALS_H
#define ADITIONALS_H



#include <GL/glut.h>
#include <vector>
#include <thread>
#include <iostream>
#include <math.h>


using namespace std;

struct molecula
{
    double x,y,z;
    float radio=0.01;
    double vx,vy;
    double ax,ay;

    molecula ()
    {
        x=0.05;
        y=0.05;
    }

    molecula(double _x, double _y)
    {
        x=_x;
        y=_y;
    }
};

struct cubo
{
  vector<molecula*> vec_cubo;
  int num_moleculas=0;

  cubo(int n)
  {
     num_moleculas=n;
      vec_cubo.resize(num_moleculas);
       for(int i=0;i<num_moleculas;i++)
       {
           vec_cubo[i]= new molecula(rand()%10 /10.0,rand()%10 /10.0);
       }
  }
};


void dibujar(molecula *mol)
{
    double x=mol->x;
    double y=mol->y;
    float radio=mol->radio;

    glPushMatrix();
        glTranslatef(x,y,0);
        glutSolidSphere(radio,20,10);
    glPopMatrix();

    glColor3f(1,1,1);
    glutWireCube(1.8);
    glFlush();

}

void limpiar()
{
    glClear(GL_COLOR_BUFFER_BIT);
}


#endif // ADITIONALS_H
