#include "aditionals.h"

double  ancho   =0.9,
        alto    =0.9;
cubo *_c=new cubo(3);
float dt=0.001;
float tfin=100;
double T=0.2;

double gasdev () {
     static bool available = false;
     static double gset;
     double fac, rsq, v1, v2;
     if (!available) {
          do {
               v1 = 2 * rand() / double(RAND_MAX) - 1.0;
               v2 = 2 * rand() / double(RAND_MAX) - 1.0;
               rsq = v1 * v1 + v2 * v2;
          } while (rsq >= 1.0 || rsq == 0.0);
          fac = sqrt(-2.0 * log(rsq) / rsq);
          gset = v1 * fac;
          available = true;
          return v2*fac;
     } else {
          available = false;
          return gset;
     }
}




void calcular_aceleracion()
{
    int num_mol=_c->num_moleculas;
    for (int i = 0; i < num_mol; i++)
        _c->vec_cubo[i]->ax=_c->vec_cubo[i]->ay= 0;
    for (int i = 0; i < num_mol; i++)
        for (int j = i + 1; j < num_mol; j++)
        {
            double dx = _c->vec_cubo[i]->x - _c->vec_cubo[j]->x;
            double dy = _c->vec_cubo[i]->y - _c->vec_cubo[j]->y;
            double dr = sqrt(dx * dx + dy * dy);
            double f = 48.0 * pow(dr, -13.0) - 24 * pow(dr, -7.0);
            _c->vec_cubo[i]->ax += f * dx / dr;
            _c->vec_cubo[i]->ay += f * dy / dr;
            _c->vec_cubo[j]->ax -= f * dx / dr;
            _c->vec_cubo[j]->ay -= f * dy / dr;
        }
}


void reescalar_velocidades()
{
    int num_mol=_c->num_moleculas;
    double vSqdSum = 0;
    for(int i=0;i<num_mol;i++)
    {
        vSqdSum += pow(_c->vec_cubo[i]->vx,2.0);
        vSqdSum += pow(_c->vec_cubo[i]->vy,2.0);
    }
    double lambda = sqrt( 3 * (num_mol-1) * T / vSqdSum );

    for (int i = 0; i < num_mol; i++)
    {
        _c->vec_cubo[i]->vx*=lambda;
        _c->vec_cubo[i]->vy*=lambda;
    }
}



void inicializar_velocidades()
{
    int num_mol=_c->num_moleculas;

    for(int i=0;i<num_mol;i++)
    {
        _c->vec_cubo[i]->vx=gasdev();
        _c->vec_cubo[i]->vy=gasdev();
    }

    double vCM[2] = {0, 0};
    for(int i=0;i<num_mol;i++)
    {
        vCM[0]+=_c->vec_cubo[i]->vx;
        vCM[1]+=_c->vec_cubo[i]->vy;
    }
    for (int i = 0; i < 2; i++)
        vCM[i] /= num_mol;
    for(int i=0;i<num_mol;i++)
    {
        _c->vec_cubo[i]->vx=vCM[0];
        _c->vec_cubo[i]->vy=vCM[1];
    }
    reescalar_velocidades();
}


void colision_limite(molecula * mol)
{
    if(mol->x > ancho )
    {
        mol->vx=-mol->vx;
        mol->x=ancho;
    }
    else if(mol->x  < -ancho)
    {
        mol->vx=-mol->vx;
        mol->x=-ancho;
    }
    else if(mol->y > alto)
    {
        mol->vy=-mol->vy;
        mol->y =alto;
    }
    else if(mol->y< -alto)
    {
        mol->vy=-mol->vy;
        mol->y= -alto;
    }
}

void velocidad_verlet(cubo *contenedor)
{
    for(int i=0;i<contenedor->num_moleculas;i++)
    {

        molecula *mol=contenedor->vec_cubo[i];
        mol->x+=mol->vx*dt+0.5*mol->ax*dt*dt;
        mol->y+=mol->vy*dt+0.5*mol->ay*dt*dt;
    }
}

void dibujar_moleculas()
{
    for(int i=0;i<_c->num_moleculas;i++)
    {
        colision_limite(_c->vec_cubo[i]);
        dibujar(_c->vec_cubo[i]);
    }
}

void correr()
{
    for(float tiempo=0;tiempo<tfin;tiempo+=dt)
    {

        //Sleep(10);
        limpiar();
        velocidad_verlet(_c);
        calcular_aceleracion();
        dibujar_moleculas();
    }

}


void recorrido()
{
  /*
    glutInitWindowPosition(700, 50);
    glutInitWindowSize(600,600);
    glutCreateWindow("MD");
    for(float tiempo=0;tiempo<tfin;tiempo+=dt)
    {
        //Sleep(10);
        dibujar(_c->vec_cubo[1]);
    }
    glutMainLoop();
    */
}



int main(int argc, char *argv[])
{

    inicializar_velocidades();
    glutInit(&argc,argv);
    glutInitWindowPosition(50, 50);
    glutInitWindowSize(600,600);
    glutCreateWindow("MD");
    thread th(recorrido);
    correr();
    th.join();
    glutMainLoop();
    return 0;
}
