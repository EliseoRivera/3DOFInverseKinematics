#ifndef ROBOT_H
#define ROBOT_H
#include "modelo3D.h"
#include<vector>
#include <cstdlib>
///Copyright (C) <2017>  <Eliseo Rivera> curso.tareas@gmail.com
class Robot
{
    public:
        Robot();
        ~Robot();
           modelo3D *base;
           modelo3D *b1;
           modelo3D *b2;
           modelo3D *b3;


void inicializar();
void initphysicsconstants();
void renderizar();
void configurarTH();
float t;
///
void Parametrica();
void initInverseKinematics();
void initPlanTrayec();
bool InverseKinematics();
void estado();
float x4,y4,THETA3,THETA30, xL ,yL;
float l1, l2;
Matrix O4,O4L, P,V; ///comentar quien es quien
void DibujarCurva();
std::vector<vector3d> curva;
bool IsMoving;
void Move();
void  comprobracion(float a, float b);  /// estas funciones se van a funciones utiles
Matrix Jf,Jfv, Jg1, Jg2, Jg3; // de tamaño 3x3 (en el plano ) xv=Jf*qv, con x las coordenadas absolutas
 void Calc_Jf(float q1, float q2, float q3);
 void Calc_Jfv(float q1v, float q2v, float q3v);  ///usar cada que se pueda notacion q para coordenadas generalizadas
  void Calc_Jg1(float q1, float q2, float q3);
 void Calc_Jg2(float q1, float q2, float q3);
 void Calc_Jg3(float q1, float q2, float q3);
 void Calc_Jg(float q1, float q2, float q3);
 void Calc_M(float q1, float q2, float q3);
 void Calc_hJtW(float q1, float q2, float q3,float q1v, float q2v, float q3v);
 Matrix hJtW,Q;
float modulo(float a, float b)
{

  return ( ((a/b) - floor(a/b))*b );
}

double menor (double a, double b, double c, double d){
    double m;


  //a=modulo(a,2*PI);
  //b=modulo(b,2*PI);
//  c=modulo(c,2*PI);
// cout<<"a  =" << a<< " , b= " << b<< " ,  c= " << c<< " , d = " <<d<<endl;
m=a;

if(fabs(d-b)<fabs(d-m)) m=b;
if (fabs(d-c)<=fabs(d-m)) m=c;

return m;
}
double menor2 (double a, double b, double c, double e,double d){
    double m;


  //a=modulo(a,2*PI);
  //b=modulo(b,2*PI);
//  c=modulo(c,2*PI);
// cout<<"a  =" << a<< " , b= " << b<< " ,  c= " << c<< " , d = " <<d<<endl;
m=a;

if(fabs(d-b)<fabs(d-m)) m=b;
if (fabs(d-c)<=fabs(d-m)) m=c;
if (fabs(d-e)<=fabs(d-m)) m=e;

return m;
}
Matrix A1,B1,C1,D1, E1;

//
void AplicarTHx(float theta, vector3d d);
void AplicarTHy(float theta, vector3d d);
void AplicarTHz(float theta, vector3d d);
Matrix THx,THy,THz,TH;

std::vector<Matrix> THList;
std::vector<vector3d> Origenes;
std::vector<modelo3D*> modelos;
/// Runge Kutta
double c3x,c3y;  //l1, l2 se declaro en la línea 29
void initRungeKutta();
void integrarRungeKutta();
void preparar();
double g, l, m,b;
unsigned int n;  //grados de libertad
Matrix Jc1,Jc2, Jc3, Jw1,Jw2,Jw3;
Matrix Y;
Matrix f1,f2,f3,f4;
double I1, I2, I3, m1, m2, m3;
double Q1, Q2, Q3,tau1,tau2,tau3;
Matrix M1, M2,M, M3;
Matrix H;
Matrix F(const Matrix &y, float t);
Matrix a(const Matrix &y, float t);
 unsigned int steps;
 double dt,t0,tf;
double t1,t2,t3,t1v,t2v,t3v;
///

Matrix X,Xv,Xvv,q,qv,qvv;
float theta1, theta2, theta3;
float z1, z2, z3;

void DefinirTHx(float theta, vector3d d);
void DefinirTHy(float theta, vector3d d);
void DefinirTHz(float theta, vector3d d);
void  Drawarrow3D( vector3d A,  vector3d B, vector3d color, double cota1,double R) ;


};

#endif // ROBOT_H
