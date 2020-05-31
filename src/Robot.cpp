#include "Robot.h"
#define PI 3.1415926535897932384626433832795
///Copyright (C) <2017>  <Eliseo Rivera> curso.tareas@gmail.com
Robot::Robot()
{
   theta1=0;theta2=0;theta3=0;
   THx.identity(4);
   THy.identity(4);
   THz.identity(4);
   int d=4;
   TH.identity(4);

      ///cinematica inversa
   O4.zero(4,1);
   O4L.zero(4,1);
   P.zero(4,1);
   V.zero(4,1);
   A1.zero(4,1);  //el tamaño 4x1 porque hay es dos posiciones iniciales y dos velocidades iniciales
   B1.zero(4,1);
   C1.zero(4,1);
   D1.zero(4,1);
 E1.zero(4,1);
   IsMoving=false;
  Jf.zero(3,3);
  Jfv.zero(3,3);
  Jg1.zero(3,3);
  Jg2.zero(3,3);
  Jg3.zero(3,3);
 X.zero(3,1); //vector de coordenzadas absolutas, lienales y angulares
  Xv.zero(3,1);
   Xvv.zero(3,1);
q.zero(3,1); //vector de coordenzadas generalizadas instantaneas
qv.zero(3,1);
qvv.zero(3,1);
hJtW.zero(3,1);
Q.zero(3,1);
   /// fin cinematica inversa
}

Robot::~Robot()
{

delete b1;
delete b2;
delete b3;

}

void Robot::inicializar(){
base=new modelo3D();
b1=new modelo3D();
b2=new modelo3D();
b3=new modelo3D();
base->leer("base.STL");
b1->leer("b1.STL");
b2->leer("b2.STL");
b3->leer("b3.STL");
base->color={0.0,0.0,0.0};
b1->color={0.89,0.74,0.003};;
b2->color={0.89,0.74,0.003};;
b3->color={0.3,0.3,0.3};
modelos.push_back(base);
modelos.push_back(b1);
modelos.push_back(b2);
modelos.push_back(b3);



}

void Robot::initphysicsconstants(){

    ///en el sistema internacional todas las constantes y variables
l1=0.5;
l2=0.3;
g=9.81;
/*
c3x=0.087; c3y=-0.04;//posicion del centroide no compuesto
m1=5, m2=2, m3=0.3;   //masas de los eslabones (eslabon ultimo compuesto)
   I1=.13; I2= 0.016, I3= 0.0012; //en kg*m^2
*/


c3x=0.197; c3y=-0.15538;//posicion del centroide  compuesto
  m1=5, m2=2, m3=1.3;   //masas de los eslabones (eslabon ultimo compuesto)
   I1=.13; I2= 0.016, I3= 1.2811e-2; //en kg*m^2


}

void Robot::configurarTH(){
theta1= 0.4;theta2=-2.1; theta3=1.05;
//theta1 = -0.2;theta2 = -0.95;theta3 = 0.55; //correcto
//theta1= 0.84105;theta2=3.8453; theta3= 1.071; //corecto en y
//theta1= 0.84105;theta2=4.1453; theta3=0.671787; //correcto en y+

//theta1= -0.640604;theta2=-7.97142; theta3=-5.74079; //correcto
//theta1 = 5.64258;theta2 = 2.44495;theta3 = 0.542396; //correcto

//theta1 = 6.07423;theta2 = 4.1453;theta3 = 0.671787; //correcto +y, +x+y revisar
DefinirTHz(0,{0,0,0}); //base
THList.push_back(THz);
DefinirTHx(0,{0,0,0}); //base
THList.push_back(THx);

DefinirTHz(theta1,{0,0,0}); //b1
THList.push_back(THz);
DefinirTHx(0,{l1,0,0}); //b1
THList.push_back(THx);

DefinirTHz( theta2,{0,0,0}); //b2
THList.push_back(THz);
DefinirTHx(0,{l2,0,0}); //b2
THList.push_back(THx);

DefinirTHz( theta3,{0,0,0}); //b3
THList.push_back(THz);
DefinirTHx(0,{0,0,0}); //b3
THList.push_back(THx);
initInverseKinematics();
initRungeKutta();

}
void Robot::initInverseKinematics()
{
    IsMoving=false;

TH.resetIdentity();
modelo3D *model;
for (int m=0;m<modelos.size();m++){

    model=modelos[m];
    TH=TH* THList[2*m]* THList[2*m+1];

}
xL=0.145;yL=-.08;  //posición local de punto de interés de trayectoria (no centroude de eslabon 3)
//xL=.087;yL=-.04;

//valores iniciales de las coordenadas absolutas del punto de interes sobre elabon final a segur trayectoria
///pero son valores iniciales a partir del control manual,por ello aún no se llama a la función Paramétrica
/// a este tiempo de la planificación de trayectorias se debe saber la velocidad inicial de punto de interés
t=0;  // es importante
O4L.entry(0,0)=xL; O4L.entry(1,0)=yL; O4L.entry(3,0)=1;
O4=TH*O4L;
P=O4;
x4=O4.entry(0,0);
y4=O4.entry(1,0);
THETA3=theta1+theta2+theta3;
THETA30=THETA3;



V.entry(0,0)=1;
V.entry(1,0)=1;
curva.resize(0);
curva.push_back(vector3d(O4.entry(0,0),O4.entry(1,0),O4.entry(2,0)));
//estado();
initPlanTrayec();
}
void Robot::initPlanTrayec(){

float vx0, vxf,x0,xf,vy0, vyf,y0,yf ,theta0, thetaf, vtheta0, vthetaf;

t0=0; tf =2;
x0=x4; y0=y4;
xf=x4+0.30; yf=y4+.6; //en metros

vx0=0; vxf=0;
vy0=0; vyf=0;

theta0=THETA30;
vtheta0=0;
thetaf=0.9;
vthetaf=0;
Matrix Sx(4,1), Sy(4,1), e1(4,1), e2(4,1), T(4,4), St(4,1);
T.entry(0,0)=t0*t0*t0;
T.entry(0,1)=t0*t0;
T.entry(0,2)=t0;
T.entry(0,3)=1;
T.entry(1,0)=3*t0*t0;
T.entry(1,1)=2*t0;
T.entry(1,2)=1;
T.entry(1,3)=0;

T.entry(2,0)=tf*tf*tf;
T.entry(2,1)=tf*tf;
T.entry(2,2)=tf;
T.entry(2,3)=1;
T.entry(3,0)=3*tf*tf;
T.entry(3,1)=2*tf;
T.entry(3,2)=1;
T.entry(3,3)=0;
//cout<<" T ="<<endl;
//T.mostrar();
Sx.entry(0,0)=x0;
Sx.entry(1,0)=vx0;
Sx.entry(2,0)=xf;
Sx.entry(3,0)=vxf;

Sy.entry(0,0)=y0;
Sy.entry(1,0)=vy0;
Sy.entry(2,0)=yf;
Sy.entry(3,0)=vyf;
St.entry(0,0)=theta0;
St.entry(1,0)=vtheta0;
St.entry(2,0)=thetaf;
St.entry(3,0)=vthetaf;

e1=T.inversa()*Sx;
e2=T.inversa()*Sy;
E1=T.inversa()*St;
/*
cout<<"e1 ="<<endl;
e1. mostrar();
cout<<"e2 ="<<endl;
e2.mostrar();
cout<<"E1 ="<<endl;
E1.mostrar();
*/
//cout<<"a1 ="<<endl;
A1.entry(0,0)=e1.entry(0,0);   //se puede reducir el tamaño de A1, debe ser de 2x1
A1.entry(1,0)=e2.entry(0,0);
//A1.mostrar();
//cout<<"b1 ="<<endl;
B1.entry(0,0)=e1.entry(1,0);
B1.entry(1,0)=e2.entry(1,0);
//B1.mostrar();
//cout<<"c1 ="<<endl;
C1.entry(0,0)=e1.entry(2,0);
C1.entry(1,0)=e2.entry(2,0);
//C1.mostrar();
//cout<<"d1 ="<<endl;
D1.entry(0,0)=e1.entry(3,0);
D1.entry(1,0)=e2.entry(3,0);
//D1.mostrar();
//cout<<"suma "<<endl;
//(A1+B1+D1).mostrar(); ///comentar

};
bool Robot::InverseKinematics(){

     cout<<"***************************inicio inverse kinematics****************"<<endl;
float a, b,c;
a=xL*cos(THETA3)-yL*sin(THETA3)-x4;
b=xL*sin(THETA3)+yL*cos(THETA3)-y4;
c=sqrt(a*a+b*b);
float df;
//cout<<"a ="<<a<<" , b = "<<b<<endl;
float ep=0.99;
float t1,t2,t3,t4=0;
if (a>0&&b<0){
    df=atan(fabs(b/a));
     float q1=(l2*l2-l1*l1-a*a-b*b)/(2*c*l1);
    if (fabs(q1)>ep) {IsMoving=false;cout<<"configuración no alcanzable, retorne usando control manual "<<endl; t=0;return false; };

    t1=acos(q1)-df;
    t2=2*PI-acos(q1)-df;
    t3=2*PI+acos(q1)-df;

     //   cout<<"t1  =" << t1<< " , t2 = " << t2<< " ,  t3= " << t3<< " , theta1 = " <<theta1<<endl;
    theta1=menor(t1,t2,t3,theta1);
    float t_=atan2(l1*sin(theta1)+b,l1*cos(theta1)+a)-theta1;
    t1=t_;
    t2=PI+t_;
    t3=2*PI+t_;
    t4=-PI+t_;
   //        cout<<"t1  =" << t1<< " , t2 = " << t2<< " ,  t3= " << t3<< " ,t4= " << t4<< "  theta2 = " <<theta2<<endl;
    theta2=menor2(t1,t2,t3,t4,theta2);
    theta3=THETA3-(theta1+theta2);


};
if (a<0&&b>0){
    df=atan(fabs(b/a));
     float q1=(l2*l2-l1*l1-a*a-b*b)/(-2*c*l1);
    if (fabs(q1)>ep) {cout<<"configuración no alcanzable, retorne usando control manual "<<endl; t=0;return false;};

    t1=acos(q1)-df;
    t2=2*PI-acos(q1)-df;
    t3=2*PI+acos(q1)-df;

  //  cout<<"t1  =" << t1<< " , t2 = " << t2<< " ,  t3= " << t3<< " , theta1 = " <<theta1<<endl;
    theta1=menor(t1,t2,t3,theta1);
    float t_=atan2(l1*sin(theta1)+b,l1*cos(theta1)+a)-theta1;
    t1=t_;
    t2=PI+t_;
    t3=2*PI+t_;
    t4=-PI+t_;
   //        cout<<"t1  =" << t1<< " , t2 = " << t2<< " ,  t3= " << t3<< " ,t4= " << t4<< "  theta2 = " <<theta2<<endl;
    theta2=menor2(t1,t2,t3,t4,theta2);
    theta3=THETA3-(theta1+theta2);

};
if (a>0&&b>0){
    df=atan(fabs(a/b));
     float q1=(l2*l2-l1*l1-a*a-b*b)/(2.0*c*l1);
    if (fabs(q1)>ep) {      cout<<"configuración no alcanzable, retorne usando control manual "<<endl; t=0 ;return false;};

    t1=asin(q1)-df;
    t2=PI-asin(q1)-df;
    t3=2*PI+asin(q1)-df;
      //  cout<<"t1  =" << t1<< " , t2 = " << t2<< " ,  t3= " << t3<< " , theta1 = " <<theta1<<endl;
    theta1=menor(t1,t2,t3,theta1);
     float t_=atan2(l1*sin(theta1)+b,l1*cos(theta1)+a)-theta1;
    t1=t_;
    t2=PI+t_;
    t3=2*PI+t_;
    t4=-PI+t_;
 //         cout<<"t1  =" << t1<< " , t2 = " << t2<< " ,  t3= " << t3<< " ,t4= " << t4<< "  theta2 = " <<theta2<<endl;
    theta2=menor2(t1,t2,t3,t4,theta2);
    theta3=THETA3-(theta1+theta2);


};
if (a<0&&b<0){
    df=atan(fabs(a/b));
    float q1=(l2*l2-l1*l1-a*a-b*b)/(-2.0*c*l1);
    if (fabs(q1)>ep) {cout<<"configuración no alcanzable, retorne usando control manual "<<endl; t=0; return false;};

   t1=asin(q1)-df;
   t2=PI-asin(q1)-df;
   t3=2*PI+asin(q1)-df;

  //     cout<<"t1  =" << t1<< " , t2 = " << t2<< " ,  t3= " << t3<< " ,t4= " << t4<< "  theta1 = " <<theta1<<endl;

    theta1=menor(t1,t2,t3,theta1);
    float t_=atan2(l1*sin(theta1)+b,l1*cos(theta1)+a)-theta1;
    t1=t_;
    t2=PI+t_;
    t3=2*PI+t_;
    t4=-PI+t_;
  //         cout<<"t1  =" << t1<< " , t2 = " << t2<< " ,  t3= " << t3<< " ,t4= " << t4<< "  theta2 = " <<theta2<<endl;
    theta2=menor2(t1,t2,t3,t4,theta2);
    theta3=THETA3-(theta1+theta2);


   // comprobracion(a,b);



};

 cout<<"***************************fin inverse kinematics****************"<<endl;
return true;
};
void Robot::Parametrica(){
 cout<<"//////////////////////////inicio parametrica-------------------------------"<<endl;
O4=t*t*t*A1+t*t*B1+t*C1+D1;  //de tamaño 4x1
x4=O4.entry(0,0);
y4=O4.entry(1,0);
THETA3=t*t*t*E1.entry(0,0)+t*t*E1.entry(1,0)+t*E1.entry(2,0)+E1.entry(3,0);
curva.push_back(vector3d(O4.entry(0,0),O4.entry(1,0),O4.entry(2,0)));




 cout<<"//////////////////////////inicio parametrica-------------------------------"<<endl;
}

void Robot::Move(){

if (IsMoving==true){
Parametrica();  /// en esta función se conocen X,XV,XVV
      bool state=InverseKinematics();
      //parametrica se llama despues
if (state==true) {

//se actualizan estados
/// si todo sale bien, antes de returnar se actualizan los valores de las coordenadas generalizadas



 cout<<"t ="<< t<<endl;
cout<<"theta1 ="<<theta1<<endl;
cout<<"theta2 ="<<theta2<<endl;
cout<<"theta3 ="<<theta3<<endl;
cout<<"0theta3 ="<<theta1+theta2+theta3<<endl;





DefinirTHz(theta1, vector3d (0,0,0));
THList[2]=THz;
DefinirTHz(theta2, vector3d (0,0,0));
THList[4]=THz;
DefinirTHz(theta3, vector3d (0,0,0));
THList[6]=THz;

///velocidades generalizadas

X.entry(0,0)=x4;
X.entry(1,0)=y4;
X.entry(2,0)=THETA3;

Xv.entry(0,0)=((3*t*t)*A1+(2*t)*B1+C1).entry(0,0);
Xv.entry(1,0)=((3*t*t)*A1+(2*t)*B1+C1).entry(1,0);
Xv.entry(2,0)=3*t*t*E1.entry(0,0)+2*t*E1.entry(1,0)+E1.entry(2,0);

Xvv.entry(0,0)=((6*t)*A1+(2)*B1).entry(0,0);
Xvv.entry(1,0)=((6*t)*A1+(2)*B1).entry(1,0);
Xvv.entry(2,0)=6*t*E1.entry(0,0)+2*E1.entry(1,0);


q.entry(0,0)=theta1;
q.entry(1,0)=theta2;
q.entry(2,0)=theta3;
cout<<"velocidades generalizadas "<<endl;

Calc_Jf(theta1,theta2,theta3);
qv=Jf.inversa()*Xv;
qv.mostrar();


Calc_Jfv(theta1,theta2,theta3);
qvv=Jf.inversa()*(Xvv-Jfv*qv);
cout<<"aceleraciones generalizadas "<<endl;
qvv.mostrar();


///Calculo de pares

Calc_M(theta1,theta2,theta3);

cout<<"M =  "<<endl;
M.mostrar();
Calc_hJtW(theta1,theta2,theta3,qv.entry(0,0),qv.entry(1,0),qv.entry(2,0) );
cout<<"h-Jt*w =  "<<endl;
hJtW.mostrar();

Q=M*qvv+hJtW;

cout<<"fuerzas generalizadas "<<endl;
Q.mostrar();

cout<<"velocidades absolutas "<<endl;
Xv.mostrar();
cout<<"aceleraciones absolutas "<<endl;
Xvv.mostrar();
t=  t+0.04; // si todo ha salido bien aumenta el tiempo, de lo contrario no tiene caso



/*
cout<<"/////////////////////////////transformations"<<endl;
 cout<<"01T"<<endl;
THList[2].mostrar();
 cout<<"12T"<<endl;
THList[4].mostrar();
 cout<<"23T"<<endl;
THList[6].mostrar();

 cout<<"02T"<<endl;
 (THList[2]*THList[4]).mostrar();
 cout<<"03T"<<endl;
 (THList[2]*THList[4]*THList[6]).mostrar();

 */
// si el tiempo transcurre entonces se recorre la parametrica conocida y se puede calcular la velocidad
}
else IsMoving=false;

}
if (t>tf){   IsMoving=false; initInverseKinematics();   };

}
void Robot::Calc_Jf(float t1, float t2, float t3)
{
    float x=xL; float y=yL;
Jf.entry(0,0)=- x*sin(t1 + t2 + t3) - l2*sin(t1 + t2) - l1*sin(t1) - y*cos(t1 + t2 + t3);
Jf.entry(0,1)=- x*sin(t1 + t2 + t3) - l2*sin(t1 + t2) - y*cos(t1 + t2 + t3);
Jf.entry(0,2)=- x*sin(t1 + t2 + t3) - y*cos(t1 + t2 + t3);

Jf.entry(1,0)=l2*cos(t1 + t2) - y*sin(t1 + t2 + t3) + l1*cos(t1) + x*cos(t1 + t2 + t3);
Jf.entry(1,1)= l2*cos(t1 + t2) - y*sin(t1 + t2 + t3) + x*cos(t1 + t2 + t3);
Jf.entry(1,2)=x*cos(t1 + t2 + t3) - y*sin(t1 + t2 + t3);

Jf.entry(2,0)=1;
Jf.entry(2,1)=1;
Jf.entry(2,2)=1;

}
void Robot::Calc_M(float t1, float t2, float t3)
{
       ////////////
M.entry(0,0)=I1+I2+I3+(c3x*c3x)*m3+(c3y*c3y)*m3+(l1*l1)*m1*2.5E-1+(l1*l1)*m2+(l1*l1)*m3+(l2*l2)*m2*2.5E-1+(l2*l2)*m3+c3x*l1*m3*cos(t2+t3)*2.0-c3y*l1*m3*sin(t2+t3)*2.0+c3x*l2*m3*cos(t3)*2.0+l1*l2*m2*cos(t2)+l1*l2*m3*cos(t2)*2.0-c3y*l2*m3*sin(t3)*2.0;
M.entry(0,1)=I2+I3+(c3x*c3x)*m3+(c3y*c3y)*m3+(l2*l2)*m2*2.5E-1+(l2*l2)*m3+c3x*l1*m3*cos(t2+t3)-c3y*l1*m3*sin(t2+t3)*1.0+c3x*l2*m3*cos(t3)*2.0+l1*l2*m2*cos(t2)*5.0E-1+l1*l2*m3*cos(t2)-c3y*l2*m3*sin(t3)*2.0;
M.entry(0,2)=I3+(c3x*c3x)*m3+(c3y*c3y)*m3+c3x*l1*m3*cos(t2+t3)-c3y*l1*m3*sin(t2+t3)*1.0+c3x*l2*m3*cos(t3)-c3y*l2*m3*sin(t3)*1.0;


M.entry(1,0)=I2+I3+(c3x*c3x)*m3+(c3y*c3y)*m3+(l2*l2)*m2*2.5E-1+(l2*l2)*m3+c3x*l1*m3*cos(t2+t3)-c3y*l1*m3*sin(t2+t3)*1.0+c3x*l2*m3*cos(t3)*2.0+l1*l2*m2*cos(t2)*5.0E-1+l1*l2*m3*cos(t2)-c3y*l2*m3*sin(t3)*2.0;
M.entry(1,1)=I2+I3+(c3x*c3x)*m3+(c3y*c3y)*m3+(l2*l2)*m2*(1.0/4.0)+(l2*l2)*m3+c3x*l2*m3*cos(t3)*2.0-c3y*l2*m3*sin(t3)*2.0;
M.entry(1,2)=I3+(c3x*c3x)*m3+(c3y*c3y)*m3+c3x*l2*m3*cos(t3)-c3y*l2*m3*sin(t3)*1.0;


M.entry(2,0)= I3+(c3x*c3x)*m3+(c3y*c3y)*m3+c3x*l1*m3*cos(t2+t3)-c3y*l1*m3*sin(t2+t3)*1.0+c3x*l2*m3*cos(t3)-c3y*l2*m3*sin(t3)*1.0;
M.entry(2,1)= I3+(c3x*c3x)*m3+(c3y*c3y)*m3+c3x*l2*m3*cos(t3)-c3y*l2*m3*sin(t3)*1.0;
M.entry(2,2)=I3+m3*(c3x*c3x+c3y*c3y);


}
void Robot::Calc_hJtW(float t1, float t2, float t3, float t1v, float t2v, float t3v){
hJtW.entry(0,0)=g*m3*(l2*cos(t1+t2)+l1*cos(t1)+c3x*cos(t1+t2+t3)-c3y*sin(t1+t2+t3)*1.0)+g*m2*(l2*cos(t1+t2)*(1.0/2.0)+l1*cos(t1))+g*l1*m1*cos(t1)*(1.0/2.0)-c3y*l1*m3*(t2v*t2v)*cos(t2+t3)*1.0-c3y*l1*m3*(t3v*t3v)*cos(t2+t3)*1.0-c3x*l1*m3*(t2v*t2v)*sin(t2+t3)*1.0-c3x*l1*m3*(t3v*t3v)*sin(t2+t3)*1.0-c3y*l2*m3*(t3v*t3v)*cos(t3)*1.0-c3x*l2*m3*(t3v*t3v)*sin(t3)*1.0-l1*l2*m2*(t2v*t2v)*sin(t2)*5.0E-1-l1*l2*m3*(t2v*t2v)*sin(t2)*1.0-c3y*l1*m3*t1v*t2v*cos(t2+t3)*2.0-c3y*l1*m3*t1v*t3v*cos(t2+t3)*2.0-c3y*l1*m3*t2v*t3v*cos(t2+t3)*2.0-c3x*l1*m3*t1v*t2v*sin(t2+t3)*2.0-c3x*l1*m3*t1v*t3v*sin(t2+t3)*2.0-c3x*l1*m3*t2v*t3v*sin(t2+t3)*2.0-c3y*l2*m3*t1v*t3v*cos(t3)*2.0-c3y*l2*m3*t2v*t3v*cos(t3)*2.0-c3x*l2*m3*t1v*t3v*sin(t3)*2.0-c3x*l2*m3*t2v*t3v*sin(t3)*2.0-l1*l2*m2*t1v*t2v*sin(t2)*1.0-l1*l2*m3*t1v*t2v*sin(t2)*2.0;

hJtW.entry(1,0)=g*m3*(l2*cos(t1+t2)+c3x*cos(t1+t2+t3)-c3y*sin(t1+t2+t3)*1.0)+g*l2*m2*cos(t1+t2)*(1.0/2.0)+c3y*l1*m3*(t1v*t1v)*cos(t2+t3)+c3x*l1*m3*(t1v*t1v)*sin(t2+t3)-c3y*l2*m3*(t3v*t3v)*cos(t3)*1.0-c3x*l2*m3*(t3v*t3v)*sin(t3)*1.0+l1*l2*m2*(t1v*t1v)*sin(t2)*5.0E-1+l1*l2*m3*(t1v*t1v)*sin(t2)-c3y*l2*m3*t1v*t3v*cos(t3)*2.0-c3y*l2*m3*t2v*t3v*cos(t3)*2.0-c3x*l2*m3*t1v*t3v*sin(t3)*2.0-c3x*l2*m3*t2v*t3v*sin(t3)*2.0;

hJtW.entry(2,0)=t2v*(t1v*(c3y*l2*m3*cos(t3)*2.0+c3x*l2*m3*sin(t3)*2.0)*(1.0/2.0)+t2v*(c3y*l2*m3*cos(t3)*2.0+c3x*l2*m3*sin(t3)*2.0)*(1.0/2.0))+t1v*(t2v*(c3y*l2*m3*cos(t3)*2.0+c3x*l2*m3*sin(t3)*2.0)*(1.0/2.0)+t1v*(c3y*l1*m3*cos(t2+t3)*2.0+c3x*l1*m3*sin(t2+t3)*2.0+c3y*l2*m3*cos(t3)*2.0+c3x*l2*m3*sin(t3)*2.0)*(1.0/2.0))+g*m3*(c3x*cos(t1+t2+t3)-c3y*sin(t1+t2+t3)*1.0);

}
void Robot::Calc_Jg(float t1, float t2, float t3)
{
    float x=xL; float y=yL;
Jg1.entry(0,0)= y*sin(t1 + t2 + t3) - l2*cos(t1 + t2) - l1*cos(t1) - x*cos(t1 + t2 + t3);
Jg1.entry(0,1)= y*sin(t1 + t2 + t3) - l2*cos(t1 + t2) - x*cos(t1 + t2 + t3);
Jg1.entry(0,2)=y*sin(t1 + t2 + t3) - x*cos(t1 + t2 + t3);

Jg1.entry(1,0)=- x*sin(t1 + t2 + t3) - l2*sin(t1 + t2) - l1*sin(t1) - y*cos(t1 + t2 + t3);
Jg1.entry(1,1)=- x*sin(t1 + t2 + t3) - l2*sin(t1 + t2) - y*cos(t1 + t2 + t3);
Jg1.entry(1,2)=- x*sin(t1 + t2 + t3) - y*cos(t1 + t2 + t3);

Jg1.entry(2,0)=0;
Jg1.entry(2,1)=0;
Jg1.entry(2,2)=0;
cout<<" Jg1---"<<endl;
Jg1.mostrar();
///----------------------------------------------

Jg2.entry(0,0)=  y*sin(t1 + t2 + t3) - l2*cos(t1 + t2) - x*cos(t1 + t2 + t3);
Jg2.entry(0,1)=y*sin(t1 + t2 + t3) - l2*cos(t1 + t2) - x*cos(t1 + t2 + t3);
Jg2.entry(0,2)=y*sin(t1 + t2 + t3) - x*cos(t1 + t2 + t3);

Jg2.entry(1,0)= - x*sin(t1 + t2 + t3) - l2*sin(t1 + t2) - y*cos(t1 + t2 + t3);
Jg2.entry(1,1)=- x*sin(t1 + t2 + t3) - l2*sin(t1 + t2) - y*cos(t1 + t2 + t3);
Jg2.entry(1,2)=- x*sin(t1 + t2 + t3) - y*cos(t1 + t2 + t3);

Jg2.entry(2,0)=0;
Jg2.entry(2,1)=0;
Jg2.entry(2,2)=0;

///-*----------------------------------------------
Jg3.entry(0,0)= y*sin(t1 + t2 + t3) - x*cos(t1 + t2 + t3);
Jg3.entry(0,1)= y*sin(t1 + t2 + t3) - x*cos(t1 + t2 + t3);
Jg3.entry(0,2)= y*sin(t1 + t2 + t3) - x*cos(t1 + t2 + t3);

Jg3.entry(1,0)= - x*sin(t1 + t2 + t3) - y*cos(t1 + t2 + t3);
Jg3.entry(1,1)=- x*sin(t1 + t2 + t3) - y*cos(t1 + t2 + t3);
Jg3.entry(1,2)=- x*sin(t1 + t2 + t3) - y*cos(t1 + t2 + t3);

Jg3.entry(2,0)=0;
Jg3.entry(2,1)=0;
Jg3.entry(2,2)=0;


}
void Robot::Calc_Jfv(float t1, float t2, float t3)
{
    Calc_Jg(t1,  t2,  t3);
Matrix Jfv1(Jg1*qv),Jfv2(Jg2*qv),Jfv3(Jg3*qv);
cout<<"Jfv1  *********"<<endl;
Jfv1.mostrar();
cout<<"Jfv2  *********"<<endl;
Jfv2.mostrar();
cout<<"Jfv3  *********"<<endl;
Jfv3.mostrar();


Jfv.entry(0,0)=Jfv1.entry(0,0);
Jfv.entry(1,0)=Jfv1.entry(1,0);
Jfv.entry(2,0)=Jfv1.entry(2,0);

Jfv.entry(0,1)=Jfv2.entry(0,0);
Jfv.entry(1,1)=Jfv2.entry(1,0);
Jfv.entry(2,1)=Jfv2.entry(2,0);

Jfv.entry(0,2)=Jfv3.entry(0,0);
Jfv.entry(1,2)=Jfv3.entry(1,0);
Jfv.entry(2,2)=Jfv3.entry(2,0);

}
void Robot::estado()
{
cout<<"******************************inicio estado****************************************************"<<endl;
cout<<"tiempo = "<<t<<endl;

cout<<"radianes----------"<<endl;
cout<<"theta1 = "<<theta1<<endl;
cout<<"theta2 = "<<theta2<<endl;
cout<<"theta3 = "<<theta3<<endl;
cout<<"0theta3 = "<<theta1+theta2+theta3<<endl; //gloval
/*
cout<<"grados------------"<<endl;
cout<<"theta1 = "<<theta1*180/PI<<endl;
cout<<"theta2 = "<<theta2*180/PI<<endl;
cout<<"theta3 = "<<theta3*180/PI<<endl;
cout<<"0theta3 = "<<(theta1+theta2+theta3)*180/PI<<endl;
*/
//calcular velocidades Jf

//calcular aceleraciones

//THList[2].mostrar();
//THList[4].mostrar();
//THList[6].mostrar();

cout<<"THETA3 = "<<THETA3<<endl;
//Matrix O34(4,1);
//O34.entry(0,0)=xL;
//O34.entry(1,0)=yL;
//O34.entry(3,0)=1;
//O34=TH*O34;  //posicion de O4 respecto del sistema de referencia global
//vector3d O4(O34.entry(0,0),O34.entry(1,0),0);


//O4.mostrar();
cout<<"-----------------------------fin estado----------------------------------------------------"<<endl;

}
void Robot::comprobracion(float a, float b){
float c1, c2; //a, b son variables locales
c1=l1*cos(theta1)+l2*cos(theta1+theta2)+a;
c2=l1*sin(theta1)+l2*sin(theta1+theta2)+b;
cout<<"-----------------------------------COMPROBACION-----------------------------------------------"<<endl;
cout<<" c1 =  "<<c1<<endl;
cout<<" c2 =  "<<c2<<endl;
cout<<"theta1 ="<<theta1<<endl;
cout<<"theta2 ="<<theta2<<endl;
cout<<"theta3 ="<<theta3<<endl;

cout<<"0theta3 = "<<theta1+theta2+theta3<<endl;
cout<<"0THETA3 = "<<THETA3<<endl;

cout<<"-----------------------------------FFFFFFIN COMPROBACION-----------------------------------------------"<<endl;
}
void Robot::initRungeKutta(){
initphysicsconstants();
 //El signo se ingreso en los vectores (pesos) wi
n=3;  //grados de libertad
Jc1.zero(3,n),Jc2.zero(3,n),Jc3.zero(3,n), Jw1.zero(3,n),Jw2.zero(3,n),Jw3.zero(3,n), Q.zero(n,1);
Y.zero(2*n,1);
f1.zero(2*n,1),f2.zero(2*n,1),f3.zero(2*n,1),f4.zero(2*n,1);
M1.zero(n,n), M2.zero(n,n), M3.zero(n,n),M.zero(n,n);
H.zero(n,1);
///////////////////////////////////////

  // I1=.13; I2= 0.016, I3= 0.0012;  //momentos de inercia de los eslabones (efector final no compuesto)
 //  m1=5, m2=2, m3=0.3;


t1=theta1;
t2=theta2;
t3=theta3;
t1v=0;t2v=0;t3v=0;

Y.entry(0,0)=t1;
Y.entry(1,0)=t2;
Y.entry(2,0)=t3;
Y.entry(3,0)=t1v;
Y.entry(4,0)=t2v;
Y.entry(5,0)=t3v;


 //  M.mostrar();

   tau1=0;  //hay que verificar la conservación del momento angular
   tau2=0;
   tau3=0;



t0=0;
dt=0.005;
t=t0;


};

Matrix Robot::a(const Matrix &y, float t){

t1=y.entry(0,0);
t2=y.entry(1,0);
t3=y.entry(2,0);

t1v=y.entry(3,0);
t2v=y.entry(4,0);
t3v=y.entry(5,0);


Calc_M(t1 ,t2 ,t3);

H.entry(0,0)=tau1-g*l2*m2*cos(t1+t2)*5.0E-1-g*l2*m3*cos(t1+t2)*1.0-g*l1*m1*cos(t1)*5.0E-1-g*l1*m2*cos(t1)*1.0-g*l1*m3*cos(t1)*1.0-c3x*g*m3*cos(t1+t2+t3)*1.0+c3y*g*m3*sin(t1+t2+t3)+c3y*l1*m3*(t2v*t2v)*cos(t2+t3)+c3y*l1*m3*(t3v*t3v)*cos(t2+t3)+c3x*l1*m3*(t2v*t2v)*sin(t2+t3)+c3x*l1*m3*(t3v*t3v)*sin(t2+t3)+c3y*l2*m3*(t3v*t3v)*cos(t3)+c3x*l2*m3*(t3v*t3v)*sin(t3)+l1*l2*m2*(t2v*t2v)*sin(t2)*5.0E-1+l1*l2*m3*(t2v*t2v)*sin(t2)+c3y*l1*m3*t1v*t2v*cos(t2+t3)*2.0+c3y*l1*m3*t1v*t3v*cos(t2+t3)*2.0+c3y*l1*m3*t2v*t3v*cos(t2+t3)*2.0+c3x*l1*m3*t1v*t2v*sin(t2+t3)*2.0+c3x*l1*m3*t1v*t3v*sin(t2+t3)*2.0+c3x*l1*m3*t2v*t3v*sin(t2+t3)*2.0+c3y*l2*m3*t1v*t3v*cos(t3)*2.0+c3y*l2*m3*t2v*t3v*cos(t3)*2.0+c3x*l2*m3*t1v*t3v*sin(t3)*2.0+c3x*l2*m3*t2v*t3v*sin(t3)*2.0+l1*l2*m2*t1v*t2v*sin(t2)+l1*l2*m3*t1v*t2v*sin(t2)*2.0;


H.entry(1,0)=  tau2-g*l2*m2*cos(t1+t2)*5.0E-1-g*l2*m3*cos(t1+t2)*1.0-c3x*g*m3*cos(t1+t2+t3)*1.0+c3y*g*m3*sin(t1+t2+t3)-c3y*l1*m3*(t1v*t1v)*cos(t2+t3)*1.0-c3x*l1*m3*(t1v*t1v)*sin(t2+t3)*1.0+c3y*l2*m3*(t3v*t3v)*cos(t3)+c3x*l2*m3*(t3v*t3v)*sin(t3)-l1*l2*m2*(t1v*t1v)*sin(t2)*5.0E-1-l1*l2*m3*(t1v*t1v)*sin(t2)*1.0+c3y*l2*m3*t1v*t3v*cos(t3)*2.0+c3y*l2*m3*t2v*t3v*cos(t3)*2.0+c3x*l2*m3*t1v*t3v*sin(t3)*2.0+c3x*l2*m3*t2v*t3v*sin(t3)*2.0;


H.entry(2,0)=tau3-c3x*g*m3*cos(t1+t2+t3)+c3y*g*m3*sin(t1+t2+t3)-c3y*l1*m3*(t1v*t1v)*cos(t2+t3)-c3x*l1*m3*(t1v*t1v)*sin(t2+t3)-c3y*l2*m3*(t1v*t1v)*cos(t3)-c3y*l2*m3*(t2v*t2v)*cos(t3)-c3x*l2*m3*(t1v*t1v)*sin(t3)-c3x*l2*m3*(t2v*t2v)*sin(t3)-c3y*l2*m3*t1v*t2v*cos(t3)*2.0-c3x*l2*m3*t1v*t2v*sin(t3)*2.0;







return M.inversa()*H;
}

Matrix Robot::F(const Matrix &y, float t){
    //n is the number of degree of freedom
Matrix f(2*n,1);
Matrix a_(a(y,t));

for (int i=0;i<n;i++){
f.entry(i,0)=y.entry(n+i,0);

}

for (int i=0;i<n;i++){
f.entry(i+n,0)=a_.entry(i,0);

}

return f;

}

void Robot::preparar(){
//Y.mostrar();
f1=dt*F(Y,t);
f2=dt*F(Y+0.5*f1,t+0.5*dt);
f3=dt*F(Y+0.5*f2,t+0.5*dt);
f4=dt*F(Y+f3,t+dt);
Y=Y+(0.16666666666666666666666666666667)*(f1+2*f2+2*f3+f4);
t=t+dt;

t1=Y.entry(0,0);
t2=Y.entry(1,0);
t3=Y.entry(2,0);
DefinirTHz(t1, vector3d (0,0,0));
THList[2]=THz;
DefinirTHz(t2, vector3d (0,0,0));
THList[4]=THz;
DefinirTHz(t3, vector3d (0,0,0));
THList[6]=THz;



}

void Robot::renderizar(){


TH.resetIdentity();

modelo3D *model;

for (int m=0;m<modelos.size();m++){

    model=modelos[m];
    TH=TH* THList[2*m]* THList[2*m+1];


vector3d ux,uy,uz,O;
ux={1,0,0};
uy={0,1,0};
uz={0,0,1};

Matrix ux4(ux,1),uy4(uy,1),uz4(uz,1),O4(O,1);


ux4=TH*ux4-TH*O4;
uy4=TH*uy4-TH*O4;
uz4=TH*uz4-TH*O4;
O4=TH*O4;


ux={ux4.aij[0][0],ux4.aij[1][0],ux4.aij[2][0]};
uy={uy4.aij[0][0],uy4.aij[1][0],uy4.aij[2][0]};
uz={uz4.aij[0][0],uz4.aij[1][0],uz4.aij[2][0]};
O={O4.aij[0][0],O4.aij[1][0],O4.aij[2][0]};

//if (m<2){
         Drawarrow3D(O,O+.1*ux,{1,0.1,0.2},0.03,0.0051);
         Drawarrow3D(O,O+.1*uy,{.1,1,0.2},0.03,0.0051);
        // Drawarrow3D(O,O+.2*uz,{0.1,0.2,1},0.03,0.01);
       //  }
         glColor4f(model->color.x,model->color.y,model->color.z, 0.5);

glEnable(GL_BLEND);
 glBegin(GL_TRIANGLES);

  glFrontFace(GL_FRONT_AND_BACK);
    for (int i=0;i<model->ntriangles;i++){

vector3d v1=model->triangulos[i].vertices[0];   //posiciones locales
vector3d v2=model->triangulos[i].vertices[1];
vector3d v3=model->triangulos[i].vertices[2];
Matrix v14(v1,1),v24(v2,1),v34(v3,1);

v14=TH*v14;
v24=TH*v24;
v34=TH*v34;
v1={v14.entry(0,0),v14.entry(1,0),v14.entry(2,0)};
v2={v24.entry(0,0),v24.entry(1,0),v24.entry(2,0)};
v3={v34.entry(0,0),v34.entry(1,0),v34.entry(2,0)};



Matrix N(4,1),d14(4,1),d24(4,1);
d14=v24-v14;
d24=v34-v14;
vector3d d1,d2,n;
d1={d14.entry(0,0),d14.entry(1,0),d14.entry(2,0)};
d2={d24.entry(0,0),d24.entry(1,0),d24.entry(2,0)};
n=d1*d2;  ///devuelve el producto vectorial
n.normalize();



        glNormal3f(n.x,n.y,n.z);
        glVertex3f(v1.x,v1.y,v1.z);
        glVertex3f(v2.x,v2.y,v2.z);
        glVertex3f(v3.x,v3.y,v3.z);
    }
glEnd();
// }
 glDisable(GL_BLEND);


///DIBUJAR EJES


//}
}



}

void Robot::DefinirTHx(float dtheta, vector3d d){

THx.aij[0][0]=1;
THx.aij[0][1]=0;
THx.aij[0][2]=0;
THx.aij[0][3]=d.x;

THx.aij[1][0]=0;
THx.aij[1][1]=cos(dtheta);
THx.aij[1][2]=-sin(dtheta);
THx.aij[1][3]=d.y;

THx.aij[2][0]=0;
THx.aij[2][1]=sin(dtheta);
THx.aij[2][2]=cos(dtheta);
THx.aij[2][3]=d.z;

THx.aij[3][0]=0;
THx.aij[3][1]=0;
THx.aij[3][2]=0;
THx.aij[3][3]=1;

}
void Robot::DefinirTHy(float dtheta, vector3d d){


THy.aij[0][0]=cos(dtheta);
THy.aij[0][1]=0;
THy.aij[0][2]=sin(dtheta);

THy.aij[1][0]=0;
THy.aij[1][1]=1;
THy.aij[1][2]=0;

THy.aij[2][0]=-sin(dtheta);
THy.aij[2][1]=0;
THy.aij[2][2]=cos(dtheta);

THy.aij[3][0]=0;
THy.aij[3][1]=0;
THy.aij[3][2]=0;
THy.aij[3][3]=1;

THy.aij[0][3]=d.x;
THy.aij[1][3]=d.y;
THy.aij[2][3]=d.z;
}
void Robot::DefinirTHz(float dtheta, vector3d d){

THz.aij[0][0]=cos(dtheta);
THz.aij[0][1]=-sin(dtheta);
THz.aij[0][2]=0;
THz.aij[0][3]=d.x;

THz.aij[1][0]=sin(dtheta);
THz.aij[1][1]=cos(dtheta);
THz.aij[1][2]=0;
THz.aij[1][3]=d.y;

THz.aij[2][0]=0;
THz.aij[2][1]=0;
THz.aij[2][2]=1;
THz.aij[2][3]=d.z;

THz.aij[3][0]=0;
THz.aij[3][1]=0;
THz.aij[3][2]=0;
THz.aij[3][3]=1;

}

void Robot::AplicarTHx(float theta, vector3d d){
theta=theta*PI/180.0;

DefinirTHx(theta,d);

}
void Robot::AplicarTHy(float theta, vector3d d){
theta=theta*PI/180.0;
DefinirTHy(theta,d);

}
void Robot::AplicarTHz(float theta, vector3d d){
theta=theta*PI/180.0;
DefinirTHz(theta,d);

}

void Robot::Drawarrow3D( vector3d A,  vector3d B, vector3d color, double cota1,double R)
{

double color1,color2,color3,a,b,c,d,e;



color1=color.x;//abs(color1/255);
color2=color.y;//abs(color2/255);
color3=color.z;//abs(color3/255);

glColor3f( color1,color2, color3);

vector3d n=B-A,np,vertex[10],normallight;
n.normalize();
if(n.z!=0)np={1,1,(-1/n.z)*(n.x+n.y)};
else if(n.y!=0)np={1,(-1/n.y)*(n.x+n.z),1};
else np={(-1/n.x)*(n.y+n.z),1,1};

np.normalize();
vertex[0]=R*np;
vertex[2]=R*(n*np).normalize();
vertex[1]=R*((0.5)*(vertex[2]-vertex[0])+vertex[0]).normalize();
vertex[4]=R*(n*vertex[2]).normalize();
vertex[3]=R*((0.5)*(vertex[4]-vertex[2])+vertex[2]).normalize();
vertex[6]=R*(n*vertex[4]).normalize();
vertex[5]=R*((0.5)*(vertex[6]-vertex[4])+vertex[4]).normalize();
vertex[7]=R*((0.5)*(vertex[0]-vertex[6])+vertex[6]).normalize();
vertex[8]=vertex[0];
vertex[9]=vertex[1];
int nx=8;
double d_thetha,fraccion=0.1,radioflecha=R+.7*R;
d_thetha=2.0f*PI/nx;


  ///tubos
 glBegin( GL_TRIANGLE_STRIP );

         for(int i=0;i<9;i++)
               {

normallight=n*(vertex[i-1]-vertex[i+1]);
normallight.normalize();
glNormal3f(normallight.x, normallight.y, normallight.z);
                 glVertex3f(vertex[i].x+A.x,vertex[i].y+A.y,vertex[i].z+A.z);

glVertex3f(vertex[i].x+B.x-fraccion*(B.x-A.x),vertex[i].y+B.y-fraccion*(B.y-A.y),vertex[i].z+B.z-fraccion*(B.z-A.z));

    // top face

                }

glEnd();



//flecha
vertex[0]=radioflecha*np;
vertex[2]=radioflecha*(n*np).normalize();
vertex[1]=radioflecha*((0.5)*(vertex[2]-vertex[0])+vertex[0]).normalize();
vertex[4]=radioflecha*(n*vertex[2]).normalize();
vertex[3]=radioflecha*((0.5)*(vertex[4]-vertex[2])+vertex[2]).normalize();
vertex[6]=radioflecha*(n*vertex[4]).normalize();
vertex[5]=radioflecha*((0.5)*(vertex[6]-vertex[4])+vertex[4]).normalize();
vertex[7]=radioflecha*((0.5)*(vertex[0]-vertex[6])+vertex[6]).normalize();
vertex[8]=vertex[0];
vertex[9]=vertex[1];
vector3d Ap(B-fraccion*(B-A));



 glBegin( GL_TRIANGLE_STRIP );  //flecha

         for(int i=0;i<9;i++)
               {

normallight=n*(vertex[i-1]-vertex[i+1]);
normallight.normalize();
glNormal3f(normallight.x, normallight.y, normallight.z);
                 glVertex3f(vertex[i].x+Ap.x,vertex[i].y+Ap.y,vertex[i].z+Ap.z);


glNormal3f(n.x, n.y, n.z);
glVertex3f(Ap.x+fraccion*(B-A).x,Ap.y+fraccion*(B-A).y,Ap.z+fraccion*(B-A).z);

    // top face

                }

glEnd();


}
void Robot::DibujarCurva(){

glBegin( GL_LINES);
for (int i=1;i<curva.size();i++){
    glVertex3f(curva[i-1].x,curva[i-1].y,curva[i-1].z);
     glVertex3f(curva[i].x,curva[i].y,curva[i].z);
    glColor3f(0,0,0);
}
glEnd();
return;
}
