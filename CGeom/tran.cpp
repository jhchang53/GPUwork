/*
   tranrot.cpp
*/
#include <stdio.h>
#include <stdlib.h>
#include "TRcoord.h"

void drawPoint(double u[], double sxy, double sz)
{
  printf("## point\n");
  printf("%.3lf %.3lf %.3lf  %.3lf %.3lf %.3lf\n",u[0],u[1],u[2], sxy,0.0,0.0);
  printf("%.3lf %.3lf %.3lf  %.3lf %.3lf %.3lf\n",u[0],u[1],u[2],-sxy,0.0,0.0);
  printf("%.3lf %.3lf %.3lf  %.3lf %.3lf %.3lf\n",u[0],u[1],u[2],0.0, sxy,0.0);
  printf("%.3lf %.3lf %.3lf  %.3lf %.3lf %.3lf\n",u[0],u[1],u[2],0.0,-sxy,0.0);
  printf("%.3lf %.3lf %.3lf  %.3lf %.3lf %.3lf\n",u[0],u[1],u[2],0.0,0.0, sz);
  printf("%.3lf %.3lf %.3lf  %.3lf %.3lf %.3lf\n",u[0],u[1],u[2],0.0,0.0,-sz);
};


void drawBox(double u0[], double u1[], double u2[], double u3[],
  double u4[], double u5[], double u6[], double u7[])
{
  /* draw box using 4 lower point and 4 upper points  */
  printf("## UV vox\n");
  printf("%.3lf %.3lf %.3lf ",u0[0],u0[1],u0[2]);
  printf("  %.3lf %.3lf %.3lf\n",u1[0]-u0[0],u1[1]-u0[1],u1[2]-u0[2]);
  printf("%.3lf %.3lf %.3lf ",u2[0],u2[1],u2[2]);
  printf("  %.3lf %.3lf %.3lf\n",u3[0]-u2[0],u3[1]-u2[1],u3[2]-u2[2]);
  //  horizontal
  printf("%.3lf %.3lf %.3lf ",u0[0],u0[1],u0[2]);
  printf("  %.3lf %.3lf %.3lf\n",u2[0]-u0[0],u2[1]-u0[1],u2[2]-u0[2]);
  printf("%.3lf %.3lf %.3lf ",u1[0],u1[1],u1[2]);
  printf("  %.3lf %.3lf %.3lf\n",u3[0]-u1[0],u3[1]-u1[1],u3[2]-u1[2]);

  printf("%.3lf %.3lf %.3lf ",u4[0],u4[1],u4[2]);
  printf("  %.3lf %.3lf %.3lf\n",u5[0]-u4[0],u5[1]-u4[1],u5[2]-u4[2]);
  printf("%.3lf %.3lf %.3lf ",u6[0],u6[1],u6[2]);
  printf("  %.3lf %.3lf %.3lf\n",u7[0]-u6[0],u7[1]-u6[1],u7[2]-u6[2]);
  //  horizontal
  printf("%.3lf %.3lf %.3lf ",u4[0],u4[1],u4[2]);
  printf("  %.3lf %.3lf %.3lf\n",u6[0]-u4[0],u6[1]-u4[1],u6[2]-u4[2]);
  printf("%.3lf %.3lf %.3lf ",u5[0],u5[1],u5[2]);
  printf("  %.3lf %.3lf %.3lf\n",u7[0]-u5[0],u7[1]-u5[1],u7[2]-u5[2]);
  //  vertical pillar
  printf("%.3lf %.3lf %.3lf ",u0[0],u0[1],u0[2]);
  printf("  %.3lf %.3lf %.3lf\n",u4[0]-u0[0],u4[1]-u0[1],u4[2]-u0[2]);
  printf("%.3lf %.3lf %.3lf ",u1[0],u1[1],u1[2]);
  printf("  %.3lf %.3lf %.3lf\n",u5[0]-u1[0],u5[1]-u1[1],u5[2]-u1[2]);
  printf("%.3lf %.3lf %.3lf ",u2[0],u2[1],u2[2]);
  printf("  %.3lf %.3lf %.3lf\n",u6[0]-u2[0],u6[1]-u2[1],u6[2]-u2[2]);
  printf("%.3lf %.3lf %.3lf ",u3[0],u3[1],u3[2]);
  printf("  %.3lf %.3lf %.3lf\n",u7[0]-u3[0],u7[1]-u3[1],u7[2]-u3[2]);
};	// drawBox

int main()
{
  TRcoord *tr = new TRcoord();
#define GEOM
#if defined(GEOM)
  double Zb = 20.0;
  double theta = 1.0;
  double phi = 0.0;
  double psi = 0.0;
#else
  double Zb = 20.0;
  double theta = 15.0;
  double phi = 0.0;
  double psi = 30.0;
#endif
  double target[] = {5.0,5.0,15.0};
  tr->set(target,Zb,theta,phi,psi);
  printf("## gnuplot load file\n");
  printf("$S << _EOD\n");
  double u0[] = {0.0,0.0,0.0};
  double u1[] = {0.0,10.0,0.0};
  double u2[] = {10.0,0.0,0.0};
  double u3[] = {10.0,10.0,0.0};
  double u4[] = {0.0,0.0,15.0};
  double u5[] = {0.0,10.0,15.0};
  double u6[] = {10.0,0.0,15.0};
  double u7[] = {10.0,10.0,15.0};
  drawBox(u0,u1,u2,u3, u4,u5,u6,u7);

  double s0[] = {0.0,0.0,1.0};
  double src[3];
  tr->S2U(s0,src);

  printf("##   P  vec\n");
  double q0[] = {1.0,1.0,1.0};
  double q1[] = {-1.0,1.0,1.0};
  double q2[] = {1.0,-1.0,1.0};
  double q3[] = {-1.0,-1.0,1.0};
  double q4[] = {1.0,1.0,5.0};
  double q5[] = {-1.0,1.0,5.0};
  double q6[] = {1.0,-1.0,5.0};
  double q7[] = {-1.0,-1.0,5.0};

  double p0[3],p1[3],p2[3],p3[3];
  double p4[3],p5[3],p6[3],p7[3];
  tr->S2U(q0,p0);
  tr->S2U(q1,p1);
  tr->S2U(q2,p2);
  tr->S2U(q3,p3);
  tr->S2U(q4,p4);
  tr->S2U(q5,p5);
  tr->S2U(q6,p6);
  tr->S2U(q7,p7);
  drawBox(p0,p1,p2,p3, p4,p5,p6,p7);
  printf("_EOD\n");
  /*  line  */
  printf("$L << _EOD\n");
  printf("#  line from source to target\n");
  printf("%.2lf %.2lf %.2lf  ",src[0],src[1],src[2]);
  printf("%.2lf %.2lf %.2lf\n",target[0]-src[0],target[1]-src[1],target[2]-src[2]);
  printf("_EOD\n");
  printf("set title 'Zb=%.1lf theta=%.1lf phi=%.1lf psi=%.1lf\n",Zb,theta,phi,psi);
  printf("set grid\n");
  printf("splot $S using 1:2:3:4:5:6 w vectors title ''");
  printf(",\\\n");
  printf(" $L using 1:2:3:4:5:6 w vectors title 'target line'\n");
  delete tr;
};

