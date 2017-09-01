#include<math.h>
#include<stdio.h>
#include"funcs.h"

extern float RR;

float omega_m,omega_l,omega_k,vhh,cbyho,LL;
int curv;

// Omega_m0, Omega_k0, H0/100 are needed as input
void initialize(float omega_m1,float omega_k1,float hh)
{
  vhh=hh;
  omega_m=omega_m1;
  omega_k=omega_k1;
  omega_l=1.-omega_m-omega_k;
  cbyho=3.e3/hh;    /*coH0 =c/H_0 in terms of Mpc*/

  curv=0; // - 1 if open 
  if(fabs(omega_k)>1.e-8)
    {
      LL=cbyho/sqrt(fabs(omega_k)); // radius of curvature 

      if(omega_k>0.) curv=-1; // open
      else curv=1; // closed
    }
}

/*----------------Returns E= H/H_0------------*/ 

float func_E(float x)
{
  float tt;
  tt=omega_m*powf(x,-3)+ omega_k*powf(x,-2.) + omega_l;
  return(sqrt(tt));
}

// used for ceta

float func( float x)
{
  return(powf(func_E(x)*x*x,-1.));
}

//calculate ceta for given scale factor

float ceta(float x)
{
  float func(float);
  float y;
  y=simp(func,10000,x,1.);
  y=y*cbyho;       
  return(y);
}

// which is conformal distance r for given scale factor  
// r=d_A *(1+z) where d_A is angular diamter distance

float rnu(float x)
{
  float ceta(float x);
    
    if(curv==0)
      {
	return(ceta(x));
      }
    else
      {
	if(curv==-1)	return(LL*sinh(ceta(x)/LL));
	else	return(LL*sin(ceta(x)/LL));

      }
}

// nu derivative of r at given scale factor

float rnup(float x)
{
  float term;

  term=cbyho/(x*x*func_E(x)*1420.);

  if(curv==-1)
    {
      term=term*cosh(ceta(x)/LL);
    }

  if(curv==1)
    {
      term=term*cos(ceta(x)/LL);
    }
  return(term);
}




/*gm(float x) returns growing mode of matter perturbation at scale factor x.This is not normalized*/

float gm(float x)
{ 
  float func_gm(float x);
  float ff;
  ff=func_E(x)*simp(func_gm,1000,1.e-6,x);
  return(ff);
} 
float func_gm(float x)
{
  float ss;
  ss=x*func_E(x);
  return(powf(ss,-3));
}

// dlogD/dloga
float ff(float x)
{
  float Ev,dEv,y;
  
  Ev=func_E(x);
  dEv=-0.5*(3.*omega_m*pow(x,-4.)+2.*omega_k*pow(x,-3))/Ev;
  y=x*dEv/Ev+1./(pow(x*Ev,2.)*gm(x));
  return(y);
}
