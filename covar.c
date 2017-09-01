// this programm calculates the visibility correlation 
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "funcs.h"

// global variables
float pi,cc=3*1.e02,betav1,betav2,norm; // cc representing the velocity of light
float U1,U2,b,d,rc,lambdac; // parameters of the central frequency
float delnuc=0.0625; // Delnu_c channel width in MHz

float a1,b1,c1,a2,b2,c2,fc,fc1,fc2,twopibyrc; // parameters for defining the func_2e(float,float,float)
float Tfac,bHI;
// done 

main()
{
  float sinc(float); // for testing

  float hh,omegam0,omegab0,omegak0,ns,sigma_8_present,T;
  int nb; // baseline number
  //int PHASE; // parameters for the input file

  float xHI; // parameters for the power spectrum
  float aac; // scale factor corresponding to central frequency

  FILE *fp; 
  fp=fopen("covar.inp","r");
  fscanf(fp,"%f%f%f%f%f%f",&omegam0,&omegab0,&omegak0,&hh,&sigma_8_present,&ns);
  fscanf(fp,"%d",&nb); // number of baselines for Phase I
  fclose(fp);
  
  initialize(omegam0,omegak0,hh);

  float tcmb=2.728; // CMBR temperature 
  TFset_parameters(omegam0*hh*hh,omegab0/omegam0,tcmb);

  norm=1.;
  norm=simp(sigma_func,100000,0.00001,3.5);
  norm=pow(sigma_8_present,2.)/norm;

  pi=4.0*atan(1.0); // value of pi

  aac=326.5/1420.; // scale factor corresponding to central frequency
  lambdac=0.21/aac; // wavelength in m
  rc=rnu(aac); // distance of observation for the central frequency
  b=30.0;// antenna width
  d=6.*1.92; // antenna length
  // done

  twopibyrc=2*pi/rc;
  bHI=2.0;
  xHI=1.e-3/omegab0;
  Tfac=4.0*(omegab0*hh*hh/0.02)*(0.7/hh)*xHI*bHI;// Tbar*xhI *bHI - z dependent factors are multiplies in function Tf()

  // defining the arguments for vis_cor
  float vis_cor(int,float,float);
 
  float yy[40][256][256]; // for the visibility corelation
  float nn[40][256],ss[40][256];
  int dd; // index for the separation of frequency
  
  // calculating the visibility corelation
  int ii,jj,kk; // two index for frequencies
  float nu1,nu2; // two differnt frequencies

  printf("%e\n",(rnu(318.5/1420)-rnu(334.5/1420))/rnu(318.5/1420));

  // initialization
   for(ii=0;ii<40;++ii)
    {
      for(jj=0;jj<128;++jj)
	{
	  nn[ii][jj]=0.0;
	  ss[ii][jj]=0.0;
	}
    } // done
   
   // visibility correlation 
   for(kk=1;kk<nb+1;++kk) // baseline loop
     {
       for(ii=0;ii<1;++ii) // 1st channel loop
	 {
	   nu1=326.5-delnuc/2.+(ii-63)*delnuc; // defining the first frequency	
	   for(jj=0;jj<1;++jj) // 2nd channel loop
	     {    
	       nu2=326.5-delnuc/2.+(jj-63)*delnuc; // defining the second frequency	       
	       yy[kk][ii][jj]=vis_cor(kk,nu1,nu2); // visibility corelation        
	     }  
	 }       
     } // visibility correlation - done 
   
     
     FILE **file;
     file=malloc(sizeof(FILE *) *50);
     char baseline[50];
     
     for(kk=5;kk<6;++kk)
       {
	 sprintf(baseline,"offdg_%d%d.out",kk,kk+1);
	 file[kk]=fopen(baseline,"w");
	 for(ii=0;ii<128;++ii)
	   {
	     for(jj=0;jj<128;++jj)
	       {		 
		 fprintf(file[kk],"%e\t",yy[kk][ii][jj]); // printing the visibility corelation
	       }	     
	     fprintf(file[kk],"\n");  // line break in file
	   }	 
	 fclose(file[kk]); // closing the file[kk]	 
       }  // done
 
} // end of the main programm


// Visibility Correlation - vis_cor(int,float,float)

float vis_cor(int nn,float nu_1, float nu_2)
{
  
  float func_2e(float qx,float qy,float qz);
  float simp_3d(float((*func_2e)(float,float,float)),int n1,float a1,float b1,int n2,float a2,float b2,int n3,float c1, float c2);
  float FF(float aa);

  // declaration of the variables
  float r1,r2,lr1,lr2,rp1,rp2,lambda1,lambda2; // defining the parameters for the two frequencies
  float nu1,nu2,aa1,aa2,AA;// scale factor for the two frequencies
  float yy; // visbility correlation
  // declaration done
  nu1=nu_1;
  nu2=nu_2;

  aa1=nu1/1420.; // scale factor for nu1
  aa2=nu2/1420.; // scale factor for nu2

  lambda1=cc/nu1; // lambda1
  lambda2=cc/nu2; // lambda2
  r1=rnu(aa1); // distance(in Mpc) for nu_1
  r2=rnu(aa2); // distance(in Mpc) for nu_2
  rp1=rnup(aa1);
  rp2=rnup(aa2);

  lr1=rc/(lambda1*r1);
  lr2=rc/(lambda2*r2);

  a1=b*lr1;
  b1=lambdac*lr1;
  c1=d*lr1;
  a2=b*lr2;
  b2=lambdac*lr2;
  c2=d*lr2;
  fc=twopibyrc*(r1-r2);
  fc1=twopibyrc*rp1*delnuc/2.;
  fc2=twopibyrc*rp2*delnuc/2.;

  betav1=ff(aa1)/bHI;
  betav2=ff(aa2)/bHI;

  U1=nn*d/lambdac; // defining the global variable U1
  U2=(nn)*d/lambdac; // defining the global variable U2
 
  AA=8.*pi*pi*pow(2*1.38/(b*d),2.)*pow(rc,-3.)*FF(aa1)*FF(aa2); // constant multiplication factor

  float u1,u2,v2,w2; // defining the limits of integration
  float fcnu;

  if(nu1<=nu2)
    {
      // off-diagonal case
      u1=nn*c1;
      u2=(nn+1)*c2;
      // end

      //diagonal case
      u1=(nn-1)*c1;
      u2=(nn+1)*c1;
      // end

      v2=a1;
      fcnu=fc2;
    }
  
  if(nu1>nu2)
    {
      if((nn-1)*c1<nn*c2)
	{		
	  u1=nn*c2;
	}		
      if((nn-1)*c1>nn*c2) 
	{
	  u1=(nn-1)*c1;
	}
      
      if((nn+1)*c1<(nn+2)*c2)
	{		
	  u2=(nn+1)*c1;		      
	}
      if((nn+1)*c1>(nn+2)*c2)
	{	      
	  u2=(nn+2)*c2;
	}
      
      v2=a2;
      fcnu=fc1;
    }
  
  // defining the number of steps for qz
  int NN=8,Nstepl=800,Nstep; 
  float sfc=10;

  w2=sfc/fcnu; 
  // qz done

  Nstep=NN*w2/(2.*pi/fc);
  Nstep=(Nstep > Nstepl)? Nstep: Nstepl;
  // Nstep done

  if(u1<u2)
    {
      yy=AA*simp_3d(func_2e,20,u1,u2,20,0.0,v2,Nstep,0.0,w2);
    }
  else
    {
      yy=0;
    }
    return(yy); // returning the value  
}


// functions - integrate
float func_2e(float qx,float qy,float qz) 
{

  float sinc(float);
  float yy;
  float kpar,kk,mu;
  kk=sqrt(qx*qx+qy*qy+qz*qz);

  if(kk>1.e-8)
    {
      mu=qz/kk;            
      yy=Pk(twopibyrc*kk)*(1+betav1*mu*mu)*(1+betav2*mu*mu)*(1-(fabs(U1*b1-qx))/c1)*(1-qy/a1)*(1-(fabs(U2*b2-qx))/c2)*(1-qy/a2)*cos(fc*qz)*sinc(fc1*qz)*sinc(fc2*qz);
    }
  else
    {
      yy=0;
    }
  return(yy);
}

float sinc(float qx) // defining the sinc function
{
  float yy; // value of the sinc(qx)

  // input condition
  if(fabs(qx)<=1.0e-7)
    {
      yy=1.;
    }
  else
    {
      yy=sin(qx)/qx;
    } // done
  return(yy); // return - yy
}

float FF(float aa)
{
  // returns Tbar * xHI * bHI * D(aa)/D(1.0)
  float yy;
  yy=Tfac*pow(aa,-2.)/func_E(aa); // Tbar*xHI*bHI
  yy*=gm(aa)/gm(1.0); // multiply by growing mode 
  return(yy);
} // done
