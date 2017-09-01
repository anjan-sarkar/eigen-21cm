static float (*funcval_2d)(float,float);
static float xval;;
static float a2val,a3val,b2val,b3val;
static int  n2val,n3val;
float simp(float (*func)(float),int,float,float);

float f2(float y)
{
  return(funcval_2d(xval,y));
}

float f1(float x)
{
  float f2(float);
  xval=x;
  return(simp(f2,n2val,a2val,b2val));
}

float  simp_2d(float (*func_2d)(float,float),int n1,float a1,float b1, int n2,float a2,float b2)
{       
  float f1(float);

  funcval_2d=func_2d;
  n2val=n2;
  a2val=a2;
  b2val=b2;

  return(simp(f1,n1,a1,b1));
}

  
