#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define NN 4992
#define NB 39
#define NR 10000
#define NC 128

// generates one realization of HI visibilities from V2

main()

{
  int ii,jj,kk;
  int rr,pp,qq;
  // double x[100][640],y[100][640],v2[640][640],a,b;
  double a,b;
  double **x,**y,**v2;
  double **vec;
  double ***mm,***ss;

  mm= (double ***)malloc(NB* sizeof(double **));
  for (ii=0; ii<NB; ii++)
    {
      mm[ii] = (double **)calloc(NC, sizeof(double *));
      for(jj=0;jj<NC;++jj)
	{
	  mm[ii][jj]=(double *)calloc(NC, sizeof(double ));
	}
    }
  
  ss= (double ***)malloc(NB* sizeof(double **));
  for (ii=0; ii<NB; ii++)
    {
      ss[ii] = (double **)calloc(NC, sizeof(double *));
      for(jj=0;jj<NC;++jj)
	{
	  ss[ii][jj]=(double *)calloc(NC, sizeof(double ));
	}
    }

  x = (double **)malloc(NR * sizeof(double *));
  for (ii=0; ii<NR; ii++)
    x[ii] = (double *)calloc(NN,sizeof(double));

  y = (double **)malloc(NR * sizeof(double *));
  for (ii=0; ii<NR; ii++)
    y[ii] = (double *)calloc(NN,sizeof(double));
  
  v2 = (double **)malloc(NN * sizeof(double *));
  for (ii=0; ii<NN; ii++)
    v2[ii] = (double *)calloc(NN,sizeof(double));
  
  vec = (double **)malloc(NN * sizeof(double *));
  for (ii=0; ii<NN; ii++)
    vec[ii] = (double *)calloc(NN,sizeof(double));
 
  double xn[NN],yn[NN],vsq[NN],Nsample=1.e1;
  double corr[NN],num[NN];
  double val[NN];
  FILE *fp,*fp1;

  // for eigenvalues 
  gsl_matrix *V2,*evec;
  gsl_eigen_symmv_workspace *w;
  gsl_vector *eval;

  V2=gsl_matrix_alloc(NN,NN);
  evec=gsl_matrix_alloc(NN,NN);
  eval=gsl_vector_alloc(NN);
  w=gsl_eigen_symmv_alloc(NN);
  // done 

  // for random number generation
  gsl_rng *r;
  unsigned long int seed;
  double gr,sigma=1.0;

  r= gsl_rng_alloc(gsl_rng_cmrg);
  seed=24264;
  gsl_rng_set (r, seed);
  // done 

  // read data covariance
  fp=fopen("vcorr_matrix.out","r");
  
  for(ii=0;ii<NN;++ii)
    for(jj=0;jj<NN;++jj)
      {
	fscanf(fp,"%lf",&v2[ii][jj]);   
	gsl_matrix_set(V2,ii,jj,v2[ii][jj]);
    }
  fclose(fp);
  // done 

  // print data covariance 
 fp=fopen("vcorr_matrix.plot","w");
 for(ii=0;ii<NN;++ii)  
   {
     for(jj=0;jj<NN;++jj)
       {
	 fprintf(fp,"%d\t%d\t%le\n",ii,jj,v2[ii][jj]);
       }
     fprintf(fp,"\n");
   }
 fclose(fp);
 // done 

 //find eigenvalues and vectors 
 ii=gsl_eigen_symmv(V2,eval,evec,w);
 printf("%d\n",ii);
 if(ii==0)
   {
     printf("code succesfully executed\n");
   }
 else
   {
     printf("Code failure!!\n");
   } 
 // done
 
 // print eigehvalues and eigenvectors
 fp=fopen("vcorr_matrix.eigval","w"); 
 fp1=fopen("vcorr_matrix.eigvec","w"); 

  for(ii=0;ii<NN;++ii)
    {
      fprintf(fp,"%e\n", gsl_vector_get(eval,ii));

      if(gsl_vector_get(eval,ii)<0.0)
	{
	  printf("eigval negative!!!!\n");
	}

      for(jj=0;jj<NN;++jj)
	fprintf(fp1,"%e\t",gsl_matrix_get(evec,ii,jj));
      fprintf(fp1,"\n");
      
    }  
  fclose(fp);
  fclose(fp1);*/
  // done 
  
 
  fp=fopen("vcorr_matrix.eigval","r"); 
  fp1=fopen("vcorr_matrix.eigvec","r"); 
 
  // read eigenvalues and eigenvectors
  for(ii=0;ii<NN;++ii) // eigenvalue loop
    {
      fscanf(fp,"%lf",&val[ii]);    
      for(jj=0;jj<NN;++jj) // component loop
	{
	 fscanf(fp1,"%lf",&vec[ii][jj]);
	}
    }
  // done
  fclose(fp);
  fclose(fp1); 
  // closing the file

  int dd; 

  //time span
  clock_t start,end;
  double cpu_time_used;
  start=clock();

  char ff[200];

  for(rr=0;rr<NR;++rr) // no of realizations
    {  

      sprintf(ff,"vis.%d.out",rr+1);
      fp=fopen(ff,"w");   

      for(pp=0;pp<NN;++pp) // eigenvalue loop
	{
	  a=pow(val[pp]/2.,0.5)*gsl_ran_gaussian(r,sigma);
	  b=pow(val[pp]/2.,0.5)*gsl_ran_gaussian(r,sigma);

	  for(qq=0;qq<NN;++qq) // component loop
	    {
	      x[rr][qq]=x[rr][qq]+a*vec[qq][pp];
	      y[rr][qq]=y[rr][qq]+b*vec[qq][pp];
	    }
	}

      for(ii=0;ii<NB-1;++ii)
	{
	  for(jj=0;jj<NC;++jj)
	    {
	      for(kk=0;kk<NC;++kk)
		{
		  // visibility covariance - diag		  
		  mm[ii][jj][kk]=mm[ii][jj][kk]+(x[rr][NC*ii+jj]*x[rr][NC*ii+kk]+y[rr][NC*ii+jj]*y[rr][NC*ii+kk]);
		  ss[ii][jj][kk]=ss[ii][jj][kk]+pow(x[rr][NC*ii+jj]*x[rr][NC*ii+kk]+y[rr][NC*ii+jj]*y[rr][NC*ii+kk],2.);
		  // done

		  // visibility covariance - off-diag
		  mm[ii][jj][kk]=mm[ii][jj][kk]+(x[rr][NC*ii+jj]*x[rr][NC*(ii+1)+kk]+y[rr][NC*ii+jj]*y[rr][NC*(ii+1)+kk]);
		  ss[ii][jj][kk]=ss[ii][jj][kk]+pow(x[rr][NC*ii+jj]*x[rr][NC*(ii+1)+kk]+y[rr][NC*ii+jj]*y[rr][NC*(ii+1)+kk],2.);
		  // done

		}
	    }
	}

      // print the visibility realizations
      for(ii=0;ii<NC;++ii)
	{
	  for(jj=0;jj<NB;++jj)
	    {
	      fprintf(fp,"%e\t%d\t%e\t%d\t",x[rr][128*jj+ii],ii+1,y[rr][128*jj+ii],ii+1);
	    }
	  fprintf(fp,"\n");
	}
      fclose(fp);
    } // done

  end=clock();
  cpu_time_used=((double)(end-start))/CLOCKS_PER_SEC;
  printf("Total time taken = %.2lf sec\n",cpu_time_used);

  // visibility covariance matrix - store  
  char mn[200];
  char vn[200];
  
  for(ii=0;ii<NB;++ii)
    {
      sprintf(mn,"dg.avg_%d%d.out",ii+1,ii+1); // for visibility covariance - diag
      fp=fopen(mn,"w");
      
      sprintf(vn,"dg.var_%d%d.out",ii+1,ii+1); // for variance of the visibility covariance - diag
      fp1=fopen(vn,"w");
      
      for(jj=0;jj<NC;++jj)
      	{
	  for(kk=0;kk<NC;++kk)
	    {
	      ss[ii][jj][kk]=ss[ii][jj][kk]/NR;
	      mm[ii][jj][kk]=mm[ii][jj][kk]/NR;
	      
	      fprintf(fp,"%e\t",mm[ii][jj][kk]);
	      fprintf(fp1,"%e\t",sqrt(ss[ii][jj][kk]-pow(mm[ii][jj][kk],2.)));	      
	    }	  
	  fprintf(fp,"\n"); // line break in the file -- matrix construction
	  fprintf(fp1,"\n"); // line break in the file -- matrix construction
	} 
      fclose(fp);
      fclose(fp1);  // close the files  
    } // done
 
  // visibilities - print  
  char bn[200];
  for(kk=0;kk<NR;++kk)
    {
     sprintf(bn,"vcorr_%d.out",kk+1);
     fp=fopen(bn,"w");
     for(ii=0;ii<NC;++ii)
       {
	 for(jj=0;jj<NB;++jj)
	   {
	     fprintf(fp,"%e\t%e\t",x[kk][NC*jj+ii],y[kk][NC*jj+ii]);
	   }
	 fprintf(fp,"\n");
       }   
     fclose(fp); // close the file     
    } // done 
  
  gsl_eigen_symmv_free(w);
  gsl_rng_free (r);
  
} // end of the programm

 
