#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define NN 4992
#define NB 39

main()
{
  double rr; // reading element of the array
  int kk,ii,jj,iv,jv; // index 

  double **mm;

  mm= (double **)malloc(NN* sizeof(double *));
  for (ii=0; ii<NN; ii++)
    mm[ii] = (double *)calloc(NN, sizeof(double));
  

  // initialization done 
  FILE *fp;
  char bn[50];

  // diagonal terms   
  for(kk=1;kk<=NB;++kk)
    {      
      sprintf(bn,"dg_%d%d.out",kk,kk);
      fp=fopen(bn,"r");
      printf("%d\t%d\n",kk,kk);
      
      for(ii=0;ii<128;++ii)
	{
	  for(jj=0;jj<128;++jj)
	    {
	      
	      fscanf(fp,"%lf",&rr);
	      iv=128*(kk-1)+ii;
	      jv=128*(kk-1)+jj;
	      mm[iv][jv]=rr;
	    }
	}
      fclose(fp);
    }

  // off diagonal 
  for(kk=1;kk<NB;++kk)
    {      
      sprintf(bn,"offdg_%d%d.out",kk,kk+1);
      fp=fopen(bn,"r");
      
      for(ii=0;ii<128;++ii)
	{
	  for(jj=0;jj<128;++jj)
	    {	      
	      // upper off diagonal 
	      fscanf(fp,"%lf",&rr);
	      iv=128*(kk-1)+ii;
	      jv=128*(kk)+jj;
	      mm[iv][jv]=rr;

	      // lower diagonal 
	      mm[jv][iv]=rr;	      
	    }
	}

      fclose(fp);
    } // storing - done

  // contruction of the total matrix
  fp=fopen("vcorr_matrix.out","w");
  
  for(ii=0;ii<NN;++ii)
    {
      for(jj=0;jj<NN;++jj)
	{	  
	  fprintf(fp,"%e\t",mm[ii][jj]);
	}
      fprintf(fp,"\n");
    }  // done

  fclose(fp);

}
