float rnu(float nu);
void initialize(float omega_m1,float omega_k1,float hh);
float simp(float (*func)(float),int n,float a,float b);
float func_E(float x);
float gm(float x);
float ff(float x);
float cera(float x);
float rnup(float x);
float TFfit_onek(float k, float *tf_baryon, float *tf_cdm);
void TFset_parameters(float omega0hh, float f_baryon, float Tcmb);
float TFzerobaryon(float omega0, float hubble, float Tcmb, float k_hmpc);
float Pk(float kk); // power spectrum 
float sigma_func(float kk); // for calculating sigma 8
/*
int TFmdm_set_cosm(float omega_matter, float omega_baryon, float omega_hdm,int degen_hdm, float omega_lambda, float hubble, float redshift, float spectral);
float TFmdm_onek_mpc(float kk);
float TFmdm_onek_hmpc(float kk);
*/

