const double hbar_c = 197.327053;
const double Nc = 1.;
const double Nq = 25.83;
const double Nm = 1.714;

const double alpha_c = 5.9;
const double alpha_q = 3.1;
const double alpha_m = 3.78;
const double beta_c = -5.2;
const double beta_q = -2.1;
const double beta_m = -2.78;
const double gamma_c = 13.9;
const double gamma_q = 7.2;
const double gamma_m = 11.4;
const double delta_c = 0.96;
const double delta_q = 1.6;
const double delta_m = 1.07;

//2nd parametarization
/*
const double alpha_c = 5.0;
const double alpha_q = 1.4;
const double alpha_m = 3.78;
const double beta_c = -4.5;
const double beta_q = -0.1;
const double beta_m = -2.78;
const double gamma_c = 11.5;
const double gamma_q = 7.7;
const double gamma_m = 11.4;
const double delta_c = 1.11;
const double delta_q = 1.7;
const double delta_m = 1.07;
*/

const double mass_omega = 782.65/1000.;//GeV
const double mass_phi = 1019.456/1000.;
const double mass_deuteron = 1875.612/1000.;//GeV
double func_Gcd(double Q2){//GeV
  double Gcd,gc,Fc;
  Fc = 1.- alpha_c - beta_c + alpha_c*(mass_omega*mass_omega/((mass_omega*mass_omega)+Q2)) +  beta_c*(mass_phi*mass_phi/((mass_phi*mass_phi)+Q2));
  gc = 1./pow(1.+(gamma_c*Q2),delta_c);
  Gcd = Nc*gc*Fc;

  return Gcd;
}
  double func_Gqd(double Q2){
  double Gqd,gq,Fq;
  Fq = 1.- alpha_q - beta_q + alpha_q*(mass_omega*mass_omega/((mass_omega*mass_omega)+Q2)) +  beta_q*(mass_phi*mass_phi/((mass_phi*mass_phi)+Q2));
  gq = 1./pow(1.+(gamma_q*Q2),delta_q);
  Gqd = Nq*gq*Fq;

  return Gqd;
}

double func_Gmd(double Q2){
  double Gmd,gm,Fm;
  Fm = 1.- alpha_m - beta_m + alpha_m*(mass_omega*mass_omega/((mass_omega*mass_omega)+Q2)) +  beta_m*(mass_phi*mass_phi/((mass_phi*mass_phi)+Q2));
  gm = 1./pow(1.+(gamma_m*Q2),delta_m);
  Gmd = Nm*gm*Fm;

  return Gmd;
}

double func_error_Gcd(double Q2,double error_Q2){
  double roundF_roundQ2 = -(alpha_c*mass_omega*mass_omega/pow(mass_omega*mass_omega+Q2,2))-(beta_c*mass_phi*mass_phi/pow(mass_phi*mass_phi+Q2,2));
  double Delta_Fcd = sqrt(pow(roundF_roundQ2*error_Q2,2));
  double roundgc_roundQ2 = -delta_c/pow(1+(gamma_c*Q2),delta_c-1);
  double Delta_gc = sqrt(pow(roundgc_roundQ2*error_Q2,2));
  double gc = 1./pow(1.+(gamma_c*Q2),delta_c);
  double Fc = 1.- alpha_c - beta_c + alpha_c*(mass_omega*mass_omega/((mass_omega*mass_omega)+Q2)) +  beta_c*(mass_phi*mass_phi/((mass_phi*mass_phi)+Q2));
  double error_Gcd = Nc*sqrt((Fc*Fc*Delta_gc*Delta_gc) + (gc*gc*Delta_Fcd*Delta_Fcd));
  return error_Gcd;
}

double func_error_Gmd(double Q2,double error_Q2){
  double roundF_roundQ2 = -(alpha_m*mass_omega*mass_omega/pow(mass_omega*mass_omega+Q2,2))-(beta_m*mass_phi*mass_phi/pow(mass_phi*mass_phi+Q2,2));
  double Delta_Fmd = sqrt(pow(roundF_roundQ2*error_Q2,2));
  double roundgm_roundQ2 = -delta_m/pow(1+(gamma_m*Q2),delta_m-1);
  double Delta_gm = sqrt(pow(roundgm_roundQ2*error_Q2,2));
  double gm = 1./pow(1.+(gamma_m*Q2),delta_m);
  double Fm = 1.- alpha_m - beta_m + alpha_m*(mass_omega*mass_omega/((mass_omega*mass_omega)+Q2)) +  beta_m*(mass_phi*mass_phi/((mass_phi*mass_phi)+Q2));
  double error_Gmd = Nm*sqrt((Fm*Fm*Delta_gm*Delta_gm) + (gm*gm*Delta_Fmd*Delta_Fmd));
  return error_Gmd;
}


//for Rosenbluth of deuteron
double func_alpha_thetad(double theta ,double Q2){//theta[deg] //Q2[GeV^2]
  double eta = Q2/(4.*mass_deuteron*mass_deuteron);
  double theta_rad = theta/2.*pi/180.;
  double alpha_D = 2./3.*eta + (4./3.*eta*(1.+eta)*pow(tan(theta_rad),2.));
  return alpha_D;
}

double func_error_alphaD(double theta, double theta_error, double Q2, double Q2_error){//Q2[GeV^2]
double theta_rad = theta*pi/180.;
double theta_error_rad = theta_error*pi/180.;
double roundalpha_roundQ2 = (Q2/6./pow(mass_deuteron,4)+(1./3./pow(mass_deuteron,2)))*pow(tan(theta_rad/2.),2)+(1./6./pow(mass_deuteron,2));
double roundalpha_roundtheta = ((Q2*Q2/12./pow(mass_deuteron,4)) + (Q2/3./pow(mass_deuteron,2)))*tan(theta_rad/2.)/pow(cos(theta_rad/2.),2);
double delta_alpha = sqrt(pow(roundalpha_roundQ2*Q2_error,2) + pow(roundalpha_roundtheta*theta_error_rad,2));
return delta_alpha;
}

double func_error_crosssection_ratio_model(double Gcd2,double Gcd2_error, double Gmd2, double Gmd2_error, double alpha, double alpha_error){
  double crosssection_ratio = Gcd2 + alpha*Gmd2;
  double delta_crosssection_ratio = sqrt(pow(Gcd2_error,2) + pow(alpha*Gmd2_error,2) + pow(Gmd2*alpha_error,2));
  return delta_crosssection_ratio;
}

double func_A_D(double Q2, double Gcd, double Gmd, double Gqd){
  double eta = Q2/(4.*mass_deuteron*mass_deuteron);
  double A_D = Gcd*Gcd + 2./3.*eta*Gmd*Gmd + 8./9.*eta*eta*Gqd*Gqd;
  return A_D;
}

double func_B_D(double Q2, double Gmd){
  double eta = Q2/(4.*mass_deuteron*mass_deuteron);
  double B_D = 4./3.*eta*(1 + eta)*Gmd*Gmd;
  return B_D;
}

double func_Crosssection_deuteron(double Mott,double theta, double A_D, double B_D){
double theta_rad = theta*pi/180.;
double crosssection = Mott*(A_D + B_D*pow(tan(theta_rad/2.),2));
return crosssection;
}
