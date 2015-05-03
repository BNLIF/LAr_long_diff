#include "LAr.h"
#include "TMath.h"

LAr::LAr(){
  T_c = 150.687; // K 
  p_c = 48.630; // bar
  rho_c = 0.5356; // g/cm^3
  
  T_t = 83.8058; // K
  p_t = 0.68891; //bar = 0.1 MPa
  
  T_NBP = 87.3; // K
}

void LAr::print_critical(){
  cout << "Properties of the critical point: " << endl; 
  cout << "    Temperature: " << T_c << " K" << endl;
  cout << "    Pressure:    " << p_c << " bar or " << p_c/10. << " MPa" << endl;
  cout << "    Density:     " << rho_c << " g/cm^3" << endl; 
}

void LAr::print_triple(){
  cout << "Properties of triple point: " << endl;
  cout << "    Temperature: " << T_t << " K" << endl;
  cout << "    Pressure:    " << p_t << " bar or " << p_t/10. << " MPa" << endl;  
}

void LAr::print_boiling(){
  cout << "Normal Boiling Point" << endl;
  cout << "    Temperature: " << T_NBP << " K" << endl; 
}

Double_t LAr::ele_lifetime(Double_t con_ppb, Double_t T, Double_t E, Int_t flag){
  Double_t density = Ldensity(T);
  Double_t Mol = 1e-9 * 1000./39.948*density*con_ppb; // Mole is mole/L
 
  Double_t results;
  Double_t p1, p2, p0, q1,q2,q3;
  if (flag==1){
    //O2
    p0 = 10.5435;
    p1 = 9701.70895;
    p2 = 6020.07535;
    q1 = 958.28688;
    q2 = 415.2558;
    q3 = 1053.92648;
    results = 1e10*ks_f3(E,p0,p1,p2,q1,q2,q3);
  }else if (flag==2){
    //CO2
    p0 = 0.529678;
    p1 = 0.504343;
    p2 = 4.481952;
    q1 = 1.69546;
    results = 1e10*ks_f1(E,p0,p1,p2,q1);
  }else if (flag==3){
    //N2O
    p0 = 1.09709;
    p1 = 6.68216;
    p2 = -0.05016;
    q1 = 0.0276336;
    results = 1e11 * ks_f1(E,p0,p1,p2,q1);
  }else if (flag==4){
    //H2O
    p0 = 1.09709;
    p1 = 6.68216;
    p2 = -0.05016;
    q1 = 0.0276336;
    results = 1.97023*1e11 * ks_f1(E,p0,p1,p2,q1);
  }else if (flag==5){
    //SF6
    p0 = 1.884705;
    p1 = 0.026195;
    p2 = 0.00009571;
    q1 = 0.240995;
    results = 1e14*ks_f1(E,p0,p1,p2,q1);
  }else if (flag==6){
    //N2
    results = 2.66347e7*(1+E)/E/density*vDrift(T,E);
  }

  results = 1./(results * Mol)*1000.;

  return results;
}

Double_t LAr::BoilAr(Double_t T){
  Double_t results;
  Double_t a1 = -5.9410;
  Double_t a2 = 1.35539;
  Double_t a3 = -0.464976;
  Double_t a4 = -1.5399;
  Double_t phi = (1-T/T_c);
  results = p_c * exp(T_c/T*(a1*phi+a2*pow(phi,1.5)
			     +a3*pow(phi,2)+a4*pow(phi,4.5)));
  return results;
}
Double_t LAr::MeltAr(Double_t T){
  Double_t results;
  Double_t a1 = -7476.2665;
  Double_t a2 = 9959.0613;
  results = p_t * (1+a1*(pow(T/T_t,1.05)-1)+a2*(pow(T/T_t,1.275)-1));
  return results;
}
Double_t LAr::SubLimeAr(Double_t T){
  Double_t results;
  Double_t a1 = -11.391604;
  Double_t a2 = -0.39513;
  Double_t phi = (1-T/T_t);
  results = p_t * exp(T_t/T*(a1*phi+a2*pow(phi,2.7)));
  return results;
}


Double_t LAr::ks_f1(Double_t x, Double_t p0, Double_t p1, Double_t p2, Double_t q1){
  Double_t results;
  results = (p0+p1*x+p2*x*x)/(1+q1*x);
  return results;
}
Double_t LAr::ks_f2(Double_t x, Double_t p0, Double_t p1, Double_t p2, Double_t q1, Double_t q2){
  Double_t results;
  results = (p0+p1*x+p2*x*x)/(1+q1*x+q2*x*x);
  return results;
}
Double_t LAr::ks_f3(Double_t x, Double_t p0, Double_t p1, Double_t p2, Double_t q1, Double_t q2, Double_t q3){
  Double_t results;
  results = (p0+p1*x+p2*x*x)/(1+q1*x+q2*x*x+q3*x*x*x);
  return results;
}



Double_t LAr::recombine_Box(Double_t dEodx,Double_t T, Double_t E, Int_t flag){
  Double_t density = Ldensity(T);
  Double_t alpha = 0.93;
  Double_t beta=0.212;
  Double_t results;
  results = log(alpha+beta*dEodx/E/density)/(beta*dEodx/E/density);
  if(flag==2){
    results = 1-0.803*results;
  }
  return results;
}

Double_t LAr::recombine_Birks(Double_t dEodx, Double_t T, Double_t E, Int_t flag){
  Double_t A = 0.8;
  Double_t k = 0.0486;
  Double_t alpha = 0.803;
  Double_t results;
  Double_t density = Ldensity(T);
  
  results = A/(1. + k/E*dEodx/density);
  if (flag ==2){
    results = 1-alpha*results;
  }
  if (E<0.1 ||E>1.0) cout << "Warning: E = " << E << " kV/cm is out of range!" << endl;
  if (dEodx < 1.5 *density || dEodx > 30*density) cout << "Warning: dE/dx = " << dEodx << " MeV/cm is out of range!" << endl; 

  return results;
}

Double_t LAr::fPade(Double_t x, Double_t t0, Double_t p1, Double_t p2, Double_t q1, Double_t q2){
  Double_t results ;
  results = (t0 + p1*x + p2*x*x)/(1+q1*x+q2*x*x);
  return results;
}

Double_t LAr::fPoly(Double_t x, Double_t t0, Double_t p1, Double_t p2, Double_t p3){
  Double_t results = t0 + p1*x+p2*x*x+p3*x*x*x;
  return results;
}

Double_t LAr::Diffusion(Double_t T, Double_t E, Int_t flag){
  Double_t tT = T_t;
  Double_t tC = T_c;

  Double_t kB = 8.6173324e-5;

  Double_t kT0T = kB * tT;
  Double_t kT0C = kB * tC;

  Double_t alpha = 2.2-1.5077*Ldensity(T);
  Double_t results;
  if (flag==1){ //transverse
    results = (1-alpha)*fPoly(E*1000./tT,kT0T,0.0026289636,0.00003204097,-2.28720038e-7) + alpha * fPade(E*1000./tC,kT0C,0.110259,0.0020160892,0.0878054884,0.00020144744);
    results *= 40./56.;
  }else{ //longitduinal
    results = (1-alpha)*fPoly(E*1000./tT,kT0T,0.0012176581,6.322883559e-6,-2.79545739e-8) + alpha * fPade(E*1000./tC,kT0C,0.05057951,0.002007786,0.5480257,0.0055775997);
    results *= 16.5/20.95;
  }

  
  
  return results;
}



Double_t LAr::vD(Double_t T, Double_t E, Int_t flag){
  Double_t p1, p2, p3, p4, p5, p6,t0;
  if (flag==1){
    //WalkowiakParameterSet
    p1 = -0.01481;
    p2 = -0.0075;
    p3 = 0.141;
    p4 = 12.4;
    p5 = 1.627;
    p6 = 0.317;
    t0 = 90.371;
  }else if (flag==2){
    //ICarusFit
    p1 = -0.04640231;
    p2 = 0.0171171;
    p3 = 1.881246;
    p4 = 0.9940772;
    p5 = 0.0117183;
    p6 = 4.202141;
    t0 = 105.7491;
  }else if (flag==3){
    // ICarus + Kalinin
    p1 = -0.0443247;
    p2 = 0.02063811;
    p3 = 1.975295;
    p4 = 1.1205820;
    p5 = 0.00594563;
    p6 = 5.807014;
    t0 = 101.05398861;
  }else if (flag==4){
    //ICarus + Kalinin + Walkowiak
    p1 = -0.0464023;
    p2 = 0.0171171;
    p3 = 1.8812463;
    p4 = 0.9940772;
    p5 = 0.0117183;
    p6 = 4.202141;
    t0 = 105.7491;
  }else if (flag==5){
    //ICarus + Kalini + Walkowiak + Aprile
    p1 = -0.05215565;
    p2 = 0.0228202;
    p3 = 1.9413135;
    p4 = -1.1397996;
    p5 = 0.0068283;
    p6 = 4.9838381;
    t0 = 99.28381;
  }
  Double_t results;
  results = (1+p1*(T-t0)) * (p3*E*log(1+fabs(p4)/E)+p5*pow(E,p6)) + p2 * (T-t0);
  return results;
}

Double_t LAr::vDrift(Double_t T, Double_t E){
  Double_t xFit = 0.0938163 - 0.0052563 *(T-87.302) - 0.000146981 *pow(T-87.302,2);
  Double_t muFit = 5.183987 + 0.01447761 * (T-87.302) - 0.0034972*pow(T-87.302,2)
    -0.0005162374*pow(T-87.302,3);
  Double_t results;
  if (E<xFit){
    results = E*muFit;
  }else if (E<=0.619){
    results = vD(T,E,2);
  }else if (E>=0.699){
    results = vD(T,E,1);
  }else{
    results = (E-0.619)/0.08*vD(T,E,1) + (0.699-E)/0.08*vD(T,E,2);
  }
  return results;
}




Double_t LAr::RRLAr(Double_t lambda, Double_t T){
  Double_t c = 29979245800.;
  Double_t kB = 1.38065e-16;
  Double_t alphaCM = 0.103276;
  Double_t alphaLL = 0.104471;
  
  Double_t omega = 2*3.1415926*c/lambda/1e-7;
  Double_t n = LInrf(lambda,T);
  Double_t epsilon = 1 + 3*(-1.+n)*(1+n)*alphaCM/(alphaCM-n*n*alphaCM+(2+n*n)*alphaLL);

  Double_t results = 1./(pow(omega,4)/6./3.1415926/pow(c,4)
			 *kB * T* IsoCom(T)*pow((-1+epsilon)*(2+epsilon)/3.,2));
  return results;
}

Double_t LAr::GInrf(Double_t lambda){
  Double_t c0 = 1.2055e-2;
  Double_t a1 = 0.2075;
  Double_t b1 = 91.012;
  Double_t a2 = 0.0415;
  Double_t b2 = 87.892;
  Double_t a3 = 4.3330;
  Double_t b3 = 214.02;
  lambda = lambda / 1000.;
  Double_t results = 1 + c0 *(a1/(b1-1./lambda/lambda) 
			      + a2/(b2-1./lambda/lambda)
			      + a3/(b3-1./lambda/lambda));
  return results;
}
Double_t LAr::LInrf(Double_t lambda, Double_t T){
  Double_t nG = GInrf(lambda);
  Double_t rhoG = 1.0034*0.0017840;
  Double_t rhoL = Ldensity(T);
  Double_t results = sqrt((2+nG*nG)*rhoG + 2*(-1+nG*nG)*rhoL)/
    sqrt((2+nG*nG)*rhoG+rhoL-nG*nG*rhoL);

  if (lambda < 110) {
    cout << "Warning: wavelength " << lambda << " nm is out of range!" << endl; 
    return 0;
  }
  return results;
}


Double_t LAr::epsilon(Double_t T){
  Double_t a = 4.12568;
  Double_t rho = Ldensity(T)/39.948;
  Double_t results = (-1.-2.*a*rho)/(-1+a*rho);
  return results;
}

Double_t LAr::Enthalpy(Double_t T){
  Double_t a = 7.98304;
  Double_t b = -0.0481275;
  Double_t c = -0.0047259;
  
  Double_t H = (a+b*T)/(1+c*T);
  if (T>100 || T < 83.8) {
    cout << "Warning: Temperature " << T << " K is out of range" << endl; 
    return 0;
  }
  return H;
}

Double_t LAr::Cp(Double_t T){
  Double_t a = 2.126910;
  Double_t b = -0.02416936;
  Double_t c = 0.00014438771;
  Double_t results = a + b *T + c*T*T;
  return results;
}
Double_t LAr::Cv(Double_t T){
  Double_t a = 0.29934286;
  Double_t b = 8.953781;
  Double_t c = -21.8602116;
  Double_t results = (a+b/T)/(1.0+c/T);
  return results; 
}

Double_t LAr::SpeedofSound(Double_t T){
  Double_t a = 1327.6370386;
  Double_t b = -7.4603194;
  Double_t c = -0.002213558263;
  Double_t results = (a+b*T)/(1.0+c*T);
  return results;
}

Double_t LAr::IsoCom(Double_t T){
  Double_t rcp = Cp(T);
  Double_t rcv = Cv(T);
  Double_t rss = SpeedofSound(T);
  Double_t rhoL = Ldensity(T);

  Double_t results = 0.1 *rcp/rcv/(1000*rhoL*rss*rss);
  return results;
}

Double_t LAr::Viscosity(Double_t T){
  Double_t a = -0.390176214;
  Double_t b = -65.2768756;
  Double_t c = 1.215505389;
  Double_t d = 128.5931277;
  Double_t e = -0.618990054;
  Double_t t = T/83.8058;
  Double_t results = (b/t+d/t/t)/(a+c/t+e/t/t);
  return results;
}

Double_t LAr::vDarIonV(Double_t T, Double_t E){
  E = E * 1000.;
  Double_t results = 0.432 * E * 0.01 / Viscosity(T);
  return results;
}
Double_t LAr::vDarIon(Double_t T, Double_t E){
  E = E * 1000.;
  Double_t a = -1.6024763e-3;
  Double_t b = 1.157022e-5;
  Double_t c = 2.86982898e-7;
  Double_t results = E * 0.01 * (a + b*T + c*T*T);
  return results;
}

Double_t LAr::VPressure(Double_t T){
  Double_t a = -5.9409785;
  Double_t b = 1.3553888;
  Double_t c = -0.46497607;
  Double_t d = -1.5399043;
  Double_t t = (1-T/T_c);
  
  Double_t p = p_c * exp(T_c/T *(a*t
				 + b*pow(t,1.5)
				 + c*pow(t,2)
				 + d*pow(t,4.5)
				 ));
  return p;
}

Double_t LAr::Ldensity(Double_t T){
  Double_t a = 1.5004262;
  Double_t b = -0.31381290;
  Double_t c = 0.086461622;
  Double_t d = -0.041477525;

  Double_t t = (1-T/T_c);
  
  Double_t rho = log(rho_c) + a * pow(t,0.334)  
    + b * pow(t,2./3.) 
    + c * pow(t,7./3.)
    + d * pow(t,4.);
  rho = exp(rho);
  return rho;
}

Double_t LAr::Gdensity(Double_t T){
  Double_t a = -1.70695656;
  Double_t b = -4.02739448;
  Double_t c = 1.55177558;
  Double_t d = -2.30683228;
  
  Double_t t = (1-T/T_c);
  Double_t rho = log(rho_c) + a * pow(t,0.345)
    + b * pow(t,5/6.)
    + c * pow(t,1)
    + d * pow(t,13./3.);
  rho = exp(rho);
  return rho;
  
}


