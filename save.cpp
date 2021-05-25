#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <valarray>
#include <cmath>
#include <complex> // Pour les nombres complexes
#include "ConfigFile.h"

using namespace std;
typedef vector<complex<double> > vec_cmplx;

// Fonction résolvant le système d'equations A * solution = rhs
// où A est une matrice tridiagonale
template <class T> void triangular_solve(vector<T> const& diag,
                                         vector<T> const& lower,
                                         vector<T> const& upper,
                                         vector<T> const& rhs,
                                         vector<T>& solution){
  vector<T> new_diag = diag;
  vector<T> new_rhs = rhs;

  // Forward elimination
  for(unsigned int i(1); i<diag.size(); ++i){
    T pivot = lower[i-1] / new_diag[i-1];
    new_diag[i] -= pivot * upper[i-1];
    new_rhs[i] -= pivot * new_rhs[i-1];
  }

  solution.resize(diag.size());

  // Solve last equation
  solution[diag.size()-1] = new_rhs[diag.size()-1] / new_diag[diag.size()-1];

  // Backward substitution
  for(int i = int(diag.size()) - 2; i >= 0; --i){
    solution[i] = (new_rhs[i] - upper[i] * solution[i+1]) / new_diag[i];
  }
}


// Potentiel V(x) : // TODO
double V_fct(double const& x, double const& omega2, double const& Delta){
  return min(0.5*omega2*pow(x - Delta, 2), 0.5*omega2*pow(x + Delta, 2));
}

// Déclaration des diagnostics de la particule d'après sa fonction d'onde psi :
//  - prob calcule la probabilité de trouver la particule entre les points de maillage nL et nR,
//  - E calcule son énergie moyenne,
//  - xmoy calcule sa position moyenne,
//  - x2moy calcule sa position au carré moyenne,
//  - pmoy calcule sa quantité de mouvement moyenne,
//  - p2moy calcule sa quantité de mouvement au carré moyenne.
double prob(vec_cmplx const& psi, int nL, int nR, double dx);
double E(vec_cmplx const& psi, vec_cmplx const& diagH, vec_cmplx const& lowerH, vec_cmplx const& upperH, double const& dx);
double xmoy(vec_cmplx const& psi, const vector<double>& x, double const& dx);
double x2moy(vec_cmplx const& psi, const vector<double>& x, double const& dx);
double pmoy(vec_cmplx const& psi, double const& dx);
double p2moy(vec_cmplx const& psi, double const& dx);

// Fonction pour normaliser une fonction d'onde :
vec_cmplx normalize(vec_cmplx const& psi, double const& dx);

// Les définitions de ces fonctions sont en dessous du main.


int main(int argc,char **argv){
  complex<double> complex_i = complex<double> (0,1); // Nombre imaginaire i

  string inputPath("configuration.in"); // Fichier d'input par défaut
  if(argc>1) // Fichier d'input specifié par l'utilisateur ("./Exercice8 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les paramètres sont lus et stockés dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complémentaires ("./Exercice8 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Paramètres physiques :
  double tfin           = configFile.get<double>("tfin");
  double xL             = configFile.get<double>("xL");
  double xR             = configFile.get<double>("xR");
  double xda             = configFile.get<double>("xda");
  double xdb             = configFile.get<double>("xdb");
  double hbar           = configFile.get<double>("hbar");;
  double m              = configFile.get<double>("mass");
  double omega          = configFile.get<double>("omega");
  double Delta          = configFile.get<double>("Delta");
  double x0             = configFile.get<double>("x0");
  double k0             = 2. * M_PI * double(configFile.get<int>("n")) / (xR-xL);
  double sigma0         = configFile.get<double>("sigma_norm") * (xR-xL);
  double t_detect       = configFile.get<double>("t_detect");
  
  double omega2 = m*omega*omega;

  // Paramètres numeriques :
  double dt      = configFile.get<double>("dt");
  int Ninters    = configFile.get<int>("Ninters");
  int Npoints    = Ninters + 1;
  double dx      = (xR-xL) / Ninters;

  // Maillage :
  vector<double> x(Npoints);
  for(int i(0); i<Npoints; ++i)
    x[i] = xL + i*dx;

  vector<double> V(Npoints); // Vecteur contenant le potentiel
  for(int i(0); i<Npoints; ++i){
    V[i] = V_fct(x[i], pow(omega, 2), Delta);
  }

  // Initialisation de la fonction d'onde :
  vec_cmplx psi(Npoints);
  // TODO: initialiser le paquet d'onde, équation (4.116) du cours
  for(int i(1); i<Npoints-1; ++i){
    psi[i] = exp(complex_i*k0*x[i])*exp(-pow(x[i]-x0, 2)/(2.*pow(sigma0,2)));
  } // MODIFIER
  // Modifications des valeurs aux bords :
  psi[0] = complex<double> (0.,0.);
  psi[Npoints-1] = complex<double> (0.,0.);
  // Normalisation :
  psi = normalize(psi, dx);

  // Matrices (d: diagonale, a: sous-diagonale, c: sur-diagonale) :
  vec_cmplx dH(Npoints), aH(Ninters), cH(Ninters); // matrice Hamiltonienne
  vec_cmplx dA(Npoints), aA(Ninters), cA(Ninters); // matrice du membre de gauche de l'équation (4.99)
  vec_cmplx dB(Npoints), aB(Ninters), cB(Ninters); // matrice du membre de droite de l'équation (4.99)

  complex<double> a = complex_i*hbar*dt/(4*m*pow(dx,2)); // Coefficient complexe a, Eq.(4.100)
  vec_cmplx b(Npoints); // Coefficient complexe b, Eq.(4.100)
  for(int i(0); i<Npoints; ++i){
    b[i] = complex_i*0.5*dt*V[i]/hbar;
  }

  // TODO: calculer les éléments des matrices A, B et H.
  // Ces matrices sont stockées sous forme tridiagonale, d:diagonale, c et a: diagonales supérieures et inférieures
  for(int i(0); i<Npoints; ++i){ // Boucle sur les points de maillage
    dH[i] = complex<double> (pow(hbar,2)/(pow(dx,2)*m) + V[i], 0.); // MODIFIER
    dA[i] = 1. + 2.*a + b[i]; // MODIFIER
    dB[i] = 1. - 2.*a - b[i]; // MODIFIER
  }
  for(int i(0); i<Ninters; ++i){ // Boucle sur les intervalles
    aH[i] = cH[i] = complex<double> (-pow(hbar,2)/(2.*pow(dx,2)*m), 0.); // MODIFIER
    aA[i] = cA[i] = -a; // MODIFIER
    aB[i] = cB[i] = a; // MODIFIER
  }

  // Conditions aux limites: psi nulle aux deux bords
  // TODO: Modifier les matrices A et B pour satisfaire les conditions aux limites
  // INSERER ICI
  cA[0] = aA[0] = cB[0] = aB[0] = cA[Npoints-2] = aA[Npoints-2] = cB[Npoints-2] = aB[Npoints-2] = 0.;
  dA[0] = dB[0] = dA[Npoints-1] = dB[Npoints-1] = 1.;

  // Fichiers de sortie :
  ofstream fichier_potentiel((configFile.get<string>("output_potential")).c_str());
  fichier_potentiel.precision(15);
  for(int i(0); i<Npoints; ++i)
    fichier_potentiel << x[i] << " " << V[i] << endl;
  fichier_potentiel.close();

  ofstream fichier_psi((configFile.get<string>("output_squared_wave")).c_str());
  fichier_psi.precision(15);

  ofstream fichier_observables((configFile.get<string>("output_observables")).c_str());
  fichier_observables.precision(15);

  // Écrire position en x
  fichier_psi << 0.0 << " ";
  for(int i(0); i<Npoints-1; ++i){
    fichier_psi << x[i] << " ";
  }
  fichier_psi << x[Npoints-1] << endl;

  // Boucle temporelle :
  valarray<double> print_array=valarray<double>(0.e0,Npoints+1);
  double t,window;
  for(t=0.; t+dt/2.<tfin; t+=dt){
    // Détection de la particule
    if(round(t/dt) == round(t_detect/dt)){
      for(int i(0); i<abs(Ninters*xL/(xL-xR)); ++i){
        psi[i] = complex<double> (0.,0.);
      }
      for(int i(abs(Ninters*xL/(xL-xR))); i<abs(Ninters*(xda-xL)/(xL-xR)); ++i){
        window = pow(sin(0.5*M_PI*x[i]/xda),2.e0);
        psi[i] = polar(window*abs(psi[i]),arg(psi[i]));
      }
      for(int i(abs(Ninters*(xdb-xL)/(xL-xR))); i<Ninters; ++i){
        window = pow(0.5*M_PI*x[i]/(xR-xdb),2.0);
        psi[i] = polar(window*abs(psi[i]),arg(psi[i]));
      }
      psi = normalize(psi, dx); // Normalise psi pour que la proba totale soit 1
    }

    // Écriture de |psi|^2 :
    print_array[0] = t;
    for(int i(1); i<Npoints+1; ++i){
      print_array[i] = norm(psi[i-1]); // La fonction C++ norm prend le module au carré
    }
    for(int i(0); i<Npoints; ++i){
      fichier_psi << print_array[i] << " ";
    }
    fichier_psi << print_array[Npoints] << endl;

    // Écriture des observables :
    fichier_observables << t << " "                                                   // Temps
                        << prob(psi,0,abs(Ninters*xL/(xL-xR)),dx) << " "              // Probabilité que la particule soit en x < 0
                        << prob(psi,abs(Ninters*xL/(xL-xR)),Ninters,dx) << " "        // Probabilité que la particule soit en x > 0
                        << prob(psi,0,Ninters,dx) << " " 		   	                      // Probabilité totale
                        << E(psi,dH,aH,cH,dx) << " "                       	          // Énergie
                        << xmoy(psi,x,dx) << " "                           	          // Position moyenne
                        << x2moy(psi,x,dx) << " "                          	          // Position^2 moyenne
                        << pmoy(psi,dx) << " "                             	          // Quantité de mouvement moyenne
                        << p2moy(psi,dx) << " "                           	          // (Quantité de mouvement)^2 moyenne
                        << sqrt(x2moy(psi,x,dx)-xmoy(psi,x,dx)*xmoy(psi,x,dx))*\
			                  sqrt(p2moy(psi,dx)-pmoy(psi,dx)*pmoy(psi,dx)) << " "          // Heisenberg index          
                        << sqrt(x2moy(psi,x,dx)-xmoy(psi,x,dx)*xmoy(psi,x,dx)) << " " // Incertitude en x
                        << sqrt(p2moy(psi,dx)-pmoy(psi,dx)*pmoy(psi,dx)) << endl;     // Incertitude en p

    // Calcul du membre de droite :
    vec_cmplx psi_tmp(Npoints,0.);

    // Multiplication psi_tmp = B * psi :
    for(int i(0); i<Npoints; ++i)
      psi_tmp[i] = dB[i] * psi[i];
    for(int i(0); i<Ninters; ++i){
      psi_tmp[i] += cB[i] * psi[i+1];
      psi_tmp[i+1] += aB[i] * psi[i];
    }

    // Résolution de A * psi = psi_tmp :
    triangular_solve(dA, aA, cA, psi_tmp, psi);
  } // Fin de la boucle temporelle

    // Écrire |psi|^2
    print_array[0] = t;
    for(int i(0); i<Npoints; ++i){
      print_array[i+1] = norm(psi[i]);
    }
    for(int i(0); i<Npoints; ++i){
      fichier_psi << print_array[i] << " ";
    }
    fichier_psi << print_array[Npoints] << endl;

    // Écriture des observables :
    fichier_observables << t << " "                                                   // Temps
                        << prob(psi,0,abs(Ninters*xL/(xL-xR)),dx) << " "              // Probabilité que la particule soit en x < 0
                        << prob(psi,abs(Ninters*xL/(xL-xR)),Ninters,dx) << " "        // Probabilité que la particule soit en x > 0
                        << prob(psi,0,Ninters,dx) << " " 		   	                      // Probabilité totale
                        << E(psi,dH,aH,cH,dx) << " "                       	          // Énergie
                        << xmoy(psi,x,dx) << " "                           	          // Position moyenne
                        << x2moy(psi,x,dx) << " "                          	          // Position^2 moyenne
                        << pmoy(psi,dx) << " "                             	          // Quantité de mouvement moyenne
                        << p2moy(psi,dx) << " "                           	          // (Quantité de mouvement)^2 moyenne
                        << sqrt(x2moy(psi,x,dx)-xmoy(psi,x,dx)*xmoy(psi,x,dx))*\
			                  sqrt(p2moy(psi,dx)-pmoy(psi,dx)*pmoy(psi,dx)) << " "          // Heisenberg index
                        << sqrt(x2moy(psi,x,dx)-xmoy(psi,x,dx)*xmoy(psi,x,dx)) << " " // Incertitude en x
                        << sqrt(p2moy(psi,dx)-pmoy(psi,dx)*pmoy(psi,dx)) << endl;     // Incertitude en p

  fichier_observables.close();
  fichier_psi.close();
}


double prob(vec_cmplx const& psi, int nL, int nR, double dx){
  //TODO: calculer la probabilité de trouver la particule entre les points nL et nR
  //Utiliser la formule des trapèzes pour l'intégration
  double resultat(0.);
  for(int i(nL); i<nR; ++i) {
    resultat += (norm(psi[i]) + norm(psi[i+1]))*0.5*dx;
  }
  return resultat;
}

//TODO: Calculer les valeurs moyennes des observables E, x, p, x^2, p^2
//Utiliser la formule des trapèzes pour l'intégration
double E(vec_cmplx const& psi, vec_cmplx const& diagH, vec_cmplx const& lowerH, vec_cmplx const& upperH, double const& dx){
  vec_cmplx psi_tmp(psi.size()); // Vecteur pour stocker H*psi
  complex<double> resultat(0,0); // Initialiser

  // H(psi): produit de la matrice H et du vecteur psi
  for(size_t i(0); i<psi.size(); ++i){
    psi_tmp[i] = diagH[i]*psi[i]; 
  }
  for(size_t i(0); i<psi.size()-1; ++i){
    psi_tmp[i] += upperH[i]*psi[i+1];
    psi_tmp[i+1] += lowerH[i]*psi[i];
  }

  // Intégrale de psi*H(psi) dx
  for(size_t i(0); i<psi.size()-1; ++i){
    resultat += (conj(psi[i+1])*psi_tmp[i+1] + conj(psi[i])*psi_tmp[i])*0.5*dx;
  }

  return real(resultat);
}


double xmoy(vec_cmplx const& psi, const vector<double>& x, double const& dx){
  double resultat(0.);
  for(size_t i(0); i<psi.size()-1; ++i) {
    resultat += (norm(psi[i])*x[i] + norm(psi[i+1])*x[i+1])*0.5*dx;
  }
  return resultat;
}

double x2moy(vec_cmplx const& psi, const vector<double>& x, double const& dx){
  double resultat(0.);
  for(size_t i(0); i<psi.size()-1; ++i) {
    resultat += norm(psi[i])*pow(x[i],2) + norm(psi[i+1])*pow(x[i+1],2)*0.5*dx;
  }
  return resultat;
}


double pmoy(vec_cmplx const& psi, double const& dx){
  complex<double> complex_i = complex<double> (0,1); // Nombre imaginaire i
  unsigned int N(psi.size()-1);
  // Utiliser la définition de p = -i hbar d/dx
  // Utiliser les différences finies centrées pour d/dx
  // Utiliser la formule des trapèzes pour l'intégration sur x
  double resultat(0.);
  double hbar(1.0);

  for(size_t i(1); i<psi.size()-2; ++i) {
    resultat -= real(0.5*hbar*complex_i*(0.5*conj(psi[i])*(psi[i+1] - psi[i-1]) + 0.5*conj(psi[i+1])*(psi[i+2] - psi[i])));
  }
  resultat -= real(0.5*hbar*complex_i*(conj(psi[0])*(psi[1] - psi[0]) + 0.5*conj(psi[1])*(psi[2] - psi[0])));  // Premier élément
  resultat -= real(0.5*hbar*complex_i*(0.5*conj(psi[N-1])*(psi[N] - psi[N-2]) + conj(psi[N])*(psi[N] - psi[N-1])));  // Dernier élément

  return resultat;
}


double p2moy(vec_cmplx const& psi, double const& dx){
  double resultat(0.);
  // Utiliser la définition de p^2 = hbar^2 d^2/dx2
  // Utiliser les différences finies centrées pour d^2/dx^2
  // Utiliser la formule des trapèzes pour l'intégration sur x
  // Ignorer la contribution du premier et du dernier point de maillage
  complex<double> complex_i = complex<double> (0,1); // Nombre imaginaire i
  double hbar(1.0);

  for(size_t i(1); i<psi.size()-2; ++i) {
    resultat -= real((pow(hbar,2)/(2.*dx))*conj(psi[i])*(psi[i+1] - 2.*psi[i] + psi[i-1]) + conj(psi[i+1])*(psi[i+2] - 2.*psi[i+1] + psi[i]));
  }

  return resultat;
}

vec_cmplx normalize(vec_cmplx const& psi, double const& dx){
  vec_cmplx psi_norm(psi.size());
  double norm = sqrt(prob(psi,0,psi.size()-1,dx));
  for(unsigned int i(0); i<psi.size(); ++i)
    psi_norm[i] = psi[i]/norm;
  return psi_norm;
}