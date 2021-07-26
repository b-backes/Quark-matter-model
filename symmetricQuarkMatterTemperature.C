/*

Equation of State solver for the DDQM model with fixed finite temperature, considering symmetric quark matter.
By Betânia Camille Backes, 2020
Should be executed with ROOT, see https://root.cern/
If the software is installed, open ROOT, compile and execute with
.X symmetricQuarkMatterTemperature.C

*/

#include <stdio.h>
#include <math.h>
using namespace TMath;

const double uCurrentMass = 5.0; // quarks masses in MeV
const double dCurrentMass = 10.0;
const double sCurrentMass = 80.0;

const double hbarc = 197.33;  // factors used to convert from orders of MeV to orders of fm^-1
const double hbarc2 = 197.33*197.33;
const double hbarc3 = 197.33*197.33*197.33;

const double piSquared = M_PI*M_PI; // for simplicity

const double kf0 = 0.0; // lower limit used while solving the integrals
const double kfmax = 1.e3; // higher limit used while solving the integrals

#define epsilon 1e-7 // used in the bissection method below
#define delta 1e-5 // used in the numerical differentiaion, the result will usually be better than with a smaller value

//**********************************************************//

// The particle densities
double density(double mass, double chemPot, double temperature, double degeneracy){

    TF1 particles("particles","1./(1+exp((sqrt(x*x + [0]) - [1])/[2]))",kf0,kfmax);
    particles.SetParameter(0,mass*mass);
    particles.SetParameter(1,chemPot);
    particles.SetParameter(2,temperature);

    TF1 antiparticles("antiparticles","1./(1+exp((sqrt(x*x + [0]) + [1])/[2]))",kf0,kfmax);
    antiparticles.SetParameter(0,mass*mass);
    antiparticles.SetParameter(1,chemPot);
    antiparticles.SetParameter(2,temperature);

    TF1 densityIntegral("densityIntegral","x*x*(particles - antiparticles)",kf0,kfmax);

    double result = densityIntegral.Integral(kf0,kfmax)*degeneracy/(2*piSquared);

    return result;
}

// The stability conditions for the chemical potentials and densities, rewritten so that the chemical potentials
// can be obtained utilising a method for finding the zeroes of a function. See the Appendix in https://arxiv.org/abs/2007.04494 

double findChemPot(double chemPot, double baryonDensity, double C, double D, double temperature, double degeneracy){

    double massU = uCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);
    double massD = dCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);
    double massS = sCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);

    double baryonNumberConservation = (density(massU,chemPot,temperature,degeneracy) + 
            density(massD,chemPot,temperature,degeneracy) + density(massS,chemPot,temperature,degeneracy))/3.;
    
    return baryonNumberConservation - baryonDensity;

}

// The bissection method below uses the findChemPot function in order to calculate the chemical potential

double chemicalPotential(double baryonDensity,double C, double D, double temperature){

    double degeneracy = 6.;

    double xleft = 0.0, xright = kfmax, xmed, prod;
        
    if(isnan(findChemPot(xright,baryonDensity,C,D,temperature,degeneracy))){
        while(isnan(findChemPot(xright,baryonDensity,C,D,temperature,degeneracy))){
            xright -= 1000;
            if(xright<0.0){
                xright = 0.0;
                break;
            }
        }
    }

    if(isnan(findChemPot(xleft,baryonDensity,C,D,temperature,degeneracy))){
        while(isnan(findChemPot(xleft,baryonDensity,C,D,temperature,degeneracy))){
            xleft += 1000;
            if(xleft > xright){
                xleft = 0.0;
                break;
            }
        }
    }        
        
    while(fabs(xright - xleft) > epsilon){

        xmed = (xright + xleft)/2;

            /*if(isnan(zeroFunc(xmed,C,D,baryonDensity))){
                while(isnan(zeroFunc(xmed,C,D,baryonDensity))){
                    xmed-=1000;
                }
            }*/

        prod = findChemPot(xmed,baryonDensity,C,D,temperature,degeneracy)*findChemPot(xright,baryonDensity,C,D,temperature,degeneracy);

        if(prod<0){
            xleft = xmed;
        }
        else if(prod>0){
            xright = xmed;
        }
        else{
            break;
        }

    }

    return xmed;
}

// The contribution of each particle to the thermodynamic potential of a free system
double omegaZeroI(double mass, double chemPot, double temperature, double degeneracy){

    TF1 particles("particles","log(1 + exp(-(sqrt(x*x + [0]) - [1])/[2]))",kf0,kfmax);
    particles.SetParameter(0,mass*mass);
    particles.SetParameter(1,chemPot);
    particles.SetParameter(2,temperature);

    TF1 antiparticles("antiparticles","log(1 + exp(-(sqrt(x*x + [0]) + [1])/[2]))",kf0,kfmax);
    antiparticles.SetParameter(0,mass*mass);
    antiparticles.SetParameter(1,chemPot);
    antiparticles.SetParameter(2,temperature);

    TF1 omegaIntegral("omegaIntegral","x*x*(particles + antiparticles)",kf0,kfmax);
    
    double integral = omegaIntegral.Integral(kf0,kfmax);
    double result = -degeneracy*temperature*integral/(2*piSquared);

    return result;

}

// The thermodynamic potential of a free system
double omegaZero(double massU, double massD, double massS, double chemPot, double temperature, double degeneracy){

    double uOmegaZero = omegaZeroI(massU,chemPot,temperature,degeneracy);
    double dOmegaZero = omegaZeroI(massD,chemPot,temperature,degeneracy);
    double sOmegaZero = omegaZeroI(massS,chemPot,temperature,degeneracy);

    return uOmegaZero + dOmegaZero + sOmegaZero;

}

void symmetricQuarkMatterTemperature(){

    // The output file
    FILE *arq;
    arq = fopen("c017d148-quarks-t25.dat","w");

    double temperature = 25.0; // the temperature, in MeV

    // The free parameters
    double C = 0.170;
    double D = 148.*148.;

    for(double baryonDensity = 0.01*hbarc3; baryonDensity <= 1.5*hbarc3; baryonDensity += 0.005*hbarc3){

        // The mass of the three lightest quarks for a specific baryonic density
        double massU = uCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);
        double massD = dCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);
        double massS = sCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);

        double degeneracy = 6.; // only quarks are being considered
     
        // The effective chemical potential, equal to all quarks
        double chemPot = chemicalPotential(baryonDensity,C,D,temperature);

        // The particle densities
        double uDensity = density(massU,chemPot,temperature,degeneracy);
        double dDensity = density(massD,chemPot,temperature,degeneracy);
        double sDensity = density(massS,chemPot,temperature,degeneracy);

        // The thermodynamic quantities are given below
        double omegaZeroTotal = omegaZero(massU,massD,massS,chemPot,temperature,degeneracy);

        double freeEnergy = omegaZeroTotal + chemPot*(uDensity + dDensity + sDensity);

        double dOmegaZerodT = (omegaZero(massU,massD,massS,chemPot,temperature+delta,degeneracy) - omegaZero(massU,massD,massS,chemPot,temperature,degeneracy))/delta;

        double energyDensity = freeEnergy - temperature*dOmegaZerodT;

        double dOmegaZerodmu = (omegaZero(massU+delta,massD,massS,chemPot,temperature,degeneracy) - omegaZero(massU,massD,massS,chemPot,temperature,degeneracy))/delta;
        double dOmegaZerodmd = (omegaZero(massU,massD+delta,massS,chemPot,temperature,degeneracy) - omegaZero(massU,massD,massS,chemPot,temperature,degeneracy))/delta;
        double dOmegaZerodms = (omegaZero(massU,massD,massS+delta,chemPot,temperature,degeneracy) - omegaZero(massU,massD,massS,chemPot,temperature,degeneracy))/delta;
        double dOmegaZerodmI = dOmegaZerodmu + dOmegaZerodmd + dOmegaZerodms;

        long double dmIdnb = C/3/pow(baryonDensity,2./3.) - D/3/pow(baryonDensity,4./3.);

        double pressure = - omegaZeroTotal + baryonDensity*dmIdnb*dOmegaZerodmI;
        
        //The real chemical potentials are given below, different from the effective ones calculated above
        double realChemPot = chemPot + dmIdnb*dOmegaZerodmI/3.;

        double baryonDensityInFermi = baryonDensity/hbarc3;
        energyDensity = energyDensity/hbarc3;
        freeEnergy = freeEnergy/hbarc3;
        pressure = pressure/hbarc3;

        // Writes the results in the output file
        fprintf(arq,"%le  %le  %le  %le  %le  %le  %le\n",freeEnergy,energyDensity,pressure,baryonDensityInFermi,freeEnergy/baryonDensityInFermi,energyDensity/baryonDensityInFermi,realChemPot);

    }

    fclose(arq);

}