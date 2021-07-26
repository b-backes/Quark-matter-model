/*

Equation of State solver for the DDQM model with zero temperature, considering symmetric quark matter.
By Betânia Camille Backes, 2020
Compile as:
gcc symmetricQuarkMatterZeroT.c -lm -Wall

*/

#include <stdio.h>
#include <math.h>

const double uCurrentMass = 5.0; // quarks masses in MeV
const double dCurrentMass = 10.0;
const double sCurrentMass = 80.0;

const double hbarc = 197.33;  // factors used to convert from orders of MeV to orders of fm^-1
const double hbarc2 = 197.33*197.33;
const double hbarc3 = 197.33*197.33*197.33;

#define epsilon 1e-7  // used in the bissection method below
#define delta 1e-5 // used in the numerical differentiaion, the result will usually be better than with a smaller value

//**************************************************************************************//

// The stability conditions for the chemical potentials and densities, rewritten so that the up quark density
// can be obtained utilising a method for finding the zeroes of a function. See the Appendix in https://arxiv.org/abs/2007.04494 

long double stabilityConditions(long double uDensity, double C, double D, double baryonDensity){

    double massU = uCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);
    double massD = dCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);
    double massS = sCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);

    long double part1 = pow(M_PI*M_PI*uDensity,2./3.) + massU*massU - massS*massS;
    long double sDensity1 = pow(part1,3./2.)/(M_PI*M_PI);

    long double part2 = pow(M_PI*M_PI*uDensity,2./3.) + massU*massU - massD*massD;
    long double part3 = pow(part2,3./2.)/(M_PI*M_PI);
    long double part4 = 3*baryonDensity - uDensity;
    long double sDensity2 = part4 - part3;

    return sDensity2 - sDensity1;

}

// The bissection method below uses the stabilityConditions function in order to calculate the up quark density

long double findUDensity(double baryonDensity, double C, double D){

    long double xleft = 0.0, xright = baryonDensity*4.0, xmed, prod;

    if(isnan(stabilityConditions(xright,C,D,baryonDensity))){
        while(isnan(stabilityConditions(xright,C,D,baryonDensity))){
            xright -= 1000;
            if(xright<0.0){
                xright = 0.0;
                break;
            }
        }
    }

    if(isnan(stabilityConditions(xleft,C,D,baryonDensity))){
        while(isnan(stabilityConditions(xleft,C,D,baryonDensity))){
            xleft += 1000;
            if(xleft > xright){
                xleft = 0.0;
                break;
            }
        }
    }

    while(fabs(xright - xleft) > epsilon*hbarc3){

        xmed = (xright + xleft)/2;


        if(isnan(stabilityConditions(xmed,C,D,baryonDensity))){
            while(isnan(stabilityConditions(xmed,C,D,baryonDensity))){
                xmed-=1000;
            }
        }

        prod = stabilityConditions(xmed,C,D,baryonDensity)*stabilityConditions(xright,C,D,baryonDensity);

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

// Thermodynamic potential of a free Fermi gas
long double omegaZero(double EffChemPot, double mass){

    long double omegaZero = -1./4/M_PI/M_PI * (EffChemPot*sqrt(EffChemPot*EffChemPot - mass*mass)*(sqrt(EffChemPot*EffChemPot - 
    mass*mass)*sqrt(EffChemPot*EffChemPot - mass*mass) - 3*mass*mass/2) + 3*pow(mass,4)/2*log((EffChemPot+sqrt(EffChemPot*EffChemPot 
    - mass*mass))/mass));
  
    return omegaZero;

}

int main(void){

    // The output file
    FILE *arq;
    arq = fopen("c0d1585-quarks-t0.dat","w");

    // The free parameters
    double C = 0.0;
    double D = 158.5*158.5;

    for(long double baryonDensity = 0.01*hbarc3; baryonDensity <= 1.5*hbarc3; baryonDensity += 0.005*hbarc3){
        
        // The mass for the three lightest quarks for a specific baryonic density
        long double massU = uCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);
        long double massD = dCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);
        long double massS = sCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);   

        // The paricle densities
        long double uDensity = findUDensity(baryonDensity,C,D);
        long double dDensity = pow(pow(M_PI*M_PI*uDensity,2./3.) + massU*massU - massD*massD,3./2.)/(M_PI*M_PI);
        long double sDensity = 3*baryonDensity - uDensity - dDensity;

        // The Fermi momentum and effective chemicalpotential of all particles        
        long double uFermiMom = cbrt(M_PI*M_PI*uDensity);
        long double dFermiMom = cbrt(M_PI*M_PI*dDensity);
        long double sFermiMom = cbrt(M_PI*M_PI*sDensity);

        long double uEffChemPot = sqrt(uFermiMom*uFermiMom + massU*massU);
        long double dEffChemPot = sqrt(dFermiMom*dFermiMom + massD*massD);
        long double sEffChemPot = sqrt(sFermiMom*sFermiMom + massS*massS);
        
        // All the terms utilised to compute the energy density and the pressure

        long double uOmegaZero = omegaZero(uEffChemPot,massU);
        long double dOmegaZero = omegaZero(dEffChemPot,massD);
        long double sOmegaZero = omegaZero(sEffChemPot,massS);
        long double totalOmegaZero = uOmegaZero + dOmegaZero + sOmegaZero; 

        long double dmIdnb = C/3/pow(baryonDensity,2./3.) - D/3/pow(baryonDensity,4./3.);

        long double dOmegaZerodEffChemPotU = (omegaZero(uEffChemPot+delta,massU) - omegaZero(uEffChemPot,massU))/delta;
        long double dOmegaZerodEffChemPotD = (omegaZero(dEffChemPot+delta,massD) - omegaZero(dEffChemPot,massD))/delta;
        long double dOmegaZerodEffChemPotS = (omegaZero(sEffChemPot+delta,massS) - omegaZero(sEffChemPot,massS))/delta;

        long double energyTerm2 = uEffChemPot*dOmegaZerodEffChemPotU + dEffChemPot*dOmegaZerodEffChemPotD + sEffChemPot*dOmegaZerodEffChemPotS;

        long double dOmegaZerodmIU = (omegaZero(uEffChemPot,massU+delta) - omegaZero(uEffChemPot,massU))/delta;
        long double dOmegaZerodmID = (omegaZero(dEffChemPot,massD+delta) - omegaZero(dEffChemPot,massD))/delta;
        long double dOmegaZerodmIS = (omegaZero(sEffChemPot,massS+delta) - omegaZero(sEffChemPot,massS))/delta;
        long double dOmegaZerodmI = dOmegaZerodmIU + dOmegaZerodmID + dOmegaZerodmIS;

        // If needed, the real chemical potentials are given below
        /*long double uRealChemPot = uEffChemPot + dmIdnb*dOmegaZerodmI/3.;
        long double dRealChemPot = dEffChemPot + dmIdnb*dOmegaZerodmI/3.;
        long double sRealChemPot = sEffChemPot + dmIdnb*dOmegaZerodmI/3.;*/

        // Assembles the previous terms to write the energy density and pressure
        long double energyDensity, pressure;

        long double baryonDensityInFermi = baryonDensity/hbarc3;
        energyDensity = (totalOmegaZero - energyTerm2)/hbarc3;
        pressure = (- totalOmegaZero + baryonDensity*dmIdnb*dOmegaZerodmI)/hbarc3;

        // Writes the results in the output file
        fprintf(arq,"%Le    %Le    %Le    %Le\n",energyDensity,pressure,baryonDensityInFermi,energyDensity/baryonDensityInFermi);

    }

    // Closes the output file
    fclose(arq);

    return 0;

}
