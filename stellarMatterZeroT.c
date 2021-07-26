/*

Equation of State solver for the DDQM model at zero temperature, considering stellar matter.
By Betânia Camille Backes, 2020
Compile as:
gcc stellarMatterZeroT.c -lm -Wall

*/

#include <stdio.h>
#include <math.h>

const double uCurrentMass = 5.0; // quarks masses in MeV
const double dCurrentMass = 10.0;
const double sCurrentMass = 80.0;
const double massE = 0.511;

const double hbarc = 197.33;  // factors used to convert from orders of MeV to orders of fm^-1
const double hbarc2 = 197.33*197.33;
const double hbarc3 = 197.33*197.33*197.33;

#define epsilon 1e-7 // used in the bissection method below
#define delta 1e-5 // used in the numerical differentiaion, the result will usually be better than with a smaller value

// The thermodynamic potential of a free system
long double omegaZero(double EffChemPot, double mass){

    long double omegaZero = -1./4/M_PI/M_PI * (EffChemPot*sqrt(EffChemPot*EffChemPot - mass*mass)*(sqrt(EffChemPot*EffChemPot - 
    mass*mass)*sqrt(EffChemPot*EffChemPot - mass*mass) - 3*mass*mass/2) + 3*pow(mass,4)/2*log((EffChemPot+sqrt(EffChemPot*EffChemPot 
    - mass*mass))/mass));
  
    return omegaZero;

}

// The stability conditions for the chemical potentials and densities, rewritten so that the strange quark density
// can be obtained utilising a method for finding the zeroes of a function. See the Appendix in https://arxiv.org/abs/2007.04494

long double stabilityConditions(double sDensity, double C, double D, double baryonDensity){

    double massU = uCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);
    double massD = dCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);
    double massS = sCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);

    long double part1 = pow(M_PI*M_PI*sDensity,2./3.);
    long double part2 = pow(part1 + massS*massS - massD*massD,3./2.)/(M_PI*M_PI);

    long double eDensity1 = 2*baryonDensity - sDensity - part2;

    long double part3 = pow(part1 + massS*massS - massD*massD,3./2.)/(M_PI*M_PI);
    long double uDensity = 3*baryonDensity - sDensity - part3;
    long double part6 = M_PI*M_PI*uDensity;
    long double part4 = pow(part6,2./3.);
    long double part5 = pow(sqrt(part4 + massU*massU) - sqrt(part1 + massS*massS),2);
    long double eDensity2 = pow(part5 - massE*massE,3./2.)/(M_PI*M_PI);

    return eDensity1 - eDensity2;

}

// The bissection method below uses the stabilityConditions function in order to calculate the strange quark density

long double findDensity(double baryonDensity, double C, double D){

    long double xleft = 0.0, xright = baryonDensity, xmed, prod;

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

    while(fabs(xright - xleft) > epsilon){

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

int main(void){

    // The output file
    FILE *arq;
    arq = fopen("eos_C05_D13575-test.dat","w");

    // The free parameters
    long double C = 0.5;
    long double D = 135.75*135.75;

    for(long double baryonDensity = 0.01*hbarc3; baryonDensity <= 1.5*hbarc3; baryonDensity += 0.005*hbarc3){

        // The mass of the three lightest quarks for a specific baryonic density
        long double massU = uCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);
        long double massD = dCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);
        long double massS = sCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);

        // The particle densities
        long double sDensity = findDensity(baryonDensity,C,D);
        long double dDensity = pow(pow(M_PI*M_PI*sDensity,2./3.) + massS*massS - massD*massD,3./2.)/(M_PI*M_PI);
        long double uDensity = 3*baryonDensity - dDensity - sDensity;
        long double eDensity = 2.*uDensity/3. - dDensity/3. - sDensity/3.;

        // The Fermi momentum and effective chemical potential of all particles.
        // For the quarks, the chemical potential is the effective ones. For leptons, the real ones.
        long double sFermiMom = cbrt(M_PI*M_PI*sDensity);
        long double uFermiMom = cbrt(M_PI*M_PI*uDensity);
        long double dFermiMom = cbrt(M_PI*M_PI*dDensity);
        long double eFermiMom = cbrt(3*M_PI*M_PI*eDensity);

        long double sEffChemPot = sqrt(sFermiMom*sFermiMom + massS*massS);
        long double uEffChemPot = sqrt(uFermiMom*uFermiMom + massU*massU);
        long double dEffChemPot = sqrt(dFermiMom*dFermiMom + massD*massD);
        long double eChemPot = sqrt(eFermiMom*eFermiMom + massE*massE);

        // All the terms utilised to compute the energy density and the pressure

        long double uOmegaZero = omegaZero(uEffChemPot,massU);
        long double dOmegaZero = omegaZero(dEffChemPot,massD);
        long double sOmegaZero = omegaZero(sEffChemPot,massS);
        // The 3.0 factor below is due to the degeneracy of the leptons being 2 instead of 6
        long double eOmegaZero = omegaZero(eChemPot,massE)/3.0;
        long double totalOmegaZero = uOmegaZero + dOmegaZero + sOmegaZero + eOmegaZero;

        long double dmIdnb = C/3/pow(baryonDensity,2./3.) - D/3/pow(baryonDensity,4./3.);

        long double dOmegaZerodEffChemPotU = (omegaZero(uEffChemPot+delta,massU) - omegaZero(uEffChemPot,massU))/delta;
        long double dOmegaZerodEffChemPotD = (omegaZero(dEffChemPot+delta,massD) - omegaZero(dEffChemPot,massD))/delta;
        long double dOmegaZerodEffChemPotS = (omegaZero(sEffChemPot+delta,massS) - omegaZero(sEffChemPot,massS))/delta;

        long double energyTerm2 = uEffChemPot*dOmegaZerodEffChemPotU + dEffChemPot*dOmegaZerodEffChemPotD + sEffChemPot*dOmegaZerodEffChemPotS - eChemPot*eDensity;

        long double dOmegaZerodmIU = (omegaZero(uEffChemPot,massU+delta) - omegaZero(uEffChemPot,massU))/delta;
        long double dOmegaZerodmID = (omegaZero(dEffChemPot,massD+delta) - omegaZero(dEffChemPot,massD))/delta;
        long double dOmegaZerodmIS = (omegaZero(sEffChemPot,massS+delta) - omegaZero(sEffChemPot,massS))/delta;
        long double dOmegaZerodmI = dOmegaZerodmIU + dOmegaZerodmID + dOmegaZerodmIS;

        // Assembles the previous terms to write the energy density and pressure

        long double energyDensity, pressure;

        long double baryonDensityInFermi = baryonDensity/hbarc3;
        energyDensity = (totalOmegaZero - energyTerm2)/hbarc3;
        pressure = (- totalOmegaZero + baryonDensity*dmIdnb*dOmegaZerodmI)/hbarc3;
        long double baryonicChemPot = (energyDensity + pressure)/baryonDensityInFermi;
        
        // Writes the results when using the EoS as input for the TOV equations
        //fprintf(arq,"%Le    %Le    %Le    %Le\n",energyDensity,pressure,baryonDensityInFermi,energyDensity/baryonDensityInFermi);

        // Writes the results in the output file
        //if(!isnan(pressure) && pressure >= 0.0){
            fprintf(arq,"%Le    %Le    %Le   %Le\n",baryonDensityInFermi,energyDensity/hbarc,pressure/hbarc,baryonicChemPot);
        //}

    }

    // Should be printed when writing a file to be used as input for the TOV equations
    // fprintf(arq,"-1.,-1.,-1.,");

    // Closes the output file
    fclose(arq);

    return 0;

}