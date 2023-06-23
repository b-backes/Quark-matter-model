/*
Equation of State solver for the DDQM model at zero temperature, considering two-flavor quark matter in chemical equilibrium.
By Betânia Camille Backes, 2020
Should be executed with ROOT, see https://root.cern/
If the software is installed, open ROOT, compile and execute with
.X stellarMatter2fZeroT.C
*/

#include <stdio.h>
#include <math.h>
using namespace TMath;

const double uCurrentMass = 5.0; // Masses of particles (MeV)
const double dCurrentMass = 10.0;
const double massE = 0.511;

const double hbarc = 197.33;   // This one and the two below are used to convert from orders of MeV to orders of fm^-1
const double hbarc2 = 197.33*197.33;
const double hbarc3 = 197.33*197.33*197.33;

const double kf0 = 0.0;   //  Lower limit of the integral to calculate the energy density/pressure for the electrons

#define epsilon 1e-7   // Used to calculate numerical derivates

// Here we define the free particles thermodynamic potential density as a funcion

long double omegaZero(double EffChemPot, double mass){

    long double omegaZero = -1./4/M_PI/M_PI * (EffChemPot*sqrt(EffChemPot*EffChemPot - mass*mass)*(sqrt(EffChemPot*EffChemPot - 
    mass*mass)*sqrt(EffChemPot*EffChemPot - mass*mass) - 3*mass*mass/2) + 3*pow(mass,4)/2*log((EffChemPot+sqrt(EffChemPot*EffChemPot 
    - mass*mass))/mass));
  
    return omegaZero;

}

// Here we define a function that finds every particles densities we need

long double zeroFunc(double dDensity, double C, double D, double baryonDensity){

    double massU = uCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);
    double massD = dCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);

    long double eDensity1 = 2*baryonDensity - dDensity;

    long double part1 = sqrt(pow(M_PI*M_PI*dDensity,2./3.) + massD*massD);
    long double part2 = sqrt(pow(M_PI*M_PI*(3*baryonDensity - dDensity),2./3.) + massU*massU);
    long double part3 = (part1 - part2)*(part1 - part2) - massE*massE;

    long double eDensity2 = pow(part3,3./2.)/(3*M_PI*M_PI);

    return eDensity1 - eDensity2;

}

/*  Here we have the main part of our code. We start by creating a file for our results and making all
    the necessary physical calculations for a large range of the C and D parameters.  */

void stellarMatter2fZeroT(){

    FILE *arq;
    arq = fopen("stellar2f.dat","w");

    long double C = 0.727273;
    long double D = 128.080808*128.080808;

    for(long double baryonDensity = 0.01*hbarc3; baryonDensity <= 1.5*hbarc3; baryonDensity += 0.001*hbarc3){

        /* Here we calculate the mass for all three quarks for a specific baryonic density */

        long double massU = uCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);
        long double massD = dCurrentMass + D/cbrt(baryonDensity) + C*cbrt(baryonDensity);

        /* Since we were having a problem where the masses could be negative for some densities, here we will 
        determine a flag that will tell us that this is the problem we're having. This will be useful when defining stability,
        since even if a model has a minimum energy lower than 930 MeV, it shouldn't be stable if the masses are negative. */

        if(massU <= 0.0){
            //break;
        }

        /* Now we calculate the down quark density using the bissection method. We use the funtion defined before the main,
        that contains the analytical expressions obtained in order for us to determine dDensity numerically. */

        long double xleft = 0.0, xright = 2*baryonDensity, xmed, prod;

        if(isnan(zeroFunc(xright,C,D,baryonDensity))){
            while(isnan(zeroFunc(xright,C,D,baryonDensity))){
                xright -= 1000;
                if(xright<0.0){
                    xright = 0.0;
                    break;
                }
            }
        }

        if(isnan(zeroFunc(xleft,C,D,baryonDensity))){
            while(isnan(zeroFunc(xleft,C,D,baryonDensity))){
                xleft += 1000;
                if(xleft > xright){
                    xleft = 0.0;
                    break;
                }
            }
        }

        while(fabs(xright - xleft) > epsilon){

            xmed = (xright + xleft)/2;

            if(isnan(zeroFunc(xmed,C,D,baryonDensity))){
                while(isnan(zeroFunc(xmed,C,D,baryonDensity))){
                    xmed-=1000;
                }
            }

            prod = zeroFunc(xmed,C,D,baryonDensity)*zeroFunc(xright,C,D,baryonDensity);

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

        /* Here we calculate the other densities, that were analytically evaluated as functions of sDensity */
        long double dDensity = xmed;
        long double uDensity = 3.*baryonDensity - dDensity;
        long double eDensity = 2.*baryonDensity - dDensity;

        long double uFermiMom = cbrt(M_PI*M_PI*uDensity);
        long double dFermiMom = cbrt(M_PI*M_PI*dDensity);
        long double eFermiMom = cbrt(3*M_PI*M_PI*eDensity);

        long double uEffChemPot = sqrt(uFermiMom*uFermiMom + massU*massU);
        long double dEffChemPot = sqrt(dFermiMom*dFermiMom + massD*massD);
        long double eChemPot = sqrt(eFermiMom*eFermiMom + massE*massE);

        long double uOmegaZero = omegaZero(uEffChemPot,massU);
        long double dOmegaZero = omegaZero(dEffChemPot,massD);
        long double totalOmegaZero = uOmegaZero + dOmegaZero;

        long double dmIdnb = C/3/pow(baryonDensity,2./3.) - D/3/pow(baryonDensity,4./3.);

        double delta = 1e-5;

        long double dOmegaZerodEffChemPotU = (omegaZero(uEffChemPot+delta,massU) - omegaZero(uEffChemPot,massU))/delta;
        long double dOmegaZerodEffChemPotD = (omegaZero(dEffChemPot+delta,massD) - omegaZero(dEffChemPot,massD))/delta;

        long double energyTerm2 = uEffChemPot*dOmegaZerodEffChemPotU + dEffChemPot*dOmegaZerodEffChemPotD;

        long double dOmegaZerodmIU = (omegaZero(uEffChemPot,massU+delta) - omegaZero(uEffChemPot,massU))/delta;
        long double dOmegaZerodmID = (omegaZero(dEffChemPot,massD+delta) - omegaZero(dEffChemPot,massD))/delta;
        long double dOmegaZerodmI = dOmegaZerodmIU + dOmegaZerodmID;

        /* The contribution of the electrons to the system's energy density and pressure are the ones from a free particle.
        We express and calculate those quantities below.*/
        
        TF1 kernelElectron("kernelElectron","x*x*sqrt(x*x + [0])",kf0,eFermiMom);
        kernelElectron.SetParameter(0,massE*massE);
        long double energyDensElectron = kernelElectron.Integral(kf0,eFermiMom)/(Pi()*Pi());

        TF1 kernelPressureElectron("kernelPressureElectron","pow(x,4)/sqrt(x*x + [0])",kf0,eFermiMom);
        kernelPressureElectron.SetParameter(0,massE*massE);
        long double pressureElectron = kernelPressureElectron.Integral(kf0,eFermiMom)/(3*Pi()*Pi());

        if(eDensity < 0.0){
            energyDensElectron = 0.0;
            pressureElectron = 0.0;
        }

        // Here we sum everything up

        long double energyDensity, pressure;

        long double baryonDensityInFermi = baryonDensity/hbarc3;
        energyDensity = (totalOmegaZero - energyTerm2 + energyDensElectron)/hbarc3;
        pressure = (- totalOmegaZero + baryonDensity*dmIdnb*dOmegaZerodmI + pressureElectron)/hbarc3;
        
        // Here we print everything in our output file
        
        fprintf(arq,"%Le    %Le    %Le    %Le\n",energyDensity,pressure,baryonDensityInFermi,energyDensity/baryonDensityInFermi);


    }

    fclose(arq);

}
