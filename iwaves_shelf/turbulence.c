/* Note: 
 * set kappa_min
 * need to define phys->TP, phys->TB, phys->TD
 * modified to be GLS + CA stability function
 * turbmodel=1 k-kl; 2 gen; 3 k-epsilon; 4 k-omega
 * File: turbulence.c
 * Author: Bing Wang
 * Institution: Stanford University
 * --------------------------------
 * Contains the Mellor-Yamada level 2.5 turbulence model.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "math.h"
#include "phys.h"
#include "grid.h"
#include "phys.h"
#include "util.h"
#include "turbulence.h"
#include "boundaries.h"
#include "scalars.h"

// Local function
static void StabilityFunctions(REAL *Sm, REAL *Sh, REAL Gh, REAL Gm);

/*
 * Function: my25
 * Usage: my25(grid,phys,prop,phys->qT,phys->lT,phys->Cn_q,phys->Cn_l,phys->nu_tv,phys->kappa_tv);
 * -----------------------------------------------------------------------------------------------
 * Computes the eddy viscosity and scalar diffusivity via the MY25 closure.  Advection of the turbulent
 * quantities q^2 and q^2l is included with the use of UpdateScalars
 *
 */
void my25(gridT *grid, physT *phys, propT *prop, REAL **wnew, REAL **q, REAL **l, REAL **Cn_q, REAL **Cn_l, REAL **nuT, REAL **kappaT, MPI_Comm comm, int myproc) {
  int i, ib, j, iptr, jptr, k, nc1, nc2, ne, kk, nf;
  REAL thetaQ=1, CdAvgT, CdAvgB, tauAvgT, *dudz, *dvdz, *drdz, z, *NN, *Gh, Gm;
  REAL Sm, Sh;
  REAL u_taub, u_taus, z0s, z0b; // hardwire holder
  REAL Ps, Pb, Diss, cpsi3;
  REAL L_min, Lcrit, Lold;
  REAL **Qold,**P, **B, **D;
  REAL ftop, fbot;
  REAL gls_m, gls_n, gls_p, cpsi1, cpsi2, cpsi3minus, cpsi3plus, sig_kpsi, sig_k, sig_psi, alpha, E2, E4;
  REAL nu_min=1e-6; // nu_min is set to be same as the molecular viscosity

  // minimum values, slightly different from Warner(2005)
  const REAL psi_min    = 1.e-12; 
  const REAL k_min      = 1.e-8;  // 7.6e-6; 
  const REAL eps_min =   1.e-10;  // used to determine L_min
  const REAL cm0     =   0.5270; //0.5477;  // Warner uses 0.5544
  
  // Define Gh smoothing parameters, Same as Warner's paper, a constant different from Burchard's paper
  const REAL Ghmax = 0.0233; //0.0233 KC; 0.0329 CA; 
  const REAL Ghmin = -0.28;  //-0.2809;
  const REAL Ghcrit = 0.03;  //0.02;

  switch (prop->turbmodel) {
	    /* k-kl */
	  case 1:
	    gls_m      =  1.0;   // power on k
	    gls_n      =  1.0;  // power on l
	    gls_p      =  0.0;   // power on cu0
	    cpsi1      =  0.9;   // c1 for psi equation, shear production
	    cpsi2      =  0.5;  // c2 for psi equation, dissipation
	    cpsi3minus =  2.38; //2.53;  // MY2.5 c3-=0.9, c3- for psi equation, buoyancy production, when unstably stratified
	    cpsi3plus  =  1.0;   // MY2.5 c3+=0.9, c3+ for psi equation, buoyance production, when stably stratified
	    //sig_kpsi   =  0.8;   // sigma_k,psi schmidt number for k equation of k-epsilon model 
	    sig_k      =  2.44;   // sigma_k, schmidt number for k equation of GOTM
	    sig_psi    =  2.44;  // sigma_psi, schmidt number for psi equation of GOTM
	    alpha      =  0.0;
	    E2         =  1.33;  // constant in the wall function
	    E4         =  0.25;   // constant in the wall function
	    break;

	    /* gen */
	  case 2:
	    gls_m      =  1.0;   // power on k
	    gls_n      = -0.67;  // power on l
	    gls_p      =  2.0;   // power on cu0
	    cpsi1      =  1.0;   // c1 for psi equation, shear production
	    cpsi2      =  1.22;  // c2 for psi equation, dissipation
	    cpsi3minus = 0.05; //0.10;  // c3+ for psi equation, buoyancy production, when unstably stratified
	    cpsi3plus  =  1.0;   // c3- for psi equation, buoyance production, when stably stratified
	    //sig_kpsi   =  0.8;   // not used, sigma_k,psi schmidt number for k equation of k-epsilon model 
	    sig_k      =  0.8;   // sigma_k, schmidt number for k equation of GOTM
	    sig_psi    =  1.07;  // sigma_psi, schmidt number for psi equation of GOTM
	    //gen_d      = -1.2;   // not used 
	    //gen_alpha  = -2.0;   // not used
	    //gen_l      =  0.2;   // not used
	    //alpha      =  0.0;   // buoyancy deduction in shear production ? not sure
	    E2      =  0.0;   // wall function
	    E4      = 0.0; 
	    break;

	    /* k-epsilon */
	  case 3:
	    gls_m      =  1.5;   // power on k
	    gls_n      =  -1.0;  // power on l
	    gls_p      =  3.0;   // power on cu0
	    cpsi1      =  1.44;   // c1 for psi equation, shear production
	    cpsi2      =  1.92;  // c2 for psi equation, dissipation
	    cpsi3minus =  -0.63;//-0.41;  // c3+ for psi equation, buoyancy production, when unstably stratified, 0.52 or 0.41 
	    cpsi3plus  =  1.0;   // c3- for psi equation, buoyance production, when stably stratified
	    //sig_kpsi   =  0.8;   // sigma_k,psi schmidt number for k equation of k-epsilon model 
	    sig_k      =  1.;   // sigma_k, schmidt number for k equation of GOTM
	    sig_psi    =  1.3;  // sigma_psi, schmidt number for psi equation of GOTM
	    //alpha      =  0.0;  // effect of buoyancy on shear production (see Gross's code or GOTM, not used here)
	    E2         =  0.0;  // constant in the wall function
	    E4         = 0.0;   // constant in the wall function
	    break;

	    /* k-omega */
	  case 4:
	    gls_m      =  0.5;   // power on k
	    gls_n      =  -1.0;  // power on l
	    gls_p      =  -1.0;   // power on cu0
	    cpsi1      =  0.555;   // c1 for psi equation, shear production
	    cpsi2      =  0.833;  // c2 for psi equation, dissipation
	    cpsi3minus = -0.64; //-0.58;  // c3- for psi equation, buoyancy production, when unstably stratified
	    cpsi3plus  =  1.0;   // c3+ for psi equation, buoyance production, when stably stratified
	    //sig_kpsi   =  0.8;   // sigma_k,psi schmidt number for k equation of k-epsilon model 
	    sig_k      =  2.;   // sigma_k, schmidt number for k equation of GOTM
	    sig_psi    =  2.;  // sigma_psi, schmidt number for psi equation of GOTM
	    //alpha      =  0.0;
	    E2         =  0.;  // constant in the wall function
	    E4         =  0.;   // constant in the wall function
    	    break;
	  }
  // Initialise minimum length scale (from Gross' code)
  L_min = pow(cm0,3)*pow(k_min,1.5)/eps_min;
  //nu_min = sqrt(2.*k_min)*L_min;

  // load array points
  NN = dudz = phys->a;
  dvdz = phys->b;
  drdz = phys->c;
  Gh = phys->d;
  Qold = phys->vold;
  P = phys->TP;
  B = phys->TB;
  D = phys->TD;

  for(i=0;i<grid->Nc;i++)
    for(k=0;k<grid->Nk[i];k++)
      P[i][k]=B[i][k]=D[i][k]=.0;
      //phys->T[i][k]=P[i][k]=B[i][k]=D[i][k]=.0;

  /********* (1) solve the tke equation *********/  
  for(i=0;i<grid->Nc;i++) {
    if(grid->ctop[i]<grid->Nk[i]) {
   // bottom and top shear stress quantities for boundary conditions
    CdAvgT = 0.0;
    CdAvgB = 0.0;
    tauAvgT = 0.0;
    //for(nf=0;nf<NFACES;nf++) {
      //ne=grid->face[i*NFACES+nf];
    for(nf=0;nf<grid->nfaces[i];nf++) {
      ne=grid->face[i*grid->nfaces[i]+nf];
      //CdAvgT += phys->CdT[ne]/3;
      //CdAvgB += phys->CdB[ne]/3;
      CdAvgT+=phys->CdT[ne]/grid->nfaces[i];
      CdAvgB+=phys->CdB[ne]/grid->nfaces[i];
      tauAvgT += fabs(phys->tau_T[ne])/grid->nfaces[i];
    }
    // modify wind stress here
    //tauAvgT = fabs(prop->tau_T); // bing, for constant windstress,this won't not affected by the normals
    //tauAvgT = 0.0;
    u_taub = sqrt(CdAvgB*(pow(phys->uc[i][grid->Nk[i]-1],2)+pow(phys->vc[i][grid->Nk[i]-1],2)));
    u_taus = sqrt(CdAvgT*(pow(phys->uc[i][grid->ctop[i]],2)+pow(phys->vc[i][grid->ctop[i]],2))+tauAvgT);


    if(grid->ctop[i]<grid->Nk[i]-1) {  // there are two or more layers
      // dudz, dvdz, and drdz store gradients at k-1/2
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
	dudz[k] = 2.0*(phys->uc[i][k-1]-phys->uc[i][k])/(grid->dzz[i][k-1]+grid->dzz[i][k]);
	dvdz[k] = 2.0*(phys->vc[i][k-1]-phys->vc[i][k])/(grid->dzz[i][k-1]+grid->dzz[i][k]);
	drdz[k] = 2.0*prop->beta*(phys->s[i][k-1]-phys->s[i][k])/(grid->dzz[i][k-1]+grid->dzz[i][k]);
      }
      if(u_taus<SMALL){ // shear-free surface
	dudz[grid->ctop[i]] = dudz[grid->ctop[i]+1];
	dvdz[grid->ctop[i]] = dvdz[grid->ctop[i]+1];
      }
      else {            // stressed surface, linear (x1) or logarithmic approximation (x3) 
	dudz[grid->ctop[i]] = 1.*dudz[grid->ctop[i]+1];  // using a factor of 3, it will give 2x of the linear gradient after taking average
	dvdz[grid->ctop[i]] = 1.*dvdz[grid->ctop[i]+1];
      }
      drdz[grid->ctop[i]]=drdz[grid->ctop[i]+1];

      
      if(u_taub<SMALL){ // shear-free surface
	dudz[grid->Nk[i]] = dudz[grid->Nk[i]-1];
	dvdz[grid->Nk[i]] = dvdz[grid->Nk[i]-1];
      }
      else {            // stressed surface, linear (x1) or logarithmic approximation (x3) 
	dudz[grid->Nk[i]] = 5.*dudz[grid->Nk[i]-1];  // when take average, factor 3. here will give 2x of the gradient
	dvdz[grid->Nk[i]] = 5.*dvdz[grid->Nk[i]-1];
      }
      drdz[grid->Nk[i]]=drdz[grid->Nk[i]-1];
    } else if(grid->ctop[i]==grid->Nk[i]-1) { // one wet layer, use logarithmic approximation, also treat depth as dz not dzz (stability)
      dudz[grid->ctop[i]] = u_taub/KAPPA_VK/grid->dz[grid->ctop[i]]
	*fabs(phys->uc[i][grid->ctop[i]])/(SMALL+sqrt(pow(phys->uc[i][grid->ctop[i]],2)+pow(phys->vc[i][grid->ctop[i]],2)));
      dvdz[grid->ctop[i]] = u_taub/KAPPA_VK/grid->dz[grid->ctop[i]]
	*fabs(phys->vc[i][grid->ctop[i]])/(SMALL+sqrt(pow(phys->uc[i][grid->ctop[i]],2)+pow(phys->vc[i][grid->ctop[i]],2)));
      drdz[grid->ctop[i]] = 0.;
    
      dudz[grid->Nk[i]] = 1.*dudz[grid->ctop[i]];
      dvdz[grid->Nk[i]] = 1.*dvdz[grid->ctop[i]];
      drdz[grid->Nk[i]] = 0.;

    }

    //temp, correct fs oscillation
    if(grid->ctop[i]!=grid->ctopold[i] && grid->Nk[i]>grid->ctop[i]+2){
      dudz[grid->ctop[i]] *=0.2; //dudz[grid->ctop[i]+1]=dudz[grid->ctop[i]+2];
      dvdz[grid->ctop[i]] *=0.2; //dvdz[grid->ctop[i]+1]=dvdz[grid->ctop[i]+2];
    }
    // q will store k=q^2/2, set minimal value and Qold will store the old value of k
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      q[i][k]*=0.5*q[i][k];
      l[i][k] = Max(l[i][k],L_min);
      q[i][k] = Max(q[i][k],k_min);
      Qold[i][k] = q[i][k];
      kappaT[i][k] = Max(kappaT[i][k],nu_min*kappaT[i][k]/nuT[i][k]); // kappaT/nuT remains unchanged
      nuT[i][k] = Max(nuT[i][k],nu_min);
    }
    
    // uold will store src1, the dissipation coefficient for k, which is the (espilon+Pbminus)/k term
    // wtmp will store src2, the production for k, which is the (Ps+Pbplus) term
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      //if(i==4852)
      //printf("dudz[%d]=%e, dvdz=%e, nu=%e, kappa=%e\n",k,dudz[k],dvdz[k],nuT[i][k],kappaT[i][k]);
      Ps = nuT[i][k]*(pow(0.5*(dudz[k]+dudz[k+1]),2)+pow(0.5*(dvdz[k]+dvdz[k+1]),2));  // use viscosity, minimum already set
      Pb = GRAV*kappaT[i][k]*0.5*(drdz[k]+drdz[k+1]); // use mass diffusivity, minimum already set
      Diss = Max(eps_min, pow(cm0,3)*pow(q[i][k],1.5)/l[i][k]);
	    
      // store these values in temporary variables for the lengthscale equation
      P[i][k] = Ps;
      B[i][k] = Pb;
      D[i][k] = Diss;




      if(Ps+Pb>0) {
	phys->wtmp[i][k] = Ps+Pb;	
	phys->uold[i][k] = Diss/q[i][k]; 
      } else {
	phys->wtmp[i][k] = Ps;
	phys->uold[i][k] = (-Pb+Diss)/q[i][k];
      }
    }

    for(k=0; k<grid->ctop[i]; k++) {
	phys->wtmp[i][k] = phys->wtmp[i][grid->ctop[i]];
	phys->uold[i][k] = phys->uold[i][grid->ctop[i]];
      }

    // kappaT will store the diffusion coefficient for k, mass diffusivity is overwritten
    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      kappaT[i][k] = phys->nu_tv[i][k]/sig_k;
    
      
    // htmp will store the value at the top boundary for k
    // hold will store it at the bottom boundary
    // no boundary flux for k
    phys->htmp[i]= 0;
    phys->hold[i]= 0; 
    } 
  
    // Specify turbulence at boundaries for use in updatescalars.  Assume that all incoming turbulence is zero and let outgoing
    // turbulence flow outward.
    for(jptr=grid->edgedist[2];jptr<grid->edgedist[5];jptr++) {
      j = grid->edgep[jptr];
      ib = grid->grad[2*j];
    
      for(k=grid->ctop[ib];k<grid->Nk[ib];k++) 
	phys->boundary_tmp[jptr-grid->edgedist[2]][k]=q[ib][k];
    }
  }


  UpdateScalars(grid,phys,prop,wnew,q,phys->boundary_tmp,Cn_q,0.0,0.0,kappaT,thetaQ,phys->uold,phys->wtmp,phys->htmp,phys->hold,0,0, comm, myproc,0,prop->TVDturb);


  /********* solve the GLS euqation *********/
  for(i=0;i<grid->Nc;i++) {
    if(grid->ctop[i]<grid->Nk[i]) {
    // bottom and top shear stress quantities for boundary conditions
    CdAvgT=0.0;
    CdAvgB=0.0;
    tauAvgT=0.0;
    //for(nf=0;nf<NFACES;nf++) {
    //  ne=grid->face[i*NFACES+nf];
    for(nf=0;nf<grid->nfaces[i];nf++) {
      ne=grid->face[i*grid->nfaces[i]+nf];
      CdAvgT+=phys->CdT[ne]/grid->nfaces[i];
      CdAvgB+=phys->CdB[ne]/grid->nfaces[i];
      //tauAvgT += fabs(phys->tau_T[ne])/3;
    }
    tauAvgT=fabs(prop->tau_T); // bing, for uniform windstress,this won't not affected by the normals
    u_taub = sqrt(CdAvgB*(pow(phys->uc[i][grid->Nk[i]-1],2)+pow(phys->vc[i][grid->Nk[i]-1],2)));
    u_taus = sqrt(CdAvgT*(pow(phys->uc[i][grid->ctop[i]],2)+pow(phys->vc[i][grid->ctop[i]],2))+tauAvgT);
    z0b = 0.5*grid->dz[grid->Nk[i]-1];  // use a factor 0.1 to 0.5, for stability, determined from tests
    z0s = 0.5*grid->dz[grid->ctop[i]];  // stability use a factor 0.1 to 0.5, for stability, determined from tests

    z = phys->h[i];
    // kappaT now stores the diffusivity for gls
    // l will store gls, old value of l is in Lold
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      z-=grid->dzz[i][k]/2; // z is used in the wall function
      kappaT[i][k] = nuT[i][k]/sig_psi;
      Lold = l[i][k];      // used in the wall function
      l[i][k] = pow(cm0,gls_p)*pow(Qold[i][k],gls_m)*pow(l[i][k],gls_n); 
      l[i][k] = Max(l[i][k],psi_min);  
 

      // choose cpsi3 based on stable or unstable density gradient
      if(B[i][k]>0)
	cpsi3 = cpsi3plus;
      else
	cpsi3 = cpsi3minus;

      // uold will store src1 for gls, 
      // wtmp will store src2 for gls, 
      Ps = cpsi1*P[i][k]*l[i][k]/Qold[i][k];
      Pb = cpsi3*B[i][k]*l[i][k]/Qold[i][k];
      Diss = cpsi2*D[i][k]*l[i][k]/Qold[i][k]*
	(1+E2*pow(Lold/KAPPA_VK/(z+grid->dv[i]),2)+E4*pow(Lold/KAPPA_VK/(z-phys->h[i]),2));

      if(Ps+Pb>0) {
	phys->wtmp[i][k] = Ps+Pb;	
	phys->uold[i][k] = Diss/l[i][k]; //2.0*q[i][k]/B1/(l[i][k]+SMALL);
      } else {
	phys->wtmp[i][k] = Ps;
	phys->uold[i][k] = (-Pb+Diss)/l[i][k];
      }
      z-=grid->dzz[i][k]/2;
    }

    // Top and bottom boundary fluxes see Warner 2005 eqn 53
     k = grid->ctop[i];
     phys->htmp[i] = (-1.)*gls_n*kappaT[i][k]*pow(cm0,gls_p)*pow(KAPPA_VK,gls_n)*
       pow(Qold[i][k],gls_m)*pow(Max(z0s,0.5*grid->dzz[i][k]),gls_n-1);
     
     k = grid->Nk[i]-1;
     phys->hold[i] = gls_n*kappaT[i][k]*pow(cm0,gls_p)*pow(KAPPA_VK,gls_n)*
       pow(Qold[i][k],gls_m)*pow(Max(z0b,0.5*grid->dzz[i][k]),gls_n-1);
     
  }

  // Specify turbulence at boundaries for use in updatescalars.  
  // Assume that all incoming turbulence is zero and let outgoing
  // turbulence flow outward.
    for(jptr=grid->edgedist[2];jptr<grid->edgedist[5];jptr++) {
      j = grid->edgep[jptr];
      ib = grid->grad[2*j];
      
      for(k=grid->ctop[ib];k<grid->Nk[ib];k++) 
	phys->boundary_tmp[jptr-grid->edgedist[2]][k]=l[ib][k];
    }  
  }
  UpdateScalars(grid,phys,prop,wnew,l,phys->boundary_tmp,Cn_l,0.0,0.0,kappaT,thetaQ,phys->uold,phys->wtmp,
		phys->htmp,phys->hold,0,0, comm, myproc,0,prop->TVDturb);

  //i=4852;k=0;
  //printf("k=%e, l=%e, P=%e, B=%e, E=%e\n",q[i][k],l[i][k],P[i][k],B[i][k],D[i][k]);


  /********* set the limits for psi and k and extract q and l *********/
  for(i=0;i<grid->Nc;i++) {
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      // recover length scale
      if(l[i][k]>psi_min && q[i][k]>k_min)
	l[i][k] = pow((l[i][k]*pow(cm0,-gls_p)*pow(q[i][k],-gls_m)),1/gls_n);
      else {
	q[i][k] = Max(q[i][k], k_min);
	l[i][k] = L_min;
      }
      q[i][k] = sqrt(2*q[i][k]);

      // length limitor based on geometry to stablize wetting and drying
      // TEMP
      //l[i][k] = Min(l[i][k],grid->dv[i]+phys->h[i]);
      //if(k==grid->Nk[i]-1 || k==grid->ctop[i])
      // l[i][k] = Min(l[i][k],KAPPA_VK*grid->dz[k]); 
    }
  }

  // Send/Recv k and gls data to neighboring processors before it computes nu_tv and kappa_tv
  ISendRecvCellData3D(q,grid,myproc,comm);
  ISendRecvCellData3D(l,grid,myproc,comm);

  for(i=0;i<grid->Nc;i++) {
    // compute NN=N^2 for the use in lengthscale limit and stability functions
    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
      //this value is at the top face of layer k
      NN[k]=-2.0*GRAV*prop->beta*(phys->s[i][k-1]-phys->s[i][k])/(grid->dzz[i][k-1]+grid->dzz[i][k]);
      //if(NN[k]<0) NN[k]=0;
    }
    NN[grid->Nk[i]] = NN[grid->Nk[i]-1];
    NN[grid->ctop[i]] = NN[grid->ctop[i]+1];
   
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      // length scale limitor
      if(NN[k]+NN[k+1] > SMALL){
	Lcrit = sqrt(0.56*q[i][k]*q[i][k]/(NN[k]+NN[k+1])); 
	//Lcrit = Max(Lcrit, pow(q[i][k],3)/prop->nu/(NN[k]+NN[k+1])*2/45); // internal wave limit
	l[i][k] = Min(l[i][k], Lcrit);
      }

      // comput Gh with smoothing function G, following Warner
      Gh[k] = Min(Ghmax, -0.5*(NN[k]+NN[k+1])*pow(l[i][k]/q[i][k],2));
      Gh[k] = Min(Gh[k], Gh[k]-pow((Gh[k]-Ghcrit),2)/(Gh[k]+Ghmax-2.*Ghcrit)); // different from Burchard and peterson 1999 by a constant
      Gh[k] = Max(Gh[k], Ghmin);   
      Gm = P[i][k]/(phys->nu_tv[i][k])*l[i][k]*l[i][k]/Qold[i][k]/2;

      //Apply stability functions 
      StabilityFunctions(&Sm,&Sh,Gh[k],Gm);

      // set minimal, disabled
      nuT[i][k]=Sm*q[i][k]*l[i][k];
      kappaT[i][k]=Sh*q[i][k]*l[i][k];

      //Set an upper bound
      // nuT[i][k]=Min(Sm*q[i][k]*l[i][k], 1.);
      // kappaT[i][k]=Min(Sh*q[i][k]*l[i][k], 1.);

      // large nut is calculated for very poorly resolved offshore locations, snohomish 
      if(nuT[i][k]>1 || nuT[i][k]!=nuT[i][k]) {
        nuT[i][k]=1;
        kappaT[i][k]=Sh/Sm*nuT[i][k];	
      }


    }

    /*
      //check top layer
     if(grid->dzz[i][grid->ctop[i]]<0.5*grid->dz[grid->ctop[i]]	 && (grid->ctop[i]!=grid->ctopold[i]) ){
       nc1=-1;
       //       printf("YYY i=%d u=%e %e %e, nut=%e nut+1=%e\n",i,phys->uc[i][grid->ctop[i]],phys->uc[i][grid->ctop[i]+1],
       //      phys->uc[i][grid->ctop[i]+2], nuT[i][grid->ctop[i]],nuT[i][grid->ctop[i]+1]);

       for(nf=0;nf<NFACES;nf++) 
	   if(grid->neigh[i*NFACES+nf]>-1 && grid->ctop[grid->neigh[i*NFACES+nf]]<=grid->ctop[i])
	     nc1=grid->neigh[i*NFACES+nf];      
       
       if(nc1>-1) {
	 for(k=grid->ctop[i];k<grid->Nk[i];k++)
	   q[i][k]=l[i][k]=kappaT[i][k]=nuT[i][k]=0;
       
	 for(  k=grid->ctop[i];k<grid->Nk[i];k++){
	   q[i][k]=q[nc1][k];
	   l[i][k]=l[nc1][k];
	   nuT[i][k]=nuT[nc1][k];
	   kappaT[i][k]=kappaT[nc1][k];
	 }
       }
       
     
     }*/

    for(k=0;k<grid->ctop[i];k++) 
      P[i][k]=B[i][k]=D[i][k]=nuT[i][k]=kappaT[i][k]=l[i][k]=q[i][k]=0.0;

  }

}

/*
 * Function: Stability Functions
 * Usage:  StabilityFunctions(&Sm,&Sh,b[k],A1,A2,B1,B2,C1);
 * --------------------------------------------------------
 * Computes the Stability functions of Blumberg et al. (1992) and 
 * places them into Sm and Sh.
 *
 */
static void StabilityFunctions(REAL *Sm, REAL *Sh, REAL Gh, REAL Gm) {

 // Define the stability function constants according 
  // Kantha & Clayson 1994 quasi-equilibrium closure
  /*
  const REAL A1 =   0.92;
  const REAL A2 =   0.74;
  const REAL B1 =   16.6;
  const REAL B2 =   10.1;
  const REAL C2 =   0.7;
  const REAL C3 =   0.2;
  *Sh = A2*(1-6*A1/B1)/(1-3*A2*Gh*(6*A1+B2*(1-C3)));
  *Sm = (pow(B1,-1.0/3.0)+(*Sh)*Gh*(18*A1*A1+9*A1*A2*(1-C2)))/
    (1-9*A1*A2*Gh);
*/  
  /*
  // Specification of constants Blumberg 1992
  const A1 = 0.92;
  const A2 = 0.74;
  const B1 = 16.6;
  const B2 = 10.1;
  const C1 = 0.08;
  *Sm = (pow(B1,-1.0/3.0)-A1*A2*Gh*((B2-3*A2)*(1-6*A1/B1)-3*C1*(B2+6*A1)))/
    ((1-3*A2*Gh*(6*A1+B2))*(1-9*A1*A2*Gh));
  *Sh = A2*(1-6*A1/B1)/(1-3*A2*Gh*(6*A1+B2));
  */

  // Canuto et al (2001)
  // this returns c*Sm as defined in Warner(2005)
  const REAL cmu0=0.5270;
  const REAL L1=0.107, L2=0.0032, L3=0.0864, L4=0.12, L5=11.9,
    L6=0.4, L7=0.0, L8=0.48;
  REAL s0 = 1.5*L1*L5*L5;
  REAL s1 = -L4*(L6+L7)+2*L4*L5*(L1-L2/3-L3)+1.5*L1*L5*L8;
  REAL s2 = -0.375*L1*(L6*L6-L7*L7);
  REAL s4 = 2.0*L5;
  REAL s5 = 2.0*L4;
  REAL s6 = (2.0/3.0)*L5*(3*L3*L3-L2*L2)-0.5*L5*L1*(3*L3-L2)+0.75*L1*(L6-L7);
  REAL b0 = 3.0*L5*L5;
  REAL b1 = L5*(7*L4+3*L8);
  REAL b2 = L5*L5*(3*L3*L3-L2*L2)-0.75*(L6*L6-L7*L7);
  REAL b3 = L4*(4*L4+3*L8);
  REAL b5 = 0.25*(L2*L2-3*L3*L3)*(L6*L6-L7*L7);
  REAL b4 = L4*(L2*L6-3*L3*L7-L5*(L2*L2-L3*L3))+L5*L8*(3*L3*L3-L2*L2);
  REAL f6= 8.0*pow(cmu0,-6);
  REAL cff;

  Gm = Min((b0/f6-b1*Gh+b3*f6*Gh*Gh)/(b2-b4*f6*Gh),Gm);
  cff = b0-b1*f6*Gh+b2*f6*Gm+b3*f6*f6*Gh*Gh-b4*f6*f6*Gh*Gm+b5*f6*f6*Gm*Gm;
  *Sm = Max(0.0, sqrt(2.)*pow(cmu0,-3)*(s0-s1*f6*Gh+s2*f6*Gm)/cff);
  *Sh = Max(0.0, sqrt(2.)*pow(cmu0,-3)*(s4-s5*f6*Gh+s6*f6*Gm)/cff); 
    
  
  // a different way of calculating CA stability functions, reference Burchard
  /*  
     Gh *= -pow(cmu0,-6)*2;
     Gm *= pow(cmu0,-6)*2;
     *Sm =  1/sqrt(2.)*pow(cmu0,-3)*(0.1070+0.0174*Gh-0.00012*Gm)/
     (1+0.2555*Gh+0.02872*Gm+0.008677*Gh*Gh+0.005222*Gh*Gm-0.0000337*Gm*Gm);
  
     *Sh =  1/sqrt(2.)*pow(cmu0,-3)*(0.1120+0.004519*Gh+0.00088*Gm)/
     (1+0.2555*Gh+0.02872*Gm+0.008677*Gh*Gh+0.005222*Gh*Gm-0.0000337*Gm*Gm);
  */
}
