//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file binary_disk.cpp
//  \brief Initializes accretion disk around binary in cartesian cooordinates(for now)

// C++ headers
#include <iostream>   // endl
#include <fstream>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>      // sqrt
#include <algorithm>  // min
#include <cstdlib>    // srand
#include <cfloat>     // FLT_MIN

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../bvals/bvals.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro_diffusion/hydro_diffusion.hpp"
#include "../utils/units.hpp"

//user boundary conditions
static void DiodeInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
     FaceField &b, Real time, Real dt, 
     int il, int iu, int jl, int ju, int kl, int ku, int ngh);
static void DiodeOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
     FaceField &b, Real time, Real dt, 
     int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void solve_u(const Real t, const Real dto, const Real e, Real *ua);
void compute_f(Real ua, const Real e, Real *f_ptr);
void compute_star_loc(const Real dto, const Real t, const Real e,
            const Real mf, Real *f_ptr, Real *u_ptr, Real *r_ptr);

static void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k);
static Real DenProfileCyl(const Real rad, const Real phi, const Real z);
static Real PoverR(const Real den);
static void VelProfileCyl(const Real den, const Real rad, const Real phi, const Real z,
                          Real &v1, Real &v2, Real &v3, const int flag);
static Real RampFunc(const Real rad, const Real phi, const Real z,
            const Real v1, const Real v2, const Real v3);
static Real TempProfileCyl(const Real den);
static Real SqSoundSpeed(const Real rad, const Real phi, const Real z);

//user defined src term
void Cooling(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
             AthenaArray<Real> &cons);
void Binary(MeshBlock *pmb, const Real time, const Real dt,
            const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
            const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
            AthenaArray<Real> &cons_scalar);

//user defined hst output
static Real hst_accm(MeshBlock *pmb, int iout);

// problem parameters which are useful to make global to this file
static Real semia,ecc,qrat,mu,incli,argp;
static Real gm0, r0, rho0, dslope, p0_over_r0, pslope, gamma_gas;
static Real dfloor,pfloor;
static Real rsoft,rsink,rin,rout,rbuf1,rbuf2;
static Real tsink; // mass removing time scale
static Real tbin;  // binary period 2pi/(GM/a^3)^1/2
static Real beta_th,alpha;
static Real sinkpower;
// parameters for computing binary orbit
static int itermax=100;
static Real dueps=1e-6;
//instantaneous binary separation, ecc and true anormaly in its orbit plane
static Real rsep,uanorm,fanorm;
//real-time binary postion in cartesian
static Real x1s,x2s,x3s,x1p,x2p,x3p;

//num of shells for history output
static int ihst1d;
static int nvar=14; //output five major variables
static Real hst1d_tstart;
static int ncbd,ncsd1,ncsd2;



// debug for binary orbit
static FILE * pf;
static Real tstart;


//Units
static Real AU = 1.495978707e13; //cm
static Real muH = 2.3;
static Real lunit = 50.*AU;
static Real dunit = muH*Constants::mH;
static Real vunit = 1.; //Constants::kms;
static Units my_unit(dunit, lunit, vunit);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // Get parameters for binary potential
  // gm0:      total mass of binary
  // qrat:     mass ratio secondary/primary (<=1 always)
  // mu:       mass fraction of the primary defined as 1/(1+qrat)
  // semia:    semi-major axis(also the initial separation)
  // incli:    inclination(angle between disk and binary orbit planes)
  // ecc:      eccentrcity
  // phi0_2:   initial phase angle of the secondary (0 as default)
  // phi0_1:   initial phase angle of the primary(PI as default)
  // argp:     argument of periapse, angle between the axis of periapse and line of nodes
  // lon:      angle between line of nodes and reference X-axis.Set to be zero always.
  gm0  = pin->GetOrAddReal("problem","GMb",1.0);
  qrat = pin->GetOrAddReal("problem","ratio",1.0);
  mu = 1.0/(1.0+qrat);
  ecc  = pin->GetOrAddReal("problem","ecc",0.0);
  incli  = pin->GetOrAddReal("problem","inclination",0.0);
  argp   = pin->GetOrAddReal("problem","periapse",0.0);
  semia  = pin->GetOrAddReal("problem","semi_major",1.0);
  tbin   = 2.0*PI/sqrt(gm0/semia/SQR(semia));
  //phi0_2 = pin->GetOrAddReal("problem","phi0",0.0);
  //phi0_1 = phi0_2+PI;
  rsoft = pin->GetOrAddReal("problem","rsoft",0.05);
  rsink = pin->GetOrAddReal("problem","rsink",0.05);
  rsoft *= semia;
  rsink *= semia;
  beta_th = pin->GetOrAddReal("problem","beta_th",0.0);
  sinkpower = pin->GetOrAddReal("problem","sink_power",1.5);

  // need to initialize binary position or from restart point
  Real dto = 0.01;
  Real torb = fmod(time,tbin);
  rsep= semia; uanorm=0.0; fanorm=0.0; //separation, ecc and true anormaly
  compute_star_loc(dto,torb,ecc,mu,&fanorm,&uanorm,&rsep);
  std::cout << "[init_cond]: dt,t,e,mu,f,u,r = "<< dto << " " << torb << " " << ecc <<" "
  << mu << " " << fanorm << " " << uanorm << " " << rsep << " " << std::endl;
  Real rcos = rsep*cos(fanorm+argp);
  Real rsin = rsep*sin(fanorm+argp);
  Real rsincos = rsin*cos(incli);
  Real rsinsin = rsin*sin(incli);
  x1p = - (1.0-mu)*rcos;
  x2p = - (1.0-mu)*rsincos;
  x3p = (1.0-mu)*rsinsin;
  x1s = mu*rcos;
  x2s = mu*rsincos;
  x3s =-mu*rsinsin;

  // Get parameters for circumbinary disk
  r0 = semia; // set length unit
  rin = pin->GetOrAddReal("problem","rin",0.0);
  rout = pin->GetOrAddReal("problem","rout",5.0);
  tsink = pin->GetOrAddReal("problem","tsink",100.0);
  rbuf1 = pin->GetOrAddReal("problem","rbuf1",8.0);
  rbuf2 = pin->GetOrAddReal("problem","rbuf2",10.0);

  // Get parameters for initial density and velocity
  rho0 = pin->GetReal("problem","rho0");
  dslope = pin->GetOrAddReal("problem","dslope",0.0);

  // Get parameters of initial pressure and cooling parameters
  if(NON_BAROTROPIC_EOS){
    p0_over_r0 = pin->GetOrAddReal("problem","p0_over_r0",0.0025);
    pslope = pin->GetOrAddReal("problem","pslope",0.0);
    gamma_gas = pin->GetReal("hydro","gamma");
  }else{
    p0_over_r0=SQR(pin->GetReal("hydro","iso_sound_speed"));
  }
  //dfloor=pin->GetOrAddReal("hydro","dfloor",(1024*(FLT_MIN)));
  dfloor=pin->GetOrAddReal("hydro","dfloor",(1024*(FLT_MIN)));
  pfloor=pin->GetOrAddReal("hydro","pfloor",(1024*(FLT_MIN)));

  int inu = pin->GetOrAddReal("problem","inu",1);
  if (inu != 0) alpha = pin->GetOrAddReal("problem","nuiso",0.1);
  else          alpha = pin->GetOrAddReal("problem","nuiso",2.5e-4)/p0_over_r0;


  // Enroll user-defined physical source terms
  if (qrat != 0.0) EnrollUserExplicitSourceFunction(Binary);
  // enroll boundary value function pointers
  EnrollUserBoundaryFunction(BoundaryFace::inner_x3, DiodeInnerX3);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x3, DiodeOuterX3);

  //Store accretion rate
  AllocateUserHistoryOutput(2);
  EnrollUserHistoryOutput(0, hst_accm, "accm1");
  EnrollUserHistoryOutput(1, hst_accm, "accm2");


  // debug binary orbit
  pf = fopen ("binary_orbit.tab","w");
  tstart = torb; //first dump

  return;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
        //Create memory to store accretion rate
        AllocateRealUserMeshBlockDataField(2);
        ruser_meshblock_data[0].NewAthenaArray(2);
        return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Initializes Keplerian accretion disk.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real rad, phi, z;
  Real v1, v2, v3;

  //  Initialize density and momenta
  for(int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie; ++i) {
      GetCylCoord(pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
      // 1) pure power law with truncation radii rin rout
      //if (rad >= rin && rad <= rout)
      //if (rad >= rin)
      //  phydro->u(IDN,k,j,i) = DenProfileCyl(rad,phi,z);
      //else
      //  phydro->u(IDN,k,j,i) = dfloor;
      // 2) power law profile with exponential taper interior to rin
      Real den0 = DenProfileCyl(rad,phi,z)*exp(-pow((rad/rin),-2.0));
      den0 = (den0 > dfloor) ?  den0 : dfloor;
      phydro->u(IDN,k,j,i) = den0;
      
      VelProfileCyl(den0, rad,phi,z,v1,v2,v3,1);
      
      phydro->u(IM1,k,j,i) = den0*v1;
      phydro->u(IM2,k,j,i) = den0*v2;
      phydro->u(IM3,k,j,i) = den0*v3;
      if (NON_BAROTROPIC_EOS){
        Real p_over_r = PoverR(den0);
        //Real press =  ( p_over_r*phydro->u(IDN,k,j,i)>pfloor ) ? p_over_r*phydro->u(IDN,k,j,i) : pfloor;
        phydro->u(IPR,k,j,i) = p_over_r*phydro->u(IDN,k,j,i);
        phydro->u(IEN,k,j,i) = p_over_r*phydro->u(IDN,k,j,i)/(gamma_gas - 1.0);
        phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                   + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
	      }
    }
   }
  }

  return;
}

//----------------------------------------------------------------------------------------
//!\f transform to cylindrical coordinate

static void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k)
{
  rad=sqrt(SQR(pco->x1v(i))+SQR(pco->x2v(j)));
  phi=atan2(pco->x2v(j),pco->x1v(i));
  z=pco->x3v(k);
  return;
}

//----------------------------------------------------------------------------------------
//! \f  computes density in cylindrical coordinates

static Real DenProfileCyl(const Real rad, const Real phi, const Real z)
{
  Real den;
  Real p_over_r = p0_over_r0;
  if (NON_BAROTROPIC_EOS) p_over_r = SqSoundSpeed(rad, phi, z);
  Real denmid = rho0*pow(rad/r0,dslope);
  Real dentem = denmid*exp(gm0/p_over_r*(1./sqrt(SQR(rad)+SQR(rsoft)+SQR(z))-1./sqrt(SQR(rad)+SQR(rsoft))));
  den = dentem;
  return std::max(den,dfloor);
}

//----------------------------------------------------------------------------------------
//! \f  computes pressure/density in cylindrical coordinates

static Real PoverR(const Real den0)
{
  Real poverr;
  // 1) simple power law with index pslope
  //poverr = p0_over_r0*pow(rad/r0, pslope);
  // 2) P/rho = Cs(r,phi,z)^2/gamma
  //          = (h/r)^2*(GM1/r1+GM2/r2)/gamma
  //Real x1 = rad*cos(phi);
  //Real x2 = rad*sin(phi);
  //Real rad1 = sqrt(SQR(x1-x1p)+SQR(x2-x2p)+SQR(z-x3p)+SQR(rsoft));
  //Real rad2 = sqrt(SQR(x1-x1s)+SQR(x2-x2s)+SQR(z-x3s)+SQR(rsoft));
  //poverr = p0_over_r0*gm0*(mu/rad1+(1.0-mu)/rad2)/gamma_gas; 
  
  //insert temperature contribution
  poverr = TempProfileCyl(den0);
  return poverr;
}


//----------------------------------------------------------------------------------------
//! \f  computes rotational velocity in cylindrical coordinates

static void VelProfileCyl(const Real den, const Real rad, const Real phi, const Real z,
                          Real &v1, Real &v2, Real &v3, const int flag)
{
  Real p_over_r = p0_over_r0;
  if (NON_BAROTROPIC_EOS) p_over_r = PoverR(den);
  Real rad1 = sqrt(SQR(rad)+SQR(rsoft));
  Real rad2 = sqrt(SQR(rad)+SQR(rsoft)+SQR(z));
  Real qrpl = 0.75*SQR(semia/rad1)*qrat/SQR(1.0+qrat);
  Real vel = sqrt((dslope+pslope)*p_over_r + gm0*SQR(rad/rad1)/rad1*(1.0+qrpl)
            +pslope*gm0*(1.0/rad1 -1.0/rad2));

  v1 = -vel*sin(phi);
  v2 = vel*cos(phi);
  v3 = 0.0;
  // correct with radial drift due to accretion
  if (NON_BAROTROPIC_EOS && flag) {
    vel = 1.5*p_over_r*gamma_gas*alpha/sqrt(gm0/sqrt(SQR(rad)+SQR(rsoft)));
    v1 -= vel*cos(phi);
    v2 -= vel*sin(phi);
  }

  if (rad < rsoft) {
    v1 = 0.0;
    v2 = 0.0;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \f  computes ramp function R/\tau
// R(x) = a*x^2 + b*x; where a=1/(r2^2-r1r2), and b=-a*r1
// tau  = 2PI/Omega(r1)

static Real RampFunc(const Real rad, const Real phi, const Real z,
 const Real v1, const Real v2, const Real v3)
{
  Real ramp,tau;
  if (rad >= rbuf2) {
    ramp = 10.0;
  } else if (rad >= rbuf1) {
    Real fac = 1.0/(SQR(rbuf2)-rbuf1*rbuf2);
    ramp = fac*SQR(rad)-fac*rbuf1*rad;
  } else {
    ramp = 0.0;
  }
  tau = 2.0*PI/sqrt(gm0/rbuf1/SQR(rbuf1));
  return ramp/tau;
}

//----------------------------------------------------------------------------------------
//! \f  computes temperature profile
// [T]=Energy
//need to fix scale here

static Real TempProfileCyl(const Real den0)
{
  Real temp;
  Real AU = 1.495978707e13; //cm
  Real Mbin = 3.; //M_sol
  Real semiabin = 50.; //AU
  Real time = sqrt(pow(semiabin*AU,3)/(Mbin*Constants::Msun*Constants::G));
  Real densref = Mbin*Constants::Msun/pow(semiabin*AU,3); //g/cm^3
  Real T0 = 15.; //K

  //Real den0 = DenProfileCyl(rad,phi,z)*exp(-pow((rad/rin),-2.0));
  //den0 = (den0 > dfloor) ?  den0 : dfloor;
  //Real densphys = den0*densref;
  Real densphys = (den0 > dfloor) ?  den0*densref : dfloor;

  if (densphys<1.e-12){
    temp = T0+1.5*1.e13*densphys;
  } else if (1.e-12 <= densphys && densphys < 1.e-11){
    temp = (T0+15.)*pow(densphys*1.e12, 0.6);
  } else if (1.e-11 < densphys && densphys <= 3.e-9){
    temp = pow(10.,0.6)*(T0+15.)*pow(densphys*1.e11, 0.44);
  }

  //temp *= 1./my_unit.Temperature * (my_unit.Pressure/my_unit.Density);
  temp *= Constants::kB/(Constants::mH*2.3);
  
  //back to code units

  temp *= pow(time/(semiabin*AU),2);

  return temp;
}

//----------------------------------------------------------------------------------------
//! \f  computes initial sound speed in cylindrical coordinates

static Real SqSoundSpeed(const Real rad, const Real phi, const Real z)
{
  Real poverr;
  // 1) simple power law with index pslope
  //poverr = p0_over_r0*pow(rad/r0, pslope);
  // 2) P/rho = Cs(r,phi,z)^2/gamma
  //          = (h/r)^2*(GM1/r1+GM2/r2)/gamma
  Real x1 = rad*cos(phi);
  Real x2 = rad*sin(phi);
  Real rad1 = sqrt(SQR(x1-x1p)+SQR(x2-x2p)+SQR(z-x3p)+SQR(rsoft));
  Real rad2 = sqrt(SQR(x1-x1s)+SQR(x2-x2s)+SQR(z-x3s)+SQR(rsoft));
  poverr = p0_over_r0*gm0*(mu/rad1+(1.0-mu)/rad2)/gamma_gas; 
  
  //insert temperature contribution
  //poverr = TempProfileCyl(rad, phi, z);
  return poverr;
}

void MeshBlock::UserWorkInLoop(void)
{

  // Clear hist data after previous loop
  if ((ihst1d==1) && (pmy_mesh->time>=hst1d_tstart)) {
    for(int n=0; n<nvar; n++) {
      for(int i=0; i<ncbd;  i++) ruser_meshblock_data[4](n,i)=0;
      for(int i=0; i<ncsd1; i++) ruser_meshblock_data[5](n,i)=0;
      for(int i=0; i<ncsd2; i++) ruser_meshblock_data[6](n,i)=0;
    }
  }
  // estimate the velocity of the binary needed for calc the hist data
  Real sini = sin(incli);
  Real cosi = cos(incli);
  Real sinf = sin(fanorm+argp);
  Real cosf = cos(fanorm+argp);
  Real v0 = sqrt(gm0/semia/(1.0-SQR(ecc)));
  Real ve = v0*(ecc-cosf);
  Real v1p = (1.0-mu)*v0*sinf;
  Real v2p = (1.0-mu)*ve*cosi;
  Real v3p = -(1.0-mu)*ve*sini;
  Real v1s = -mu*v0*sinf;
  Real v2s = -mu*ve*cosi;
  Real v3s = mu*ve*sini;

  // DEBUG: binary orbit
  if (gid == 0 && pmy_mesh->time >= tstart && pmy_mesh->time <= (2.0*PI)) {
    tstart += 0.01*tbin;
    fprintf(pf,"%19.15e %19.15e %19.15e %19.15e %19.15e %19.15e %19.15e \
      %19.15e %19.15e %19.15e %19.15e %19.15e %19.15e %19.15e %19.15e \
      %19.15e\n",pmy_mesh->time, fanorm,uanorm,
      rsep,x1s,x2s,x3s,x1p,x2p,x3p,v1s,v2s,v3s,v1p,v2p,v3p);
  }
  //// Prepare index bounds
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js;
  int ju = je;
  if (block_size.nx2 > 1) { //2D
    jl -= (NGHOST);
    ju += (NGHOST);
  }
  int kl = ks;
  int ku = ke;
  if (block_size.nx3 > 1) { //3D
    kl -= (NGHOST);
    ku += (NGHOST);
  }
  //int il = is;
  //int iu = ie;
  //int jl = js;
  //int ju = je;
  //int kl = ks;
  //int ku = ke;
  for(int k=kl; k<=ku; k++) {
    for(int j=jl; j<=ju; j++) {
      for(int i=il; i<=iu; i++) {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        Real radp = sqrt(SQR(x1-x1p)+SQR(x2-x2p)+SQR(x3-x3p));
        Real rads = sqrt(SQR(x1-x1s)+SQR(x2-x2s)+SQR(x3-x3s));
        Real rad, phi, z;
        GetCylCoord(pcoord,rad,phi,z,i,j,k);

        Real& u_d  = phydro->u(IDN,k,j,i);
        Real& u_m1 = phydro->u(IM1,k,j,i);
        Real& u_m2 = phydro->u(IM2,k,j,i);
        Real& u_m3 = phydro->u(IM3,k,j,i);
        Real& u_e  = phydro->u(IEN,k,j,i);

        Real& w_d  = phydro->w(IDN,k,j,i);
        Real& w_vx = phydro->w(IVX,k,j,i);
        Real& w_vy = phydro->w(IVY,k,j,i);
        Real& w_vz = phydro->w(IVZ,k,j,i);
        Real& w_p  = phydro->w(IEN,k,j,i);

        Real dt = pmy_mesh->dt;

        // apply sink cells condition within rsink( ususally >= rsoft)
        if ((radp <= rsink)||(rads <=rsink)) {
          Real u_d0 = u_d;
          Real trm = std::max(dt,tsink*pow(std::min(radp,rads)/rsink, sinkpower));
          // steeper than alpha model
          //Real trm = std::max(dt,tsink*pow(std::min(radp,rads)/rsink, 3.0));
          u_d -= dt*u_d/trm;
          // apply density floor, without changing momentum or energy
          u_d = (u_d > dfloor) ?  u_d : dfloor;
          w_d = u_d;
          // adjust the momentum
          Real fac = u_d/u_d0;
          u_m1 *= fac;
          u_m2 *= fac;
          u_m3 *= fac;
          Real di = 1.0/u_d;
          w_vx = u_m1*di;
          w_vy = u_m2*di;
          w_vz = u_m3*di;
          
          // store the accretion rate
          if ((i<=ie && i>=is) && (j<=je && j>=js) && (k<=ke && k>=ks)){
            Real vol = pcoord->GetCellVolume(k,j,i);
            Real& accm1 = ruser_meshblock_data[0](0);
            Real& accm2 = ruser_meshblock_data[0](1);
            if (radp <= rsink) accm1 += (u_d0-u_d)*vol;
            else accm2 += (u_d0-u_d)*vol;
            //if (radp <= rsink) accm1 += (u_d0-u_d)*vol/dt;
            //else accm2 += (u_d0-u_d)*vol/dt;
          }

        }
        // apply wave-killing zone within [rbuf1,rbuf2] to quench m=4 mode
        if (rad >= rbuf1) {
          Real v1, v2, v3;
          Real den0 = DenProfileCyl(rad,phi,z);
          VelProfileCyl(den0,rad,phi,z,v1,v2,v3,0);
          Real ramp = RampFunc(rad,phi,z,v1,v2,v3);
          u_d  -= dt*(u_d-den0)*ramp;
          u_m1 -= dt*(u_m1-den0*v1)*ramp;
          u_m2 -= dt*(u_m2-den0*v2)*ramp;
          u_m3 -= dt*(u_m3-den0*v3)*ramp;
          Real di = 1.0/u_d;
          w_vx = u_m1*di;
          w_vy = u_m2*di;
          w_vz = u_m3*di;

          //if (NON_BAROTROPIC_EOS){
          //  Real p_over_r = PoverR(rad,phi,z);
          //  Real pres0 = p_over_r*den0;
          //  w_p -= dt*(w_p-pres0)*ramp;
          //  w_p = (w_p > pfloor) ?  w_p : pfloor;
          //  Real ke = 0.5*di*(SQR(u_m1) + SQR(u_m2) + SQR(u_m3));
          //  u_e = w_p/(gamma_gas-1.0)+ke;
          //}
        }
        // apply velocity cap
        u_d = (u_d > dfloor) ?  u_d : dfloor;
        w_d = u_d;
        Real vcap = 50.0*sqrt(p0_over_r0);
        w_vz = (fabs(w_vz) <= vcap) ? w_vz : 0.0;
        u_m3 = u_d*w_vz;
        w_vx = (fabs(w_vx) <= vcap) ? w_vx : 0.0;
        u_m1 = u_d*w_vx;
        w_vy = (fabs(w_vy) <= vcap) ? w_vy : 0.0;
        u_m2 = u_d*w_vy;

        // apply extremely short cooling
        if (NON_BAROTROPIC_EOS){
          Real pres0 = u_d*PoverR(u_d);
          w_p = pres0;
          w_p = (w_p > pfloor) ?  w_p : pfloor;
          Real di = 1.0/u_d;
          Real gmi = 1.0/(gamma_gas-1.0);
          Real ke = 0.5*di*(SQR(u_m1) + SQR(u_m2) + SQR(u_m3));
          u_e = w_p*gmi+ke;
        }

      }
    }
  }
  if ((ihst1d==1) && (pmy_mesh->time>=hst1d_tstart)) {
    for(int k=ks; k<=ke; k++) {
      for(int j=js; j<=je; j++) {
        for(int i=is; i<=ie; i++) {
          // calculate the shell integrated binary torque,tilt angle, and precession angle
          //CalculateTorque();
          Real x1 = pcoord->x1v(i);
          Real x2 = pcoord->x2v(j);
          Real x3 = pcoord->x3v(k);
          Real radp = sqrt(SQR(x1-x1p)+SQR(x2-x2p)+SQR(x3-x3p));
          Real rads = sqrt(SQR(x1-x1s)+SQR(x2-x2s)+SQR(x3-x3s));
          Real rad, phi, z;
          GetCylCoord(pcoord,rad,phi,z,i,j,k);

          Real& u_d  = phydro->u(IDN,k,j,i);
          Real& u_m1 = phydro->u(IM1,k,j,i);
          Real& u_m2 = phydro->u(IM2,k,j,i);
          Real& u_m3 = phydro->u(IM3,k,j,i);
          Real& u_e  = phydro->u(IEN,k,j,i);

          Real& w_d  = phydro->w(IDN,k,j,i);
          Real& w_vx = phydro->w(IVX,k,j,i);
          Real& w_vy = phydro->w(IVY,k,j,i);
          Real& w_vz = phydro->w(IVZ,k,j,i);
          Real& w_p  = phydro->w(IEN,k,j,i);

          Real dt = pmy_mesh->dt;

          AthenaArray<Real>& rcbd = ruser_meshblock_data[1];
          AthenaArray<Real>& rcsd1 = ruser_meshblock_data[2];
          AthenaArray<Real>& rcsd2 = ruser_meshblock_data[3];
          AthenaArray<Real>& data_cbd = ruser_meshblock_data[4];
          AthenaArray<Real>& data_csd1 = ruser_meshblock_data[5];
          AthenaArray<Real>& data_csd2 = ruser_meshblock_data[6];
          Real rsph = sqrt(SQR(x1)+SQR(x2)+SQR(x3));

          // define a pointer to hydro diffusion class
          HydroDiffusion *phdif = &phydro->hdif;
          // viscous stress tensor T_{ij}. For example,
          //      x1flx(IM1,k,j,i) = T_xx (k,j,i)
          //      x1flx(IM2,k,j,i) = T_xy (k,j,i)
          //      x1flx(IM3,k,j,i) = T_xz (k,j,i)
          //      x2flx(IM1,k,j,i) = T_yx (k,j,i)
          AthenaArray<Real> &x1flx=phdif->visflx[X1DIR]; //x1flx([IM1,IM2,IM3],k,j,i)
          AthenaArray<Real> &x2flx=phdif->visflx[X2DIR]; //x2flx([IM1,IM2,IM3],k,j,i)
          AthenaArray<Real> &x3flx=phdif->visflx[X3DIR]; //x3flx([IM1,IM2,IM3],k,j,i)

          for (int n=0; n<ncbd; n++) {
            Real rc = rcbd(0,n);
            Real rl = rc-0.5*rcbd(1,n);
            Real rr = rc+0.5*rcbd(1,n);
            if (rsph>=rl && rsph<rr) {
              Real a1 = (1.0-mu)*rsep;
              Real a2 = mu*rsep;
              Real phiprime = atan2(x2*cos(incli)-x3*sin(incli),x1);
              Real sinphi = sin(fanorm+argp-phiprime);
              Real vol = pcoord->GetCellVolume(k,j,i);
              Real rsphp = pow(SQR(x1-x1p)+SQR(x2-x2p)+SQR(x3-x3p)+SQR(rsoft),1.5); //r1^3
              Real rsphs = pow(SQR(x1-x1s)+SQR(x2-x2s)+SQR(x3-x3s)+SQR(rsoft),1.5); //r2^3
              Real rsph_p = sqrt(SQR(x1-x1p)+SQR(x2-x2p)+SQR(x3-x3p)+SQR(rsoft)); //r1
              Real rsph_s = sqrt(SQR(x1-x1s)+SQR(x2-x2s)+SQR(x3-x3s)+SQR(rsoft)); //r2
              // Z-TORQUE AVERAGE OVER SHELL WIDTH
              Real dphidx = 0.5*((x1-x1p)/rsphp+(x1-x1s)/rsphs);
              Real dphidy = 0.5*((x2-x2p)/rsphp+(x2-x2s)/rsphs);
              Real dphidz = 0.5*((x3-x3p)/rsphp+(x3-x3s)/rsphs);
              data_cbd(0,n) += u_d*vol/rcbd(1,n)*(x1*dphidy-x2*dphidx);
              // MDOT AVERAGE OVER SHELL WIDTH
              Real vrad = (x1*w_vx+x2*w_vy+x3*w_vz)/rsph;
              data_cbd(1,n) += u_d*vrad*vol/rcbd(1,n);
              // SHELL INTEGRATED ANGULAR MOMENTUM \rho r x p
              Real rxp1 = x2*w_vz-x3*w_vy;
              Real rxp2 = x1*w_vz-x3*w_vx;
              Real rxp3 = x1*w_vy-x2*w_vx;
              data_cbd(2,n) += u_d*rxp1*vol;
              data_cbd(3,n) += u_d*rxp2*vol;
              data_cbd(4,n) += u_d*rxp3*vol;
              // SHELL INTEGRATED SURFACE DENSITY (sigma*r*dphi)
              data_cbd(5,n) += u_d*vol/rcbd(1,n);
              // X, Y GRAV TORQUE
              data_cbd(6,n) += u_d*vol/rcbd(1,n)*(x2*dphidz-x3*dphidy);
              data_cbd(7,n) += u_d*vol/rcbd(1,n)*(x3*dphidx-x1*dphidz);
              // TORQUE FROM ADVECTION (X,Y,Z)
              data_cbd(8,n) += u_d*vol/rcbd(1,n)/rsph*(x3*w_vy-x2*w_vz)*(x1*w_vx+x2*w_vy+x3*w_vz);
              data_cbd(9,n) += u_d*vol/rcbd(1,n)/rsph*(x1*w_vz-x3*w_vx)*(x1*w_vx+x2*w_vy+x3*w_vz);
              data_cbd(10,n) += u_d*vol/rcbd(1,n)/rsph*(x2*w_vx-x1*w_vy)*(x1*w_vx+x2*w_vy+x3*w_vz);
              // VISCOUS TORQUE (X,Y,Z)
              Real TdotRhat_x = (x1flx(IM1,k,j,i)*x1+x1flx(IM2,k,j,i)*x2+x1flx(IM3,k,j,i)*x3)/rsph;
              Real TdotRhat_y = (x2flx(IM1,k,j,i)*x1+x2flx(IM2,k,j,i)*x2+x2flx(IM3,k,j,i)*x3)/rsph;
              Real TdotRhat_z = (x3flx(IM1,k,j,i)*x1+x3flx(IM2,k,j,i)*x2+x3flx(IM3,k,j,i)*x3)/rsph;
              data_cbd(11,n) += vol/rcbd(1,n)*(x3*TdotRhat_y-x2*TdotRhat_z);
              data_cbd(12,n) += vol/rcbd(1,n)*(x1*TdotRhat_z-x3*TdotRhat_x);
              data_cbd(13,n) += vol/rcbd(1,n)*(x2*TdotRhat_x-x1*TdotRhat_y);
            }
          }

          for (int n=0; n<ncsd1; n++) {
            Real rc = rcsd1(0,n);
            Real rl = rc-0.5*rcsd1(1,n);
            Real rr = rc+0.5*rcsd1(1,n);
            if (radp>=rl && radp<rr) {
              Real a1 = (1.0-mu)*rsep;
              Real a2 = mu*rsep;
              Real phiprime = atan2(x2*cos(incli)-x3*sin(incli),x1);
              Real sinphi = sin(fanorm+argp-phiprime);
              Real vol = pcoord->GetCellVolume(k,j,i);
              Real rsphp = pow(SQR(x1-x1p)+SQR(x2-x2p)+SQR(x3-x3p)+SQR(rsoft),1.5); //r1^3
              Real rsphs = pow(SQR(x1-x1s)+SQR(x2-x2s)+SQR(x3-x3s)+SQR(rsoft),1.5); //r2^3
              // z-torque average over shell width
              data_csd1(0,n) += u_d*rad*gm0*sinphi*(-mu*a1/rsphp+(1.0-mu)*a2/rsphs)*vol/rcsd1(1,n);
              // mdot average over shell width
              Real rrp = sqrt(SQR(x1-x1p)+SQR(x2-x2p)+SQR(x3-x3p));//\vec{r}-\vec{rp}
              Real vrad = ((x1-x1p)*(w_vx-v1p)+(x2-x2p)*(w_vy-v2p)+(x3-x3p)*(w_vz-v3p))/rrp;
              data_csd1(1,n) += u_d*vrad*vol/rcsd1(1,n);
              // shell integrated angular momentum \rho r x p
              Real rxp1 = (x2-x2p)*(w_vz-v3p)-(x3-x3p)*(w_vy-v2p);
              Real rxp2 = (x1-x1p)*(w_vz-v3p)-(x3-x3p)*(w_vx-v1p);
              Real rxp3 = (x1-x1p)*(w_vy-v2p)-(x2-x2p)*(w_vx-v1p);
              data_csd1(2,n) += u_d*rxp1*vol;
              data_csd1(3,n) += u_d*rxp2*vol;
              data_csd1(4,n) += u_d*rxp3*vol;
              // shell integrated surface density
              data_csd1(5,n) += u_d*vol/rcsd1(1,n);
              // others
              data_csd1(6,n) += 1.0;
              data_csd1(7,n) += 1.0;
              data_csd1(8,n) += 1.0;
              data_csd1(9,n) += 1.0;
              data_csd1(10,n) += 1.0;
              data_csd1(11,n) += 1.0;
              data_csd1(12,n) += 1.0;
              data_csd1(13,n) += 1.0;
            }
          }

          for (int n=0; n<ncsd2; n++) {
            Real rc = rcsd2(0,n);
            Real rl = rc-0.5*rcsd2(1,n);
            Real rr = rc+0.5*rcsd2(1,n);
            if (rads>=rl && rads<rr) {
              Real a1 = (1.0-mu)*rsep;
              Real a2 = mu*rsep;
              Real phiprime = atan2(x2*cos(incli)-x3*sin(incli),x1);
              Real sinphi = sin(fanorm+argp-phiprime);
              Real vol = pcoord->GetCellVolume(k,j,i);
              Real rsphp = pow(SQR(x1-x1p)+SQR(x2-x2p)+SQR(x3-x3p)+SQR(rsoft),1.5); //r1^3
              Real rsphs = pow(SQR(x1-x1s)+SQR(x2-x2s)+SQR(x3-x3s)+SQR(rsoft),1.5); //r2^3
              // z-torque average over shell width
              data_csd2(0,n) += u_d*rad*gm0*sinphi*(-mu*a1/rsphp+(1.0-mu)*a2/rsphs)*vol/rcsd2(1,n);
              // mdot average over shell width
              Real rrs = sqrt(SQR(x1-x1s)+SQR(x2-x2s)+SQR(x3-x3s));//\vec{r}-\vec{rs}
              Real vrad = ((x1-x1s)*(w_vx-v1s)+(x2-x2s)*(w_vy-v2s)+(x3-x3s)*(w_vz-v3s))/rrs;
              data_csd2(1,n) += u_d*vrad*vol/rcsd2(1,n);
              // shell integrated angular momentum \rho r x p
              Real rxp1 = (x2-x2s)*(w_vz-v3s)-(x3-x3s)*(w_vy-v2s);
              Real rxp2 = (x1-x1s)*(w_vz-v3s)-(x3-x3s)*(w_vx-v1s);
              Real rxp3 = (x1-x1s)*(w_vy-v2s)-(x2-x2s)*(w_vx-v1s);
              data_csd2(2,n) += u_d*rxp1*vol;
              data_csd2(3,n) += u_d*rxp2*vol;
              data_csd2(4,n) += u_d*rxp3*vol;
              // shell integrated surface density
              data_csd2(5,n) += u_d*vol/rcsd2(1,n);
              // others
              data_csd2(6,n) += 1.0;
              data_csd2(7,n) += 1.0;
              data_csd2(8,n) += 1.0;
              data_csd2(9,n) += 1.0;
              data_csd2(10,n) += 1.0;
              data_csd2(11,n) += 1.0;
              data_csd2(12,n) += 1.0;
              data_csd2(13,n) += 1.0;
            }
          }

        }
      }
    } //end second loop - calc hst1d data
  }

} //end UserWorkInLoop


void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  // DEBUG: binary orbit
  if (time > (2.0*PI)) fclose(pf);
  return;
}

//void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
//{
//  for(int k=ks; k<=ke; k++) {
//   for(int j=js; j<=je; j++) {
//      for(int i=is; i<=ie; i++) {
//        //user_out_var(0,k,j,i) = phydro->w(IPR,k,j,i)/phydro->w(IDN,k,j,i)*my_unit.Temperature;
//        user_out_var(0,k,j,i) = phydro->w(IPR,k,j,i);
//      }
//    }
//  }
//}


//----------------------------------------------------------------------------------------
//!\fn void compute_star_loc(const Real dto, const Real t, const Real e,
//                            const Real mf, Real *f_ptr, Real *u_ptr, Real *r_ptr)
// \brief computes the binary orbit based on its orbit elements (a,e,i) and time t.
// PURPOSE: Functions to compute the locations of binary star members
//   t:            time in units where 2PI is the orbit period
//   ecc:          eccentricity
//   mu:           mass fraction of the primary (>= 0.5)
//   x1p,x2p,x3p:  position of the primary in cartesian coord
//   x1s,x2s,x3s:  position of the secondary
//   dt:           time step since the last evaluation (need not be accurate)
//   uanorm:       the eccentric anomaly
//   fanorm:       the true anomaly
//   rsep:         the real-time separation

void solve_u(const Real t, const Real dto, const Real e, Real *u_ptr)
{
  Real du=1.0;
  Real ua = *u_ptr;
  ua += dto/(1.0-e*cos(ua));

  for (int i=1; i<=itermax; i++){
    du = - (ua-e*sin(ua)-t)/(1.0-e*cos(ua));
    ua += du;
    if (fabs(du) < dueps) {
      *u_ptr = ua;
      return;
    }
  }

  std::cout << "[solve_u error]: exceed the limit ITERMAX = " << itermax << std::endl;
  exit(1);
}

void compute_f(Real ua, const Real e, Real *f_ptr)
{
  Real asp = sqrt((1.0+e)/(1.0-e));
  ua = fmod(ua,2.0*PI);
  *f_ptr = 2.0*atan2(asp*sin(0.5*ua),cos(0.5*ua));
}

void compute_star_loc(const Real dto, const Real t, const Real e, const Real mf, Real *f_ptr, Real *u_ptr, Real *r_ptr)
{
    solve_u(t, dto, e, u_ptr);
    compute_f(*u_ptr, e, f_ptr);
    *r_ptr = (1.0 - SQR(e))/(1.0 + e*cos(*f_ptr));
}


//----------------------------------------------------------------------------------------
//! \fn void Binary(MeshBlock *pmb, const Real time, const Real dt,
//       const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
//  \brief Add binary force as a user defined source term
//
void Binary(MeshBlock *pmb, const Real time, const Real dt,
            const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
            const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
            AthenaArray<Real> &cons_scalar)
{
  Coordinates *pco = pmb->pcoord;
  Real src[NHYDRO];

  Real dto = dt;
  //Real torb = fmod(time,2.0*PI);
  // the following scheme works for VL2 integrator only
  Real torb = time;
  if (dt < pmb->pmy_mesh->dt) torb = fmod((time+dt),2.0*PI);//step1
  //else torb = fmod((time+0.5*dt),2.0*PI);                   //step2
  else torb = fmod((time),2.0*PI);                   //step2

  //compute the binary position in its orbital plane
  //then project it to disk plane
  compute_star_loc(dto,torb,ecc,mu,&fanorm,&uanorm,&rsep);
  rsep *=semia;
  Real rcos = rsep*cos(fanorm+argp);
  Real rsin = rsep*sin(fanorm+argp);
  Real rsincos = rsin*cos(incli);
  Real rsinsin = rsin*sin(incli);
  x1p = - (1.0-mu)*rcos;
  x2p = - (1.0-mu)*rsincos;
  x3p = (1.0-mu)*rsinsin;
  x1s = mu*rcos;
  x2s = mu*rsincos;
  x3s =-mu*rsinsin;

  //compute the binary potential according to binary
  //location
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real den = prim(IDN,k,j,i);
        Real x1  = pmb->pcoord->x1v(i);
        Real x2  = pmb->pcoord->x2v(j);
        Real x3  = pmb->pcoord->x3v(k);
        Real radp = pow(SQR(x1-x1p)+SQR(x2-x2p)+SQR(x3-x3p)+SQR(rsoft),1.5);
        //Real rads = pow(SQR(x1-x1s)+SQR(x2-x2s)+SQR(x3-x3s)+SQR(rsoft*qrat),1.5);
        Real rads = pow(SQR(x1-x1s)+SQR(x2-x2s)+SQR(x3-x3s)+SQR(rsoft),1.5);
        Real srcp = dt*den*gm0*mu/radp;
        Real srcs = dt*den*gm0*(1.0-mu)/rads;
        cons(IM1,k,j,i) -= srcp*(x1-x1p)+srcs*(x1-x1s);
        cons(IM2,k,j,i) -= srcp*(x2-x2p)+srcs*(x2-x2s);
        cons(IM3,k,j,i) -= srcp*(x3-x3p)+srcs*(x3-x3s);
        if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) -=
          srcp*((x1-x1p)*prim(IVX,k,j,i)+(x2-x2p)*prim(IVY,k,j,i)+(x3-x3p)*prim(IVZ,k,j,i))+
          srcs*((x1-x1s)*prim(IVX,k,j,i)+(x2-x2s)*prim(IVY,k,j,i)+(x3-x3s)*prim(IVZ,k,j,i));
  }}}

  // we delay the super-short cooling in UserWorkInLoop.
  //if (NON_BAROTROPIC_EOS && beta_th != 0.0)
  //  Cooling(pmb,time,dt,prim,bcc,cons);

}

// uncomment the following if the cooling is physical and evaluate at each substep
//
void Cooling(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  // apply extremely short cooling
  //Real gam = pmb->peos->GetGamma();
  //Real gam1 = pmb->peos->GetGamma()-1.0;
  Real gam1 = gamma_gas-1.0;
  Real rad, phi, z;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        GetCylCoord(pmb->pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
        Real eth = cons(IEN,k,j,i)-0.5*(SQR(cons(IM1,k,j,i))+SQR(cons(IM2,k,j,i))+
                   SQR(cons(IM3,k,j,i)))/cons(IDN,k,j,i);
        Real eth0 = cons(IDN,k,j,i)*PoverR(cons(IDN,k,j,i))/gam1;
        Real tcool= std::max(beta_th*2.0*PI/sqrt(gm0/SQR(rad)/rad),dt);
        cons(IEN,k,j,i) -= (eth-eth0)*dt/tcool;
      }
    }
  }
  return;
}

static Real hst_accm(MeshBlock *pmb, int iout)
{
        //std::cout <<"gid = "<< pmb->gid << " accm1= " << pmb->accm1 << " accm2= " << pmb->accm2 << std::endl;
        if (iout == 0) return pmb->ruser_meshblock_data[0](0);
        else return pmb->ruser_meshblock_data[0](1);
        //if (iout == 0) return accm1;
        //else return accm2;
}

void DiodeInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
     FaceField &b, Real time, Real dt, 
     int il, int iu, int jl, int ju, int kl, int ku, int ngh)
{
  int is = il + ngh;
  int ie = iu - ngh;
  int js = jl + ngh;
  int je = ju - ngh;
  int ks = kl + ngh;
  // copy hydro variables into ghost zones
  for (int k=1; k<=ngh; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        a(IDN,ks-k,j,i) = a(IDN,ks,j,i);
        a(IVX,ks-k,j,i) = a(IVX,ks,j,i);
        a(IVY,ks-k,j,i) = a(IVY,ks,j,i);
        a(IVZ,ks-k,j,i) = std::min(a(IVZ,ks,j,i),0.0);
        if(NON_BAROTROPIC_EOS) {
          a(IPR,ks-k,j,i) = a(IPR,ks,j,i);
        }
      }
  }}
  return;
}

void DiodeOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
     FaceField &b, Real time, Real dt, 
     int il, int iu, int jl, int ju, int kl, int ku, int ngh)
{
  int is = il + ngh;
  int ie = iu - ngh;
  int js = jl + ngh;
  int je = ju - ngh;
  int ks = kl + ngh;
  int ke = ku - ngh;
  // copy hydro variables into ghost zones
  for (int k=1; k<=ngh; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        a(IDN,ke+k,j,i) = a(IDN,ke,j,i);
        a(IVX,ke+k,j,i) = a(IVX,ke,j,i);
        a(IVY,ke+k,j,i) = a(IVY,ke,j,i);
        a(IVZ,ke+k,j,i) = std::max(a(IVZ,ke,j,i),0.0);
        if(NON_BAROTROPIC_EOS) {
          a(IPR,ke+k,j,i) = a(IPR,ke,j,i);
        }
      }
  }}
  return;
}


//// calc torque torqcbd as a function of rcbd
//// [4] AthenaArray(4, ncbd) torq,mdot,beta,pres
//// [5] AthenaArray(4, ncsd1) torq,mdot,beta,pres
//// [6] AthenaArray(4, ncsd2) torq,mdot,beta,pres
//void CalculateTorque()
//{
//
//}
