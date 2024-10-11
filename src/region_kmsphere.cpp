/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "region_kmsphere.h"

#include "error.h"
#include "input.h"
#include "update.h"
#include "variable.h"
#include "compute.h"
#include "modify.h"
#include "math_const.h"
#include "force.h"
#include "math_special.h"

#include <cmath>

using namespace LAMMPS_NS;
using MathConst::MY_PI;
using MathSpecial::powint;

enum { CONSTANT,VARIABLE, AUTO, MIXING };

/* ---------------------------------------------------------------------- */

RegKMSphere::RegKMSphere(LAMMPS *lmp, int narg, char **arg) :
    Region(lmp, narg, arg), deltastr(nullptr)
{
  options(narg - 20, &arg[20]);
  xc = xscale * utils::numeric(FLERR, arg[2], false, lmp);
  yc = yscale * utils::numeric(FLERR, arg[3], false, lmp);
  zc = zscale * utils::numeric(FLERR, arg[4], false, lmp);
  radius = xscale * utils::numeric(FLERR, arg[5], false, lmp);
  vradius = xscale * utils::numeric(FLERR, arg[6], false, lmp);
  CB= xscale*utils::numeric(FLERR, arg[7], false, lmp);
  rho=utils::numeric(FLERR, arg[8], false, lmp);
  Pinfty=utils::numeric(FLERR, arg[9], false, lmp);
  sigma=utils::numeric(FLERR, arg[10], false, lmp);
  miu=utils::numeric(FLERR, arg[11], false, lmp);
  cutoff=utils::numeric(FLERR, arg[12], false, lmp);
  radius_in=radius-cutoff;
  freq=utils::numeric(FLERR, arg[13], false, lmp);
  w=2.0*MY_PI*freq;
  PA=utils::numeric(FLERR, arg[14], false, lmp)*101325;
  tstart=utils::numeric(FLERR, arg[15], false, lmp);
  scalefactor=utils::numeric(FLERR, arg[16], false, lmp);
  Tinfty=utils::numeric(FLERR, arg[17], false, lmp)*scalefactor;
  alpha1=utils::numeric(FLERR, arg[18], false, lmp);
  // delta variable checking
  if (utils::strmatch(arg[19], "^v_")) {
    deltastr = utils::strdup(arg[19] + 2);
    SL_delta=0.0;
    SL_deltaold=0.0;
    deltastyle=VARIABLE;
    varshape=1;
    }
  else if (utils::strmatch(arg[19], "^m_"))
  {
    deltastr = utils::strdup(arg[19] + 2);
    SL_delta=0.0;
    SL_deltaold=0.0;
    deltastyle=MIXING;
    varshape=1;
  }
  else {
      SL_delta=utils::numeric(FLERR, arg[19], false, lmp);
      SL_deltaold=SL_delta;
      deltastyle=AUTO;
    }
    utils::logmesg(lmp,"REGION KMSPHERE. deltasytle is {} \n",deltastyle);

  varshape=1; // enable the shape_update()
  if (varshape)
  {
    variable_check();
    if (deltastyle == VARIABLE || deltastyle==MIXING)
    delta_update();
  }


  //
  // error check

  if (radius < 0.0) error->all(FLERR, "Illegal region sphere radius: {} \n", radius);

  // extent of sphere
  // for variable radius, uses initial radius and origin for variable center

  if (interior) {
    bboxflag = 1;
    extent_xlo = xc - radius;
    extent_xhi = xc + radius;
    extent_ylo = yc - radius;
    extent_yhi = yc + radius;
    extent_zlo = zc - radius;
    extent_zhi = zc + radius;
  } else
    bboxflag = 0;

  cmax = 1;
  contact = new Contact[cmax];
  tmax = 1;

  utils::logmesg(lmp,"REGION KMSPHERE, Initial radius is {}\n",radius);
}

/* ---------------------------------------------------------------------- */

RegKMSphere::~RegKMSphere()
{
  delete[] deltastr;
  delete[] contact;
}

/* ---------------------------------------------------------------------- */

void RegKMSphere::init()
{
  Region::init();
  SL_radius=radius;
  SL_lastradius=radius;
  SL_vradius=vradius;
  pressure=0.0;
  pressure_old=0.0;
  if (varshape) variable_check();
  
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegKMSphere::inside(double x, double y, double z)
{
  double delx = x - xc;
  double dely = y - yc;
  double delz = z - zc;
  double r = sqrt(delx * delx + dely * dely + delz * delz);

  if (r <= radius) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   one contact if 0 <= x < cutoff from inner surface of sphere
   no contact if outside (possible if called from union/intersect)
   delxyz = vector from nearest point on sphere to x
   special case: no contact if x is at center of sphere
------------------------------------------------------------------------- */

int RegKMSphere::surface_interior(double *x, double cutoff)
{
  double delx = x[0] - xc;
  double dely = x[1] - yc;
  double delz = x[2] - zc;
  double r = sqrt(delx * delx + dely * dely + delz * delz);
  if (r > radius || r == 0.0) return 0;

  double delta = radius - r;
  if (delta < cutoff) {
    contact[0].r = delta;
    contact[0].delx = delx * (1.0 - radius / r);
    contact[0].dely = dely * (1.0 - radius / r);
    contact[0].delz = delz * (1.0 - radius / r);
    contact[0].radius = -radius;
    contact[0].iwall = 0;
    contact[0].varflag = 1;
    return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   one contact if 0 <= x < cutoff from outer surface of sphere
   no contact if inside (possible if called from union/intersect)
   delxyz = vector from nearest point on sphere to x
------------------------------------------------------------------------- */

int RegKMSphere::surface_exterior(double *x, double cutoff)
{
  double delx = x[0] - xc;
  double dely = x[1] - yc;
  double delz = x[2] - zc;
  double r = sqrt(delx * delx + dely * dely + delz * delz);
  if (r < radius) return 0;

  double delta = r - radius;
  if (delta < cutoff) {
    contact[0].r = delta;
    contact[0].delx = delx * (1.0 - radius / r);
    contact[0].dely = dely * (1.0 - radius / r);
    contact[0].delz = delz * (1.0 - radius / r);
    contact[0].radius = radius;
    contact[0].iwall = 0;
    contact[0].varflag = 1;
    return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   change region shape via variable evaluation
------------------------------------------------------------------------- */

void RegKMSphere::shape_update()
{
  
  if (update->ntimestep>1)
  {
  //radius=radius+vradius*update->dt;
  radius_in=radius-cutoff;
  //pressure+=0.1;
  //auto compute_press=modify->get_compute_by_id(utils::strdup(press_id + 2));
  //if (!compute_press)
  //        error->all(FLERR, "Could not find compute ID in KMSPHERE: {}",press_id);
  //press=compute_press->scalar;
  volume=4.0*MY_PI/3.0*(radius*radius*radius-radius_in*radius_in*radius_in);
  pressure=-stress/(3.0*volume)*101325;
  

  // invoke the RK4 update and log down thw old pressure
  SL_lastradius=radius;
  
  if (deltastyle==VARIABLE)
  {
    delta_update();
    SLRK4(pressure,pressure_old);
  }
  else if (deltastyle==AUTO)
  {
    SL_deltaold=SL_delta;
    SLRK4(pressure,pressure_old,SL_Tblliquid,SL_Tblliquidold);

    //if(update->ntimestep%1000==0)
    //utils::logmesg(lmp,"Tbl is {}, Tblold is {} \n",SL_Tblgas/scalefactor,SL_Tblgasold/scalefactor);
    //SLRK4(pressure,pressure_old,SL_Tblgas,SL_Tblgasold);
  }
  else if (deltastyle==MIXING)
  {
    if(update->ntimestep<100000) 
    {
      delta_update();
      SLRK4(pressure,pressure_old);
    }
    else {
       SL_deltaold=SL_delta;
       SLRK4(pressure,pressure_old,SL_Tblliquid,SL_Tblliquidold);
    }
  }
  
  pressure_old=pressure;
  SL_radius=radius;
  SL_vradius=vradius;
  
  }
  else delta_update();

  if(update->ntimestep%1==0)
  {
    //utils::logmesg(lmp,"REGION KMSPHERE, THE Press IN STEP {} IS {}, the RADIUS IS {}, the vradius is {}, number of atoms is {}\n",update->ntimestep,pressure/101325,radius,vradius,numatoms);
    //utils::logmesg(lmp,"REGION KMSPHERE, In step {}, Current delta is {}, deltaold is {}\n",update->ntimestep,SL_delta,SL_deltaold);
    
  }

}

/* ----------------------------------------------------------------------
   error check on existence of variable
------------------------------------------------------------------------- */

void RegKMSphere::variable_check()
{
  if (deltastyle == VARIABLE || deltastyle==MIXING) {
    deltavar=input->variable->find(deltastr);
    if (deltavar < 0) error->all(FLERR, "Variable {} for heat bath wall does not exist", deltastr);
    if (!input->variable->equalstyle(deltavar))
      error->all(FLERR, "Variable {} for  heat bath wall is invalid style", deltastr);
  }

}

/* ----------------------------------------------------------------------
   Set values needed to calculate velocity due to shape changes.
   These values do not depend on the contact, so this function is
   called once per timestep by fix/wall/gran/region.

------------------------------------------------------------------------- */

void RegKMSphere::set_velocity_shape()
{
  xcenter[0] = xc;
  xcenter[1] = yc;
  xcenter[2] = zc;
  forward_transform(xcenter[0], xcenter[1], xcenter[2]);
  if (update->ntimestep > 0)
    rprev = prev[4];
  else
    rprev = radius;
  prev[4] = radius;
}

/* ----------------------------------------------------------------------
   add velocity due to shape change to wall velocity
------------------------------------------------------------------------- */

void RegKMSphere::velocity_contact_shape(double *vwall, double *xc)
{
  double delx, dely, delz;    // Displacement of contact point in x,y,z

  delx = (xc[0] - xcenter[0]) * (1 - rprev / radius);
  dely = (xc[1] - xcenter[1]) * (1 - rprev / radius);
  delz = (xc[2] - xcenter[2]) * (1 - rprev / radius);

  vwall[0] += delx / update->dt;
  vwall[1] += dely / update->dt;
  vwall[2] += delz / update->dt;
}


void RegKMSphere::SLRK4(const double pressure,const double pressure_old)
{
  // y1: position
  // y2: velocity
  // dy1/dt=y2
  // dy2/dt=ddR
  
  double h=(update->dt*1e-15/force->femtosecond);
  double current_t=(tstart+update->ntimestep*update->dt)*1e-15/force->femtosecond;
  double diffpress=(pressure-pressure_old)/h;
  //utils::logmesg(lmp,"SLRK4, Input radius is {}, vradius is {}, diffpress is {}\n",radius,vradius,diffpress);
  double k1y1=vradius*1e5;   // standard unit
  double k1y2=calf2(current_t,radius*1e-10,vradius*1e5,pressure,pressure_old);

  // h*k1y2 -> standard unit
  double k2y1=vradius*1e5+(h*k1y2/2.0);
  double k2y2=calf2(current_t+h/2.0,radius*1e-10+h/2.0*k1y1,vradius*1e5+h/2.0*k1y2,pressure,pressure_old);

  double k3y1=vradius*1e5+(h*k2y2/2.0);
  double k3y2=calf2(current_t+h/2.0,radius*1e-10+h/2.0*k2y1,vradius*1e5+h/2.0*k2y2,pressure,pressure_old);

  double k4y1=vradius*1e5+(h*k3y1);
  double k4y2=calf2(current_t+h,radius*1e-10+h*k3y1,vradius*1e5+h*k3y2,pressure,pressure_old);
  
  radius=radius+h/6.0*(k1y1+2.0*k2y1+2.0*k3y1+k4y1)*1e10;
  vradius=vradius+h/6.0*(k1y2+2.0*k2y2+2.0*k3y2+k4y2)*1e-5;
  //utils::logmesg(lmp,"SLRK4, OUTPUT radius is {}, vradius is {}\n",radius,vradius);
  //double debugvradius=1.0/6.0*(k1y1+2.0*k2y1+2.0*k3y1+k4y1)*1e-5;
  //utils::logmesg(lmp,"debugvradius is {},calculated vradius is {} \n",debugvradius,vradius);
  //utils::logmesg(lmp,"k1y1 is {}, k2y1 is {}, k3y1 is {}, k4y1 is {}, calculated vradius is {} \n",k1y1,k2y1,k3y1,k4y1,vradius);
  //utils::logmesg(lmp,"k1y2 is {}, k2y2 is {}, k3y2 is {}, k4y2 is {} \n",k1y2,k2y2,k3y2,k4y2);

}

double RegKMSphere::calf2(double update_t, double update_radius, double update_vradius, double pressure, double pressure_old)
{
  // the units of the input should be standard
  double h=(update->dt*1e-15/force->femtosecond);
  double diffpress=(pressure-pressure_old)/h;
  double f2=-((3.0*powint(update_vradius,2)*(update_vradius/(3.0*CB) - 1.0))/2.0 - ((update_vradius/CB + 1)*(Pinfty - pressure - PA*sin(w*(update_t + update_radius/CB)) \
  + (2.0*sigma)/update_radius + (4.0*update_vradius*miu)/update_radius) - (update_radius*(diffpress + (2.0*update_vradius*sigma)/powint(update_radius,2) + (4.0*powint(update_vradius,2)*miu)/powint(update_radius,2) \
  + PA*w*cos(w*(update_t + update_radius/CB))*(update_vradius/CB + 1)))/CB)/rho)\
  /(update_radius*(update_vradius/CB - 1.0) - (4.0*miu)/(CB*rho));

  //return standard output
  return f2;
}

void RegKMSphere::SLRK4(const double pressure,const double pressure_old, const double Tbl, const double Tbl_old)
{
  // y1: position
  // y2: velocity
  // dy1/dt=y2
  // dy2/dt=ddR

  double h=(update->dt*1e-15/force->femtosecond);
  double current_t=(tstart+update->ntimestep*update->dt)*1e-15/force->femtosecond;
  double diffpress=(pressure-pressure_old)/h;
  double diffTbl=(Tbl-Tbl_old)/h;
  double k1y1=vradius*1e5;   // standard unit
  double k1y2=calf2(current_t,radius*1e-10,vradius*1e5,pressure,pressure_old);
  double k1y3=calf3(current_t,radius*1e-10,vradius*1e5,SL_delta*1e-10,Tbl,Tbl_old);
  
  // h*k1y2 -> standard unit
  double k2y1=vradius*1e5+(h*k1y2/2.0);
  double k2y2=calf2(current_t+h/2.0,radius*1e-10+h/2.0*k1y1,vradius*1e5+h/2.0*k1y2,pressure,pressure_old);
  double k2y3=calf3(current_t+h/2.0,radius*1e-10+h/2.0*k1y1,vradius*1e5+h/2.0*k1y2,SL_delta*1e-10+h/2.0*k1y3,Tbl,Tbl_old);

  double k3y1=vradius*1e5+(h*k2y2/2.0);
  double k3y2=calf2(current_t+h/2.0,radius*1e-10+h/2.0*k2y1,vradius*1e5+h/2.0*k2y2,pressure,pressure_old);
  double k3y3=calf3(current_t+h/2.0,radius*1e-10+h/2.0*k2y1,vradius*1e5+h/2.0*k2y2,SL_delta*1e-10+h/2.0*k2y3,Tbl,Tbl_old);

  double k4y1=vradius*1e5+(h*k3y1);
  double k4y2=calf2(current_t+h,radius*1e-10+h*k3y1,vradius*1e5+h*k3y2,pressure,pressure_old);
  double k4y3=calf3(current_t+h,radius*1e-10+h*k3y1,vradius*1e5+h*k3y2,SL_delta*1e-10+h*k3y3,Tbl,Tbl_old);
  
  radius=radius+h/6.0*(k1y1+2.0*k2y1+2.0*k3y1+k4y1)*1e10;
  vradius=vradius+h/6.0*(k1y2+2.0*k2y2+2.0*k3y2+k4y2)*1e-5;
  SL_delta=SL_delta+h/6.0*(k1y3+2.0*k2y3+2.0*k3y3+k4y3)*1e10;
  //double debugvradius=1.0/6.0*(k1y1+2.0*k2y1+2.0*k3y1+k4y1)*1e-5;
  //utils::logmesg(lmp,"debugvradius is {},calculated vradius is {} \n",debugvradius,vradius);
  //utils::logmesg(lmp,"k1y1 is {}, k2y1 is {}, k3y1 is {}, k4y1 is {}, calculated vradius is {} \n",k1y1,k2y1,k3y1,k4y1,vradius);
  //utils::logmesg(lmp,"k1y2 is {}, k2y2 is {}, k3y2 is {}, k4y2 is {} \n",k1y2,k2y2,k3y2,k4y2);
  //if(update->ntimestep%1000==0)
  //utils::logmesg(lmp,"k1y3 is {}, k2y3 is {}, k3y3 is {}, k4y3 is {}, Tbl is {}, Tblold is {}, delta is {} \n",k1y3,k2y3,k3y3,k4y3,Tbl,Tbl_old,SL_delta);

}

double RegKMSphere::calf3(double update_t, double update_radius, double update_vradius, const double update_delta, const double Tbl, const double Tbl_old)
{
  // the units of the input should be standard
  double h=(update->dt*1e-15/force->femtosecond);
  double diffpress=(pressure-pressure_old)/h;
  double diffTbl=(Tbl-Tbl_old)/h;
  double denominator=1.0+update_delta/update_radius+3.0/10.0*powint(update_delta/update_radius,2);
  double numerator=6*alpha1/update_delta-(2.0*update_delta/update_radius+1.0/2.0*powint((update_delta/update_radius),2))*update_vradius-update_delta*(1.0+update_delta/(2.0*update_radius)+1.0/10.0*powint((update_delta/update_radius),2))*1.0/(Tbl-Tinfty)*diffTbl;
  double f3=numerator/denominator;
  //if(update->ntimestep%1000==0)
  utils::logmesg(lmp,"numerator is {}, denominator is {}, Tbl is {}, Tblold is {}, delta is {}, deltaold is {} \n",numerator,denominator,Tbl,Tbl_old,SL_delta,SL_deltaold);
  //return standard output
  return f3;
}

void RegKMSphere::delta_update()
{
  //utils::logmesg(lmp,"Deltastr is {}\n",deltastr);
  //utils::logmesg(lmp,"Delta update begins! Delta style is Constant {} Variable {} Auto {}\n",deltastyle==CONSTANT,deltastyle==VARIABLE,deltastyle==AUTO);
  if (deltastyle == VARIABLE || deltastyle==MIXING) 
  {
    SL_deltaold=SL_delta;
    SL_delta= input->variable->compute_equal(deltavar);
    //utils::logmesg(lmp,"In step {}, Current delta is {} \n",update->ntimestep,SL_delta);
  }
  if (SL_delta < 0.0) error->one(FLERR, "Variable delta evaluation in heat bath boundary gave bad value");
  
}