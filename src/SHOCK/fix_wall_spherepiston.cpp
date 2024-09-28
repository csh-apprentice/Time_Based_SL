// clang-format off
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

#include "fix_wall_spherepiston.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "lattice.h"
#include "math_const.h"
#include "random_mars.h"
#include "update.h"
#include "cassert"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixWallSpherePiston::FixWallSpherePiston(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), randomt(nullptr), gfactor1(nullptr), gfactor2(nullptr)
{
  force_reneighbor = 1;
  next_reneighbor = -1;

  if (narg < 4) error->all(FLERR,"Illegal fix wall/piston command");

  randomt = nullptr;
  gfactor1 = gfactor2 = nullptr;
  scaleflag = 1;
 
  roughdist = 0.0;
  x0=y0=z0=0.0;
  radius_flag=origin_flag=0;
  vr=0.0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"origin") == 0)
     {
      origin_flag=1;
      x0=utils::numeric(FLERR,arg[iarg+1],false,lmp);
      y0=utils::numeric(FLERR,arg[iarg+2],false,lmp);
      z0=utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg+=4;
     }
  else if (strcmp(arg[iarg],"radius") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall/piston command");
      radius_flag=1;
      radius = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    }
    
    else if (strcmp(arg[iarg],"vel") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall/piston command");
      vr = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall/piston command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix wall/piston command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix wall/piston command");
  }

  // setup scaling
  const double xscale = (scaleflag) ? domain->lattice->xlattice : 1.0;
  const double yscale = (scaleflag) ? domain->lattice->ylattice : 1.0;
  const double zscale = (scaleflag) ? domain->lattice->zlattice : 1.0;
  x0 *= xscale;
  y0*=yscale;
  z0 *= zscale;
  //roughdist *= zscale;

}

/* ---------------------------------------------------------------------- */

FixWallSpherePiston::~FixWallSpherePiston()
{
  delete[] gfactor2;
  delete[] gfactor1;
  delete randomt;
}

/* ---------------------------------------------------------------------- */

int FixWallSpherePiston::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallSpherePiston::initial_integrate(int /*vflag*/)
{
  next_reneighbor = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixWallSpherePiston::post_integrate()
{
  double zlo;

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double t    = (update->ntimestep - update->beginstep) * update->dt;
  double tott = (update->endstep - update->beginstep) * update->dt;
  double tt = t * t;
  double ttt = tt * t;
  double tttt = tt * tt;
  double t0p5 = sqrt(t/tott);
  double t1p5 = t0p5*t0p5*t0p5;
  double t2p5 = t1p5*t0p5*t0p5;
  float point2origin2_old=0.0;
  float point2origin_old=0.0;
  float point2origin2_new=0.0;
  float point2origin_new=0.0;
  //float point2origin_sphere=0.0;
  float nx,ny,nz;
  nx=ny=nz=0.0;
  float pos_scale=0.0;
  float normal_scale=0.0;
  float sphere_scale=0.0;

  float v_hor=0.0;
  float v_ver=0.0;



  if (radius_flag) { radius = radius+vr * update->dt; }
  

  if ((update->ntimestep % 100 == 0) && (comm->me == 0))
    utils::logmesg(lmp,"SHOCK: step {} t {} origin {} {} {} radius {}\n",
                   update->ntimestep, t, x0,y0,z0,radius);

  float radius2=radius*radius;
  // VIRIAL PRESSURE CONTRIBUTION?
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) { 
      roughoff = 0.1;
      point2origin2_old=(x[i][0]-x0)*(x[i][0]-x0)+(x[i][1]-y0)*(x[i][1]-y0)+(x[i][2]-z0)*(x[i][2]-z0);
      if (radius_flag &&point2origin2_old>radius2 ) {
        point2origin_old=sqrt(point2origin2_old);
        point2origin_new= 2.0 * (radius - roughoff) - point2origin_old;
        //pos_scale=point2origin2_new/point2origin2_old;
        normal_scale=1/point2origin_old;
        sphere_scale=radius/(point2origin_old+roughoff);  //<1
        //sphere_scale=0.9;
        float x2norm=(x[i][0]-x0)*(x[i][0]-x0);
        float y2norm=(x[i][1]-y0)*(x[i][1]-y0);
        float z2norm=(x[i][2]-z0)*(x[i][2]-z0);
        nx=sqrt(x2norm/(x2norm+y2norm+z2norm));
        ny=sqrt(y2norm/(x2norm+y2norm+z2norm));
        nz=sqrt(z2norm/(x2norm+y2norm+z2norm));

        float oldx=x[i][0];
        float oldy=x[i][1];
        float oldz=x[i][2];
        x[i][0] =2.0*((x[i][0]-x0)*sphere_scale)-(x[i][0]-x0)+x0;
        x[i][1] =2.0*((x[i][1]-y0)*sphere_scale)-(x[i][1]-y0)+y0;
        x[i][2] =2.0*((x[i][2]-z0)*sphere_scale)-(x[i][2]-z0)+z0;

        //x[i][0] =x[i][0]*0.9;//2.0*x[i][0]*sphere_scale-x[i][0];
        //x[i][1] =x[i][1]*0.9;//2.0*x[i][1]*sphere_scale-x[i][1];
        //x[i][2] =x[i][2]*0.9;//2.0*x[i][2]*sphere_scale-x[i][2];
        
        float new_radius=sqrt(x[i][0]*x[i][0]+x[i][1]*x[i][1]+x[i][2]*x[i][2]);
        float v_norm=sqrt(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2]);
        if(new_radius>radius)
        {
          utils::logmesg(lmp,"RADIUS OUT OF BOUNDARY!!step {} new_radius {} radius {} sphere scale {} oldx {} oldy {} oldz {} v_norm {}\n",update->ntimestep,new_radius,radius,sphere_scale,oldx,oldy,oldz,v_norm);
        }


        //x[i][0] =(rand() % (1000) / (float)(1000))*radius/2;
        //x[i][1] =(rand() % (1000) / (float)(1000))*radius/2;
        //x[i][2] =(rand() % (1000) / (float)(1000))*radius/2;

        v_ver=nx*v[i][0]+ny*v[i][1]+nz*v[i][2];

        //v[i][0] =0.0;//v[i][0]-2.0*nx*v_ver;
        //v[i][1]=0.0;//v[i][1]-2.0*ny*v_ver;
        //v[i][2]=0.0;//v[i][2]-2.0*nz*v_ver;
        
        if(v_ver>v_norm+0.1)
        {
          utils::logmesg(lmp,"Speed Calculation ERROR!! nx {} ny {} nz {} v_norm {} v_ver {}\n",nx,ny,nz,v_norm,v_ver);
        }

        v[i][0] =v[i][0]-2.0*nx*v_ver;
        v[i][1]=v[i][1]-2.0*ny*v_ver;
        v[i][2]=v[i][2]-2.0*nz*v_ver;

        
      }
    }
  }

}
