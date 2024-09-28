/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
#ifdef FIX_CLASS
// clang-format off
FixStyle(wall/spherepiston,FixWallSpherePiston);
// clang-format on
#else

#ifndef LMP_FIX_WALL_SPHEREPISTON_H
#define LMP_FIX_WALL_SPHEREPISTON_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallSpherePiston : public Fix {
 public:
  FixWallSpherePiston(class LAMMPS *, int, char **);
  ~FixWallSpherePiston() override;
  int setmask() override;
  void post_integrate() override;
  void initial_integrate(int) override;

 private:
  int radius_flag, origin_flag;
  int scaleflag;
  double roughdist, roughoff, x0, y0, z0, vx, vy, vz, maxvx, maxvy, maxvz, paccelx, paccely,
      paccelz, angfreq,vr,radius;
  class RanMars *randomt;
  double *gfactor1, *gfactor2;
};

}    // namespace LAMMPS_NS

#endif
#endif
