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

#ifndef LMP_FIX_WALL_H
#define LMP_FIX_WALL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWall : public Fix {
 public:
  int nwall;
  int wallwhich[8];
  double coord0[8];
  double coordradius;
  //double coordorigin[3];
  int xflag;    // 1 if any wall position is a variable
  int xstyle[8];
  int xindex[8];
  char *xstr[8];

  FixWall(class LAMMPS *, int, char **);
  ~FixWall() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void pre_force(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double compute_scalar() override;
  double compute_vector(int) override;

  virtual void precompute(int) = 0;
  virtual void wall_particle(int, int, double) = 0;

 protected:
  double epsilon[8], sigma[8], alpha[8], cutoff[8];
  double ewall[9], ewall_all[9];
  double xscale, yscale, zscale;
  int estyle[8], sstyle[8], astyle[8], wstyle[8];
  int eindex[8], sindex[8];
  char *estr[8], *sstr[8], *astr[8], *lstr[8], *fstr[8], *kstr[8];
  int varflag;    // 1 if any wall position,epsilon,sigma is a variable
  int eflag;      // per-wall flag for energy summation
  int ilevel_respa;
  int fldflag;
};

}    // namespace LAMMPS_NS

#endif
