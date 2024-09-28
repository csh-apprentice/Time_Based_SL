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

#ifdef REGION_CLASS
// clang-format off
RegionStyle(kmsphere,RegKMSphere);
// clang-format on
#else

#ifndef LMP_REGION_KMSPHERE_H
#define LMP_REGION_KMSPHERE_H

#include "region.h"

namespace LAMMPS_NS {

class RegKMSphere : public Region {
 public:
  RegKMSphere(class LAMMPS *, int, char **);
  ~RegKMSphere() override;
  void init() override;
  int inside(double, double, double) override;
  int surface_interior(double *, double) override;
  int surface_exterior(double *, double) override;
  void shape_update() override;
  void set_velocity_shape() override;
  void velocity_contact_shape(double *, double *) override;

 private:
  double xc, yc, zc;

  int xstyle, xvar;
  int ystyle, yvar;
  int zstyle, zvar;
  int rstyle, svar;
  char *xstr, *ystr, *zstr, *sstr;

  double radius;
  double radius_in;
  double vradius;  // dR
  double CB; //the velocity of sound in mediin 
  double rho;  // the density of the medium
  double Pinfty;  //the medium pressure
  double sigma; //surface tension
  double  miu;  //dynamics viscosity of liquid
  //auto  compute_press;
  double stresssum; // sum of the stress
  double pressure_old;  //current pressure
  double cutoff; // cutoff in the fix/wall
  double volume; //volume of the outer shell

  double freq; // frequency of the sound
  double w;
  double PA; // amplitude of the sound

  double tstart; // the starting time in fs
  //double 
  void variable_check();
  void SLRK4(double pressure, double pressure_old);
  double calf2(double update_t, double update_radius, double update_vradius, double pressure, double pressure_old);
};

}    // namespace LAMMPS_NS

#endif
#endif
