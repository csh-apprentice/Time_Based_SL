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

/* ----------------------------------------------------------------------
   Contributing author: Ray Shan (SNL)
------------------------------------------------------------------------- */

#include "pair_coul_long_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;


#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairCoulLongKokkos<DeviceType>::PairCoulLongKokkos(LAMMPS *lmp):PairCoulLong(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairCoulLongKokkos<DeviceType>::~PairCoulLongKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->destroy_kokkos(k_cutsq,cutsq);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairCoulLongKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);
  k_cutsq.template sync<DeviceType>();
  k_params.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  qqrd2e = force->qqrd2e;
  newton_pair = force->newton_pair;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];
  special_coul[0] = force->special_coul[0];
  special_coul[1] = force->special_coul[1];
  special_coul[2] = force->special_coul[2];
  special_coul[3] = force->special_coul[3];

  // loop over neighbors of my atoms

  copymode = 1;

  EV_FLOAT ev;
  if (ncoultablebits)
    ev = pair_compute<PairCoulLongKokkos<DeviceType>,CoulLongTable<1> >
      (this,(NeighListKokkos<DeviceType>*)list);
  else
    ev = pair_compute<PairCoulLongKokkos<DeviceType>,CoulLongTable<0> >
      (this,(NeighListKokkos<DeviceType>*)list);


  if (eflag) {
    eng_vdwl += ev.evdwl;
    eng_coul += ev.ecoul;
  }
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  copymode = 0;
}

/* ----------------------------------------------------------------------
   compute coulomb pair force between atoms i and j
   ---------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS,  class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairCoulLongKokkos<DeviceType>::
compute_fcoul(const F_FLOAT& rsq, const int& /*i*/, const int&j,
              const int& /*itype*/, const int& /*jtype*/,
              const F_FLOAT& factor_coul, const F_FLOAT& qtmp) const {
  if (Specialisation::DoTable && rsq > tabinnersq) {
    union_int_float_t rsq_lookup;
    rsq_lookup.f = rsq;
    const int itable = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;
    const F_FLOAT fraction = (rsq_lookup.f - d_rtable[itable]) * d_drtable[itable];
    const F_FLOAT table = d_ftable[itable] + fraction*d_dftable[itable];
    F_FLOAT forcecoul = qtmp*q[j] * table;
    if (factor_coul < 1.0) {
      const F_FLOAT table = d_ctable[itable] + fraction*d_dctable[itable];
      const F_FLOAT prefactor = qtmp*q[j] * table;
      forcecoul -= (1.0-factor_coul)*prefactor;
    }
    return forcecoul/rsq;
  } else {
    const F_FLOAT r = sqrt(rsq);
    const F_FLOAT grij = g_ewald * r;
    const F_FLOAT expm2 = exp(-grij*grij);
    const F_FLOAT t = 1.0 / (1.0 + EWALD_P*grij);
    const F_FLOAT rinv = 1.0/r;
    const F_FLOAT erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
    const F_FLOAT prefactor = qqrd2e * qtmp*q[j]*rinv;
    F_FLOAT forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
    if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;

    return forcecoul*rinv*rinv;
  }
}

/* ----------------------------------------------------------------------
   compute coulomb pair potential energy between atoms i and j
   ---------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairCoulLongKokkos<DeviceType>::
compute_ecoul(const F_FLOAT& rsq, const int& /*i*/, const int&j,
              const int& /*itype*/, const int& /*jtype*/,
              const F_FLOAT& factor_coul, const F_FLOAT& qtmp) const {
  if (Specialisation::DoTable && rsq > tabinnersq) {
    union_int_float_t rsq_lookup;
    rsq_lookup.f = rsq;
    const int itable = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;
    const F_FLOAT fraction = (rsq_lookup.f - d_rtable[itable]) * d_drtable[itable];
    const F_FLOAT table = d_etable[itable] + fraction*d_detable[itable];
    F_FLOAT ecoul = qtmp*q[j] * table;
    if (factor_coul < 1.0) {
      const F_FLOAT table = d_ctable[itable] + fraction*d_dctable[itable];
      const F_FLOAT prefactor = qtmp*q[j] * table;
      ecoul -= (1.0-factor_coul)*prefactor;
    }
    return ecoul;
  } else {
    const F_FLOAT r = sqrt(rsq);
    const F_FLOAT grij = g_ewald * r;
    const F_FLOAT expm2 = exp(-grij*grij);
    const F_FLOAT t = 1.0 / (1.0 + EWALD_P*grij);
    const F_FLOAT erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
    const F_FLOAT prefactor = qqrd2e * qtmp*q[j]/r;
    F_FLOAT ecoul = prefactor * erfc;
    if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
    return ecoul;
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairCoulLongKokkos<DeviceType>::allocate()
{
  PairCoulLong::allocate();

  int n = atom->ntypes;

  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  d_cut_coulsq = typename AT::t_ffloat_2d("pair:cut_coulsq",n+1,n+1);

  k_params = Kokkos::DualView<params_coul**,Kokkos::LayoutRight,DeviceType>("PairCoulLong::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

template<class DeviceType>
void PairCoulLongKokkos<DeviceType>::init_tables(double cut_coul, double *cut_respa)
{
  Pair::init_tables(cut_coul,cut_respa);

  typedef typename ArrayTypes<DeviceType>::t_ffloat_1d table_type;
  typedef typename ArrayTypes<LMPHostType>::t_ffloat_1d host_table_type;

  int ntable = 1;
  for (int i = 0; i < ncoultablebits; i++) ntable *= 2;


  // Copy rtable and drtable
  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);
  for (int i = 0; i < ntable; i++) {
    h_table(i) = rtable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_rtable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);
  for (int i = 0; i < ntable; i++) {
    h_table(i) = drtable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_drtable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  // Copy ftable and dftable
  for (int i = 0; i < ntable; i++) {
    h_table(i) = ftable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_ftable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  for (int i = 0; i < ntable; i++) {
    h_table(i) = dftable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_dftable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  // Copy ctable and dctable
  for (int i = 0; i < ntable; i++) {
    h_table(i) = ctable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_ctable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  for (int i = 0; i < ntable; i++) {
    h_table(i) = dctable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_dctable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  // Copy etable and detable
  for (int i = 0; i < ntable; i++) {
    h_table(i) = etable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_etable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  for (int i = 0; i < ntable; i++) {
    h_table(i) = detable[i];
  }
  Kokkos::deep_copy(d_table,h_table);
  d_detable = d_table;
  }
}


/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

template<class DeviceType>
void PairCoulLongKokkos<DeviceType>::settings(int narg, char **arg)
{
  if (narg > 2) error->all(FLERR,"Illegal pair_style command");

  PairCoulLong::settings(narg,arg);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairCoulLongKokkos<DeviceType>::init_style()
{
  PairCoulLong::init_style();

  Kokkos::deep_copy(d_cut_coulsq,cut_coulsq);

  // error if rRESPA with inner levels

  if (update->whichflag == 1 && utils::strmatch(update->integrate_style,"^respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    if (respa)
      error->all(FLERR,"Cannot use Kokkos pair style with rRESPA inner/middle");
  }

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairCoulLongKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairCoulLong::init_one(i,j);

  k_params.h_view(i,j).cut_coulsq = cut_coulsq;

  k_params.h_view(j,i) = k_params.h_view(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.h_view(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = cutone*cutone;
    m_cut_coulsq[j][i] = m_cut_coulsq[i][j] = cut_coulsq;
  }

  k_cutsq.h_view(i,j) = cutone*cutone;
  k_cutsq.template modify<LMPHostType>();
  k_params.template modify<LMPHostType>();

  return cutone;
}

namespace LAMMPS_NS {
template class PairCoulLongKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairCoulLongKokkos<LMPHostType>;
#endif
}
