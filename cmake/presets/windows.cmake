set(WIN_PACKAGES
  AMOEBA
  ASPHERE
  AWPMD
  BOCS
  BODY
  BPM
  BROWNIAN
  CG-DNA
  CG-SPICA
  CLASS2
  COLLOID
  COLVARS
  CORESHELL
  DIELECTRIC
  DIFFRACTION
  DIPOLE
  DPD-BASIC
  DPD-MESO
  DPD-REACT
  DPD-SMOOTH
  DRUDE
  EFF
  ELECTRODE
  EXTRA-COMPUTE
  EXTRA-DUMP
  EXTRA-FIX
  EXTRA-MOLECULE
  EXTRA-PAIR
  FEP
  GRANULAR
  INTERLAYER
  KSPACE
  LEPTON
  MACHDYN
  MANIFOLD
  MANYBODY
  MC
  MEAM
  MESONT
  MISC
  ML-IAP
  ML-POD
  ML-SNAP
  MOFFF
  MOLECULE
  MOLFILE
  OPENMP
  OPT
  ORIENT
  PERI
  PHONON
  POEMS
  PLUGIN
  PTM
  QEQ
  QTB
  REACTION
  REAXFF
  REPLICA
  RIGID
  SHOCK
  SMTBQ
  SPH
  SPIN
  SRD
  TALLY
  UEF
  YAFF)

foreach(PKG ${WIN_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()
