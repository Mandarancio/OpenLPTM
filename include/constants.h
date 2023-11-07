#ifndef _CONSTANTS_H__
#define _CONSTANTS_H__

#define _CNST_F64 static const double

namespace lptm {
// J / (kg * K)
namespace SpecificHeat {
// metals
_CNST_F64 Ag = 233;
_CNST_F64 Al = 897;
_CNST_F64 Au = 129;
_CNST_F64 Be = 1820;
_CNST_F64 Bi = 123;
_CNST_F64 Cd = 231;
_CNST_F64 Cr = 449;
_CNST_F64 Cu = 385;
_CNST_F64 Fe = 444;
_CNST_F64 Hg = 139.5;
_CNST_F64 Li = 3580;
_CNST_F64 Mg = 1020;
_CNST_F64 Na = 1230;
_CNST_F64 Pb = 129;
_CNST_F64 Sn = 227;
_CNST_F64 Ti = 523;
_CNST_F64 U = 116;
_CNST_F64 W = 134;
_CNST_F64 Zn = 387;
// gas
_CNST_F64 Ar = 520.3;
_CNST_F64 CH4 = 2191;
_CNST_F64 CO2 = 839;
_CNST_F64 H2 = 14300;
_CNST_F64 He = 5193.2;
_CNST_F64 N = 1040;
_CNST_F64 O2 = 918;
// liquid
_CNST_F64 H2O = 4183;
} // namespace SpecificHeat

// W / (m * K)
namespace ThermalConductivity {
_CNST_F64 Ag = 406;
_CNST_F64 Al = 237;
_CNST_F64 Au = 315;
_CNST_F64 Al2O3 = 30;
_CNST_F64 Cu = 401;
_CNST_F64 Fe = 83.5;
_CNST_F64 Mn = 7.81;
_CNST_F64 W = 197;
// Liquid
_CNST_F64 H2O = 0.5918;
} // namespace ThermalConductivity
} // namespace lptm
#endif
