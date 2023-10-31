#ifndef _EGUN_CONSTANTS_H__
#define _EGUN_CONSTANTS_H__
#include "utils.h"

namespace Density {
const double Tungsten = 1.93e4;
const double ImprTungsten = 0.8 * Tungsten;
const double Alumina = 1.93e4;
const double Molybdenum = 1.022e4;
}; // namespace Density

namespace SpecificHeat {
const double Alumina = 134.4;
const double Tungsten = 134.4;
const double ImprTungsten = 142;
const double Molybdenum = 251;
}; // namespace SpecificHeat

namespace Conductivity {
const double Alumina = 170;
const double Tungsten = 170;
const double ImprTungsten = 170;
const double Molybdenum = 138;
const double Oil = 100;
}; // namespace Conductivity

namespace Emissivity {
const double Tungsten = 0.5;
const double ImprTungsten = 0.6;
const double Alumina = 0.5;
const double Molybdenum = 0.4;
const double Oil = 1;
}; // namespace Emissivity

namespace def {
const unsigned int n_filamets = 8;
const unsigned int n_lower_shields = 4;
const unsigned int n_upper_shields = 3;
const unsigned int n_em_segments = 3;
const unsigned int n_inner_bodies = 10;
const unsigned int n_outer_bodies = 11;

const double fil_masses[2] = {
    Density::Tungsten * 1.03e-7,
    Density::Tungsten * 1.10e-7,
};

const double alu_masses[2] = {
    Density::Alumina * 8.190e-7,
    Density::Alumina * 8.758e-7,
};

const double alumina_fialemnt[2] = {
    radiationResistance(Emissivity::Tungsten, 592.36e-6, Emissivity::Alumina,
                        961.3e-6),
    radiationResistance(Emissivity::Tungsten, 644.40e-6, Emissivity::Alumina,
                        1027.9e-6),
};

const double ls_masses[n_lower_shields] = {
    Density::Molybdenum * 568.1971e-9,
    Density::Molybdenum * 568.1971e-9,
    Density::Molybdenum * 568.1971e-9,
    Density::Molybdenum * 568.1971e-9,
};

const double us_masses[n_upper_shields] = {
    Density::Molybdenum * 3089e-9,
    Density::Molybdenum * 3615.5426e-9,
    Density::Molybdenum * 4816.7602e-9,
};

const double ns_mass = Density::Molybdenum * 76004.12e-9;
const double pr_mass = Density::Molybdenum * 73254.25e-9;

const double em_masses[n_em_segments] = {
    Density::ImprTungsten * 1,
    Density::ImprTungsten * 1,
    Density::ImprTungsten * 1,
};

const double oil_nosecone = radiationResistance(Emissivity::Molybdenum, 3.8232e-2, Emissivity::Oil, 1);
const double oil_prolongator = radiationResistance(Emissivity::Molybdenum, 4.18842e-2, Emissivity::Oil, 1); 
}; // namespace def
#endif
