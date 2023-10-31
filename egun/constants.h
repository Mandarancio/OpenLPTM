#ifndef _EGUN_CONSTANTS_H__
#define _EGUN_CONSTANTS_H__
#include "utils.h"

namespace Density {
const double Tungsten = 1.93e4;
const double Alumina = 1.93e4;
}; // namespace Density

namespace SpecificHeat {
const double Alumina = 134.4;
const double Tungsten = 134.4;
}; // namespace SpecificHeat

namespace Conductivity {
const double Alumina = 170;
const double Tungsten = 170;
}; // namespace Conductivity

namespace Emissivity {
const double Tungsten = 0.5;
const double Alumina = 0.5;
}; // namespace Emissivity

namespace def {
const unsigned int n_filamets = 8;

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
}; // namespace def
#endif
