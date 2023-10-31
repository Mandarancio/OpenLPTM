#ifndef _EGUN_UTILS_H__
#define _EGUN_UTILS_H__

inline double radiationResistance(double e_b1, double A_b1, double e_b2,
                                  double A_b2) {
  double energy_resistance =
      (1 - e_b1) / (e_b1 * A_b1) + 1 / A_b1 + (1 - e_b2) / (e_b2 * A_b2);
  return 1. / energy_resistance;
}

#endif
