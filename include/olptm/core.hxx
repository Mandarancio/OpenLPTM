#ifndef _LTM_H__
#define _LTM_H__

// TODO: constant are not constant and depends from the temperature!!!

#include <eigen3/Eigen/Core>
#include <stdint.h>
#include <vector>

// Power 4
#define _P4_(x) x *x *x *x

namespace lptm {

typedef float f32;
typedef double f64;
typedef uint32_t u32;
typedef int32_t i32;

// Stefan-Boltzman constant
static const f64 SIGMA = 5.6704e-8;

// simple conduction function
inline f64 __conduction_f(f64 T1, f64 T2, f64 R) { return (T2 - T1) / R; }
// simple radiation function
inline f64 __radiation_f(f64 T1, f64 T2, f64 alpha) {
  return alpha * (_P4_(T2) - _P4_(T1));
}

struct body_t {
  f64 (*inv_capacity)(f64 T); // capacity in  deg C / J # remember W = J / s
  f64 inv_mass;               // body mass in 1 / kg
  f64 T0;                     // intial temperature in deg C
  u32 id;                     // index in the system
};

struct relation_t {
  u32 id1;   // body 1 index
  u32 id2;   // body 2 index
  f64 exc_k; // exchange constant
  /**
   * Function used to compute heat exchange.
   * Parameters:
   * ----------
   *  - T1: temperature body 1 (unit: deg C)
   *  - T2: temperature body 2 (unit: deg C)
   *  - exc_k: exchange constant (unit: deg C / W)
   * Returns:
   * --------
   * Heat exchanged in W
   **/
  f64 (*fn)(f64 T1, f64 T2, f64 exc_k);
};

// Simple system descriptor
struct system_t {
  Eigen::VectorX<f64> temperatures;   // unit deg C
  Eigen::VectorX<f64> heats;          // unit W
  Eigen::VectorX<f64> inv_capacities; // unit deg C / J
  std::vector<body_t> bodies;         // list of bodies
  std::vector<relation_t> relations;  // list of heat exchange relations
};

inline void add_body(system_t &sys, body_t &body) {
  body.id = sys.bodies.size();
  sys.bodies.push_back(body);
}

inline void add_exchange(system_t &sys, relation_t rel) {
  sys.relations.push_back(rel);
}

inline void evaluate(system_t &sys, f64 dt) {
  sys.heats.setZero();
  // compute all heat exchange
  for (relation_t r : sys.relations) {
    u32 x = r.id1;
    u32 y = r.id2;
    f64 heat = r.fn(sys.temperatures[x], sys.temperatures[y], r.exc_k);
    sys.heats[x] += heat;
    sys.heats[y] -= heat;
  }
  // integrate results
  sys.temperatures += sys.heats.cwiseProduct(sys.inv_capacities) * dt;
  // update thermal capacity
  for (body_t b : sys.bodies) {
    sys.inv_capacities[b.id] =
        b.inv_mass * b.inv_capacity(sys.temperatures[b.id]);
  }
}

/**
 * Create a body with a constant specific heat
 * Unit:
 * -----
 *  - mass: kg
 *  - specific_heat: J / (kg * C)
 **/
inline body_t create(const f64 mass, f64 (*specific_heat)(f64 T),
                     const f64 T0) {
  return {specific_heat, 1.0 / mass, T0};
}

// create a body at constant temperature (infinite capcaity)
inline body_t create(f64 T) {
  return {[](f64) -> f64 { return 0; }, 0, T};
}

// create a conduction heat exchange relationship between 2 bodies
inline relation_t conduction(body_t b1, body_t b2, f64 Re1, f64 Re2) {
  f64 Req = Re1 + Re2;
  return {
      b1.id,
      b2.id,
      Req,
      __conduction_f,
  };
}

// create a radiation heat exchange relationship between 2 bodies
inline relation_t radiation(body_t b1, body_t b2, f64 correction_factor,
                            f64 view_factor) {
  f64 alpha = correction_factor * view_factor * SIGMA;
  return {
      b1.id,
      b2.id,
      alpha,
      __radiation_f,
  };
}

/**
 * Equivalent resistance for a surface contact
 * Units:
 * ------
 * - surface_area: m
 * - thermal_conductivity: W / (m * deg C)
 * - return: W / deg C
 **/
inline f64 Req(f64 surface_area, f64 thermal_conductity) {
  return surface_area / thermal_conductity;
}

// initialize system
inline void initialize(system_t &sys) {
  sys.heats = Eigen::VectorX<f64>(sys.bodies.size());
  sys.temperatures = Eigen::VectorX<f64>(sys.bodies.size());
  sys.inv_capacities = Eigen::VectorX<f64>(sys.bodies.size());
  for (u32 i = 0; i < sys.bodies.size(); i++) {
    sys.temperatures[i] = sys.bodies[i].T0;
    sys.inv_capacities[i] =
        sys.bodies[i].inv_capacity(sys.bodies[i].T0) * sys.bodies[i].inv_mass;
  }
}

// check if syste is initialized
inline bool is_initialized(system_t sys) {
  return sys.bodies.size() == sys.temperatures.size();
}
} // namespace lptm

#endif
