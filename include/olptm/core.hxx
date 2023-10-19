#ifndef _LTM_H__
#define _LTM_H__

// TODO: constant are not constant and depends from the temperature!!!

#include <eigen3/Eigen/Core>
#include <ratio>
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
  std::vector<relation_t> exchanges;  // list of heat exchange relations
};

template <typename T> inline void insert(Eigen::VectorX<T> &v, T val) {
  u32 size = v.size();
  Eigen::VectorX<T> tmp = v;
  v.resize(size + 1);
  v.head(size) = tmp;
  v[size] = val;
}

inline void add_body(system_t &sys, body_t &body) {
  u32 id = sys.bodies.size();
  body.id = id;

  sys.bodies.push_back(body);
  insert<f64>(sys.temperatures, body.T0);
  insert<f64>(sys.heats, 0);
  insert<f64>(sys.inv_capacities, body.inv_capacity(body.T0) * body.inv_mass);
}

inline void add_exchange(system_t &sys, relation_t rel) {
  sys.exchanges.push_back(rel);
}

inline void evaluate(system_t &sys, f64 dt) {
  sys.heats.setZero();
  // compute all heat exchange
  for (relation_t r : sys.exchanges) {
    u32 x = r.id1;
    u32 y = r.id2;
    f64 heat = r.fn(sys.temperatures[x], sys.temperatures[y], r.exc_k);
    sys.heats[x] += heat;
    sys.heats[y] -= heat;
  }
  // integrate results
  sys.temperatures += sys.heats.cwiseProduct(sys.inv_capacities) * dt;
}

inline void evaluate(system_t &sys, const f64 DT, f64 &dt, f64 min_dt = 0.001,
                     f64 max_dt = 1.0) {
  sys.heats.setZero();
  // compute all heat exchange
  for (relation_t r : sys.exchanges) {
    u32 x = r.id1;
    u32 y = r.id2;
    f64 heat = r.fn(sys.temperatures[x], sys.temperatures[y], r.exc_k);
    sys.heats[x] += heat;
    sys.heats[y] -= heat;
  }
  Eigen::VectorX<f64> Q = sys.heats.cwiseProduct(sys.inv_capacities);
  f64 max_Q = Q.cwiseAbs().maxCoeff();
  if (max_Q == 0) { // constant case
    dt = max_dt;
    return;
  }
  f64 rdt = DT / max_Q;
  dt = (rdt < min_dt ? min_dt : (rdt > max_dt ? max_dt : rdt));
  // integrate results
  sys.temperatures += Q * dt;
}

inline void update_parameters(system_t &sys) {
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

} // namespace lptm

#endif
