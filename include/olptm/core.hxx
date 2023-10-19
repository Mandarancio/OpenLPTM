#ifndef _LTM_H__
#define _LTM_H__

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
typedef uint64_t u64;
typedef int32_t i32;

// variable value depandant on temperature
typedef std::function<f64(f64 T)> fnT_t;

// Stefan-Boltzman constant
static const f64 SIGMA = 5.6704e-8;

// simple conduction function
inline f64 __conduction_f(f64 T1, f64 T2, f64 R) { return (T2 - T1) / R; }
// simple radiation function
inline f64 __radiation_f(f64 T1, f64 T2, f64 alpha) {
  return alpha * (_P4_(T2) - _P4_(T1));
}

struct body_t {
  fnT_t inv_specific_heat_fn;    // capacity in  deg C / J # remember W = J / s
  fnT_t thermal_conductivity_fn; // thermal conductivy in function of T
  f64 inv_mass;                  // body mass in 1 / kg
  f64 T0;                        // intial temperature in deg C
  f64 thermal_conductivy;        // current thermal conductivy
  u32 id;                        // index in the system
};

struct exchange_t {
  body_t *body_1; // body 1 index
  body_t *body_2; // body 2 index
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
  std::function<f64(f64 T1, f64 T2, body_t *b1, body_t *b2)> fn;
};

// Simple system descriptor
struct system_t {
  Eigen::VectorX<f64> temperatures;   // unit deg C
  Eigen::VectorX<f64> heats;          // unit W
  Eigen::VectorX<f64> inv_capacities; // unit deg C / J
  std::vector<body_t> bodies;         // list of bodies
  std::vector<exchange_t> exchanges;  // list of heat exchange relations
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
  insert<f64>(sys.inv_capacities,
              body.inv_specific_heat_fn(body.T0) * body.inv_mass);
}

inline void add_exchange(system_t &sys, exchange_t rel) {
  sys.exchanges.push_back(rel);
}

// compute the overall thermal exchange of the system.
inline void evaluate_exchanges(system_t &sys) {
  // set heats to 0
  sys.heats.setZero();
  // compute all heat exchange
  for (exchange_t r : sys.exchanges) {
    u32 x = r.body_1->id;
    u32 y = r.body_2->id;
    f64 heat =
        r.fn(sys.temperatures[x], sys.temperatures[y], r.body_1, r.body_2);
    sys.heats[x] += heat;
    sys.heats[y] -= heat;
  }
}

// evaluate the system with a constant time step
inline void evaluate(system_t &sys, f64 dt) {
  evaluate_exchanges(sys);
  // integrate results
  sys.temperatures += sys.heats.cwiseProduct(sys.inv_capacities) * dt;
}

/**
 * Evaluate the system with a dynamic time step driven by the desire maximum
 *change of temperature.
 **/
inline void evaluate(system_t &sys, const f64 dT, f64 &dt, f64 min_dt = 0.001,
                     f64 max_dt = 1.0) {
  evaluate_exchanges(sys);
  Eigen::VectorX<f64> Q = sys.heats.cwiseProduct(sys.inv_capacities);
  f64 max_Q = Q.cwiseAbs().maxCoeff();
  if (max_Q == 0) { // constant case
    dt = max_dt;
    return;
  }
  f64 rdt = dT / max_Q;
  dt = (rdt < min_dt ? min_dt : (rdt > max_dt ? max_dt : rdt));
  // integrate results
  sys.temperatures += Q * dt;
}

inline void update_parameters(system_t &sys) {
  for (body_t b : sys.bodies) {
    // update thermal capacity
    sys.inv_capacities[b.id] =
        b.inv_mass * b.inv_specific_heat_fn(sys.temperatures[b.id]);
    // TODO update thermal conductivity
    b.thermal_conductivy = b.thermal_conductivity_fn(sys.temperatures[b.id]);
  }
}

/**
 * Create a body with a constant specific heat
 * Unit:
 * -----
 *  - mass: kg
 *  - specific_heat: J / (kg * C)
 **/
inline body_t body(const f64 mass, fnT_t specific_heat, fnT_t conductivity,
                   const f64 T0) {
  return {specific_heat, conductivity, 1.0 / mass, T0, conductivity(T0)};
}

// create a body at constant temperature (infinite capcaity)
inline body_t constant_temperature_body(f64 T, f64 conductivity) {
  return {[](f64) -> f64 { return 0; },
          [conductivity](f64) -> f64 { return 1.0 / conductivity; }, 0, T,
          conductivity};
}

// create a conduction heat exchange relationship between 2 bodies
inline exchange_t conduction(body_t &b1, body_t &b2, f64 surf_area_b1,
                             f64 surf_area_b2) {
  return {
      &b1,
      &b2,
      [surf_area_b1, surf_area_b2](f64 T1, f64 T2, body_t *b1,
                                   body_t *b2) -> f64 {
        return (T2 - T1) / (surf_area_b1 / b1->thermal_conductivy +
                            surf_area_b2 / b2->thermal_conductivy);
      },
  };
}

// create a radiation heat exchange relationship between 2 bodies
inline exchange_t radiation(body_t &b1, body_t &b2, f64 correction_factor,
                            f64 view_factor) {
  const f64 alpha = correction_factor * view_factor * SIGMA;
  return {
      &b1,
      &b2,
      [alpha](f64 T1, f64 T2, body_t *b1, body_t *b2) -> f64 {
        return alpha * (_P4_(T2) - _P4_(T1));
      },
  };
}

} // namespace lptm

#endif
