#ifndef _LTM_H_
#define _LTM_H_
#include <stdint.h>
#include <stdio.h>

typedef uint8_t u8;
typedef int8_t i8;
typedef uint32_t u32;
typedef int32_t i32;
typedef uint64_t u64;
typedef int64_t i64;
typedef double f64;
typedef float f32;

/**
 * Thermal object abstract body.
 **/
typedef struct ltm_body_t {
  char label[256]; // name of the object
  void *self;      // internal payload
  void (*update)(struct ltm_body_t *,
                 f64); // update temperature using current heat
  void (*add_heat)(struct ltm_body_t *, f64); // add delta heat to the body
  f64 (*temperature)(struct ltm_body_t *);    // get current temperature
  f64 (*temperature4)(
      struct ltm_body_t *); // get current temperature at the fourth power
  f64 (*heat)(struct ltm_body_t *);     // get current heat
  void (*destroy)(struct ltm_body_t *); // body destructor
} ltm_body_t;

/**
 * Generic constructor
 **/
ltm_body_t *ltm_body(const char *label, void *payload,
                     void (*update)(struct ltm_body_t *, f64),
                     void (*add_heat)(struct ltm_body_t *, f64),
                     f64 (*temperature)(struct ltm_body_t *),
                     f64 (*temperature4)(struct ltm_body_t *),
                     f64 (*heat)(struct ltm_body_t *),
                     void (*destroy)(struct ltm_body_t *));

/**
 * Constructor for a dynamic body with a given initial temperature, mass and
 * specific heat
 **/
ltm_body_t *ltm_dyn_body(const char *label, f64 T0, f64 mass, f64 spec_heat);

/**
 * Constructor for a dynamic body with a given initial temperature and
 * equivalent thermal capacity
 **/
ltm_body_t *ltm_capacity_body(const char *label, f64 T0, f64 capacity);

/**
 * Constructor for a constant temperature body with a given temperature
 **/
ltm_body_t *ltm_const_body(const char *label, f64 T0);

/**
 * Generic thermal exchange structure.
 **/
typedef struct ltm_exchange_t {
  void *self;                               // internal data
  void (*eval)(struct ltm_exchange_t *);    // evaluate exchange
  void (*destroy)(struct ltm_exchange_t *); // destructor
} ltm_exchange_t;

/**
 * Constructor for a conduction thermal exchange between to given bodies and
 * with a given equivalent thermal resistance.
 **/
ltm_exchange_t *ltm_conduction(ltm_body_t *a, ltm_body_t *b,
                               f64 thermal_resistance);
/**
 * Constructor for a radiation thermal exchange between to given bodies and
 * with a given equivalent thermal resistance.
 **/
ltm_exchange_t *ltm_radiation(ltm_body_t *a, ltm_body_t *b,
                              f64 thermal_resistance);

/**
 * System structure
 **/
typedef struct {
  char label[256];            // name of the system
  ltm_body_t **bodies;        // list of bodies
  ltm_exchange_t **exchanges; // list of exchanges
  u32 n_bodies;               // number of bodies
  u32 n_exchanges;            // number of exchanges
} ltm_system_t;

/**
 * System constructor
 **/
ltm_system_t *ltm_system(const char *name);

/**
 * System destructor
 **/
void ltm_sys_destroy(ltm_system_t *sys);
/**
 * Add body to a system
 **/
void ltm_sys_add_body(ltm_system_t *sys, ltm_body_t *body);
/**
 * Add an exchange to a system
 **/
void ltm_sys_add_exchange(ltm_system_t *sys, ltm_exchange_t *exch);
/**
 * Evaluate system for a given time interval dt
 **/
void ltm_sys_evaluate(ltm_system_t *sys, f64 dt);

/**********/
/* OUTPUT */
/**********/
/**
 * print state of a system
 **/
void ltm_sys_print(ltm_system_t *sys);
/**
 * write CSV header for a system
 **/
void ltm_sys_csv_header(ltm_system_t *sys, FILE *fptr);
/**
 * write system state as a CSV line
 **/
void ltm_sys_to_csv(ltm_system_t *sys, FILE *fptr);
/**
 * Initialize binary file with system information
 **/
FILE *ltm_sys_bin_file(ltm_system_t *sys, const char *fname);
/**
 * write system state to binary file
 **/
void ltm_sys_to_bin(ltm_system_t *sys, FILE *fptr);
#endif
