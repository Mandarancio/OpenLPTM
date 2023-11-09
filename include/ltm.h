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

typedef struct {
  char label[256];
  void *self;
  void (*update)(void *, f64);
  void (*add_heat)(void *, f64);
  f64 (*temperature)(void *);
  f64 (*heat)(void *);
  void (*destroy)(void *);
} ltm_body_t;

ltm_body_t *ltm_body(const char *label, void *payload,
                     void (*update)(void *, f64), void (*add_heat)(void *, f64),
                     f64 (*temperature)(void *), f64 (*heat)(void *),
                     void (*destroy)(void *));

ltm_body_t *ltm_dyn_body(const char *label, f64 T0, f64 mass, f64 spec_heat);

ltm_body_t *ltm_const_body(const char *label, f64 T0);

typedef struct {
  void *self;
  void (*eval)(void *);
  void (*destroy)(void *);
} ltm_exchange_t;

ltm_exchange_t *ltm_conduction(ltm_body_t *a, ltm_body_t *b,
                               f64 thermal_resistance);
ltm_exchange_t *ltm_radiation(ltm_body_t *a, ltm_body_t *b,
                              f64 thermal_resistance);



typedef struct {
  char label[256];
  ltm_body_t **bodies;
  ltm_exchange_t **exchanges;
  u32 n_bodies;
  u32 n_exchanges;
} ltm_system_t;

ltm_system_t *ltm_system(const char *name);
void ltm_sys_destroy(ltm_system_t *sys);
void ltm_sys_add_body(ltm_system_t *sys, ltm_body_t *body);
void ltm_sys_add_exchange(ltm_system_t *sys, ltm_exchange_t *exch);
void ltm_sys_evaluate(ltm_system_t *sys, f64 dt);
void ltm_sys_print(ltm_system_t *sys);
void ltm_sys_csv_header(ltm_system_t *sys, FILE *fptr);
void ltm_sys_to_csv(ltm_system_t *sys, FILE *fptr) ;
FILE *ltm_sys_bin_file(ltm_system_t *sys, const char *fname);
void ltm_sys_to_bin(ltm_system_t *sys, FILE *fptr);
#endif
