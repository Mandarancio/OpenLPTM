#include "ltm.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

// Stefan-Boltzman constant (W m-2 K-4)
#define sigmaSB 5.670374419e-8

typedef struct {
  f64 T;
  f64 T4;
  f64 H;
  f64 iC;
} _ltm_body_t;

void _b_updt_(struct ltm_body_t *self, f64 dt) {
  _ltm_body_t *b = self->self;
  b->T += b->H * b->iC * dt;
  b->T4 = b->T * b->T * b->T * b->T;
  b->H = 0;
};

void _b_const_(struct ltm_body_t *self, f64 dt) {
  ((_ltm_body_t *)self->self)->H = 0;
}

void _b_aheat_(struct ltm_body_t *self, f64 dH) {
  _ltm_body_t *b = self->self;
  b->H += dH;
}

f64 _b_temp_(struct ltm_body_t *self) { return ((_ltm_body_t *)self->self)->T; }

f64 _b_temp4_(struct ltm_body_t *self) {
  return ((_ltm_body_t *)self->self)->T4;
}

f64 _b_heat_(struct ltm_body_t *self) { return ((_ltm_body_t *)self->self)->H; }

void _b_free_(struct ltm_body_t *self) {
  _ltm_body_t *b = (_ltm_body_t *)self->self;
  free(b);
}

ltm_body_t *ltm_body(const char *label, f64 em, void *payload,
                     void (*update)(struct ltm_body_t *, f64),
                     void (*add_heat)(struct ltm_body_t *, f64),
                     f64 (*temperature)(struct ltm_body_t *),
                     f64 (*temperature4)(struct ltm_body_t *),
                     f64 (*heat)(struct ltm_body_t *),
                     void (*destroy)(struct ltm_body_t *)) {

  ltm_body_t *body = malloc(sizeof(ltm_body_t));
  strcpy(body->label, label);
  body->emissivity = em;
  body->self = payload;
  body->update = update;
  body->add_heat = add_heat;
  body->temperature = temperature;
  body->temperature4 = temperature4;
  body->heat = heat;
  body->destroy = destroy;
  return body;
}

ltm_body_t *ltm_dyn_body(const char *label, f64 em, f64 T0, f64 mass,
                         f64 spec_heat) {
  _ltm_body_t *payload = (_ltm_body_t *)malloc(sizeof(_ltm_body_t));
  payload->T = T0;
  payload->T4 = T0 * T0 * T0 * T0;
  payload->H = 0;
  payload->iC = 1. / (mass * spec_heat);
  return ltm_body(label, em, payload, _b_updt_, _b_aheat_, _b_temp_, _b_temp4_,
                  _b_heat_, _b_free_);
}

ltm_body_t *ltm_capacity_body(const char *label, f64 em, f64 T0, f64 capacity) {
  _ltm_body_t *payload = (_ltm_body_t *)malloc(sizeof(_ltm_body_t));
  payload->T = T0;
  payload->T4 = T0 * T0 * T0 * T0;
  payload->H = 0;
  payload->iC = 1. / capacity;
  return ltm_body(label, em, payload, _b_updt_, _b_aheat_, _b_temp_, _b_temp4_,
                  _b_heat_, _b_free_);
}

ltm_body_t *ltm_const_body(const char *label, f64 em, f64 T0) {
  _ltm_body_t *payload = (_ltm_body_t *)malloc(sizeof(_ltm_body_t));
  payload->T = T0;
  payload->T4 = T0 * T0 * T0 * T0;
  payload->H = 0;
  payload->iC = 0;
  return ltm_body(label, em, payload, _b_const_, _b_aheat_, _b_temp_, _b_temp4_,
                  _b_heat_, _b_free_);
}

typedef struct {
  ltm_body_t *a;
  ltm_body_t *b;
  f64 R;
} _ltm_binary_exc_t;

void _be_free_(struct ltm_exchange_t *self) {
  _ltm_binary_exc_t *e = (_ltm_binary_exc_t *)self->self;
  free(e);
}

void _be_cond_(struct ltm_exchange_t *self) {
  _ltm_binary_exc_t *e = (_ltm_binary_exc_t *)self->self;
  f64 dh = (e->a->temperature(e->a) - e->b->temperature(e->b)) / e->R;
  e->a->add_heat(e->a, -dh);
  e->b->add_heat(e->b, dh);
}

void _be_rad_(struct ltm_exchange_t *self) {
  _ltm_binary_exc_t *e = (_ltm_binary_exc_t *)self->self;
  f64 dh =
      sigmaSB * (e->a->temperature4(e->a) - e->b->temperature4(e->b)) / e->R;
  e->a->add_heat(e->a, -dh);
  e->b->add_heat(e->b, dh);
}

ltm_exchange_t *ltm_conduction(const char *label, ltm_body_t *a, ltm_body_t *b,
                               f64 thermal_resistance) {
  _ltm_binary_exc_t *self = malloc(sizeof(_ltm_binary_exc_t));
  self->a = a;
  self->b = b;
  self->R = thermal_resistance;
  ltm_exchange_t *exc = malloc(sizeof(ltm_exchange_t));
  strcpy(exc->label, label);
  exc->self = self;
  exc->eval = _be_cond_;
  exc->destroy = _be_free_;
  return exc;
}

ltm_exchange_t *ltm_radiation(const char *label, ltm_body_t *a, ltm_body_t *b,
                              f64 thermal_resistance) {
  _ltm_binary_exc_t *self = malloc(sizeof(_ltm_binary_exc_t));
  self->a = a;
  self->b = b;
  self->R = thermal_resistance;
  ltm_exchange_t *exc = malloc(sizeof(ltm_exchange_t));
  strcpy(exc->label, label);
  exc->self = self;
  exc->eval = _be_rad_;
  exc->destroy = _be_free_;
  return exc;
}

ltm_system_t *ltm_system(const char *name) {
  ltm_system_t *sys = malloc(sizeof(ltm_system_t));
  strcpy(sys->label, name);
  sys->bodies = NULL;
  sys->exchanges = NULL;
  sys->n_bodies = 0u;
  sys->n_exchanges = 0u;
  return sys;
}

void ltm_sys_destroy(ltm_system_t *sys) {
  u32 i;
  for (i = 0; i < sys->n_bodies; i++) {
    ltm_body_t *b = sys->bodies[i];
    b->destroy(b);
    free(b);
  }
  for (i = 0; i < sys->n_exchanges; i++) {
    ltm_exchange_t *e = sys->exchanges[i];
    e->destroy(e);
    free(e);
  }
  free(sys->bodies);
  free(sys->exchanges);
  sys->bodies = NULL;
  sys->exchanges = NULL;
}

void *_extend_(void *arr, u32 n, size_t size) {
  if (n == 0) {
    return malloc(size);
  }
  return realloc(arr, (n + 1) * size);
}

void ltm_sys_add_body(ltm_system_t *sys, ltm_body_t *body) {
  sys->bodies = _extend_(sys->bodies, sys->n_bodies, sizeof(ltm_body_t *));
  sys->bodies[sys->n_bodies++] = body;
}

void ltm_sys_add_exchange(ltm_system_t *sys, ltm_exchange_t *exch) {
  sys->exchanges =
      _extend_(sys->exchanges, sys->n_exchanges, sizeof(ltm_exchange_t *));
  sys->exchanges[sys->n_exchanges++] = exch;
}

ltm_body_t *ltm_sys_get_body(ltm_system_t *sys, const char *label) {
  u32 i;
  for (i = 0; i < sys->n_bodies; i++) {
    if (strcmp(sys->bodies[i]->label, label) == 0)
      return sys->bodies[i];
  }
  return 0;
}

ltm_exchange_t *ltm_sys_get_exchange(ltm_system_t *sys, const char *label) {
  u32 i;
  for (i = 0; i < sys->n_exchanges; i++) {
    if (strcmp(sys->exchanges[i]->label, label) == 0)
      return sys->exchanges[i];
  }
  return 0;
}

void ltm_sys_evaluate(ltm_system_t *sys, f64 dt) {
  u32 i;
  for (i = 0; i < sys->n_exchanges; i++) {
    ltm_exchange_t *exch = sys->exchanges[i];
    exch->eval(exch);
  }
  for (i = 0; i < sys->n_bodies; i++) {
    ltm_body_t *body = sys->bodies[i];
    body->update(body, dt);
  }
}

void ltm_sys_print(ltm_system_t *sys) {
  u32 i;
  for (i = 0; i < sys->n_bodies; i++) {
    ltm_body_t *body = sys->bodies[i];
    printf(" - %s: %f\n", body->label, body->temperature(body));
  }
}

void ltm_sys_csv_header(ltm_system_t *sys, FILE *fptr) {
  u32 i;
  for (i = 0; i < sys->n_bodies; i++) {
    ltm_body_t *body = sys->bodies[i];
    fprintf(fptr, "\"%s\",", body->label);
  }
  fprintf(fptr, "\n");
}

void ltm_sys_to_csv(ltm_system_t *sys, FILE *fptr) {
  u32 i;
  for (i = 0; i < sys->n_bodies; i++) {
    ltm_body_t *body = sys->bodies[i];
    fprintf(fptr, "%f,", body->temperature(body));
  }
  fprintf(fptr, "\n");
}

FILE *ltm_sys_bin_file(ltm_system_t *sys, const char *fname) {
  FILE *f = fopen(fname, "wb");
  fwrite(&sys->n_bodies, sizeof(u32), 1, f);
  u32 i;
  for (i = 0; i < sys->n_bodies; i++) {
    ltm_body_t *b = sys->bodies[i];
    u32 len = strlen(b->label);
    fwrite(&len, sizeof(u32), 1, f);
    fwrite(b->label, len, 1, f);
  }
  return f;
}
void ltm_sys_to_bin(ltm_system_t *sys, FILE *fptr) {
  u8 end[2] = {0, 0};
  u32 i;
  for (i = 0; i < sys->n_bodies; i++) {
    ltm_body_t *body = sys->bodies[i];
    f64 t = body->temperature(body);
    fwrite(&t, sizeof(f64), 1, fptr);
  }
  fwrite(&end, 2, 1, fptr);
}
