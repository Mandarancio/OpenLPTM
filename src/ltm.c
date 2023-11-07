#include "ltm.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  f64 T;
  f64 H;
  f64 iC;
} _ltm_body_t;

void _b_updt_(void *self, f64 dt) {
  _ltm_body_t *b = self;
  b->T += b->H * b->iC * dt;
  b->H = 0;
};

void _b_const_(void *self, f64 dt) { ((_ltm_body_t *)self)->H = 0; }

void _b_aheat_(void *self, f64 dH) {
  _ltm_body_t *b = self;
  b->H += dH;
}

f64 _b_temp_(void *self) { return ((_ltm_body_t *)self)->T; }

f64 _b_heat_(void *self) { return ((_ltm_body_t *)self)->H; }

void _b_free_(void *self) {
  _ltm_body_t *b = (_ltm_body_t *)self;
  free(b);
}

ltm_body_t *ltm_body(const char *label, void *payload,
                     void (*update)(void *, f64), void (*add_heat)(void *, f64),
                     f64 (*temperature)(void *), f64 (*heat)(void *),
                     void (*destroy)(void *)) {

  ltm_body_t *body = malloc(sizeof(ltm_body_t));
  body->label = label;
  body->self = payload;
  body->update = update;
  body->add_heat = add_heat;
  body->temperature = temperature;
  body->heat = heat;
  body->destroy = destroy;
  return body;
}

ltm_body_t *ltm_dyn_body(const char *label, f64 T0, f64 mass, f64 spec_heat) {
  _ltm_body_t *payload = (_ltm_body_t *)malloc(sizeof(_ltm_body_t));
  payload->T = T0;
  payload->H = 0;
  payload->iC = 1. / (mass * spec_heat);
  return ltm_body(label, payload, _b_updt_, _b_aheat_, _b_temp_, _b_heat_,
                  _b_free_);
}

ltm_body_t *ltm_const_body(const char *label, f64 T0) {
  _ltm_body_t *payload = (_ltm_body_t *)malloc(sizeof(_ltm_body_t));
  payload->T = T0;
  payload->H = 0;
  payload->iC = 0;
  return ltm_body(label, payload, _b_const_, _b_aheat_, _b_temp_, _b_heat_,
                  _b_free_);
}

typedef struct {
  ltm_body_t *a;
  ltm_body_t *b;
  f64 R;
} _ltm_binary_exc_t;

void _be_free_(void *self) {
  _ltm_binary_exc_t *e = (_ltm_binary_exc_t *)self;
  free(e);
}

void _be_cond_(void *self) {
  _ltm_binary_exc_t *e = (_ltm_binary_exc_t *)self;
  f64 dh =
      e->R * (e->a->temperature(e->a->self) - e->b->temperature(e->b->self));
  e->a->add_heat(e->a->self, -dh);
  e->b->add_heat(e->b->self, dh);
}

void _be_rad_(void *self) {
  _ltm_binary_exc_t *e = (_ltm_binary_exc_t *)self;
  f64 dh = e->R * (pow(e->a->temperature(e->a->self), 4) -
                   pow(e->b->temperature(e->b->self), 4));
  e->a->add_heat(e->a->self, -dh);
  e->b->add_heat(e->b->self, dh);
}

ltm_exchange_t *ltm_conduction(ltm_body_t *a, ltm_body_t *b,
                               f64 thermal_resistance) {
  _ltm_binary_exc_t *self = malloc(sizeof(_ltm_binary_exc_t));
  self->a = a;
  self->b = b;
  self->R = thermal_resistance;
  ltm_exchange_t *exc = malloc(sizeof(ltm_exchange_t));
  exc->self = self;
  exc->eval = _be_cond_;
  exc->destroy = _be_free_;
  return exc;
}

ltm_exchange_t *ltm_radiation(ltm_body_t *a, ltm_body_t *b,
                              f64 thermal_resistance) {
  _ltm_binary_exc_t *self = malloc(sizeof(_ltm_binary_exc_t));
  self->a = a;
  self->b = b;
  self->R = thermal_resistance;
  ltm_exchange_t *exc = malloc(sizeof(ltm_exchange_t));
  exc->self = self;
  exc->eval = _be_rad_;
  exc->destroy = _be_free_;
  return exc;
}
void ltm_list_destroy(ltm_node_l *n) {
  if (n != NULL && n->next != NULL) {
    ltm_list_destroy(n->next);
    n->next = NULL;
  }
  if (n != NULL)
    free(n);
  n = NULL;
}

void _lst_push_(ltm_node_l *n, ltm_node_l *m) {
  if (n->next == NULL) {
    n->next = m;
  } else {
    _lst_push_(n->next, m);
  }
}

ltm_node_l *ltm_list_append(ltm_node_l *head, void *data) {
  ltm_node_l *n = malloc(sizeof(ltm_node_l));
  n->data = data;
  n->next = NULL;
  if (head == NULL) {
    return n;
  } else {
    _lst_push_(head, n);
  }
  return head;
}

ltm_system_t *ltm_system(const char *name) {
  ltm_system_t *sys = malloc(sizeof(ltm_system_t));
  sys->label = name;
  sys->bodies = NULL;
  sys->exchanges = NULL;
  sys->n_bodies = 0u;
  sys->n_exchanges = 0u;
  return sys;
}

void ltm_sys_destroy(ltm_system_t *sys) {
  ltm_node_l *n;
  for (n = sys->bodies; n != NULL; n = n->next) {
    ltm_body_t *b = n->data;
    b->destroy(b->self);
    free(b);
  }
  for (n = sys->exchanges; n != NULL; n = n->next) {
    ltm_exchange_t *e = n->data;
    e->destroy(e->self);
    free(e);
  }
  ltm_list_destroy(sys->bodies);
  ltm_list_destroy(sys->exchanges);
  sys->bodies = NULL;
  sys->exchanges = NULL;
}

void ltm_sys_add_body(ltm_system_t *sys, ltm_body_t *body) {
  sys->bodies = ltm_list_append(sys->bodies, body);
  sys->n_bodies++;
}

void ltm_sys_add_exchange(ltm_system_t *sys, ltm_exchange_t *exch) {
  sys->exchanges = ltm_list_append(sys->exchanges, exch);
  sys->n_exchanges++;
}

void ltm_sys_evaluate(ltm_system_t *sys, f64 dt) {
  ltm_node_l *it = sys->exchanges;
  for (; it != NULL; it = it->next) {
    ltm_exchange_t *exch = it->data;
    exch->eval(exch->self);
  }
  it = sys->bodies;
  for (; it != NULL; it = it->next) {
    ltm_body_t *body = it->data;
    body->update(body->self, dt);
  }
}

void ltm_sys_print(ltm_system_t *sys) {
  ltm_node_l *it = sys->bodies;
  for (; it != NULL; it = it->next) {
    ltm_body_t *body = it->data;
    printf(" - %s: %f\n", body->label, body->temperature(body->self));
  }
}

void ltm_sys_csv_header(ltm_system_t *sys, FILE *fptr) {
  ltm_node_l *it = sys->bodies;
  for (; it != NULL; it = it->next) {
    ltm_body_t *body = it->data;
    fprintf(fptr, "'%s',", body->label);
  }
  fprintf(fptr, "\n");
}

void ltm_sys_to_csv(ltm_system_t *sys, FILE *fptr) {
  ltm_node_l *it = sys->bodies;
  for (; it != NULL; it = it->next) {
    ltm_body_t *body = it->data;
    fprintf(fptr, "%f,", body->temperature(body->self));
  }
  fprintf(fptr, "\n");
}

FILE *ltm_sys_bin_file(ltm_system_t *sys, const char *fname) {
  FILE *f = fopen(fname, "wb");
  fwrite(&sys->n_bodies, sizeof(u32), 1, f);
  ltm_node_l *it;
  for (it = sys->bodies; it != NULL; it = it->next) {
    ltm_body_t *b = it->data;
    u32 len = strlen(b->label);
    fwrite(&len, sizeof(u32), 1, f);
    fwrite(b->label, len, 1, f);
  }
  return f;
}
void ltm_sys_to_bin(ltm_system_t *sys, FILE *fptr) {
  u8 end[2] = {0, 0};
  ltm_node_l *it = sys->bodies;
  for (; it != NULL; it = it->next) {
    ltm_body_t *body = it->data;
    f64 t = body->temperature(body->self);
    fwrite(&t, sizeof(f64), 1, fptr);
  }
  fwrite(&end, 2, 1, fptr);
}
