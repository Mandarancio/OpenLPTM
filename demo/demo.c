#include "ltm.h"
#include <stdio.h>
#include <string.h>
#include <time.h>

typedef enum {
  LTM_F_NONE,
  LTM_F_CSV,
  LTM_F_BIN,
} LTM_F_MODE;

int main(int argc, char **argv) {
  union {
    FILE *ptr;
  } file;
  LTM_F_MODE fmode = argc == 1 ? LTM_F_NONE : LTM_F_CSV;
  if (argc > 2) {
    if (strcmp(argv[1], "-b") == 0 || strcmp(argv[1], "--bin") == 0) {
      fmode = LTM_F_BIN;
    }
  }

  ltm_body_t *a = ltm_dyn_body("Body A", 270, 0.5, 400);
  ltm_body_t *b = ltm_dyn_body("Body B", 600, 0.5, 800.);

  ltm_exchange_t *exch = ltm_conduction(a, b, 1.4);

  ltm_system_t *sys = ltm_system("demo");
  ltm_sys_add_body(sys, a);
  ltm_sys_add_body(sys, b);
  ltm_sys_add_exchange(sys, exch);
  int N = 1000000;
  if (fmode == LTM_F_BIN) {
    file.ptr = ltm_sys_bin_file(sys, argv[argc - 1]);
  } else if (fmode == LTM_F_CSV) {
    file.ptr = fopen(argv[argc - 1], "w");
    ltm_sys_csv_header(sys, file.ptr);
  }
  int i = 0;
  clock_t start = clock();
  for (i = 0; i < N - 1; i++) {
    ltm_sys_evaluate(sys, 0.001);
    if (fmode == LTM_F_BIN)
      ltm_sys_to_bin(sys, file.ptr);
    else if (fmode == LTM_F_CSV)
      ltm_sys_to_csv(sys, file.ptr);
  }
  clock_t end = clock();
  f64 duration = (f64)(end - start) * 1000. / (f64)CLOCKS_PER_SEC;
  printf("duration: %f ms\n", duration);
  if (fmode != LTM_F_NONE)
    fclose(file.ptr);

  ltm_sys_destroy(sys);
  return 0;
}
