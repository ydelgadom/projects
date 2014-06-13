#ifndef RANLXD_H
#define RANLXD_H

#ifdef __cplusplus
extern "C" {
#endif

  extern void ranlxd(double r[],int n);
  extern void rlxd_init(int level,int seed);
  extern void rlxd_get(int state[]);
  extern void rlxd_reset(int state[]);

#ifdef __cplusplus
}
#endif

#endif
