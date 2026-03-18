#ifndef NTT_H
#define NTT_H

#include <stdint.h>
#include "params.h"

#define zetas KYBER_NAMESPACE(zetas)
extern const int16_t zetas[128];

#define ntt KYBER_NAMESPACE(ntt)
void ntt(int16_t poly[KYBER_N]);

#define invntt KYBER_NAMESPACE(invntt)
void invntt(int16_t poly[KYBER_N]);

#define basemul KYBER_NAMESPACE(basemul)
#if (KYBER_N == 64) || (KYBER_N == 128)
void basemul(int16_t r[1], const int16_t a[1], const int16_t b[1]);
#elif KYBER_N == 256
void basemul(int16_t r[2], const int16_t a[2], const int16_t b[2], int16_t zeta);
#else
  #error "KYBER_N not supported"
#endif

#endif
