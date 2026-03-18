#include <stddef.h>
#include <stdint.h>
#include <immintrin.h>
#include <string.h>
#include "align.h"
#include "params.h"
#include "indcpa.h"
#include "polyvec.h"
#include "poly.h"
#include "ntt.h"
#include "cbd.h"
#include "rejsample.h"
#include "symmetric.h"
#include "randombytes.h"

/* NTT(y+2) in Montgomery domain = NTT(2285y + 2285) */
#ifdef PRECOMPUTE_TWIST
  poly precomputed_twist = { .coeffs = {
	1241, 1241, 1241, 1241, 1241, 1241, 1241, 1241, 
	1241, 1241, 1241, 1241, 1241, 1241, 1241, 1241, 
	2285, 2285, 2285, 2285, 2285, 2285, 2285, 2285, 
	2285, 2285, 2285, 2285, 2285, 2285, 2285, 2285, 
	1241, 1241, 1241, 1241, 1241, 1241, 1241, 1241, 
	1241, 1241, 1241, 1241, 1241, 1241, 1241, 1241, 
	2285, 2285, 2285, 2285, 2285, 2285, 2285, 2285, 
	2285, 2285, 2285, 2285, 2285, 2285, 2285, 2285, 
	1241, 1241, 1241, 1241, 1241, 1241, 1241, 1241, 
	1241, 1241, 1241, 1241, 1241, 1241, 1241, 1241, 
	2285, 2285, 2285, 2285, 2285, 2285, 2285, 2285, 
	2285, 2285, 2285, 2285, 2285, 2285, 2285, 2285, 
	1241, 1241, 1241, 1241, 1241, 1241, 1241, 1241, 
	1241, 1241, 1241, 1241, 1241, 1241, 1241, 1241, 
	2285, 2285, 2285, 2285, 2285, 2285, 2285, 2285, 
	2285, 2285, 2285, 2285, 2285, 2285, 2285, 2285, 
	1241, 1241, 1241, 1241, 1241, 1241, 1241, 1241, 
	1241, 1241, 1241, 1241, 1241, 1241, 1241, 1241, 
	2285, 2285, 2285, 2285, 2285, 2285, 2285, 2285, 
	2285, 2285, 2285, 2285, 2285, 2285, 2285, 2285, 
	1241, 1241, 1241, 1241, 1241, 1241, 1241, 1241, 
	1241, 1241, 1241, 1241, 1241, 1241, 1241, 1241, 
	2285, 2285, 2285, 2285, 2285, 2285, 2285, 2285, 
	2285, 2285, 2285, 2285, 2285, 2285, 2285, 2285, 
	1241, 1241, 1241, 1241, 1241, 1241, 1241, 1241, 
	1241, 1241, 1241, 1241, 1241, 1241, 1241, 1241, 
	2285, 2285, 2285, 2285, 2285, 2285, 2285, 2285, 
	2285, 2285, 2285, 2285, 2285, 2285, 2285, 2285, 
	1241, 1241, 1241, 1241, 1241, 1241, 1241, 1241, 
	1241, 1241, 1241, 1241, 1241, 1241, 1241, 1241, 
	2285, 2285, 2285, 2285, 2285, 2285, 2285, 2285, 
	2285, 2285, 2285, 2285, 2285, 2285, 2285, 2285,
  }};
#endif // PRECOMPUTE_TWIST

/*************************************************
* Name:        pack_pk
*
* Description: Serialize the public key as concatenation of the
*              serialized vector of polynomials pk and the
*              public seed used to generate the matrix A.
*              The polynomial coefficients in pk are assumed to
*              lie in the invertal [0,q], i.e. pk must be reduced
*              by polyvec_reduce().
*
* Arguments:   uint8_t *r: pointer to the output serialized public key
*              polyvec *pk: pointer to the input public-key polyvec
*              const uint8_t *seed: pointer to the input public seed
**************************************************/
static void pack_pk(uint8_t r[KYBER_INDCPA_PUBLICKEYBYTES],
                    polyvec *pk,
                    const uint8_t seed[KYBER_SYMBYTES])
{
  polyvec_tobytes(r, pk);
  memcpy(r+KYBER_POLYVECBYTES, seed, KYBER_SYMBYTES);
}

/*************************************************
* Name:        unpack_pk
*
* Description: De-serialize public key from a byte array;
*              approximate inverse of pack_pk
*
* Arguments:   - polyvec *pk: pointer to output public-key polynomial vector
*              - uint8_t *seed: pointer to output seed to generate matrix A
*              - const uint8_t *packedpk: pointer to input serialized public key
**************************************************/
static void unpack_pk(polyvec *pk,
                      uint8_t seed[KYBER_SYMBYTES],
                      const uint8_t packedpk[KYBER_INDCPA_PUBLICKEYBYTES])
{
  polyvec_frombytes(pk, packedpk);
  memcpy(seed, packedpk+KYBER_POLYVECBYTES, KYBER_SYMBYTES);
}

/*************************************************
* Name:        pack_sk
*
* Description: Serialize the secret key.
*              The polynomial coefficients in sk are assumed to
*              lie in the invertal [0,q], i.e. sk must be reduced
*              by polyvec_reduce().
*
* Arguments:   - uint8_t *r: pointer to output serialized secret key
*              - polyvec *sk: pointer to input vector of polynomials (secret key)
**************************************************/
static void pack_sk(uint8_t r[KYBER_INDCPA_SECRETKEYBYTES], polyvec *sk)
{
  polyvec_tobytes(r, sk);
}

/*************************************************
* Name:        unpack_sk
*
* Description: De-serialize the secret key; inverse of pack_sk
*
* Arguments:   - polyvec *sk: pointer to output vector of polynomials (secret key)
*              - const uint8_t *packedsk: pointer to input serialized secret key
**************************************************/
static void unpack_sk(polyvec *sk, const uint8_t packedsk[KYBER_INDCPA_SECRETKEYBYTES])
{
  polyvec_frombytes(sk, packedsk);
}

/*************************************************
* Name:        pack_ciphertext
*
* Description: Serialize the ciphertext as concatenation of the
*              compressed and serialized vector of polynomials b
*              and the compressed and serialized polynomial v.
*              The polynomial coefficients in b and v are assumed to
*              lie in the invertal [0,q], i.e. b and v must be reduced
*              by polyvec_reduce() and poly_reduce(), respectively.
*
* Arguments:   uint8_t *r: pointer to the output serialized ciphertext
*              poly *pk: pointer to the input vector of polynomials b
*              poly *v: pointer to the input polynomial v
**************************************************/
static void pack_ciphertext(uint8_t r[KYBER_INDCPA_BYTES], polyvec *b, poly *v)
{
  polyvec_compress(r, b);
  poly_compress(r+KYBER_POLYVECCOMPRESSEDBYTES, v);
}

/*************************************************
* Name:        unpack_ciphertext
*
* Description: De-serialize and decompress ciphertext from a byte array;
*              approximate inverse of pack_ciphertext
*
* Arguments:   - polyvec *b: pointer to the output vector of polynomials b
*              - poly *v: pointer to the output polynomial v
*              - const uint8_t *c: pointer to the input serialized ciphertext
**************************************************/
static void unpack_ciphertext(polyvec *b, poly *v, const uint8_t c[KYBER_INDCPA_BYTES])
{
  polyvec_decompress(b, c);
  poly_decompress(v, c+KYBER_POLYVECCOMPRESSEDBYTES);
}

/*************************************************
* Name:        rej_uniform
*
* Description: Run rejection sampling on uniform random bytes to generate
*              uniform random integers mod q
*
* Arguments:   - int16_t *r: pointer to output array
*              - unsigned int len: requested number of 16-bit integers (uniform mod q)
*              - const uint8_t *buf: pointer to input buffer (assumed to be uniformly random bytes)
*              - unsigned int buflen: length of input buffer in bytes
*
* Returns number of sampled 16-bit integers (at most len)
**************************************************/
static unsigned int rej_uniform(int16_t *r,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen)
{
  unsigned int ctr, pos;
  uint16_t val0, val1;

  ctr = pos = 0;
  while(ctr < len && pos <= buflen - 3) {  // buflen is always at least 3
    val0 = ((buf[pos+0] >> 0) | ((uint16_t)buf[pos+1] << 8)) & 0xFFF;
    val1 = ((buf[pos+1] >> 4) | ((uint16_t)buf[pos+2] << 4)) & 0xFFF;
    pos += 3;

    if(val0 < KYBER_Q)
      r[ctr++] = val0;
    if(ctr < len && val1 < KYBER_Q)
      r[ctr++] = val1;
  }

  return ctr;
}

#define gen_a(A,B)  gen_matrix(A,B,0)
#define gen_at(A,B) gen_matrix(A,B,1)

/*************************************************
* Name:        gen_matrix
*
* Description: Deterministically generate matrix A (or the transpose of A)
*              from a seed. Entries of the matrix are polynomials that look
*              uniformly random. Performs rejection sampling on output of
*              a XOF
*
* Arguments:   - polyvec *a: pointer to ouptput matrix A
*              - const uint8_t *seed: pointer to input seed
*              - int transposed: boolean deciding whether A or A^T is generated
**************************************************/
#if KYBER_K == 2
void gen_matrix(polyvec *a, const uint8_t seed[32], int transposed)
{
  unsigned int ctr0, ctr1, ctr2, ctr3;
  ALIGNED_UINT8(REJ_UNIFORM_AVX_NBLOCKS*SHAKE128_RATE) buf[4];
  __m256i f;
  keccakx4_state state;

  f = _mm256_loadu_si256((__m256i *)seed);
  _mm256_store_si256(buf[0].vec, f);
  _mm256_store_si256(buf[1].vec, f);
  _mm256_store_si256(buf[2].vec, f);
  _mm256_store_si256(buf[3].vec, f);

  if(transposed) {
    buf[0].coeffs[32] = 0;
    buf[0].coeffs[33] = 0;
    buf[1].coeffs[32] = 0;
    buf[1].coeffs[33] = 1;
    buf[2].coeffs[32] = 1;
    buf[2].coeffs[33] = 0;
    buf[3].coeffs[32] = 1;
    buf[3].coeffs[33] = 1;
  }
  else {
    buf[0].coeffs[32] = 0;
    buf[0].coeffs[33] = 0;
    buf[1].coeffs[32] = 1;
    buf[1].coeffs[33] = 0;
    buf[2].coeffs[32] = 0;
    buf[2].coeffs[33] = 1;
    buf[3].coeffs[32] = 1;
    buf[3].coeffs[33] = 1;
  }

  shake128x4_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 34);
  shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, REJ_UNIFORM_AVX_NBLOCKS, &state);

  ctr0 = rej_uniform_avx(a[0].vec[0].coeffs, buf[0].coeffs, KYBER_N);
  ctr1 = rej_uniform_avx(a[0].vec[1].coeffs, buf[1].coeffs, KYBER_N);
  ctr2 = rej_uniform_avx(a[1].vec[0].coeffs, buf[2].coeffs, KYBER_N);
  ctr3 = rej_uniform_avx(a[1].vec[1].coeffs, buf[3].coeffs, KYBER_N);

  while(ctr0 < KYBER_N || ctr1 < KYBER_N || ctr2 < KYBER_N || ctr3 < KYBER_N) {
    shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 1, &state);

    ctr0 += rej_uniform(a[0].vec[0].coeffs + ctr0, KYBER_N - ctr0, buf[0].coeffs, SHAKE128_RATE);
    ctr1 += rej_uniform(a[0].vec[1].coeffs + ctr1, KYBER_N - ctr1, buf[1].coeffs, SHAKE128_RATE);
    ctr2 += rej_uniform(a[1].vec[0].coeffs + ctr2, KYBER_N - ctr2, buf[2].coeffs, SHAKE128_RATE);
    ctr3 += rej_uniform(a[1].vec[1].coeffs + ctr3, KYBER_N - ctr3, buf[3].coeffs, SHAKE128_RATE);
  }

  poly_nttunpack(&a[0].vec[0]);
  poly_nttunpack(&a[0].vec[1]);
  poly_nttunpack(&a[1].vec[0]);
  poly_nttunpack(&a[1].vec[1]);
}
#elif KYBER_K == 3
void gen_matrix(polyvec *a, const uint8_t seed[32], int transposed)
{
  unsigned int ctr0, ctr1, ctr2, ctr3;
  ALIGNED_UINT8(REJ_UNIFORM_AVX_NBLOCKS*SHAKE128_RATE) buf[4];
  __m256i f;
  keccakx4_state state;
  keccak_state state1x;

  f = _mm256_loadu_si256((__m256i *)seed);
  _mm256_store_si256(buf[0].vec, f);
  _mm256_store_si256(buf[1].vec, f);
  _mm256_store_si256(buf[2].vec, f);
  _mm256_store_si256(buf[3].vec, f);

  if(transposed) {
    buf[0].coeffs[32] = 0;
    buf[0].coeffs[33] = 0;
    buf[1].coeffs[32] = 0;
    buf[1].coeffs[33] = 1;
    buf[2].coeffs[32] = 0;
    buf[2].coeffs[33] = 2;
    buf[3].coeffs[32] = 1;
    buf[3].coeffs[33] = 0;
  }
  else {
    buf[0].coeffs[32] = 0;
    buf[0].coeffs[33] = 0;
    buf[1].coeffs[32] = 1;
    buf[1].coeffs[33] = 0;
    buf[2].coeffs[32] = 2;
    buf[2].coeffs[33] = 0;
    buf[3].coeffs[32] = 0;
    buf[3].coeffs[33] = 1;
  }

  shake128x4_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 34);
  shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, REJ_UNIFORM_AVX_NBLOCKS, &state);

  ctr0 = rej_uniform_avx(a[0].vec[0].coeffs, buf[0].coeffs, KYBER_N);
  ctr1 = rej_uniform_avx(a[0].vec[1].coeffs, buf[1].coeffs, KYBER_N);
  ctr2 = rej_uniform_avx(a[0].vec[2].coeffs, buf[2].coeffs, KYBER_N);
  ctr3 = rej_uniform_avx(a[1].vec[0].coeffs, buf[3].coeffs, KYBER_N);

  while(ctr0 < KYBER_N || ctr1 < KYBER_N || ctr2 < KYBER_N || ctr3 < KYBER_N) {
    shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 1, &state);

    ctr0 += rej_uniform(a[0].vec[0].coeffs + ctr0, KYBER_N - ctr0, buf[0].coeffs, SHAKE128_RATE);
    ctr1 += rej_uniform(a[0].vec[1].coeffs + ctr1, KYBER_N - ctr1, buf[1].coeffs, SHAKE128_RATE);
    ctr2 += rej_uniform(a[0].vec[2].coeffs + ctr2, KYBER_N - ctr2, buf[2].coeffs, SHAKE128_RATE);
    ctr3 += rej_uniform(a[1].vec[0].coeffs + ctr3, KYBER_N - ctr3, buf[3].coeffs, SHAKE128_RATE);
  }

  poly_nttunpack(&a[0].vec[0]);
  poly_nttunpack(&a[0].vec[1]);
  poly_nttunpack(&a[0].vec[2]);
  poly_nttunpack(&a[1].vec[0]);

  f = _mm256_loadu_si256((__m256i *)seed);
  _mm256_store_si256(buf[0].vec, f);
  _mm256_store_si256(buf[1].vec, f);
  _mm256_store_si256(buf[2].vec, f);
  _mm256_store_si256(buf[3].vec, f);

  if(transposed) {
    buf[0].coeffs[32] = 1;
    buf[0].coeffs[33] = 1;
    buf[1].coeffs[32] = 1;
    buf[1].coeffs[33] = 2;
    buf[2].coeffs[32] = 2;
    buf[2].coeffs[33] = 0;
    buf[3].coeffs[32] = 2;
    buf[3].coeffs[33] = 1;
  }
  else {
    buf[0].coeffs[32] = 1;
    buf[0].coeffs[33] = 1;
    buf[1].coeffs[32] = 2;
    buf[1].coeffs[33] = 1;
    buf[2].coeffs[32] = 0;
    buf[2].coeffs[33] = 2;
    buf[3].coeffs[32] = 1;
    buf[3].coeffs[33] = 2;
  }

  shake128x4_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 34);
  shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, REJ_UNIFORM_AVX_NBLOCKS, &state);

  ctr0 = rej_uniform_avx(a[1].vec[1].coeffs, buf[0].coeffs, KYBER_N);
  ctr1 = rej_uniform_avx(a[1].vec[2].coeffs, buf[1].coeffs, KYBER_N);
  ctr2 = rej_uniform_avx(a[2].vec[0].coeffs, buf[2].coeffs, KYBER_N);
  ctr3 = rej_uniform_avx(a[2].vec[1].coeffs, buf[3].coeffs, KYBER_N);

  while(ctr0 < KYBER_N || ctr1 < KYBER_N || ctr2 < KYBER_N || ctr3 < KYBER_N) {
    shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 1, &state);

    ctr0 += rej_uniform(a[1].vec[1].coeffs + ctr0, KYBER_N - ctr0, buf[0].coeffs, SHAKE128_RATE);
    ctr1 += rej_uniform(a[1].vec[2].coeffs + ctr1, KYBER_N - ctr1, buf[1].coeffs, SHAKE128_RATE);
    ctr2 += rej_uniform(a[2].vec[0].coeffs + ctr2, KYBER_N - ctr2, buf[2].coeffs, SHAKE128_RATE);
    ctr3 += rej_uniform(a[2].vec[1].coeffs + ctr3, KYBER_N - ctr3, buf[3].coeffs, SHAKE128_RATE);
  }

  poly_nttunpack(&a[1].vec[1]);
  poly_nttunpack(&a[1].vec[2]);
  poly_nttunpack(&a[2].vec[0]);
  poly_nttunpack(&a[2].vec[1]);

  f = _mm256_loadu_si256((__m256i *)seed);
  _mm256_store_si256(buf[0].vec, f);
  buf[0].coeffs[32] = 2;
  buf[0].coeffs[33] = 2;
  shake128_absorb_once(&state1x, buf[0].coeffs, 34);
  shake128_squeezeblocks(buf[0].coeffs, REJ_UNIFORM_AVX_NBLOCKS, &state1x);
  ctr0 = rej_uniform_avx(a[2].vec[2].coeffs, buf[0].coeffs, KYBER_N);
  while(ctr0 < KYBER_N) {
    shake128_squeezeblocks(buf[0].coeffs, 1, &state1x);
    ctr0 += rej_uniform(a[2].vec[2].coeffs + ctr0, KYBER_N - ctr0, buf[0].coeffs, SHAKE128_RATE);
  }

  poly_nttunpack(&a[2].vec[2]);
}
#elif KYBER_K == 4
void gen_matrix(polyvec *a, const uint8_t seed[32], int transposed)
{
  unsigned int i, ctr0, ctr1, ctr2, ctr3;
  ALIGNED_UINT8(REJ_UNIFORM_AVX_NBLOCKS*SHAKE128_RATE) buf[4];
  __m256i f;
  keccakx4_state state;

  for(i=0;i<4;i++) {
    f = _mm256_loadu_si256((__m256i *)seed);
    _mm256_store_si256(buf[0].vec, f);
    _mm256_store_si256(buf[1].vec, f);
    _mm256_store_si256(buf[2].vec, f);
    _mm256_store_si256(buf[3].vec, f);

    if(transposed) {
      buf[0].coeffs[32] = i;
      buf[0].coeffs[33] = 0;
      buf[1].coeffs[32] = i;
      buf[1].coeffs[33] = 1;
      buf[2].coeffs[32] = i;
      buf[2].coeffs[33] = 2;
      buf[3].coeffs[32] = i;
      buf[3].coeffs[33] = 3;
    }
    else {
      buf[0].coeffs[32] = 0;
      buf[0].coeffs[33] = i;
      buf[1].coeffs[32] = 1;
      buf[1].coeffs[33] = i;
      buf[2].coeffs[32] = 2;
      buf[2].coeffs[33] = i;
      buf[3].coeffs[32] = 3;
      buf[3].coeffs[33] = i;
    }

    shake128x4_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 34);
    shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, REJ_UNIFORM_AVX_NBLOCKS, &state);

    ctr0 = rej_uniform_avx(a[i].vec[0].coeffs, buf[0].coeffs, KYBER_N);
    ctr1 = rej_uniform_avx(a[i].vec[1].coeffs, buf[1].coeffs, KYBER_N);
    ctr2 = rej_uniform_avx(a[i].vec[2].coeffs, buf[2].coeffs, KYBER_N);
    ctr3 = rej_uniform_avx(a[i].vec[3].coeffs, buf[3].coeffs, KYBER_N);

    while(ctr0 < KYBER_N || ctr1 < KYBER_N || ctr2 < KYBER_N || ctr3 < KYBER_N) {
      shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 1, &state);

      ctr0 += rej_uniform(a[i].vec[0].coeffs + ctr0, KYBER_N - ctr0, buf[0].coeffs, SHAKE128_RATE);
      ctr1 += rej_uniform(a[i].vec[1].coeffs + ctr1, KYBER_N - ctr1, buf[1].coeffs, SHAKE128_RATE);
      ctr2 += rej_uniform(a[i].vec[2].coeffs + ctr2, KYBER_N - ctr2, buf[2].coeffs, SHAKE128_RATE);
      ctr3 += rej_uniform(a[i].vec[3].coeffs + ctr3, KYBER_N - ctr3, buf[3].coeffs, SHAKE128_RATE);
    }

    poly_nttunpack(&a[i].vec[0]);
    poly_nttunpack(&a[i].vec[1]);
    poly_nttunpack(&a[i].vec[2]);
    poly_nttunpack(&a[i].vec[3]);
  }
}
#endif

static void compute_twist(poly *r)
{
  #ifdef PRECOMPUTE_TWIST
    *r = precomputed_twist;
  #else
    // Initialize the polynomial y+2 in Montgomery domain
    r->coeffs[0] = 1241; // 2
    r->coeffs[1] = 2285; // 1
    for(unsigned int i = 2; i < KYBER_N; i++) r->coeffs[i] = 0;
    poly_ntt(r);
  #endif

  // Convenience for printing and storing as precomputed version
  //printf("twist\n");
  //for(unsigned int i = 0; i < KYBER_N/8; i++) {
  //  for(unsigned int j = 0; j < 8; j++) {
  //    printf("%d, ", r->coeffs[8*i+j]);
  //  }
  //  printf("\n");
  //}
  //printf("\n");
}

#if KYBER_K == 2
void gen_matrix_ring(polyvec *a, const uint8_t seed[KYBER_SYMBYTES])
{
  unsigned int i, j;
  unsigned int ctr0, ctr1, ctr2, ctr3;
  ALIGNED_UINT8(REJ_UNIFORM_AVX_NBLOCKS*SHAKE128_RATE) buf[4];
  __m256i f;
  keccakx4_state state;
  polyvec at;

  // Sample as 4 x 128-bit coefficients to use SHAKE better

  f = _mm256_loadu_si256((__m256i *)seed);
  _mm256_store_si256(buf[0].vec, f);
  _mm256_store_si256(buf[1].vec, f);
  _mm256_store_si256(buf[2].vec, f);
  _mm256_store_si256(buf[3].vec, f);

  buf[0].coeffs[32] = 0;
  buf[0].coeffs[33] = 0;
  buf[1].coeffs[32] = 1;
  buf[1].coeffs[33] = 0;
  buf[2].coeffs[32] = 0;
  buf[2].coeffs[33] = 1;
  buf[3].coeffs[32] = 1;
  buf[3].coeffs[33] = 1;

  shake128x4_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 34);
  shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, REJ_UNIFORM_AVX_NBLOCKS, &state);

  ctr0 = rej_uniform_avx(a[0].vec[0].coeffs, buf[0].coeffs, KYBER_N/2);
  ctr1 = rej_uniform_avx(a[1].vec[0].coeffs, buf[1].coeffs, KYBER_N/2);
  ctr2 = rej_uniform_avx(a[0].vec[0].coeffs + KYBER_N/2, buf[2].coeffs, KYBER_N/2);
  ctr3 = rej_uniform_avx(a[1].vec[0].coeffs + KYBER_N/2, buf[3].coeffs, KYBER_N/2);

  while(ctr0 < KYBER_N/2 || ctr1 < KYBER_N/2 || ctr2 < KYBER_N/2 || ctr3 < KYBER_N/2) {
    shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 1, &state);

    ctr0 += rej_uniform(a[0].vec[0].coeffs + ctr0, KYBER_N/2 - ctr0, buf[0].coeffs, SHAKE128_RATE);
    ctr1 += rej_uniform(a[1].vec[0].coeffs + ctr1, KYBER_N/2 - ctr1, buf[1].coeffs, SHAKE128_RATE);
    ctr2 += rej_uniform(a[0].vec[0].coeffs + KYBER_N/2 + ctr2, KYBER_N/2 - ctr2, buf[2].coeffs, SHAKE128_RATE);
    ctr3 += rej_uniform(a[1].vec[0].coeffs + KYBER_N/2 + ctr3, KYBER_N/2 - ctr3, buf[3].coeffs, SHAKE128_RATE);
  }

  poly_nttunpack(&a[0].vec[0]);
  poly_nttunpack(&a[1].vec[0]);

  /* Copy for the full non-twisted matrix */
  for(i=1;i<KYBER_K;i++) { // column
    for(j=i;j<KYBER_K;j++) { // row
      memcpy(&a[j].vec[i], &a[j-i].vec[0], 2*KYBER_N);
    }
  }

  // Can be pre-computed, but minor impact for AVX2
  /* Compute NTT(y+1) */
  poly twist;
  compute_twist(&twist);

  /* Generate the K-1 twists (multiplied by (y+1)) */
  for(i=1;i<KYBER_K;i++) {
    poly_basemul_montgomery(&at.vec[i], &a[i].vec[0], &twist); // Montgomery factor cancels out
  }

  /* Copy for the full twisted matrix */
  for(i=0;i<KYBER_K-1;i++) { // row
    for(j=i+1;j<KYBER_K;j++) { // column
      memcpy(&a[i].vec[j], &at.vec[KYBER_K+i-j], 2*KYBER_N);
    }
  }
}
#elif KYBER_K == 3
void gen_matrix_ring(polyvec *a, const uint8_t seed[KYBER_SYMBYTES])
{
  unsigned int i, j;
  unsigned int ctr0, ctr1, ctr2;
  ALIGNED_UINT8(REJ_UNIFORM_AVX_NBLOCKS*SHAKE128_RATE) buf[4];
  __m256i f;
  keccakx4_state state;
  polyvec at;

  f = _mm256_loadu_si256((__m256i *)seed);
  _mm256_store_si256(buf[0].vec, f);
  _mm256_store_si256(buf[1].vec, f);
  _mm256_store_si256(buf[2].vec, f);
  _mm256_store_si256(buf[3].vec, f);

  buf[0].coeffs[32] = 0;
  buf[0].coeffs[33] = 0;
  buf[1].coeffs[32] = 1;
  buf[1].coeffs[33] = 0;
  buf[2].coeffs[32] = 2;
  buf[2].coeffs[33] = 0;
  buf[3].coeffs[32] = 0;
  buf[3].coeffs[33] = 1;

  shake128x4_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 34);
  shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, REJ_UNIFORM_AVX_NBLOCKS, &state);

  ctr0 = rej_uniform_avx(a[0].vec[0].coeffs, buf[0].coeffs, KYBER_N);
  ctr1 = rej_uniform_avx(a[1].vec[0].coeffs, buf[1].coeffs, KYBER_N);
  ctr2 = rej_uniform_avx(a[2].vec[0].coeffs, buf[2].coeffs, KYBER_N);

  while(ctr0 < KYBER_N || ctr1 < KYBER_N || ctr2 < KYBER_N) {
    shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 1, &state);

    ctr0 += rej_uniform(a[0].vec[0].coeffs + ctr0, KYBER_N - ctr0, buf[0].coeffs, SHAKE128_RATE);
    ctr1 += rej_uniform(a[1].vec[0].coeffs + ctr1, KYBER_N - ctr1, buf[1].coeffs, SHAKE128_RATE);
    ctr2 += rej_uniform(a[2].vec[0].coeffs + ctr2, KYBER_N - ctr2, buf[2].coeffs, SHAKE128_RATE);
  }

  poly_nttunpack(&a[0].vec[0]);
  poly_nttunpack(&a[1].vec[0]);
  poly_nttunpack(&a[2].vec[0]);

  /* Copy for the full non-twisted matrix */
  for(i=1;i<KYBER_K;i++) { // column
    for(j=i;j<KYBER_K;j++) { // row
      memcpy(&a[j].vec[i], &a[j-i].vec[0], 2*KYBER_N);
    }
  }

  // Can be pre-computed, but minor impact for AVX2
  /* Compute NTT(y+1) */
  poly twist;
  compute_twist(&twist);

  /* Generate the K-1 twists (multiplied by (y+1)) */
  for(i=1;i<KYBER_K;i++) {
    poly_basemul_montgomery(&at.vec[i], &a[i].vec[0], &twist); // Montgomery factor cancels out
  }

  /* Copy for the full twisted matrix */
  for(i=0;i<KYBER_K-1;i++) { // row
    for(j=i+1;j<KYBER_K;j++) { // column
      memcpy(&a[i].vec[j], &at.vec[KYBER_K+i-j], 2*KYBER_N);
    }
  }
}
#elif KYBER_K == 4
void gen_matrix_ring(polyvec *a, const uint8_t seed[KYBER_SYMBYTES])
{
  unsigned int i, j;
  unsigned int ctr0, ctr1, ctr2, ctr3;
  ALIGNED_UINT8(REJ_UNIFORM_AVX_NBLOCKS*SHAKE128_RATE) buf[4];
  __m256i f;
  keccakx4_state state;
  polyvec at;

  f = _mm256_loadu_si256((__m256i *)seed);
  _mm256_store_si256(buf[0].vec, f);
  _mm256_store_si256(buf[1].vec, f);
  _mm256_store_si256(buf[2].vec, f);
  _mm256_store_si256(buf[3].vec, f);

  buf[0].coeffs[32] = 0;
  buf[0].coeffs[33] = 0;
  buf[1].coeffs[32] = 1;
  buf[1].coeffs[33] = 0;
  buf[2].coeffs[32] = 2;
  buf[2].coeffs[33] = 0;
  buf[3].coeffs[32] = 3;
  buf[3].coeffs[33] = 0;

  shake128x4_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 34);
  shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, REJ_UNIFORM_AVX_NBLOCKS, &state);

  ctr0 = rej_uniform_avx(a[0].vec[0].coeffs, buf[0].coeffs, KYBER_N);
  ctr1 = rej_uniform_avx(a[1].vec[0].coeffs, buf[1].coeffs, KYBER_N);
  ctr2 = rej_uniform_avx(a[2].vec[0].coeffs, buf[2].coeffs, KYBER_N);
  ctr3 = rej_uniform_avx(a[3].vec[0].coeffs, buf[3].coeffs, KYBER_N);

  while(ctr0 < KYBER_N || ctr1 < KYBER_N || ctr2 < KYBER_N || ctr3 < KYBER_N) {
    shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 1, &state);

    ctr0 += rej_uniform(a[0].vec[0].coeffs + ctr0, KYBER_N - ctr0, buf[0].coeffs, SHAKE128_RATE);
    ctr1 += rej_uniform(a[1].vec[0].coeffs + ctr1, KYBER_N - ctr1, buf[1].coeffs, SHAKE128_RATE);
    ctr2 += rej_uniform(a[2].vec[0].coeffs + ctr2, KYBER_N - ctr2, buf[2].coeffs, SHAKE128_RATE);
    ctr3 += rej_uniform(a[3].vec[0].coeffs + ctr3, KYBER_N - ctr3, buf[3].coeffs, SHAKE128_RATE);
  }

  poly_nttunpack(&a[0].vec[0]);
  poly_nttunpack(&a[1].vec[0]);
  poly_nttunpack(&a[2].vec[0]);
  poly_nttunpack(&a[3].vec[0]);

  /* Copy for the full non-twisted matrix */
  for(i=1;i<KYBER_K;i++) { // column
    for(j=i;j<KYBER_K;j++) { // row
      memcpy(&a[j].vec[i], &a[j-i].vec[0], 2*KYBER_N);
    }
  }

  // Can be pre-computed, but minor impact for AVX2
  /* Compute NTT(y+1) */
  poly twist;
  compute_twist(&twist);

  /* Generate the K-1 twists (multiplied by (y+1)) */
  for(i=1;i<KYBER_K;i++) {
    poly_basemul_montgomery(&at.vec[i], &a[i].vec[0], &twist); // Montgomery factor cancels out
  }

  /* Copy for the full twisted matrix */
  for(i=0;i<KYBER_K-1;i++) { // row
    for(j=i+1;j<KYBER_K;j++) { // column
      memcpy(&a[i].vec[j], &at.vec[KYBER_K+i-j], 2*KYBER_N);
    }
  }
}
#else
/* Updated matrix for ring setting */
#define GEN_MATRIX_NBLOCKS ((12*KYBER_N/8*(1 << 12)/KYBER_Q + XOF_BLOCKBYTES)/XOF_BLOCKBYTES)
/* Generate ring element a[0] of K coefficients and twisted elements by (y+1) in a[1] */
void gen_matrix_ring(polyvec *a, const uint8_t seed[KYBER_SYMBYTES])
{
  unsigned int ctr, i, j;
  unsigned int buflen;
  uint8_t buf[GEN_MATRIX_NBLOCKS*XOF_BLOCKBYTES];
  xof_state state;
  polyvec at;

  //for(i=0;i<KYBER_K;i++) {
  //  for(j=0;j<KYBER_K;j++) {
  //    for(unsigned int l = 0; l < KYBER_N; l++) a[i].vec[j].coeffs[l] = 0;
  //  }
  //}

  /* Generate K random polynomials */
  for(i=0;i<KYBER_K;i++) {
      xof_absorb(&state, seed, i, 0);
      xof_squeezeblocks(buf, GEN_MATRIX_NBLOCKS, &state);
      buflen = GEN_MATRIX_NBLOCKS*XOF_BLOCKBYTES;

      ctr = rej_uniform(a[i].vec[0].coeffs, KYBER_N, buf, buflen);

      while(ctr < KYBER_N) {
        xof_squeezeblocks(buf, 1, &state);
        buflen = XOF_BLOCKBYTES;
        ctr += rej_uniform(a[i].vec[0].coeffs + ctr, KYBER_N - ctr, buf, buflen);
      }
  }

  /* Copy for the full non-twisted matrix */
  for(i=1;i<KYBER_K;i++) { // column
    for(j=i;j<KYBER_K;j++) { // row
      memcpy(&a[j].vec[i], &a[j-i].vec[0], 2*KYBER_N);
    }
  }

  // Can be pre-computed, but minor impact for AVX2
  /* Compute NTT(y+1) */
  poly twist;
  compute_twist(&twist);

  /* Generate the K-1 twists (multiplied by (y+1)) */
  for(i=1;i<KYBER_K;i++) {
    poly_basemul_montgomery(&at.vec[i], &a[i].vec[0], &twist); // Montgomery factor cancels out
  }

  /* Copy for the full twisted matrix */
  for(i=0;i<KYBER_K-1;i++) { // row
    for(j=i+1;j<KYBER_K;j++) { // column
      memcpy(&a[i].vec[j], &at.vec[KYBER_K+i-j], 2*KYBER_N);
    }
  }
}
#endif

/* Ring multiplication r = a*b, returning the p most significant coefficients */
static void polyvec_partial_inner_prod(polyvec *r, polyvec *a, polyvec *b, int p)
{
  int i, j;
  poly t;

  // Can be pre-computed, but minor impact for AVX2
  /* Compute NTT(y+1) */
  poly twist;
  compute_twist(&twist);

  for(i = KYBER_K; i > KYBER_K-p; i--) { // Only final rows

    poly_basemul_montgomery(&r->vec[KYBER_K-i], &a->vec[KYBER_K-1], &b->vec[i % KYBER_K]);
    for(j = KYBER_K-1; j > 0; j--) {
      poly_basemul_montgomery(&t, &a->vec[(j-1) % KYBER_K], &b->vec[(i-j+KYBER_K) % KYBER_K]);
      poly_add(&r->vec[KYBER_K-i], &r->vec[KYBER_K-i], &t);
    }

    // Twist element of a for next iteration
    if (i != (KYBER_K-p+1)) {
      poly_basemul_montgomery(&a->vec[i-1], &a->vec[i-1], &twist);
    }
  }
}

/*************************************************
* Name:        indcpa_keypair_derand
*
* Description: Generates public and private key for the CPA-secure
*              public-key encryption scheme underlying Kyber
*
* Arguments:   - uint8_t *pk: pointer to output public key
*                             (of length KYBER_INDCPA_PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key
*                             (of length KYBER_INDCPA_SECRETKEYBYTES bytes)
*              - const uint8_t *coins: pointer to input randomness
*                             (of length KYBER_SYMBYTES bytes)
**************************************************/
void indcpa_keypair_derand(uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                           uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES],
                           const uint8_t coins[KYBER_SYMBYTES])
{
  unsigned int i;
  uint8_t buf[2*KYBER_SYMBYTES];
  const uint8_t *publicseed = buf;
  const uint8_t *noiseseed = buf + KYBER_SYMBYTES;
  polyvec a[KYBER_K], e, pkpv, skpv;

  memcpy(buf, coins, KYBER_SYMBYTES);
  buf[KYBER_SYMBYTES] = KYBER_K;
  hash_g(buf, buf, KYBER_SYMBYTES+1);

  gen_matrix_ring(a, publicseed);

#if KYBER_K == 2
  poly_getnoise_eta1_4x(skpv.vec+0, skpv.vec+1, e.vec+0, e.vec+1, noiseseed, 0, 1, 2, 3);
#elif KYBER_K == 3
  poly_getnoise_eta1_4x(skpv.vec+0, skpv.vec+1, skpv.vec+2, e.vec+0, noiseseed, 0, 1, 2, 3);
  poly_getnoise_eta1_4x(e.vec+1, e.vec+2, pkpv.vec+0, pkpv.vec+1, noiseseed, 4, 5, 6, 7);
#elif KYBER_K == 4
  poly_getnoise_eta1_4x(skpv.vec+0, skpv.vec+1, skpv.vec+2, skpv.vec+3, noiseseed,  0, 1, 2, 3);
  poly_getnoise_eta1_4x(e.vec+0, e.vec+1, e.vec+2, e.vec+3, noiseseed, 4, 5, 6, 7);
#endif

  polyvec_ntt(&skpv);
  polyvec_reduce(&skpv);
  polyvec_ntt(&e);

  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++) {
    polyvec_basemul_acc_montgomery(&pkpv.vec[i], &a[i], &skpv);
    poly_tomont(&pkpv.vec[i]);
  }

  polyvec_add(&pkpv, &pkpv, &e);
  polyvec_reduce(&pkpv);

  pack_sk(sk, &skpv);
  pack_pk(pk, &pkpv, publicseed);
}

/*************************************************
* Name:        indcpa_enc
*
* Description: Encryption function of the CPA-secure
*              public-key encryption scheme underlying Kyber.
*
* Arguments:   - uint8_t *c: pointer to output ciphertext
*                            (of length KYBER_INDCPA_BYTES bytes)
*              - const uint8_t *m: pointer to input message
*                                  (of length KYBER_INDCPA_MSGBYTES bytes)
*              - const uint8_t *pk: pointer to input public key
*                                   (of length KYBER_INDCPA_PUBLICKEYBYTES)
*              - const uint8_t *coins: pointer to input random coins used as seed
*                                      (of length KYBER_SYMBYTES) to deterministically
*                                      generate all randomness
**************************************************/
void indcpa_enc(uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                const uint8_t coins[KYBER_SYMBYTES])
{
  unsigned int i;
  uint8_t seed[KYBER_SYMBYTES];
  polyvec sp, pkpv, ep, at[KYBER_K], b;
  poly k, epp;
  polyvec v;

  unpack_pk(&pkpv, seed, pk);
  poly_frommsg(&k, m);
  gen_matrix_ring(at, seed);

#if KYBER_K == 2
  poly_getnoise_eta1122_4x(sp.vec+0, sp.vec+1, ep.vec+0, ep.vec+1, coins, 0, 1, 2, 3);
  poly_getnoise_eta2(&epp, coins, 4);
#elif KYBER_K == 3
  poly_getnoise_eta1_4x(sp.vec+0, sp.vec+1, sp.vec+2, ep.vec+0, coins, 0, 1, 2 ,3);
  poly_getnoise_eta1_4x(ep.vec+1, ep.vec+2, &epp, b.vec+0, coins,  4, 5, 6, 7);
#elif KYBER_K == 4
  poly_getnoise_eta1_4x(sp.vec+0, sp.vec+1, sp.vec+2, sp.vec+3, coins, 0, 1, 2, 3);
  poly_getnoise_eta1_4x(ep.vec+0, ep.vec+1, ep.vec+2, ep.vec+3, coins, 4, 5, 6, 7);
  poly_getnoise_eta2(&epp, coins, 8);
#endif

  polyvec_ntt(&sp);

  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++)
    polyvec_basemul_acc_montgomery(&b.vec[i], &at[i], &sp);

  polyvec_partial_inner_prod(&v, &pkpv, &sp, KYBER_SS_NUM_COEFFS);

  polyvec_invntt_tomont(&b);
  for(i=0;i<KYBER_SS_NUM_COEFFS;i++) {
    poly_invntt_tomont(&v.vec[i]);
  }

  polyvec_add(&b, &b, &ep);
  for(i=0;i<KYBER_SS_NUM_COEFFS;i++) {
    poly_add(&v.vec[i], &v.vec[i], &epp);
    poly_add(&v.vec[i], &v.vec[i], &k);
    poly_reduce(&v.vec[i]);
  }

  polyvec_reduce(&b);

  pack_ciphertext(c, &b, &v.vec[0]);
}

/*************************************************
* Name:        indcpa_dec
*
* Description: Decryption function of the CPA-secure
*              public-key encryption scheme underlying Kyber.
*
* Arguments:   - uint8_t *m: pointer to output decrypted message
*                            (of length KYBER_INDCPA_MSGBYTES)
*              - const uint8_t *c: pointer to input ciphertext
*                                  (of length KYBER_INDCPA_BYTES)
*              - const uint8_t *sk: pointer to input secret key
*                                   (of length KYBER_INDCPA_SECRETKEYBYTES)
**************************************************/
void indcpa_dec(uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES])
{
  unsigned int i;
  polyvec b, skpv, mp;
  polyvec v;

  unpack_ciphertext(&b, &v.vec[0], c);
  unpack_sk(&skpv, sk);

  polyvec_ntt(&b);

  //polyvec_basemul_acc_montgomery(&mp, &skpv, &b);
  //poly_invntt_tomont(&mp);
  polyvec_partial_inner_prod(&mp, &b, &skpv, KYBER_SS_NUM_COEFFS);

  for(i=0;i<KYBER_SS_NUM_COEFFS;i++) {
    poly_invntt_tomont(&mp.vec[i]);
  }

  //poly_sub(&mp, &v, &mp);
  //poly_reduce(&mp);
  for(i=0;i<KYBER_SS_NUM_COEFFS;i++) {
    poly_sub(&mp.vec[i], &v.vec[i], &mp.vec[i]);
    poly_reduce(&mp.vec[i]);
  }

  poly_tomsg(m, &mp.vec[0]);
}
