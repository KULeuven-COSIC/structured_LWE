#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "params.h"
#include "indcpa.h"
#include "polyvec.h"
#include "poly.h"
#include "ntt.h"
#include "symmetric.h"
#include "randombytes.h"

/* NTT(y+2) in Montgomery domain = NTT(2285y + 2285) */
#if KYBER_N == 64
poly precomputed_twist = { .coeffs = {
    -865, 18, -1436, 589, 689, -1536, -1073, 226,
    -52, -795, -597, -250, 959, 1523, -303, -544,
    -1572, 725, 1233, 1249, 921, 1561, 575, -1422,
    -377, -470, 79, -926, 1367, 1115, -619, -228,
    388, -1235, 1151, 1331, 970, 1512, -1258, 411,
    1348, 1134, -180, -667, 994, 1488, 290, -1137,
    843, 1639, -1127, 280, -267, -580, 516, -1363,
    -1640, 793, 176, -1023, -1411, 564, -34, -813
}};
#elif KYBER_N == 128
poly precomputed_twist = { .coeffs = {
    138, -985, -1658, 811, -1533, 686, -1245, 398,
    -10, -837, -1217, 370, -538, -309, 1346, 1136,
    1663, 819, -1501, 654, 1418, 1064, 1006, 1476,
    950, 1532, 781, -1628, -514, -333, -435, -412,
    995, 1487, -1310, 463, -929, 82, 1094, 1388,
    464, -1311, -605, -242, 639, -1486, -969, 122,
    -349, -498, -1444, 597, 369, -1216, 1590, 892,
    1659, 823, 1570, 912, 1085, 1397, 1166, 1316,
    -1271, 424, -991, 144, -1485, 638, -1478, 631,
    -766, -81, -44, -803, -224, -623, 1625, 857,
    26, -873, 1105, 1377, -870, 23, -94, -753,
    367, -1214, 1461, 1021, 54, -901, -418, -429,
    56, -903, -289, -558, -37, -810, -1294, 447,
    -269, -578, 387, -1234, 371, -1218, -1610, 763,
    1133, 1349, 933, 1549, -1092, 245, -1097, 250,
    -1130, 283, -219, -628, -566, -281, -460, -387,
}};
#elif KYBER_N == 256
poly precomputed_twist = { .coeffs = { 
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044,
    1241, -1044, 1241, -1044, 1241, -1044, 1241, -1044
}};
#else
  #error "KYBER_N not supported"
#endif

/*************************************************
* Name:        pack_pk
*
* Description: Serialize the public key as concatenation of the
*              serialized vector of polynomials pk
*              and the public seed used to generate the matrix A.
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
* Description: Serialize the secret key
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
*              and the compressed and serialized polynomial v
*
* Arguments:   uint8_t *r: pointer to the output serialized ciphertext
*              poly *pk: pointer to the input vector of polynomials b
*              poly *v: pointer to the input polynomial v
**************************************************/
static void pack_ciphertext(uint8_t r[KYBER_INDCPA_BYTES], polyvec *b, polyvec *v)
{
  polyvec_compress(r, b);
  for(unsigned int i=0;i<KYBER_SS_NUM_COEFFS;i++) {
    poly_compress(r+KYBER_POLYVECCOMPRESSEDBYTES+i*KYBER_POLYCOMPRESSEDBYTES, &v->vec[i]);
  }
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
static void unpack_ciphertext(polyvec *b, polyvec *v, const uint8_t c[KYBER_INDCPA_BYTES])
{
  polyvec_decompress(b, c);
  for(unsigned int i=0;i<KYBER_SS_NUM_COEFFS;i++) {
    poly_decompress(&v->vec[i], c+KYBER_POLYVECCOMPRESSEDBYTES+i*KYBER_POLYCOMPRESSEDBYTES);
  }
}

/*************************************************
* Name:        rej_uniform
*
* Description: Run rejection sampling on uniform random bytes to generate
*              uniform random integers mod q
*
* Arguments:   - int16_t *r: pointer to output buffer
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
  while(ctr < len && pos + 3 <= buflen) {
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
#if(XOF_BLOCKBYTES % 3)
#error "Implementation of gen_matrix assumes that XOF_BLOCKBYTES is a multiple of 3"
#endif

#define GEN_MATRIX_NBLOCKS ((12*KYBER_N/8*(1 << 12)/KYBER_Q + XOF_BLOCKBYTES)/XOF_BLOCKBYTES)
// Not static for benchmarking
void gen_matrix(polyvec *a, const uint8_t seed[KYBER_SYMBYTES], int transposed)
{
  unsigned int ctr, i, j;
  unsigned int buflen;
  uint8_t buf[GEN_MATRIX_NBLOCKS*XOF_BLOCKBYTES];
  xof_state state;

  for(i=0;i<KYBER_K;i++) {
    for(j=0;j<KYBER_K;j++) {
      if(transposed)
        xof_absorb(&state, seed, i, j);
      else
        xof_absorb(&state, seed, j, i);

      xof_squeezeblocks(buf, GEN_MATRIX_NBLOCKS, &state);
      buflen = GEN_MATRIX_NBLOCKS*XOF_BLOCKBYTES;
      ctr = rej_uniform(a[i].vec[j].coeffs, KYBER_N, buf, buflen);

      while(ctr < KYBER_N) {
        xof_squeezeblocks(buf, 1, &state);
        buflen = XOF_BLOCKBYTES;
        ctr += rej_uniform(a[i].vec[j].coeffs + ctr, KYBER_N - ctr, buf, buflen);
      }
    }
  }
}

static void compute_twist(poly *r)
{
  #ifdef PRECOMPUTE_TWIST
    *r = precomputed_twist;
  #else
    // Initialize the polynomial y+2 in Montgomery domain
    r->coeffs[0] = 1241; // 2
    r->coeffs[1] = 2285; // 1
    for(unsigned int i = 2; i < KYBER_N; i++) r->coeffs[i] = 0;
    ntt(&r->coeffs[0]);
    poly_reduce(r);
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

  for(i=0;i<KYBER_K;i++) {
    for(j=0;j<KYBER_K;j++) {
      for(unsigned int l = 0; l < KYBER_N; l++) a[i].vec[j].coeffs[l] = 0;
    }
  }

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

/* Ring multiplication r = a*b, returning the p most significant coefficients */
static void polyvec_partial_inner_prod(polyvec *r, polyvec *a, polyvec *b, int p)
{
  int i, j;
  poly t, twist;

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
  const uint8_t *noiseseed = buf+KYBER_SYMBYTES;
  uint8_t nonce = 0;
  polyvec a[KYBER_K], e, pkpv, skpv;

  memcpy(buf, coins, KYBER_SYMBYTES);
  buf[KYBER_SYMBYTES] = KYBER_K;
  hash_g(buf, buf, KYBER_SYMBYTES+1);

  gen_matrix_ring(a, publicseed);

  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta1(&skpv.vec[i], noiseseed, nonce++);
  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta1(&e.vec[i], noiseseed, nonce++);

  polyvec_ntt(&skpv);
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
  uint8_t nonce = 0;
  polyvec sp, pkpv, ep, epp, at[KYBER_K], b;
  poly k[KYBER_SS_NUM_COEFFS];
  polyvec v;

  unpack_pk(&pkpv, seed, pk);
  poly_frommsg(k, m);
  gen_matrix_ring(at, seed);

  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta1(sp.vec+i, coins, nonce++);
  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta2(ep.vec+i, coins, nonce++);
  for(i=0;i<KYBER_SS_NUM_COEFFS;i++)
    poly_getnoise_eta2(epp.vec+i, coins, nonce++);

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
    poly_add(&v.vec[i], &v.vec[i], &epp.vec[i]);
    poly_add(&v.vec[i], &v.vec[i], &k[i]);
    poly_reduce(&v.vec[i]);
  }
  polyvec_reduce(&b);

  pack_ciphertext(c, &b, &v);
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
  polyvec b, skpv;
  polyvec v; 
  polyvec mp;

  unpack_ciphertext(&b, &v, c);
  unpack_sk(&skpv, sk);

  polyvec_ntt(&b);
  polyvec_partial_inner_prod(&mp, &b, &skpv, KYBER_SS_NUM_COEFFS);

  for(i=0;i<KYBER_SS_NUM_COEFFS;i++) {
    poly_invntt_tomont(&mp.vec[i]);
  }

  for(i=0;i<KYBER_SS_NUM_COEFFS;i++) {
    poly_sub(&mp.vec[i], &v.vec[i], &mp.vec[i]);
    poly_reduce(&mp.vec[i]);
  }

  poly_tomsg(m, &mp);
}
