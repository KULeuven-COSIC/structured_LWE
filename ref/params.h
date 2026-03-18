#ifndef PARAMS_H
#define PARAMS_H

#ifndef BIT_SECURITY_TARGET
#define BIT_SECURITY_TARGET 128
#endif

/* Fixed across all parameter sets */
#define KYBER_Q 3329
#define KYBER_SYMBYTES 32 /* size in bytes of hashes, and seeds */
#define KYBER_ETA2 2
#define PRECOMPUTE_TWIST

#if   (BIT_SECURITY_TARGET == 120)
#define KYBER_N 256
#define KYBER_POLYBYTES	(3 * KYBER_N / 2)
#define KYBER_SS_NUM_COEFFS 1 /* size in number of coefficients of shared secret polynomial */
#define KYBER_K 2
#define KYBER_ETA1 3
#define KYBER_POLYVECCOMPRESSEDBYTES (10 * KYBER_K * KYBER_N / 8)
#define KYBER_POLYCOMPRESSEDBYTES    (4 * KYBER_N / 8)
#define KYBER_NAMESPACE(s) pqcrystals_kyber120_ref_##s
#elif   (BIT_SECURITY_TARGET == 128)
#define KYBER_N 64
#define KYBER_POLYBYTES	(3 * KYBER_N / 2)
#define KYBER_SS_NUM_COEFFS 2 /* size in number of coefficients of shared secret polynomial */
#define KYBER_K 9
#define KYBER_ETA1 2
#define KYBER_POLYVECCOMPRESSEDBYTES (10 * KYBER_K * KYBER_N / 8)
#define KYBER_POLYCOMPRESSEDBYTES    (4 * KYBER_N / 8)
#define KYBER_NAMESPACE(s) pqcrystals_kyber128_ref_##s
#elif   (BIT_SECURITY_TARGET == 145)
#define KYBER_N 128
#define KYBER_POLYBYTES	(3 * KYBER_N / 2)
#define KYBER_SS_NUM_COEFFS 1 /* size in number of coefficients of shared secret polynomial */
#define KYBER_K 5
#define KYBER_ETA1 2
#define KYBER_POLYVECCOMPRESSEDBYTES (10 * KYBER_K * KYBER_N / 8)
#define KYBER_POLYCOMPRESSEDBYTES    (4 * KYBER_N / 8)
#define KYBER_NAMESPACE(s) pqcrystals_kyber145_ref_##s
#elif   (BIT_SECURITY_TARGET == 163)
#define KYBER_N 64
#define KYBER_POLYBYTES	(3 * KYBER_N / 2)
#define KYBER_SS_NUM_COEFFS 3 /* size in number of coefficients of shared secret polynomial */
#define KYBER_K 11
#define KYBER_ETA1 2
#define KYBER_POLYVECCOMPRESSEDBYTES (10 * KYBER_K * KYBER_N / 8)
#define KYBER_POLYCOMPRESSEDBYTES    (4 * KYBER_N / 8)
#define KYBER_NAMESPACE(s) pqcrystals_kyber163_ref_##s
#elif (BIT_SECURITY_TARGET == 192)
#define KYBER_N 256
#define KYBER_POLYBYTES	(3 * KYBER_N / 2)
#define KYBER_SS_NUM_COEFFS 1 /* size in number of coefficients of shared secret polynomial */
#define KYBER_K 3
#define KYBER_ETA1 2
#define KYBER_POLYVECCOMPRESSEDBYTES (10 * KYBER_K * KYBER_N / 8)
#define KYBER_POLYCOMPRESSEDBYTES    (4 * KYBER_N / 8)
#define KYBER_NAMESPACE(s) pqcrystals_kyber192_ref_##s
#elif   (BIT_SECURITY_TARGET == 199)
#define KYBER_N 64
#define KYBER_POLYBYTES	(3 * KYBER_N / 2)
#define KYBER_SS_NUM_COEFFS 4 /* size in number of coefficients of shared secret polynomial */
#define KYBER_K 13
#define KYBER_ETA1 2
#define KYBER_POLYVECCOMPRESSEDBYTES (10 * KYBER_K * KYBER_N / 8)
#define KYBER_POLYCOMPRESSEDBYTES    (4 * KYBER_N / 8)
#define KYBER_NAMESPACE(s) pqcrystals_kyber199_ref_##s
#elif (BIT_SECURITY_TARGET == 224)
#define KYBER_N 128
#define KYBER_POLYBYTES	(3 * KYBER_N / 2)
#define KYBER_SS_NUM_COEFFS 2 /* size in number of coefficients of shared secret polynomial */
#define KYBER_K 7
#define KYBER_ETA1 2
#define KYBER_POLYVECCOMPRESSEDBYTES (11 * KYBER_K * KYBER_N / 8)
#define KYBER_POLYCOMPRESSEDBYTES    (5 * KYBER_N / 8)
#define KYBER_NAMESPACE(s) pqcrystals_kyber224_ref_##s
#elif   (BIT_SECURITY_TARGET == 235)
#define KYBER_N 64
#define KYBER_POLYBYTES	(3 * KYBER_N / 2)
#define KYBER_SS_NUM_COEFFS 4 /* size in number of coefficients of shared secret polynomial */
#define KYBER_K 15
#define KYBER_ETA1 2
#define KYBER_POLYVECCOMPRESSEDBYTES (11 * KYBER_K * KYBER_N / 8)
#define KYBER_POLYCOMPRESSEDBYTES    (5 * KYBER_N / 8)
#define KYBER_NAMESPACE(s) pqcrystals_kyber235_ref_##s
#elif (BIT_SECURITY_TARGET == 256)
#define KYBER_N 256
#define KYBER_POLYBYTES	(3 * KYBER_N / 2)
#define KYBER_SS_NUM_COEFFS 1 /* size in number of coefficients of shared secret polynomial */
#define KYBER_K 4
#define KYBER_ETA1 2
#define KYBER_POLYVECCOMPRESSEDBYTES (11 * KYBER_K * KYBER_N / 8)
#define KYBER_POLYCOMPRESSEDBYTES    (5 * KYBER_N / 8)
#define KYBER_NAMESPACE(s) pqcrystals_kyber256_ref_##s
#elif (BIT_SECURITY_TARGET == 512)
#define KYBER_N 256
#define KYBER_POLYBYTES	(3 * KYBER_N / 2)
#define KYBER_SS_NUM_COEFFS 2 /* size in number of coefficients of shared secret polynomial */
#define KYBER_K 8
#define KYBER_ETA1 2
#define KYBER_POLYVECCOMPRESSEDBYTES (11 * KYBER_K * KYBER_N / 8)
#define KYBER_POLYCOMPRESSEDBYTES    (5 * KYBER_N / 8)
#define KYBER_NAMESPACE(s) pqcrystals_kyber512_ref_##s
#endif

#define KYBER_SSBYTES                (KYBER_SS_NUM_COEFFS * KYBER_N/8) /* size in bytes of shared key */
#define KYBER_POLYVECBYTES	         (KYBER_K * KYBER_POLYBYTES)

#define KYBER_INDCPA_MSGBYTES       (KYBER_SSBYTES)
#define KYBER_INDCPA_PUBLICKEYBYTES (KYBER_POLYVECBYTES + KYBER_SYMBYTES)
#define KYBER_INDCPA_SECRETKEYBYTES (KYBER_POLYVECBYTES)
#define KYBER_INDCPA_BYTES          (KYBER_POLYVECCOMPRESSEDBYTES + KYBER_SS_NUM_COEFFS * KYBER_POLYCOMPRESSEDBYTES)

#define KYBER_PUBLICKEYBYTES  (KYBER_INDCPA_PUBLICKEYBYTES)
/* 32 bytes of additional space to save H(pk) */
#define KYBER_SECRETKEYBYTES  (KYBER_INDCPA_SECRETKEYBYTES + KYBER_INDCPA_PUBLICKEYBYTES + 2*KYBER_SYMBYTES)
#define KYBER_CIPHERTEXTBYTES (KYBER_INDCPA_BYTES)

#endif
