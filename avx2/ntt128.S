#include "consts.h"
.include "shuffle.inc"

.macro mul0 rh0,rh1,zl0=15,zl1=15,zh0=2,zh1=2
vpmullw		%ymm\zl0,%ymm\rh0,%ymm12
vpmullw		%ymm\zl0,%ymm\rh1,%ymm13
vpmulhw		%ymm\zh0,%ymm\rh0,%ymm\rh0
vpmulhw		%ymm\zh0,%ymm\rh1,%ymm\rh1
.endm

.macro mul1 rh2,rh3,zl0=15,zl1=15,zh0=2,zh1=2
vpmullw		%ymm\zl1,%ymm\rh2,%ymm14
vpmullw		%ymm\zl1,%ymm\rh3,%ymm15
vpmulhw		%ymm\zh1,%ymm\rh2,%ymm\rh2
vpmulhw		%ymm\zh1,%ymm\rh3,%ymm\rh3
.endm

.macro mul2 rh0,rh1,zl0=15,zl1=15,zh0=2,zh1=2
vpmullw		%ymm\zl0,%ymm\rh0,%ymm12
vpmullw		%ymm\zl1,%ymm\rh1,%ymm13
vpmulhw		%ymm\zh0,%ymm\rh0,%ymm\rh0
vpmulhw		%ymm\zh1,%ymm\rh1,%ymm\rh1
.endm

.macro mul3 rh2,rh3,zl0=15,zl1=15,zh0=2,zh1=2
vpmullw		%ymm\zl0,%ymm\rh2,%ymm14
vpmullw		%ymm\zl1,%ymm\rh3,%ymm15
vpmulhw		%ymm\zh0,%ymm\rh2,%ymm\rh2
vpmulhw		%ymm\zh1,%ymm\rh3,%ymm\rh3
.endm

.macro mul rh0,rh1,rh2,rh3,zl0=15,zl1=15,zh0=2,zh1=2
vpmullw		%ymm\zl0,%ymm\rh0,%ymm12
vpmullw		%ymm\zl0,%ymm\rh1,%ymm13

vpmullw		%ymm\zl1,%ymm\rh2,%ymm14
vpmullw		%ymm\zl1,%ymm\rh3,%ymm15

vpmulhw		%ymm\zh0,%ymm\rh0,%ymm\rh0
vpmulhw		%ymm\zh0,%ymm\rh1,%ymm\rh1

vpmulhw		%ymm\zh1,%ymm\rh2,%ymm\rh2
vpmulhw		%ymm\zh1,%ymm\rh3,%ymm\rh3
.endm

.macro reduce
vpmulhw		%ymm0,%ymm12,%ymm12
vpmulhw		%ymm0,%ymm13,%ymm13

vpmulhw		%ymm0,%ymm14,%ymm14
vpmulhw		%ymm0,%ymm15,%ymm15
.endm

.macro update rln,rl0,rl1,rl2,rl3,rh0,rh1,rh2,rh3
vpaddw		%ymm\rh0,%ymm\rl0,%ymm\rln
vpsubw		%ymm\rh0,%ymm\rl0,%ymm\rh0
vpaddw		%ymm\rh1,%ymm\rl1,%ymm\rl0

vpsubw		%ymm\rh1,%ymm\rl1,%ymm\rh1
vpaddw		%ymm\rh2,%ymm\rl2,%ymm\rl1
vpsubw		%ymm\rh2,%ymm\rl2,%ymm\rh2

vpaddw		%ymm\rh3,%ymm\rl3,%ymm\rl2
vpsubw		%ymm\rh3,%ymm\rl3,%ymm\rh3

vpsubw		%ymm12,%ymm\rln,%ymm\rln
vpaddw		%ymm12,%ymm\rh0,%ymm\rh0
vpsubw		%ymm13,%ymm\rl0,%ymm\rl0

vpaddw		%ymm13,%ymm\rh1,%ymm\rh1
vpsubw		%ymm14,%ymm\rl1,%ymm\rl1
vpaddw		%ymm14,%ymm\rh2,%ymm\rh2

vpsubw		%ymm15,%ymm\rl2,%ymm\rl2
vpaddw		%ymm15,%ymm\rh3,%ymm\rh3
.endm

.macro level0 off
vpbroadcastq	(_ZETAS_EXP+0)*2(%rsi),%ymm15
vmovdqa		(64)*2(%rdi),%ymm8
vmovdqa		(80)*2(%rdi),%ymm9
vmovdqa		(96)*2(%rdi),%ymm10
vmovdqa		(112)*2(%rdi),%ymm11
vpbroadcastq	(_ZETAS_EXP+4)*2(%rsi),%ymm2

mul		8,9,10,11

vmovdqa		( 0)*2(%rdi),%ymm4
vmovdqa		(16)*2(%rdi),%ymm5
vmovdqa		(32)*2(%rdi),%ymm6
vmovdqa		(48)*2(%rdi),%ymm7

reduce
update		3,4,5,6,7,8,9,10,11

.endm

.macro levels1t6 off
/* level 1 */

vmovdqa		(_ZETAS_EXP+224*0+16)*2(%rsi),%ymm15
vmovdqa		(_ZETAS_EXP+224*0+32)*2(%rsi),%ymm2

mul0		5,6

vmovdqa		(_ZETAS_EXP+224*1+16)*2(%rsi),%ymm15
vmovdqa		(_ZETAS_EXP+224*1+32)*2(%rsi),%ymm2

mul1        10,11

reduce
update		7,3,4,8,9,5,6,10,11

/* level 2 */
// [7, 3, 5, 6, 4, 8, 10, 11]

shuffle8	7,5,9,5
shuffle8	3,6,7,6

// [9, 7, 5, 6, 4, 8, 10, 11]

vmovdqa		(_ZETAS_EXP+224*0+48)*2(%rsi),%ymm15
vmovdqa		(_ZETAS_EXP+224*0+64)*2(%rsi),%ymm2

mul0		7,6

shuffle8	4,10,3,10
shuffle8	8,11,4,11

// [9, 7, 5, 6, 3, 4, 10, 11]

vmovdqa		(_ZETAS_EXP+224*1+48)*2(%rsi),%ymm15
vmovdqa		(_ZETAS_EXP+224*1+64)*2(%rsi),%ymm2

mul1		4,11

reduce
update		8,9,5,3,10,7,6,4,11

//vmovdqa		%ymm8, (  32)*2(%rdi)
//vmovdqa		%ymm9, (  16)*2(%rdi)
//vmovdqa		%ymm4, (  32)*2(%rdi)
//vmovdqa		%ymm11, (  48)*2(%rdi)

//vmovdqa		%ymm8, (  0)*2(%rdi)
//vmovdqa		%ymm9, ( 16)*2(%rdi)
//vmovdqa		%ymm7, ( 32)*2(%rdi)
//vmovdqa		%ymm6, ( 48)*2(%rdi)
//vmovdqa		%ymm5, ( 64)*2(%rdi)
//vmovdqa		%ymm3, ( 80)*2(%rdi)
//vmovdqa		%ymm4,( 96)*2(%rdi)
//vmovdqa		%ymm11,(112)*2(%rdi)

/* level 3 */
// [8, 9, 7, 6, 5, 3, 4, 11]

//vmovdqa		%ymm8, (  0)*2(%rdi)
//vmovdqa		%ymm9, ( 16)*2(%rdi)
//vmovdqa		%ymm7, ( 32)*2(%rdi)
//vmovdqa		%ymm6, ( 48)*2(%rdi)

shuffle4	8,7,10,7
shuffle4	9,6,8,6

// [10, 8, 7, 6, 5, 3, 4, 11]

//vmovdqa		%ymm10, (  0)*2(%rdi)
//vmovdqa		%ymm9, ( 16)*2(%rdi)
//vmovdqa		%ymm8, ( 32)*2(%rdi)
//vmovdqa		%ymm6, ( 48)*2(%rdi)

vmovdqa		(_ZETAS_EXP+224*0+80)*2(%rsi),%ymm15
vmovdqa		(_ZETAS_EXP+224*0+96)*2(%rsi),%ymm2

mul0        8,6

shuffle4	5,4,9,4
shuffle4	3,11,5,11

// [10, 8, 7, 6, 9, 5, 4, 11]

vmovdqa		(_ZETAS_EXP+224*1+80)*2(%rsi),%ymm15
vmovdqa		(_ZETAS_EXP+224*1+96)*2(%rsi),%ymm2

mul1        5,11

reduce
update		3,10,7,9,4,8,6,5,11

//vmovdqa		%ymm3, (  0)*2(%rdi)
//vmovdqa		%ymm10, ( 16)*2(%rdi)
//vmovdqa		%ymm8, ( 32)*2(%rdi)
//vmovdqa		%ymm6, ( 48)*2(%rdi)
//vmovdqa		%ymm7, ( 64)*2(%rdi)
//vmovdqa		%ymm9, ( 80)*2(%rdi)
//vmovdqa		%ymm5,( 96)*2(%rdi)
//vmovdqa		%ymm11,(112)*2(%rdi)

///* level 4 */
// [3, 10, 8, 6, 7, 9, 5, 11]

shuffle2	3,8,4,8
shuffle2	10,6,3,6

// [4, 3, 8, 6, 7, 9, 5, 11]

vmovdqa		(_ZETAS_EXP+224*0+112)*2(%rsi),%ymm15
vmovdqa		(_ZETAS_EXP+224*0+128)*2(%rsi),%ymm2

mul0    3,6

shuffle2	7,5,10,5
shuffle2	9,11,7,11

vmovdqa		(_ZETAS_EXP+224*1+112)*2(%rsi),%ymm15
vmovdqa		(_ZETAS_EXP+224*1+128)*2(%rsi),%ymm2

// [4, 3, 8, 6, 10, 7, 5, 11]

mul1    7,11

reduce
update		9,4,8,10,5,3,6,7,11

vmovdqa		%ymm9, (  0)*2(%rdi)
vmovdqa		%ymm4, ( 16)*2(%rdi)
vmovdqa		%ymm3, ( 32)*2(%rdi)
vmovdqa		%ymm6, ( 48)*2(%rdi)
vmovdqa		%ymm8, ( 64)*2(%rdi)
vmovdqa		%ymm10, ( 80)*2(%rdi)
vmovdqa		%ymm7,( 96)*2(%rdi)
vmovdqa		%ymm11,(112)*2(%rdi)

/* level 5 */
// [9, 4, 3, 6, 8, 10, 7, 11]

shuffle1	9,3,5,3
shuffle1	4,6,9,6

// [5, 9, 3, 6, 8, 10, 7, 11]

vmovdqa		(_ZETAS_EXP+224*0+144)*2(%rsi),%ymm15
vmovdqa		(_ZETAS_EXP+224*0+160)*2(%rsi),%ymm2

mul0        9,6

shuffle1	8,7,4,7
shuffle1	10,11,8,11

// [5, 9, 3, 6, 4, 8, 7, 11]

vmovdqa		(_ZETAS_EXP+224*1+144)*2(%rsi),%ymm15
vmovdqa		(_ZETAS_EXP+224*1+160)*2(%rsi),%ymm2

mul1        8,11

reduce
update		10,5,3,4,7,9,6,8,11

vmovdqa		%ymm10, (  0)*2(%rdi)
vmovdqa		%ymm5, ( 16)*2(%rdi)
vmovdqa		%ymm9, ( 32)*2(%rdi)
vmovdqa		%ymm6, ( 48)*2(%rdi)
vmovdqa		%ymm3, ( 64)*2(%rdi)
vmovdqa		%ymm4, ( 80)*2(%rdi)
vmovdqa		%ymm8,( 96)*2(%rdi)
vmovdqa		%ymm11,(112)*2(%rdi)

/* level 6 */
// [10, 5, 9, 6, 3, 4, 8, 11]

vmovdqa		(_ZETAS_EXP+224*0+176)*2(%rsi),%ymm14
vmovdqa		(_ZETAS_EXP+224*0+208)*2(%rsi),%ymm15
vmovdqa		(_ZETAS_EXP+224*0+192)*2(%rsi),%ymm8
vmovdqa		(_ZETAS_EXP+224*0+224)*2(%rsi),%ymm2

vmovdqa		%ymm14, (  0)*2(%rdi)
vmovdqa		%ymm15, ( 16)*2(%rdi)
vmovdqa		%ymm8, ( 32)*2(%rdi)
vmovdqa		%ymm2, ( 48)*2(%rdi)

mul2        5,6,14,15,8,2

vmovdqa		(_ZETAS_EXP+224*1+176)*2(%rsi),%ymm14
vmovdqa		(_ZETAS_EXP+224*1+208)*2(%rsi),%ymm15
vmovdqa		(_ZETAS_EXP+224*1+192)*2(%rsi),%ymm8
vmovdqa		(_ZETAS_EXP+224*1+224)*2(%rsi),%ymm2

vmovdqa		%ymm14, (  0)*2(%rdi)
vmovdqa		%ymm15, ( 16)*2(%rdi)
vmovdqa		%ymm8, ( 32)*2(%rdi)
vmovdqa		%ymm2, ( 48)*2(%rdi)

mul3        4,11,14,15,8,2

// [10, 5, 9, 6, 3, 4, 8, 11]

reduce
update		7,10,9,3,8,5,6,4,11

vmovdqa		%ymm7, (  0)*2(%rdi)
vmovdqa		%ymm10, ( 16)*2(%rdi)
vmovdqa		%ymm5, ( 32)*2(%rdi)
vmovdqa		%ymm6, ( 48)*2(%rdi)
vmovdqa		%ymm9, ( 64)*2(%rdi)
vmovdqa		%ymm3, ( 80)*2(%rdi)
vmovdqa		%ymm4,( 96)*2(%rdi)
vmovdqa		%ymm11,(112)*2(%rdi)

// [7, 10, 5, 6, 9, 3, 4, 11], free: 8

/* Unshuffle, only for testing */
//shuffle1 7,5,8,5
//shuffle1 10,6,7,6
//shuffle1 9,4,10,4
//shuffle1 3,11,9,11
//// [8, 7, 5, 6, 10, 9, 4, 11], free: 8
//vmovdqa		%ymm8, (  0)*2(%rdi)
//vmovdqa		%ymm7, ( 16)*2(%rdi)
//vmovdqa		%ymm5, ( 32)*2(%rdi)
//vmovdqa		%ymm6, ( 48)*2(%rdi)
//vmovdqa		%ymm10, ( 64)*2(%rdi)
//vmovdqa		%ymm9, ( 80)*2(%rdi)
//vmovdqa		%ymm4,( 96)*2(%rdi)
//vmovdqa		%ymm11,(112)*2(%rdi)

//shuffle2 8,7,3,7
//shuffle2 5,6,8,6
//shuffle2 10,9,5,9
//shuffle2 4,11,10,11
//// [3, 7, 8, 6, 5, 9, 10, 11], free: 4
//vmovdqa		%ymm3, (  0)*2(%rdi)
//vmovdqa		%ymm7, ( 16)*2(%rdi)
//vmovdqa		%ymm8, ( 32)*2(%rdi)
//vmovdqa		%ymm6, ( 48)*2(%rdi)
//vmovdqa		%ymm5, ( 64)*2(%rdi)
//vmovdqa		%ymm9, ( 80)*2(%rdi)
//vmovdqa		%ymm10,( 96)*2(%rdi)
//vmovdqa		%ymm11,(112)*2(%rdi)
//
//shuffle4 3,8,4,8
//shuffle4 7,6,3,6
//shuffle4 5,10,7,10
//shuffle4 9,11,5,11
//// [4, 3, 8, 6, 7, 5, 10, 11], free: 9
//vmovdqa		%ymm4, (  0)*2(%rdi)
//vmovdqa		%ymm3, ( 16)*2(%rdi)
//vmovdqa		%ymm8, ( 32)*2(%rdi)
//vmovdqa		%ymm6, ( 48)*2(%rdi)
//vmovdqa		%ymm7, ( 64)*2(%rdi)
//vmovdqa		%ymm5, ( 80)*2(%rdi)
//vmovdqa		%ymm10,( 96)*2(%rdi)
//vmovdqa		%ymm11,(112)*2(%rdi)
//
//shuffle8 4,3,9,3
//shuffle8 8,6,4,6
//shuffle8 7,5,8,5
//shuffle8 10,11,7,11
//// [9, 3, 4, 6, 8, 5, 7, 11], free: 10
//vmovdqa		%ymm9, (  0)*2(%rdi)
//vmovdqa		%ymm3, ( 16)*2(%rdi)
//vmovdqa		%ymm4, ( 32)*2(%rdi)
//vmovdqa		%ymm6, ( 48)*2(%rdi)
//vmovdqa		%ymm8, ( 64)*2(%rdi)
//vmovdqa		%ymm5, ( 80)*2(%rdi)
//vmovdqa		%ymm7,( 96)*2(%rdi)
//vmovdqa		%ymm11,(112)*2(%rdi)




//shuffle2 8,5,3,5
//shuffle2 7,6,8,6
//shuffle2 10,4,7,4
//shuffle2 9,11,10,11
//// [3, 8, 5, 6, 7, 10, 4, 11], free: 9
//vmovdqa		%ymm3, (  0)*2(%rdi)
//vmovdqa		%ymm8, ( 16)*2(%rdi)
//vmovdqa		%ymm5, ( 32)*2(%rdi)
//vmovdqa		%ymm6, ( 48)*2(%rdi)
//vmovdqa		%ymm7, ( 64)*2(%rdi)
//vmovdqa		%ymm10, ( 80)*2(%rdi)
//vmovdqa		%ymm4,( 96)*2(%rdi)
//vmovdqa		%ymm11,(112)*2(%rdi)
//shuffle4 3,5,9,5
//shuffle4 8,6,3,6
//shuffle4 7,4,8,4
//shuffle4 10,11,7,11
//// [9, 3, 5, 6, 8, 7, 4, 11], free: 10
//vmovdqa		%ymm9, (  0)*2(%rdi)
//vmovdqa		%ymm3, ( 16)*2(%rdi)
//vmovdqa		%ymm5, ( 32)*2(%rdi)
//vmovdqa		%ymm6, ( 48)*2(%rdi)
//vmovdqa		%ymm8, ( 64)*2(%rdi)
//vmovdqa		%ymm7, ( 80)*2(%rdi)
//vmovdqa		%ymm4,( 96)*2(%rdi)
//vmovdqa		%ymm11,(112)*2(%rdi)
//shuffle8 9,5,10,5
//shuffle8 3,6,9,6
//shuffle8 8,4,3,4
//shuffle8 7,11,8,11
//// [10, 9, 5, 6, 3, 8, 4, 11], free: 7
//vmovdqa		%ymm10, (  0)*2(%rdi)
//vmovdqa		%ymm9, ( 16)*2(%rdi)
//vmovdqa		%ymm5, ( 32)*2(%rdi)
//vmovdqa		%ymm6, ( 48)*2(%rdi)
//vmovdqa		%ymm3, ( 64)*2(%rdi)
//vmovdqa		%ymm8, ( 80)*2(%rdi)
//vmovdqa		%ymm4,( 96)*2(%rdi)
//vmovdqa		%ymm11,(112)*2(%rdi)

.endm

.text
.global cdecl(ntt_avx)
cdecl(ntt_avx):
vmovdqa		_16XQ*2(%rsi),%ymm0

level0		0

levels1t6	0

ret
