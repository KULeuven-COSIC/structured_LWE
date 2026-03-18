CC = gcc
CFLAGS += -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wpointer-arith -O3 -fomit-frame-pointer
NISTFLAGS += -Wno-unused-result -O3 -fomit-frame-pointer
RM = /bin/rm

SOURCES = kem.c indcpa.c polyvec.c poly.c ntt.c cbd.c reduce.c verify.c
SOURCESKECCAK = $(SOURCES) fips202.c symmetric-shake.c
HEADERS = params.h kem.h indcpa.h polyvec.h poly.h ntt.h cbd.h reduce.c verify.h symmetric.h
HEADERSKECCAK = $(HEADERS) fips202.h

.PHONY: all speed clean

all: test speed

test: \
  test/test_kyber120 \
  test/test_kyber128 \
  test/test_kyber145 \
  test/test_kyber163 \
  test/test_kyber192 \
  test/test_kyber199 \
  test/test_kyber224 \
  test/test_kyber235 \
  test/test_kyber256 \
  #test/test_kyber512 \

speed: \
  test/test_speed120 \
  test/test_speed128 \
  test/test_speed145 \
  test/test_speed163 \
  test/test_speed192 \
  test/test_speed199 \
  test/test_speed224 \
  test/test_speed235 \
  test/test_speed256 \
  #test/test_speed512 \

test/test_kyber120: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DBIT_SECURITY_TARGET=120 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -o $@

test/test_kyber128: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DBIT_SECURITY_TARGET=128 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -o $@

test/test_kyber145: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DBIT_SECURITY_TARGET=145 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -o $@

test/test_kyber163: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DBIT_SECURITY_TARGET=163 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -o $@

test/test_kyber192: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DBIT_SECURITY_TARGET=192 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -o $@

test/test_kyber199: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DBIT_SECURITY_TARGET=199 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -o $@

test/test_kyber224: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DBIT_SECURITY_TARGET=224 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -o $@

test/test_kyber235: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DBIT_SECURITY_TARGET=235 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -o $@

test/test_kyber256: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DBIT_SECURITY_TARGET=256 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -o $@

#test/test_kyber512: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
#	$(CC) $(CFLAGS) -DBIT_SECURITY_TARGET=512 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -o $@

test/test_speed120: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DBIT_SECURITY_TARGET=120 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -o $@

test/test_speed128: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DBIT_SECURITY_TARGET=128 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -o $@

test/test_speed145: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DBIT_SECURITY_TARGET=145 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -o $@

test/test_speed163: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DBIT_SECURITY_TARGET=163 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -o $@

test/test_speed192: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DBIT_SECURITY_TARGET=192 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -o $@

test/test_speed199: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DBIT_SECURITY_TARGET=199 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -o $@

test/test_speed224: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DBIT_SECURITY_TARGET=224 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -o $@

test/test_speed235: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DBIT_SECURITY_TARGET=235 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -o $@

test/test_speed256: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DBIT_SECURITY_TARGET=256 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -o $@

#test/test_speed512: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
#	$(CC) $(CFLAGS) -DBIT_SECURITY_TARGET=512 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -o $@

clean:
	-$(RM) -f *.gcno *.gcda *.lcov *.o *.so
	-$(RM) -f test/test_kyber120
	-$(RM) -f test/test_kyber128
	-$(RM) -f test/test_kyber145
	-$(RM) -f test/test_kyber163
	-$(RM) -f test/test_kyber192
	-$(RM) -f test/test_kyber199
	-$(RM) -f test/test_kyber224
	-$(RM) -f test/test_kyber235
	-$(RM) -f test/test_kyber256
#	-$(RM) -f test/test_kyber512
	-$(RM) -f test/test_speed120
	-$(RM) -f test/test_speed128
	-$(RM) -f test/test_speed145
	-$(RM) -f test/test_speed163
	-$(RM) -f test/test_speed192
	-$(RM) -f test/test_speed199
	-$(RM) -f test/test_speed224
	-$(RM) -f test/test_speed235
	-$(RM) -f test/test_speed256
#	-$(RM) -f test/test_speed512
	-$(RM) -rf lib/

