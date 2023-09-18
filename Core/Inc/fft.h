/*
 * fft.h
 *
 *  Created on: Jan 22, 2021
 *      Author: bpaik
 */

#ifndef FFT_H_
#define FFT_H_

#include "stm32g4xx_hal.h"
#include "math.h"

typedef struct {
  float_t re;
  float_t im;
} COMPLEX;

typedef struct {
  uint16_t log2_length;
  uint16_t length;
  uint16_t stages;
  uint16_t *index_map;
  uint16_t input_index;

  float_t *input;
  COMPLEX *temp;
  COMPLEX *out;
  COMPLEX *w;

  float_t *freqs;
  float_t *mag_out;
  float_t *phase_out;
} FFT;

#ifndef TRUE
  #define TRUE  1
#endif
#ifndef FALSE
  #define FALSE   0
#endif


static inline COMPLEX complex_sum(COMPLEX x1, COMPLEX x2)
{
  COMPLEX result;
  result.re = x1.re + x2.re;
  result.im = x1.im + x2.im;
  return result;
}

static inline COMPLEX complex_diff(COMPLEX x1, COMPLEX x2)
{
  COMPLEX result;
  result.re = x1.re - x2.re;
  result.im = x1.im - x2.im;
  return result;
}

static inline COMPLEX complex_product(COMPLEX x1, COMPLEX x2)
{
  COMPLEX result;
  result.re = (x1.re * x2.re) - (x1.im * x2.im);
  result.im = (x1.re * x2.im) + (x1.im * x2.re);
  return result;
}

static inline float_t complex_mag(COMPLEX x)
{
	return sqrtf((x.re * x.re) + (x.im * x.im));
}

static inline float_t complex_phase(COMPLEX x)
{
	return atan2f(x.im, x.re);
}

void fft_setup(FFT *fft, uint16_t log2_length, float_t fs);
uint16_t fft_input(FFT *fft, float_t input);
void fft_compute(FFT *fft);
void fft_convert_out(FFT *fft);

#endif /* FFT_H_ */
