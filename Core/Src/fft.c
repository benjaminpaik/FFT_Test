/*
 * fft.c
 *
 *  Created on: Jan 22, 2021
 *      Author: bpaik
 */

#include "stdlib.h"
#include "string.h"
#include "fft.h"

void fft_indexing(FFT *fft);
void fft_index_stage(uint16_t index_map[], uint16_t length);
void fft_butterfly(FFT *fft, uint16_t offset, uint16_t length, uint16_t w_step);

void fft_setup(FFT *fft, uint16_t log2_length, float_t fs)
{
  float_t angle, nyquist_freq = fs / 2.0F;
  uint16_t i;
  // calculate buffer lengths
  fft->log2_length = log2_length;
  fft->stages = fft->log2_length - 1;
  fft->length = (1U << fft->log2_length);
  fft->input_index = 0;
  // allocate memory
  fft->index_map = (uint16_t*) malloc(fft->length * sizeof(uint16_t));
  fft->input = (float_t*) malloc(fft->length * sizeof(float_t));
  fft->temp = (COMPLEX*) malloc(fft->length * sizeof(COMPLEX));
  fft->out = (COMPLEX*) malloc(fft->length * sizeof(COMPLEX));
  fft->w = (COMPLEX*) malloc(((fft->length / 2) * sizeof(COMPLEX)));
  fft->freqs = (float_t*) malloc((fft->length / 2) * sizeof(float_t));
  fft->mag_out = (float_t*) malloc((fft->length / 2) * sizeof(float_t));
  fft->phase_out = (float_t*) malloc((fft->length / 2) * sizeof(float_t));
  // generate index map
  fft_indexing(fft);
  // angle increment
  angle = (-2.0F * M_PI) / fft->length;

  for (i = 0; i < (fft->length / 2); i++) {
    fft->w[i].re = cosf(angle * i);
    fft->w[i].im = sinf(angle * i);
    fft->freqs[i] = (nyquist_freq * i) / (fft->length / 2);
  }
}

void fft_indexing(FFT *fft)
{
  uint16_t i, j, segments, elements, start_index;

  // initialize the index map
  for (i = 0; i < fft->length; i++) {
    fft->index_map[i] = i;
  }

  // re-assign indices
  for (i = 0; i < fft->stages; i++) {
    segments = (1U << i);
    elements = fft->length / segments;

    for (j = 0; j < segments; j++) {
      start_index = elements * j;
      fft_index_stage(fft->index_map + start_index, elements);
    }
  }
}

void fft_index_stage(uint16_t index_map[], uint16_t length)
{
  uint16_t i, input_index;
  uint16_t *even, *odd;

  // allocate memory
  even = (uint16_t*) malloc((length / 2) * sizeof(uint16_t));
  odd = (uint16_t*) malloc((length / 2) * sizeof(uint16_t));

  for (i = 0; i < (length / 2); i++) {
    input_index = 2 * i;
    even[i] = index_map[input_index];
    odd[i] = index_map[input_index + 1];
  }

  for (i = 0; i < (length / 2); i++) {
    index_map[i] = even[i];
    index_map[i + (length / 2)] = odd[i];
  }
  free(even);
  free(odd);
}

void fft_butterfly(FFT *fft, uint16_t offset, uint16_t length, uint16_t w_step)
{
  COMPLEX *temp, *output, x_scaled;
  uint16_t half_length = length / 2;
  uint16_t i;

  // set relative buffers
  temp = fft->temp + offset;
  output = fft->out + offset;
  // load the temporary buffer
  memcpy((void*) temp, (void*) output, sizeof(COMPLEX) * length);

  for (i = 0; i < half_length; i++) {
    x_scaled = complex_product(fft->w[w_step * i], temp[half_length + i]);
    output[i] = complex_sum(temp[i], x_scaled);
    output[i + half_length] = complex_diff(temp[i], x_scaled);
  }
}

uint16_t fft_input(FFT *fft, float_t input)
{
  if (fft->input_index < fft->length) {
    fft->input[fft->input_index++] = input;
    return (fft->input_index == fft->length);
  }
  else {
    return FALSE;
  }
}

void fft_compute(FFT *fft)
{
  uint16_t i, j, sections, offset;
  uint16_t b_length = 2;

  for (i = 0; i < fft->length; i++) {
    fft->out[i].re = fft->input[fft->index_map[i]];
    fft->out[i].im = 0;
  }

  for (i = 0; i <= fft->stages; i++) {
    sections = (1U << (fft->stages - i));
    for (j = 0; j < sections; j++) {
      offset = j * b_length;
      fft_butterfly(fft, offset, b_length, sections);
    }
    b_length <<= 1;
  }

  fft->input_index = 0;
}

void fft_convert_out(FFT *fft)
{
  uint16_t i = 0;

  for(i = 0; i < (fft->length / 2); i++) {
    fft->mag_out[i] = complex_mag(fft->out[i]);
    fft->phase_out[i] = complex_phase(fft->out[i]);
  }
}

