/*
Copyright (c) 2015-2016, Apple Inc. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:  

1.  Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2.  Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
    in the documentation and/or other materials provided with the distribution.

3.  Neither the name of the copyright holder(s) nor the names of any contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "lzfse_internal.h"

// Compare 64-bit unsigned integers for qsort.
static int compare_counts(const void *a, const void *b) {
  uint64_t aa = load8(a);
  uint64_t bb = load8(b);
  if (aa < bb)
    return -1;
  else if (aa > bb)
    return 1;
  return 0;
}

static inline unsigned get_state_generator(unsigned nstates) {
	return (nstates >> 1) | (nstates >> 3) | 3;
}

static inline unsigned ilog2_ceil(uint32_t n) {
	if (n <= 1)
		return 0;
	return 32 - __builtin_clz(n - 1);
}

// Initialize encoder table T[NSYMBOLS].
// NSTATES = sum FREQ[i] is the number of states (a power of 2)
// NSYMBOLS is the number of symbols.
// FREQ[NSYMBOLS] is a normalized histogram of symbol frequencies, with FREQ[i]
// >= 0.
// Some symbols may have a 0 frequency.  In that case, they should not be
// present in the data.
void fse_init_encoder_table(int num_states, int alphabet_size,
                            const uint16_t *__restrict state_counts,
                            fse_encoder_entry *__restrict sym_encinfo,
                            uint16_t *__restrict next_statesx) {
  const int log2_num_states = 31 - __builtin_clz(num_states);
  const unsigned state_generator = get_state_generator(num_states);
  const unsigned state_mask = num_states - 1;
  unsigned cumul_total;
  unsigned sym;
  unsigned state;
  unsigned count;
  unsigned max_bits;
  uint16_t cumul_state_counts[alphabet_size];
  uint8_t state_to_symbol[num_states];

  cumul_total = 0;
  for (sym = 0; sym < alphabet_size; sym++) {

    count = state_counts[sym];

    if (count == 0) /* Unused symbol? */
      continue;

    cumul_state_counts[sym] = cumul_total;

    max_bits = log2_num_states - ilog2_ceil(count) + 1;

    sym_encinfo[sym].adjusted_num_states_in_big_ranges =
      ((uint32_t)max_bits << 16) -
      ((uint32_t)count << max_bits);

    sym_encinfo[sym].next_states_begin = (int32_t)cumul_total - (int32_t)count;

    cumul_total += count;
  }

  state = 0;
  for (sym = 0; sym < alphabet_size; sym++) {
    count = state_counts[sym];
    while (count--) {
      state_to_symbol[state] = sym;
      state = (state + state_generator) & state_mask;
    }
  }

  for (state = 0; state < num_states; state++) {
    unsigned symbol = state_to_symbol[state];
    unsigned position = cumul_state_counts[symbol]++;
    next_statesx[position] = num_states + state;
  }
}

// Initialize decoder table T[NSTATES].
// NSTATES = sum FREQ[i] is the number of states (a power of 2)
// NSYMBOLS is the number of symbols.
// FREQ[NSYMBOLS] is a normalized histogram of symbol frequencies, with FREQ[i]
// >= 0.
// Some symbols may have a 0 frequency.  In that case, they should not be
// present in the data.
int fse_init_decoder_table(int nstates, int nsymbols,
                           uint16_t *__restrict freq,
                           int32_t *__restrict t) {
  fse_decoder_entry *decode_table = (fse_decoder_entry *)t;
  const int log2_num_states =  31 - __builtin_clz(nstates);
  const unsigned state_generator = get_state_generator(nstates);
  const unsigned state_mask = nstates - 1;
  unsigned state = 0;
  uint32_t total_count = 0;
  unsigned sym;

  for (sym = 0; sym < nsymbols; sym++) {
    unsigned count = freq[sym];
    if (count == 0)
      continue;
    total_count += count;
    do {
      decode_table[state].symbol = sym;
      state = (state + state_generator) & state_mask;
    } while (--count);
  }

  if (total_count != nstates)
	  return 0;

  for (state = 0; state < nstates; state++) {

    uint8_t sym = decode_table[state].symbol;
    uint32_t counter = freq[sym]++;
    unsigned num_bits = log2_num_states - (31 - __builtin_clz(counter));
    uint32_t destination_range_start = (counter << num_bits) - nstates;

    decode_table[state].k = num_bits;
    decode_table[state].delta = destination_range_start;
  }

  return 1;
}

// Initialize value decoder table T[NSTATES].
// NSTATES = sum FREQ[i] is the number of states (a power of 2)
// NSYMBOLS is the number of symbols.
// FREQ[NSYMBOLS] is a normalized histogram of symbol frequencies, with FREQ[i]
// >= 0.
// SYMBOL_VBITS[NSYMBOLS] and SYMBOLS_VBASE[NSYMBOLS] are the number of value
// bits to read and the base value for each symbol.
// Some symbols may have a 0 frequency.  In that case, they should not be
// present in the data.
void fse_init_value_decoder_table(int nstates, int nsymbols,
                                  uint16_t *__restrict freq,
                                  const uint8_t *__restrict symbol_vbits,
                                  const int32_t *__restrict symbol_vbase,
                                  fse_value_decoder_entry *__restrict t) {

  fse_value_decoder_entry *decode_table = t;
  const int log2_num_states =  31 - __builtin_clz(nstates);
  const unsigned state_generator = get_state_generator(nstates);
  const unsigned state_mask = nstates - 1;
  unsigned state = 0;
  uint32_t total_count = 0;
  unsigned sym;

  for (sym = 0; sym < nsymbols; sym++) {
    unsigned count = freq[sym];
    if (count == 0)
      continue;
    total_count += count;
    do {
      decode_table[state].total_bits = sym; // temporary
      state = (state + state_generator) & state_mask;
    } while (--count);
  }

  if (total_count != nstates)
	  return; // TODO

  for (state = 0; state < nstates; state++) {

    uint8_t sym = decode_table[state].total_bits;
    uint32_t counter = freq[sym]++;
    unsigned num_bits = log2_num_states - (31 - __builtin_clz(counter));
    uint32_t destination_range_start = (counter << num_bits) - nstates;

    decode_table[state].total_bits = symbol_vbits[sym] + num_bits;
    decode_table[state].value_bits = symbol_vbits[sym];
    decode_table[state].delta = destination_range_start;
    decode_table[state].vbase = symbol_vbase[sym];
  }
}

// Normalize a table T[NSYMBOLS] of symbols,occurrences to FREQ[NSYMBOLS].
// IMPORTANT: T will be modified (sorted) by this call.
// Return 1 if OK, and 0 on failure.
int fse_normalize_freq(int nstates, int nsymbols, fse_occurrence_entry *t,
                       uint16_t *freq) {
  // Sort T in increasing count order. Entry seen as a 64-bit value is (count <<
  // 32) + symbol_id
  qsort(t, nsymbols, sizeof(t[0]), compare_counts);

  // Get sum of T.count
  uint32_t s_count = 0;
  for (int i = 0; i < nsymbols; i++)
    s_count += t[i].count;

  // Start from low values
  uint32_t available = (uint32_t)nstates; // remaining available states
  for (int i = 0; i < nsymbols; i++) {
    uint32_t symbol = t[i].symbol;
    uint32_t count = t[i].count;

    if (count == 0) {
      freq[symbol] = 0;
      continue;
    } // no states used
    if (available <= 0 || s_count <= 0)
      return 0; // failed
    uint32_t k =
        (uint32_t)((double)count * (double)available / (double)s_count);
    if (k == 0)
      k = 1; // if count is not 0, we need at least 1 state to represent it
    if (i == nsymbols - 1)
      k = available; // force usage of all states
    freq[symbol] = (uint16_t)k;
    s_count -= count;
    available -= k;
  }
  return 1; // OK
}
