/*
 * rs — systematic shortened Reed–Solomon over GF(256), prim poly 0x11d.
 *
 * One fixed 255-byte codeword region per G4FSK data frame: k = 255 - P data
 * bytes (type | len | payload | crc | zero pad) followed by P parity bytes.
 * Shortened RS(255, 255-P) behaves exactly like the full code whose leading
 * 255-k message bytes are zero — the implicit zeros add nothing to the
 * syndromes — so encode/decode work directly on the 255 wire bytes.
 *
 * Parity tiers P in {16, 32, 64, 84} -> t = P/2 corrected bytes per block
 * (minimum distance P+1).  Validated by nsfsk_params_validate(); passing
 * anything else to these entry points is a bug (encode aborts, decode
 * returns -1).
 *
 * Position convention: block[j] carries the coefficient of x^(254-j) (first
 * wire byte = highest degree).  Syndrome S_i = r(alpha^i), i = 1..P.
 *
 * Decode runs syndromes -> Berlekamp-Massey -> Chien -> Forney, corrects in
 * place, then RE-VERIFIES the syndromes: a correction that does not zero
 * them is rejected (-1), never emitted.  No dynamic allocation.
 */
#ifndef NSFSK_RS_H
#define NSFSK_RS_H

#include <stdint.h>

#define NSFSK_RS_N     255u  /* codeword bytes (data region + parity)      */
#define NSFSK_RS_PMAX   84u  /* deepest supported tier                     */
#define NSFSK_RS_TMAX   42u  /* max corrected bytes per block              */

/* data is the k = 255 - P byte region; parity receives P bytes. */
void nsfsk_rs_encode(unsigned P, const uint8_t *data, uint8_t *parity);

/* Correct block[0..254] in place; on success copy the clean k-byte data
 * region to out_data (which may alias block) and return bytes corrected.
 * Return -1 when uncorrectable (> t errors or degenerate): the block then
 * holds garbage and MUST NOT be used. */
int nsfsk_rs_decode(unsigned P, uint8_t *block, uint8_t *out_data);

#endif /* NSFSK_RS_H */
