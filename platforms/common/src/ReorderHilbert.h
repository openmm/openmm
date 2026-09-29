/* Hilbert Curve implementation copyright 1998, Rice University.
 * Original algorithm: Doug Moore; Copyright (c) 1998-2000, Rice University.
 * Modified by Codex on 2026-09-29: exact 3D/8-bit table specialization for
 * the existing atom-reorder call, with the original generic fallback.
 */
/* LICENSE
 *
 * This software is copyrighted by Rice University.  It may be freely copied,
 * modified, and redistributed, provided that the copyright notice is
 * preserved on all copies.
 *
 * There is no warranty or other guarantee of fitness for this software,
 * it is provided solely "as is".  Bug reports or fixes may be sent
 * to the author, who may or may not act on them as he desires.
 *
 * You may include this software in a program or other software product,
 * but must display the notice:
 *
 * Hilbert Curve implementation copyright 1998, Rice University
 *
 * in any place where the end-user would see your own copyright.
 *
 * If you modify this software, you should include a notice giving the
 * name of the person performing the modification, the date of modification,
 * and the reason for such modification.
 */

#ifndef OPENMM_REORDER_HILBERT_H_
#define OPENMM_REORDER_HILBERT_H_

#include "hilbert.h"
#include <cstdint>

namespace OpenMM {

// Specialize the bundled Doug Moore/Rice University Hilbert algorithm for the
// existing nDims=3, nBits=8 call. The bundled hilbert.cpp license applies to
// that algorithm. No floating-point arithmetic, clamping or key changes.
inline int computeReorderHilbert3D8(int x, int y, int z, bool enabled) {
    if (!enabled || static_cast<unsigned int>(x) > 255 ||
            static_cast<unsigned int>(y) > 255 || static_cast<unsigned int>(z) > 255) {
        const bitmask_t coords[3] = {static_cast<bitmask_t>(x), static_cast<bitmask_t>(y), static_cast<bitmask_t>(z)};
        return static_cast<int>(hilbert_c2i(3, 8, coords));
    }
    // Spread the eight coordinate bits into Morton bit positions 0,3,...,21.
    static const std::uint32_t spread[256] = {
        0x000000u, 0x000001u, 0x000008u, 0x000009u, 0x000040u, 0x000041u, 0x000048u, 0x000049u,
        0x000200u, 0x000201u, 0x000208u, 0x000209u, 0x000240u, 0x000241u, 0x000248u, 0x000249u,
        0x001000u, 0x001001u, 0x001008u, 0x001009u, 0x001040u, 0x001041u, 0x001048u, 0x001049u,
        0x001200u, 0x001201u, 0x001208u, 0x001209u, 0x001240u, 0x001241u, 0x001248u, 0x001249u,
        0x008000u, 0x008001u, 0x008008u, 0x008009u, 0x008040u, 0x008041u, 0x008048u, 0x008049u,
        0x008200u, 0x008201u, 0x008208u, 0x008209u, 0x008240u, 0x008241u, 0x008248u, 0x008249u,
        0x009000u, 0x009001u, 0x009008u, 0x009009u, 0x009040u, 0x009041u, 0x009048u, 0x009049u,
        0x009200u, 0x009201u, 0x009208u, 0x009209u, 0x009240u, 0x009241u, 0x009248u, 0x009249u,
        0x040000u, 0x040001u, 0x040008u, 0x040009u, 0x040040u, 0x040041u, 0x040048u, 0x040049u,
        0x040200u, 0x040201u, 0x040208u, 0x040209u, 0x040240u, 0x040241u, 0x040248u, 0x040249u,
        0x041000u, 0x041001u, 0x041008u, 0x041009u, 0x041040u, 0x041041u, 0x041048u, 0x041049u,
        0x041200u, 0x041201u, 0x041208u, 0x041209u, 0x041240u, 0x041241u, 0x041248u, 0x041249u,
        0x048000u, 0x048001u, 0x048008u, 0x048009u, 0x048040u, 0x048041u, 0x048048u, 0x048049u,
        0x048200u, 0x048201u, 0x048208u, 0x048209u, 0x048240u, 0x048241u, 0x048248u, 0x048249u,
        0x049000u, 0x049001u, 0x049008u, 0x049009u, 0x049040u, 0x049041u, 0x049048u, 0x049049u,
        0x049200u, 0x049201u, 0x049208u, 0x049209u, 0x049240u, 0x049241u, 0x049248u, 0x049249u,
        0x200000u, 0x200001u, 0x200008u, 0x200009u, 0x200040u, 0x200041u, 0x200048u, 0x200049u,
        0x200200u, 0x200201u, 0x200208u, 0x200209u, 0x200240u, 0x200241u, 0x200248u, 0x200249u,
        0x201000u, 0x201001u, 0x201008u, 0x201009u, 0x201040u, 0x201041u, 0x201048u, 0x201049u,
        0x201200u, 0x201201u, 0x201208u, 0x201209u, 0x201240u, 0x201241u, 0x201248u, 0x201249u,
        0x208000u, 0x208001u, 0x208008u, 0x208009u, 0x208040u, 0x208041u, 0x208048u, 0x208049u,
        0x208200u, 0x208201u, 0x208208u, 0x208209u, 0x208240u, 0x208241u, 0x208248u, 0x208249u,
        0x209000u, 0x209001u, 0x209008u, 0x209009u, 0x209040u, 0x209041u, 0x209048u, 0x209049u,
        0x209200u, 0x209201u, 0x209208u, 0x209209u, 0x209240u, 0x209241u, 0x209248u, 0x209249u,
        0x240000u, 0x240001u, 0x240008u, 0x240009u, 0x240040u, 0x240041u, 0x240048u, 0x240049u,
        0x240200u, 0x240201u, 0x240208u, 0x240209u, 0x240240u, 0x240241u, 0x240248u, 0x240249u,
        0x241000u, 0x241001u, 0x241008u, 0x241009u, 0x241040u, 0x241041u, 0x241048u, 0x241049u,
        0x241200u, 0x241201u, 0x241208u, 0x241209u, 0x241240u, 0x241241u, 0x241248u, 0x241249u,
        0x248000u, 0x248001u, 0x248008u, 0x248009u, 0x248040u, 0x248041u, 0x248048u, 0x248049u,
        0x248200u, 0x248201u, 0x248208u, 0x248209u, 0x248240u, 0x248241u, 0x248248u, 0x248249u,
        0x249000u, 0x249001u, 0x249008u, 0x249009u, 0x249040u, 0x249041u, 0x249048u, 0x249049u,
        0x249200u, 0x249201u, 0x249208u, 0x249209u, 0x249240u, 0x249241u, 0x249248u, 0x249249u,
    };
    // State=(rotation*4+flipOrdinal), flipOrdinal maps 0,1,2,3 to 0,1,2,4.
    // Each entry packs nextState in the high bits and the output digit in 0..2.
    // Generated from rotateRight(), flipBit=1<<rotation and adjust_rotation().
    static const unsigned char transition[12][8] = {
        {40, 73, 10, 75, 44, 77, 14, 79},
        {73, 40, 75, 10, 77, 44, 79, 14},
        {10, 75, 40, 73, 14, 79, 44, 77},
        {44, 77, 14, 79, 40, 73, 10, 75},
        {80, 84, 17, 21, 50, 54, 19, 23},
        {84, 80, 21, 17, 54, 50, 23, 19},
        {17, 21, 80, 84, 19, 23, 50, 54},
        {50, 54, 19, 23, 80, 84, 17, 21},
        {24, 90, 28, 94, 57, 59, 61, 63},
        {90, 24, 94, 28, 59, 57, 63, 61},
        {28, 94, 24, 90, 61, 63, 57, 59},
        {57, 59, 61, 63, 24, 90, 28, 94},
    };
    std::uint32_t coords = spread[x] | (spread[y] << 1) | (spread[z] << 2);
    coords ^= coords >> 3;
    unsigned int state = 0;
    std::uint32_t index = 0;
    for (int shift = 21; shift >= 0; shift -= 3) {
        const unsigned int step = transition[state][(coords >> shift)&7u];
        index = (index << 3) | (step&7u);
        state = step >> 3;
    }
    // Original nthbits=(2^24-1)/7 followed by the same Gray decoding.
    index ^= 0x124924u;
    index ^= index >> 1;
    index ^= index >> 2;
    index ^= index >> 4;
    index ^= index >> 8;
    index ^= index >> 16;
    return static_cast<int>(index);
}

} // namespace OpenMM
#endif
