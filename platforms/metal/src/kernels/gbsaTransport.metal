/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * Adapted from the OpenMM Common kernels used by CUDA, HIP, and OpenCL:
 * platforms/common/src/kernels/gbsaObc.cc
 * Original code: Portions copyright (c) 2009-2026 Stanford University and
 * the Authors. Authors: Peter Eastman and the OpenMM contributors.
 * Metal code: Portions copyright (c) 2026 Chun-Chi Hung.
 * Author: Chun-Chi Hung.
 * This program is free software under the GNU Lesser General Public License,
 * version 3 or (at your option) any later version. It is distributed without
 * any warranty; see <http://www.gnu.org/licenses/> for the license.
 * -------------------------------------------------------------------------- */

/**
 * @brief Broadcast the Born-sum inputs of one lane for a diagonal tile.
 * Call before particle-validity/cutoff branches so every source lane is active.
 * The Common kernel owns AtomData1; this template changes only its transport.
 */
template <class Atom>
inline Atom metalGbsaBroadcastBorn(Atom data, uint lane) {
    Atom result = {};
    result.x = simdShuffle(data.x, lane);
    result.y = simdShuffle(data.y, lane);
    result.z = simdShuffle(data.z, lane);
    result.radius = simdShuffle(data.radius, lane);
    result.scaledRadius = simdShuffle(data.scaledRadius, lane);
    return result;
}

/** @brief Rotate a complete Born-sum particle, including its lane accumulator. */
template <class Atom>
inline Atom metalGbsaRotateBorn(Atom data, uint lane) {
    Atom result;
    result.x = simdShuffle(data.x, lane);
    result.y = simdShuffle(data.y, lane);
    result.z = simdShuffle(data.z, lane);
    result.q = simdShuffle(data.q, lane);
    result.radius = simdShuffle(data.radius, lane);
    result.scaledRadius = simdShuffle(data.scaledRadius, lane);
    result.bornSum = simdShuffle(data.bornSum, lane);
    return result;
}

/** @brief Broadcast the Force1 inputs of one lane before pair-validity tests. */
template <class Atom>
inline Atom metalGbsaBroadcastForce(Atom data, uint lane) {
    Atom result = {};
    result.x = simdShuffle(data.x, lane);
    result.y = simdShuffle(data.y, lane);
    result.z = simdShuffle(data.z, lane);
    result.q = simdShuffle(data.q, lane);
    result.bornRadius = simdShuffle(data.bornRadius, lane);
    return result;
}

/**
 * @brief Rotate coordinates and all four Force1 accumulators together.
 * After TILE_SIZE=32 rotations, the accumulated particle returns to its owner.
 */
template <class Atom>
inline Atom metalGbsaRotateForce(Atom data, uint lane) {
    Atom result;
    result.x = simdShuffle(data.x, lane);
    result.y = simdShuffle(data.y, lane);
    result.z = simdShuffle(data.z, lane);
    result.q = simdShuffle(data.q, lane);
    result.fx = simdShuffle(data.fx, lane);
    result.fy = simdShuffle(data.fy, lane);
    result.fz = simdShuffle(data.fz, lane);
    result.fw = simdShuffle(data.fw, lane);
    result.bornRadius = simdShuffle(data.bornRadius, lane);
    return result;
}
