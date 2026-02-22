#include "warpcore.hpp"

#include <simd_detector.h>

#include "dsp_scalar.hpp"

namespace warpcore {
std::unique_ptr<IWarpCore> CreateDsp() {
    return std::make_unique<WarpCoreScalar>();
    // simd_detector::is_supported(simd_detector::InstructionSet::SSE2);
    // simd_detector::is_supported(simd_detector::InstructionSet::SSE4_1);
    // simd_detector::is_supported(simd_detector::InstructionSet::AVX);
    // simd_detector::is_supported(simd_detector::InstructionSet::AVX2);
}
}
