#pragma once
#include <numbers>
#include <complex>
#include "pluginshared/simd.hpp"
#include "global.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"

namespace warpcore {
static consteval int AlignUp(int x, int align) noexcept {
    return ((x + (align - 1)) / align) * align;
}

template <simd::IsSimdFloat SimdT>
struct SvfLaneN {
    static constexpr int kSvfCoeffSize = AlignUp(global::kMaxPoles, simd::LaneSize<SimdT>);

    std::array<float, kSvfCoeffSize> g;
    std::array<float, kSvfCoeffSize> d;
    std::array<float, kSvfCoeffSize> last_g;
    std::array<float, kSvfCoeffSize> last_d;

    struct SvfState {
       SimdT s1_re_l;
       SimdT s1_re_r;
       SimdT s1_im_l;
       SimdT s1_im_r;
       SimdT s2_re_l;
       SimdT s2_re_r;
       SimdT s2_im_l;
       SimdT s2_im_r;
    };
    simd::Array256<SvfState, global::kMaxBands * global::kMaxPoles / simd::LaneSize<SimdT>> state;

    void Init() noexcept {
        g.fill(0.0f);
        d.fill(0.0f);
        last_d.fill(0.0f);
        last_g.fill(0.0f);
    }

    void Reset() noexcept {
        state.fill(SvfState{});
    }

    void SetPole(int pole_idx, float w, float Q, float fmul) noexcept {
        float r2 = 1.0f / Q;
        float g_val = std::tan(w / 2.0f) * fmul;
        g[pole_idx] = g_val;
        d[pole_idx] = 1.0f / (1.0f + r2 * g_val + g_val * g_val);
    }

    void SetFreq(float w_lowpass, int num_filters) noexcept {
        constexpr float atten = 0.5f;
        constexpr float square_epsi = (1.0f - atten * atten) / atten;
        float analog_fmul = 1.0f / std::pow(square_epsi, 0.25f / static_cast<float>(num_filters));

        int n = 2 * num_filters;
        for (int k = 1; k <= num_filters; ++k) {
            float phi =
                (2.0f * static_cast<float>(k) - 1.0f) * std::numbers::pi_v<float> / (2.0f * static_cast<float>(n));
            float Q = 1.0f / (2.0f * std::sin(phi));
            SetPole(k - 1, w_lowpass, Q, analog_fmul);
        }
    }
};

enum class FreqDistrbution {
    k0_n,
    k1_n,
    k0_2n,
    k1_2n,
};

struct ProcessorState {
    float fs{};
    int num_warps{};
    int poles{};
    float base_mix{};
    bool pitch_affect{};
    FreqDistrbution freq_distribution{};

    // complex sine generator
    float pre_osc_phase{};
    float pre_osc_phase_inc{};
    float last_pre_osc_phase_inc{};

    float post_osc_phase{};
    float post_osc_phase_inc{};
    float last_post_osc_phase_inc{};

    // complex lowpass filter
    union {
        SvfLaneN<simd::Float128> svf128;
        SvfLaneN<simd::Float256> svf256;
    };

    std::complex<float> band0_s1[global::kMaxPoles * 2]{};
    std::complex<float> band0_s2[global::kMaxPoles * 2]{};
};

struct Param {
    int bands;
    float f_low;
    float f_high;
    float filter_scale;
    int filter_order;
    float pitch_shift;
    float base_mix;
    bool pitch_affect;
    FreqDistrbution freq_distribution;
};

struct ProcessorDsp {
    void(*init)(ProcessorState& state, float fs) noexcept;
    void(*reset)(ProcessorState& state) noexcept;
    void(*update)(ProcessorState& state, const Param& p) noexcept;
    void(*process)(ProcessorState& state, float* left, float* right, int num_samples) noexcept;

    const char* name;

    bool IsValid() const noexcept {
        return init != nullptr;
    }
};

ProcessorDsp GetProcessorDsp() noexcept;
}

#pragma GCC diagnostic pop
