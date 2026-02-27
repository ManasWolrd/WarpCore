#include "dsp_state.hpp"

#include <cassert>
#include <complex>

namespace warpcore {
template <int kPoles, bool kPitchAffect>
static void ProcessInternal(warpcore::ProcessorState& state, float* left, float* right, int num_samples) noexcept {
    constexpr simd::Array256<simd::Float256, 8> kBandGainLut{
        simd::Float256{1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f},
        simd::Float256{1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
        simd::Float256{1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
        simd::Float256{1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
        simd::Float256{1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f},
        simd::Float256{1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f},
        simd::Float256{1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f},
        simd::Float256{1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f},
    };
    simd::Float256 tail_gain = kBandGainLut[state.num_warps & 7];

    float band0_dry_mix = state.base_mix;
    float band0_wet_mix = 1.0f - band0_dry_mix;
    int simd_loop_count = (state.num_warps + 7) / 8;

    if (state.num_warps <= 8) {
        tail_gain[0] = band0_wet_mix;
    }

    for (int i = 0; i < num_samples; i++) {
        // -------------------- tick complex sine generators --------------------
        state.pre_osc_phase += state.pre_osc_phase_inc;
        state.pre_osc_phase -= std::floor(state.pre_osc_phase);

        state.post_osc_phase += state.post_osc_phase_inc;
        state.post_osc_phase -= std::floor(state.post_osc_phase);

        constexpr float twopi = 2.0f * std::numbers::pi_v<float>;
        std::complex<float> pre_osc_f32 = {std::cos(state.pre_osc_phase * twopi), std::sin(state.pre_osc_phase * twopi)};
        std::complex<float> post_osc_f32 = {std::cos(state.post_osc_phase * twopi), std::sin(state.post_osc_phase * twopi)};
        if constexpr (kPitchAffect) {
            std::swap(pre_osc_f32, post_osc_f32);
        }

        const auto pre_osc_f32_0 = pre_osc_f32;
        const auto pre_osc_f32_1 = pre_osc_f32 * pre_osc_f32;
        const auto pre_osc_f32_2 = pre_osc_f32 * pre_osc_f32 * pre_osc_f32;
        const auto pre_osc_f32_3 = pre_osc_f32 * pre_osc_f32 * pre_osc_f32 * pre_osc_f32;
        const auto pre_osc_f32_4 = pre_osc_f32_0 * pre_osc_f32_3;
        const auto pre_osc_f32_5 = pre_osc_f32_1 * pre_osc_f32_3;
        const auto pre_osc_f32_6 = pre_osc_f32_2 * pre_osc_f32_3;
        const auto pre_osc_f32_7 = pre_osc_f32_3 * pre_osc_f32_3;

        const auto post_osc_f32_0 = post_osc_f32;
        const auto post_osc_f32_1 = post_osc_f32 * post_osc_f32;
        const auto post_osc_f32_2 = post_osc_f32 * post_osc_f32 * post_osc_f32;
        const auto post_osc_f32_3 = post_osc_f32 * post_osc_f32 * post_osc_f32 * post_osc_f32;
        const auto post_osc_f32_4 = post_osc_f32_0 * post_osc_f32_3;
        const auto post_osc_f32_5 = post_osc_f32_1 * post_osc_f32_3;
        const auto post_osc_f32_6 = post_osc_f32_2 * post_osc_f32_3;
        const auto post_osc_f32_7 = post_osc_f32_3 * post_osc_f32_3;

        const simd::Complex256 pre_osc{
            .re = simd::BroadcastF256(pre_osc_f32_7.real()),
            .im = simd::BroadcastF256(pre_osc_f32_7.imag()),
        };
        const simd::Complex256 post_osc{
            .re = simd::BroadcastF256(post_osc_f32_7.real()),
            .im = simd::BroadcastF256(post_osc_f32_7.imag()),
        };
        simd::Complex256 pre_osc_n_val{
            .re = {pre_osc_f32_0.real(), pre_osc_f32_1.real(), pre_osc_f32_2.real(), pre_osc_f32_3.real(),
                   pre_osc_f32_4.real(), pre_osc_f32_5.real(), pre_osc_f32_6.real(), pre_osc_f32_7.real()},
            .im = {pre_osc_f32_0.imag(), pre_osc_f32_1.imag(), pre_osc_f32_2.imag(), pre_osc_f32_3.imag(),
                   pre_osc_f32_4.imag(), pre_osc_f32_5.imag(), pre_osc_f32_6.imag(), pre_osc_f32_7.imag()},
        };
        simd::Complex256 post_osc_n_val{
            .re = {post_osc_f32_0.real(), post_osc_f32_1.real(), post_osc_f32_2.real(), post_osc_f32_3.real(),
                   post_osc_f32_4.real(), post_osc_f32_5.real(), post_osc_f32_6.real(), post_osc_f32_7.real()},
            .im = {post_osc_f32_0.imag(), post_osc_f32_1.imag(), post_osc_f32_2.imag(), post_osc_f32_3.imag(),
                   post_osc_f32_4.imag(), post_osc_f32_5.imag(), post_osc_f32_6.imag(), post_osc_f32_7.imag()},
        };

        // -------------------- process first band --------------------
        float x_left = left[i];
        float x_right = right[i];

        std::complex<float> cpx_x_left = x_left * pre_osc_f32;
        std::complex<float> cpx_x_right = x_right * pre_osc_f32;
        #pragma unroll
        for (int k = 0; k < kPoles; ++k) {
            const float gk = state.svf256.g[k];
            const float dk = state.svf256.d[k];
            auto& s1_l = state.band0_s1[2 * k];
            auto& s1_r = state.band0_s1[2 * k + 1];
            auto& s2_l = state.band0_s2[2 * k];
            auto& s2_r = state.band0_s2[2 * k + 1];

            auto bp_l = dk * (gk * (cpx_x_left - s2_l) + s1_l);
            auto bp_r = dk * (gk * (cpx_x_right - s2_r) + s1_r);
            auto v1_l = bp_l - s1_l;
            auto v1_r = bp_r - s1_r;
            auto v2_l = gk * bp_l;
            auto v2_r = gk * bp_r;
            auto lp_l = v2_l + s2_l;
            auto lp_r = v2_r + s2_r;
            
            s1_l = bp_l + v1_l;
            s1_r = bp_r + v1_r;
            s2_l = lp_l + v2_l;
            s2_r = lp_r + v2_r;
            cpx_x_left = lp_l;
            cpx_x_right = lp_r;
        }
        float y_l = (cpx_x_left * std::conj(post_osc_f32)).real() * band0_dry_mix;
        float y_r = (cpx_x_right * std::conj(post_osc_f32)).real() * band0_dry_mix;

        // -------------------- process bands --------------------
        auto* svf_state = state.svf256.state.data();
        simd::Float256 band_gain{band0_wet_mix, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};

        for (int j = 0; j < simd_loop_count - 1; ++j) {
            // std::complex<float> tmp = x * pre_osc_n_val;
            // pre_osc_n_val *= pre_osc;
            auto tmp_l = simd::BroadcastF256(x_left) * pre_osc_n_val;
            auto tmp_r = simd::BroadcastF256(x_right) * pre_osc_n_val;
            pre_osc_n_val *= pre_osc;

            #pragma unroll
            for (int k = 0; k < kPoles; ++k) {
                const float gk = state.svf256.g[k];
                const float dk = state.svf256.d[k];

                auto s1_re = svf_state->s1_re;
                auto s1_im = svf_state->s1_im;
                auto s2_re = svf_state->s2_re;
                auto s2_im = svf_state->s2_im;

                auto bp_re = dk * (gk * (tmp_l.re - s2_re) + s1_re);
                auto bp_im = dk * (gk * (tmp_l.im - s2_im) + s1_im);
                auto v1_re = bp_re - s1_re;
                auto v1_im = bp_im - s1_im;
                auto v2_re = gk * bp_re;
                auto v2_im = gk * bp_im;
                auto lp_re = v2_re + s2_re;
                auto lp_im = v2_im + s2_im;

                svf_state->s1_re = bp_re + v1_re;
                svf_state->s1_im = bp_im + v1_im;
                svf_state->s2_re = lp_re + v2_re;
                svf_state->s2_im = lp_im + v2_im;

                ++svf_state;
                
                tmp_l.re = lp_re;
                tmp_l.im = lp_im;
            }

            #pragma unroll
            for (int k = 0; k < kPoles; ++k) {
                const float gk = state.svf256.g[k];
                const float dk = state.svf256.d[k];

                auto s1_re = svf_state->s1_re;
                auto s1_im = svf_state->s1_im;
                auto s2_re = svf_state->s2_re;
                auto s2_im = svf_state->s2_im;

                auto bp_re = dk * (gk * (tmp_r.re - s2_re) + s1_re);
                auto bp_im = dk * (gk * (tmp_r.im - s2_im) + s1_im);
                auto v1_re = bp_re - s1_re;
                auto v1_im = bp_im - s1_im;
                auto v2_re = gk * bp_re;
                auto v2_im = gk * bp_im;
                auto lp_re = v2_re + s2_re;
                auto lp_im = v2_im + s2_im;

                svf_state->s1_re = bp_re + v1_re;
                svf_state->s1_im = bp_im + v1_im;
                svf_state->s2_re = lp_re + v2_re;
                svf_state->s2_im = lp_im + v2_im;

                ++svf_state;
                
                tmp_r.re = lp_re;
                tmp_r.im = lp_im;
            }

            // y += (tmp * post_osc_n_val).real();
            // post_osc_n_val *= post_osc;
            auto band_out_l = tmp_l * post_osc_n_val;
            auto band_out_r = tmp_r * post_osc_n_val;
            y_l += simd::ReduceAdd(band_out_l.re * band_gain);
            y_r += simd::ReduceAdd(band_out_r.re * band_gain);
            post_osc_n_val *= post_osc;

            band_gain = simd::BroadcastF256(1.0f);
        }

        // -------------------- here we have: 1/2/3/4 --------------------
        // std::complex<float> tmp = x * pre_osc_n_val;
        // pre_osc_n_val *= pre_osc;
        auto tmp_l = simd::BroadcastF256(x_left) * pre_osc_n_val;
        auto tmp_r = simd::BroadcastF256(x_right) * pre_osc_n_val;

        #pragma unroll
        for (int k = 0; k < kPoles; ++k) {
            const float gk = state.svf256.g[k];
            const float dk = state.svf256.d[k];

            auto s1_re = svf_state->s1_re;
            auto s1_im = svf_state->s1_im;
            auto s2_re = svf_state->s2_re;
            auto s2_im = svf_state->s2_im;

            auto bp_re = dk * (gk * (tmp_l.re - s2_re) + s1_re);
            auto bp_im = dk * (gk * (tmp_l.im - s2_im) + s1_im);
            auto v1_re = bp_re - s1_re;
            auto v1_im = bp_im - s1_im;
            auto v2_re = gk * bp_re;
            auto v2_im = gk * bp_im;
            auto lp_re = v2_re + s2_re;
            auto lp_im = v2_im + s2_im;

            svf_state->s1_re = bp_re + v1_re;
            svf_state->s1_im = bp_im + v1_im;
            svf_state->s2_re = lp_re + v2_re;
            svf_state->s2_im = lp_im + v2_im;

            ++svf_state;
            
            tmp_l.re = lp_re;
            tmp_l.im = lp_im;
        }

        #pragma unroll
        for (int k = 0; k < kPoles; ++k) {
            const float gk = state.svf256.g[k];
            const float dk = state.svf256.d[k];

            auto s1_re = svf_state->s1_re;
            auto s1_im = svf_state->s1_im;
            auto s2_re = svf_state->s2_re;
            auto s2_im = svf_state->s2_im;

            auto bp_re = dk * (gk * (tmp_r.re - s2_re) + s1_re);
            auto bp_im = dk * (gk * (tmp_r.im - s2_im) + s1_im);
            auto v1_re = bp_re - s1_re;
            auto v1_im = bp_im - s1_im;
            auto v2_re = gk * bp_re;
            auto v2_im = gk * bp_im;
            auto lp_re = v2_re + s2_re;
            auto lp_im = v2_im + s2_im;

            svf_state->s1_re = bp_re + v1_re;
            svf_state->s1_im = bp_im + v1_im;
            svf_state->s2_re = lp_re + v2_re;
            svf_state->s2_im = lp_im + v2_im;

            ++svf_state;
            
            tmp_r.re = lp_re;
            tmp_r.im = lp_im;
        }

        // y += (tmp * post_osc_n_val).real();
        // post_osc_n_val *= post_osc;
        auto band_out_l = tmp_l * post_osc_n_val;
        auto band_out_r = tmp_r * post_osc_n_val;
        y_l += simd::ReduceAdd(band_out_l.re * tail_gain);
        y_r += simd::ReduceAdd(band_out_r.re * tail_gain);

        left[i] = y_l;
        right[i] = y_r;
    }
}

template <bool kPitchAffect>
static void ProcessPoles(warpcore::ProcessorState& state, float* left, float* right, int num_samples) noexcept {
    switch (state.poles) {
        case 1:
            ProcessInternal<1, kPitchAffect>(state, left, right, num_samples);
            break;
        case 2:
            ProcessInternal<2, kPitchAffect>(state, left, right, num_samples);
            break;
        case 3:
            ProcessInternal<3, kPitchAffect>(state, left, right, num_samples);
            break;
        case 4:
            ProcessInternal<4, kPitchAffect>(state, left, right, num_samples);
            break;
        case 5:
            ProcessInternal<5, kPitchAffect>(state, left, right, num_samples);
            break;
        case 6:
            ProcessInternal<6, kPitchAffect>(state, left, right, num_samples);
            break;
        case 7:
            ProcessInternal<7, kPitchAffect>(state, left, right, num_samples);
            break;
        case 8:
            ProcessInternal<8, kPitchAffect>(state, left, right, num_samples);
            break;
        default:
            assert(false);
            break;
    }
}

// ----------------------------------------
// dsp processor
// ----------------------------------------

static void Init(warpcore::ProcessorState& state, float fs) noexcept {
    state.fs = fs;
}

static void Reset(warpcore::ProcessorState& state) noexcept {
    state.svf256.Reset();
    state.pre_osc_phase = 0.0f;
    state.post_osc_phase = 0.0f;

    std::fill_n(state.band0_s1, global::kMaxPoles * 2, std::complex<float>{});
    std::fill_n(state.band0_s2, global::kMaxPoles * 2, std::complex<float>{});
}

static void Update(warpcore::ProcessorState& state, const warpcore::Param& p) noexcept {
    state.num_warps = p.bands;
    state.base_mix = p.base_mix;
    state.pitch_affect = p.pitch_affect;

    float fhigh = p.f_high;
    if (fhigh > 40000.0f) {
        fhigh = state.fs;
    }
    fhigh = std::min(fhigh, state.fs);

    float finc = fhigh / static_cast<float>(p.bands);
    float fbase = finc / 2;
    float fshit = p.pitch_affect ? -p.pitch_shift : p.pitch_shift;
    fshit = std::exp2(fshit / 12.0f);
    state.osc_base_freq = fbase;
    state.pre_osc_phase_inc = fbase / state.fs;
    state.post_osc_phase_inc = state.pre_osc_phase_inc * fshit;

    // butterworth lowpass
    float wbase = fbase * 2 * std::numbers::pi_v<float> / state.fs;
    if (p.pitch_affect) {
        wbase *= fshit;
    }
    float filter_w = wbase * p.filter_scale;
    filter_w = std::min(filter_w, std::numbers::pi_v<float> - 0.1f);

    if (state.poles != p.filter_order) {
        state.svf256.Reset();
    }
    state.poles = p.filter_order;
    state.svf256.SetFreq(filter_w, p.filter_order);
}

static void Process(warpcore::ProcessorState& state, float* left, float* right, int num_samples) noexcept {
    if (state.pitch_affect) {
        ProcessPoles<true>(state, left, right, num_samples);
    }
    else {
        ProcessPoles<false>(state, left, right, num_samples);
    }
}

// ----------------------------------------
// export
// ----------------------------------------

#ifndef DSP_EXPORT_NAME
#error "不应该编译这个文件,在其他cpp包含这个cpp并定义DSP_EXPORT_NAME=`dsp_dispatch.cpp里的变量`"
#endif

ProcessorDsp DSP_EXPORT_NAME{
    Init,
    Reset,
    Update,
    Process
};
}
