#pragma once
#include <cmath>
#include <complex>
#include <numbers>

namespace warpcore {
class SvfTPT {
public:
    void Reset() noexcept {
        s1_ = 0;
        s2_ = 0;
    }

    /**
     * @param w [0, pi]
     * @param r2 <0:不稳定 =0:无阻尼 0~1:复共轭极点 >1:分裂成两个单极点
     */
    void SetCoeffSVF(float w, float r2) noexcept {
        g_ = std::tan(w / 2);
        R2_ = r2;
        d_ = 1 / (1 + r2 * g_ + g_ * g_);
    }

    void SetCoeffQ(float w, float Q) noexcept {
        SetCoeffSVF(w, 1 / Q);
    }

    std::complex<float> TickLowpass(std::complex<float> x) noexcept {
        auto bp = d_ * (g_ * (x - s2_) + s1_);
        auto v1 = bp - s1_;
        auto v2 = g_ * bp;
        auto lp = v2 + s2_;
        s1_ = bp + v1;
        s2_ = lp + v2;
        return lp;
    }
private:
    float R2_{};
    float g_{};
    float d_{};
    std::complex<float> s1_{};
    std::complex<float> s2_{};
};

class VicSineOsc {
public:
    void Reset() noexcept {
        u_ = 1.0f;
        v_ = 0.0f;
    }

    void KeepAmp() noexcept {
        float g = 1.0f / std::sqrt(u_ * u_ + v_ * v_);
        u_ *= g;
        v_ *= g;
    }

    float Tick() noexcept {
        float w = u_ - k1_ * v_;
        v_ = v_ + k2_ * w;
        u_ = w - k1_ * v_;
        return v_;
    }

    /**
     * @param f hz
     * @param fs sample rate
     */
    void SetFreq(float f, float fs) noexcept {
        auto omega = f / fs * std::numbers::pi_v<float> * 2.0f;
        k1_ = std::tan(omega / 2.0f);
        k2_ = 2 * k1_ / (1 + k1_ * k1_);
    }

    /**
     * @param w [-pi, pi]
     */
    void SetFreq(float w) noexcept {
        k1_ = std::tan(w / 2.0f);
        k2_ = 2 * k1_ / (1 + k1_ * k1_);
    }

    std::complex<float> GetCpx() const noexcept {
        return {u_, v_};
    }
private:
    float k1_{};
    float k2_{};
    float u_{1.0f};
    float v_{};
};

class SingleWarp {
public:
    float Tick(std::complex<float> osc, float x) noexcept {
        auto x1 = x * osc;
        auto x2 = svf_.TickLowpass(x1);
        x2 = svf2_.TickLowpass(x2);
        x2 = svf3_.TickLowpass(x2);
        x2 = svf4_.TickLowpass(x2);
        auto x3 = x2 * osc;
        return x3.real();
    }

    void SetFreq(float fcenter) noexcept {
        svf_.SetCoeffQ(fcenter, 0.50979558f);
        svf2_.SetCoeffQ(fcenter, 0.60134489f);
        svf3_.SetCoeffQ(fcenter, 0.89997622f);
        svf4_.SetCoeffQ(fcenter, 2.5629154f);
    }

    void Reset() noexcept {
        svf_.Reset();
    }
private:
    SvfTPT svf_;
    SvfTPT svf2_;
    SvfTPT svf3_;
    SvfTPT svf4_;
};

class WarpCore {
public:
    void Init(float fs) {
        fs_ = fs;
    }

    void Reset() noexcept {
        for (auto& w : warps_) {
            w.Reset();
        }
        osc_.Reset();
    }

    void Process(float* left, float* right, int num_samples) noexcept {
        for (int i = 0; i < num_samples; i++) {
            osc_.Tick();
            std::complex<float> osc_mul = osc_.GetCpx();
            osc_mul = std::conj(osc_mul);
            std::complex<float> osc = osc_mul;

            auto x = left[i];
            float y = 0;
            for (int j = 0; j < num_warps_; j++) {
                y += warps_[j].Tick(osc, x);
                osc *= osc_mul;
            }
            
            left[i] = y;
        }
        std::copy_n(left, num_samples, right);
    }

    struct Param {
        bool changed;
        int bands;
        float f_low;
        float f_high;
        float scale;
    };
    void Update(const Param& p) noexcept {
        num_warps_ = p.bands;

        float finc = p.f_high / static_cast<float>(p.bands);
        float fbase = finc / 2;
        float wbase = Freq2Omega(fbase);
        osc_.SetFreq(wbase);

        float filter_w = wbase * p.scale;
        filter_w = std::min(filter_w, std::numbers::pi_v<float> / 2 - 1e-5f);
        for (int i = 0; i < p.bands; i++) {
            warps_[i].SetFreq(filter_w);
        }
    }

    float Freq2Omega(float f) noexcept {
        return 2 * std::numbers::pi_v<float> * f / fs_;
    }
private:
    float fs_{};
    SingleWarp warps_[512];

    int num_warps_{};
    VicSineOsc osc_{};
};
}
