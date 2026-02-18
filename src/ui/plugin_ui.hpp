#pragma once
#include "pluginshared/component.hpp"
#include "pluginshared/preset_panel.hpp"

class EmptyAudioProcessor;

class PluginUi : public juce::Component {
public:
    static constexpr int kWidth = 480;
    static constexpr int kHeight = 320;

    explicit PluginUi(EmptyAudioProcessor& p);

    void resized() override;
    void paint(juce::Graphics& g) override;
    
    std::function<void(int,int)> on_want_new_size;
private:
    void TrySetSize(int width, int height) {
        if (on_want_new_size) {
            on_want_new_size(width, height); 
        }
    }
    pluginshared::PresetPanel preset_;

    ui::Dial warp_{"warp"};
    ui::Dial f_low_{"scale"};
    ui::Dial f_high_{"f_high"};
};
