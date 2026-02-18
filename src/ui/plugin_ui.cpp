#include "plugin_ui.hpp"

#include "../PluginProcessor.h"

PluginUi::PluginUi(EmptyAudioProcessor& p)
    : preset_(*p.preset_manager_) {
    addAndMakeVisible(preset_);

    auto& apvt = *p.value_tree_;
    warp_.BindParam(apvt, "warp");
    addAndMakeVisible(warp_);

    f_low_.BindParam(apvt, "scale");
    addAndMakeVisible(f_low_);

    f_high_.BindParam(apvt, "f_high");
    addAndMakeVisible(f_high_);

    setSize(400, 300);
}

void PluginUi::resized() {
    auto b = getLocalBounds();
    preset_.setBounds(b.removeFromTop(30));

    auto line = b.removeFromTop(65);
    warp_.setBounds(line.removeFromLeft(50));
    f_low_.setBounds(line.removeFromLeft(50));
    f_high_.setBounds(line.removeFromLeft(50));
}

void PluginUi::paint(juce::Graphics& g) {
    
}
