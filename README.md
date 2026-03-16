# WarpCore
尝试重现 Prosoniq PiWarp 或 Zynaptiq Wormhole 插件的效果。  
This project attempts to recreate the effect style of Prosoniq PiWarp or Zynaptiq Wormhole.

WarpCore 是一个`多段频谱反转`插件，类似于 PiWarp/Wormhole 的`时域局部频谱反转`效果。  
WarpCore is a `multi-band spectrum inversion` plugin, similar to the `local time-domain spectrum inversion` effect of PiWarp/Wormhole.

[YouTube](https://www.youtube.com/watch?v=7CM1Xm0MM6E)  
[Bilibili](https://www.bilibili.com/video/BV1UVAkzsEvP)

## 功能
## Features

- [x] 可配置极点数滤波器  
  Configurable filter pole count.

  包括：滤波器极点数量、滤波器截止频率缩放。  
  Includes: filter pole count and filter cutoff frequency scaling.

- [ ] 内置带通分频器  
  Built-in band-splitting crossover.

  包括：低分频和高分频（最后再做）。  
  Includes: low split and high split (to be done later).

- [ ] 可调整带宽  
  Adjustable bandwidth.

  包括：必须是递增函数的带宽分配。  
  Includes: bandwidth allocation that must follow a monotonic increasing function.

- [x] 共振峰移动  
  Formant shifting.

  应该是加在输出的振荡器上面。  
  This should be applied on the output oscillator stage.

## 图形界面
## GUI

![display](doc/gui.png)

> [!TIP]
> if you want a close preset to **default Piwarp**  
> set paramter  
> **Warp** = 50  
> **Freq High** = Full  
> **Scale** = 1.0  
> **Poles** = 2  
> **Pitch** = 0  
> **Base Mix** = 0  
> **Freq Mode** = music: 0 + 2n  

## macOS

```bash
sudo xattr -dr com.apple.quarantine /path/to/your/plugins/plugin_name.component
sudo xattr -dr com.apple.quarantine /path/to/your/plugins/plugin_name.vst3
sudo xattr -dr com.apple.quarantine /path/to/your/plugins/plugin_name.lv2
```

## 构建

```bash
git clone --recurse https://github.com/ManasWorld/WarpCore.git

# Windows
cmake -G "Ninja" -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang -DCMAKE_BUILD_TYPE=Release -S . -B build
cmake --build build --config Release

# Linux
sudo apt update
sudo apt-get install libx11-dev libfreetype-dev libfontconfig1-dev libasound2-dev libxrandr-dev libxinerama-dev libxcursor-dev
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -S . -B build
cmake --build build --config Release

# macOS
cmake -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" -S . -B build
cmake --build build --config Release
```
