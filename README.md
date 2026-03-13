# WarpCore
尝试重现 prosoniq piwarp 或者 zynaptiq wormhole 插件的效果  
WarpCore是一个`多段频谱反转`插件，类似于piwarp/wormhole的`时域局部频谱反转`效果  

[youtube](https://www.youtube.com/watch?v=7CM1Xm0MM6E)  
[bilibili](https://www.bilibili.com/video/BV1UVAkzsEvP)

## featrue
- [x] 可配置极点数滤波器
    包括: 滤波器极点数量 滤波器截止频率缩放
- [ ] 内置带通分频器
    包括: 低分频和高分频(最后再做)
- [ ] 可调整带宽
    包括: 必须是递增的函数的带宽分配
- [x] 共振峰移动
    应该是加在输出的振荡器上面

## GUI
![display](doc/gui.png)

## MacOS

```bash
sudo xattr -dr com.apple.quarantine /path/to/your/plugins/plugin_name.component
sudo xattr -dr com.apple.quarantine /path/to/your/plugins/plugin_name.vst3
sudo xattr -dr com.apple.quarantine /path/to/your/plugins/plugin_name.lv2
```
## build

```bash
git clone --recurse https://github.com/ManasWorld/plugin-template.git

# windows
cmake -G "Ninja" -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang -DCMAKE_BUILD_TYPE=Release -S . -B build
cmake --build build --config Release

# linux
sudo apt update
sudo apt-get install libx11-dev libfreetype-dev libfontconfig1-dev libasound2-dev libxrandr-dev libxinerama-dev libxcursor-dev
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -S . -B .build
cmake --build build --config Release

# macOS
cmake -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" -S . -B build
cmake --build build --config Release
```
