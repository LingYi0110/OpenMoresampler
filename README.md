# OpenMoresampler

OpenMoresampler 是 Moresampler 0.8.4 的非官方开源实现，与 Moresampler 没有隶属关系。本项目的许可证仅适用于项目作者有权授权的代码，不包含、也不暗示授予对 Moresampler 原版程序及其名称、商标、模型、数据或其他受保护内容的任何许可。相关权利归其原作者或相应权利人所有。使用者应根据具体用途自行确认并取得所需授权。

项目仍处于兼容验证阶段，不保证输出与 Moresampler 0.8.4 完全一致。


## 构建

请使用 GCC 编译器

### 64 位

在 PowerShell 中进入源码目录后运行：

```powershell
cmake -S . -B build -G "MinGW Makefiles" -DMORESAMPLER_ARCH=64 -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

生成的程序位于 `build\moresampler.exe`。

### 32 位

32 位构建需要真正的 32 位 MinGW 工具链，或支持 `-m32` 且已安装对应运行库的 GCC/Clang：

```powershell
cmake -S . -B build32 -G "MinGW Makefiles" -DMORESAMPLER_ARCH=32 -DCMAKE_BUILD_TYPE=Release
cmake --build build32 --parallel
```

`MORESAMPLER_ARCH` 的默认值是 `32`。只安装了常见的 64 位 MinGW 时，请显式指定 `-DMORESAMPLER_ARCH=64`。


## 使用

### 在 UTAU 宿主中使用

将 `moresampler.exe` 同时配置为 resampler 和 wavtool。程序会根据宿主传入的参数自动选择运行模式，不提供图形界面。

首次处理一个音源时会生成分析缓存。建议在批量渲染前先备份音源目录，并确认生成的 `.llsm`、`.mrq` 和音频结果可以正常使用。

兼容的 resampler 参数格式为：

```text
moresampler <输入音频> <输出音频> <音高> <辅音速度> [flags] [offset] [length] [consonant] [cutoff] [volume] [modulation] [!tempo] [pitchbend]
```

一般用户不需要手动拼接这条命令，UTAU 或 OpenUtau 会负责传参。

### 自动生成 oto.ini

> 此功能尚未经过实际测试，可能出现分析失败、参数错误或生成的 `oto.ini` 无法使用等问题。运行前请备份已有的 `oto.ini` 和音源目录。

把音源目录作为唯一参数传入：

```powershell
.\moresampler.exe "\path\to\MyVoicebank"
```

程序会询问录音格式、统一后缀等信息，然后在目标目录写入 `oto.ini`。如果目标目录已有 `oto.ini`，程序会读取并可能覆盖它。自动生成至少需要 15 个可识别的 WAV 文件，也可以从同一目录的 `index.csv` 读取文件和发音信息。

### 合成 LLSM 数据

把兼容的 LLSM/LMIX 数据文件作为唯一参数传入时，程序会在当前工作目录写出 `output.wav`：

```powershell
.\moresampler.exe ".\sample.llsm"
```


## 目录结构

| 路径 | 内容 |
| --- | --- |
| `main.c` | 运行模式检测和程序入口 |
| `resampler.*` | UTAU resampler 参数解析与合成流程 |
| `wavtool.*` | 音符拼接、包络和最终音频输出 |
| `lmix/` | LLSM 中间数据的分析、修改和合成 |
| `xchg.*` | LLSM/LMIX 数据交换格式 |
| `nebula/` | Nebula F0 估计 |
| `otogen/` | 自动原音设定生成 |
| `external/` | 随源码编译的第三方组件 |

## 已知限制

- 32位与原版存在少量数值差异
- 当前只支持 Windows。
- 尚不能保证所有 flag、配置项和缓存版本与原版一致。
- 32 位和 64 位构建可能在浮点结果、性能和旧缓存兼容性上存在差异。
- oto.ini 自动生成未经测试，可能存在错误


## 致谢

- Moresampler 及其作者 [Kanru Hua](https://github.com/Sleepwalking)
- [libllsm](https://github.com/Sleepwalking/libllsm)
- [Nebula](https://github.com/Sleepwalking/nebula)
- [SHIRO](https://github.com/Sleepwalking/SHIRO)
- [ciglet](https://github.com/Sleepwalking/ciglet)
- [cJSON](https://github.com/DaveGamble/cJSON)
- [fastapprox](https://github.com/etheory/fastapprox)
- [FastLZ](https://github.com/ariya/FastLZ)
- [libaiff](https://aifftools.sourceforge.net/libaiff/)
- [liblrhsmm](https://github.com/Sleepwalking/liblrhsmm)
- [libgvps](https://github.com/Sleepwalking/libgvps)
- [Lua 5.3.3](https://www.lua.org/versions.html#5.3)
- [mrq](https://github.com/Sleepwalking/mrq)
