# kfbf-converter — Pure-Python KFBio Fluorescent Slide Converter

## Overview

**kfbf-converter** is a fast, low-memory converter for computational pathology. It
turns proprietary **KFBio `.kfbf` fluorescent whole-slide images** into QuPath-ready
multi-channel **OME-TIFFs** and standard **BigTIFFs**. QuPath, Bio-Formats, napari,
OMERO and most other image-analysis tools can open the results.

## Key Features

- **Pyramidal OME-TIFF** (`.ome.tif`): every raw fluorescent channel in one file,
  with:
  - the physical pixel size from the slide header (e.g. **0.2528 µm/px**)
  - the channel labels (DAPI, FITC, TRITC, Cy5, …) and their display colours
- **Composite RGB BigTIFF** (`.rgb.tif`): the channels merged into one 8-bit RGB
  image, with percentile contrast stretching (0.1–99.9 %) so dim markers stay
  visible next to bright DAPI.
- **Split grayscale BigTIFFs** (`_<channel>.tif`): one file per channel,
  e.g. `slide_DAPI.tif`.
- **Pure Python**: no C++/Rust compiler, vendor SDK or DLLs needed. All
  dependencies install from pip.
- **Low RAM use**: whole slides are processed tile by tile, so a full-resolution
  image never has to fit in memory. A 20 000 × 12 000 px, 3-channel slide converts
  to every format using about 250 MB of RAM.
- **Batch mode and interactive menu**: point it at a folder to convert every
  slide, or just run it and pick a format.

Every output is a tiled, pyramidal file. It holds the full resolution plus
downsampled levels, so it opens and zooms instantly in QuPath.

---

## Quick Start (Zero Coding / Non-Technical Users)

1. **Install Python 3.10+** from [python.org](https://www.python.org/downloads/).
   During setup, tick **"Add Python to PATH"**.
2. **Install the dependencies.** Open a command prompt in this folder and run:
   ```bash
   pip install -r requirements.txt
   ```
3. **Copy both `kfbf_converter.py` and `run_converter.bat`** into the folder
   that contains your `.kfbf` files. The launcher looks for the script next to
   itself, so the two files must stay together.
4. **Double-click `run_converter.bat`** (or run `python kfbf_converter.py` in
   that folder). You will see:
   ```text
   Found 3 slide(s).
   Select output format:
     [1] OME-TIFF (Default / QuPath)
     [2] Composite RGB BigTIFF
     [3] Split Grayscale Channels
     [4] All of the above
   Enter choice [1-4] (Press Enter for 1):
   ```
5. **Press Enter** to choose the default QuPath OME-TIFF mode.

The converted files appear next to your slides. Existing outputs are skipped,
so running it again only converts new slides.

---

## Command Line Usage

**Single file**
```bash
python kfbf_converter.py slide.kfbf --mode ome
```

**Entire folder (batch)**
```bash
python kfbf_converter.py "C:\path\to\slides" --mode all -o "C:\path\to\output"
```

**Inspect a slide** (dimensions, pixel size, channels, pyramid levels), without
converting it:
```bash
python kfbf_converter.py slide.kfbf --info
```

### Output modes

| `--mode` | Creates | Best for |
|---|---|---|
| `ome` *(default)* | `slide.ome.tif` | QuPath and Bio-Formats analysis with the raw channels |
| `rgb` | `slide.rgb.tif` | Quick viewing, figures, tools that only accept RGB |
| `split` | `slide_DAPI.tif`, `slide_FITC.tif`, … | Tools that expect one channel per file |
| `all` | All of the above | |

### Common options

| Option | Default | What it does |
|---|---|---|
| `-o, --output-dir` | next to the input | Where to save the converted files |
| `--compression` | `deflate` | `deflate` or `lzw` (lossless), or `jpeg` (smaller, lossy) |
| `--overwrite` | off | Replace files that already exist |
| `--recursive` | off | Also convert slides in sub-folders |
| `--channel-names` | from slide | Rename the channels, e.g. `--channel-names DAPI CD3 CD8` |
| `--mpp` | `0.2528` | Pixel size (µm/px), used only if the slide header lacks one |
| `--tile-size` | `512` | Output tile size in pixels |
| `--min-level-size` | `256` | Smallest pyramid level to keep, in pixels |
| `-q, --quiet` | off | Hide progress bars |

Run `python kfbf_converter.py --help` for the full list.

---

## QuPath Integration Guide

1. **Drag and drop** the `.ome.tif` file into the QuPath window. If QuPath asks
   for an image type, choose **Fluorescence**.
2. Open **Brightness & Contrast** (**Shift + C**) to show or hide each channel
   and set its display range.
3. Channel names, colours and the pixel size (µm) are already set, so
   measurements and scale bars are calibrated from the start.

---

## Troubleshooting

| Problem | Fix |
|---|---|
| `Python is not installed or not in PATH` | Reinstall Python and tick **Add Python to PATH** |
| `No module named ...` | Run `pip install -r requirements.txt` |
| `No .kfbf files found` | Put the slides in the same folder you run the converter from |
| Wrong or generic channel names | Set them with `--channel-names` |
| Image looks scrambled or wrong | Run `--info` on the slide and share the output |

## Requirements

Python 3.10+ with `tifffile`, `imagecodecs`, `numpy`, `pillow` and `tqdm`
(see `requirements.txt`).

---

<details>
<summary><strong>Technical Details &amp; Python API</strong></summary>

### The `.kfbf` format (reverse-engineered)

`.kfbf` is an undocumented, proprietary format. The layout below comes from
reverse engineering and may not cover every scanner or software version. All
values are little-endian.

**File header**

| Offset | Type | Field |
|---|---|---|
| `0x00` | 4 bytes | Magic `F1 01 EE EE` |
| `0x04` | 4 bytes | ASCII `KFBF` |
| `0x0C` | float32 | Format version (e.g. 2.1) |
| `0x10` | uint32 | Tile record count N |
| `0x14` | uint32 | Image height (level 0) |
| `0x18` | uint32 | Image width (level 0) |
| `0x1C` | uint32 | Scan magnification (e.g. 40) |
| `0x20` | 8 bytes | Codec string (`JPEG`) |
| `0x44` | uint64 | Offset of the tile index table |
| `0x4C` | float32 | Pixel size in µm/px (e.g. 0.2528) |
| `0x5C` | — | Metadata block: marker `FF 01 EE EE`, then `tag:uint32, length:uint32, value` entries. Tag `0x4B` = channel count; tags `0x4C`+ = channel labels |

**Tile index**: marker `F1 04 EE EE`, followed by N records of 64 bytes each:

| Offset | Type | Field |
|---|---|---|
| `+0x00` | uint32 | Tile X in **pixels** within its level (`col = x // 256`) |
| `+0x04` | uint32 | Tile Y in **pixels** within its level (`row = y // 256`) |
| `+0x08` / `+0x0C` | uint32 | Unknown (shown by `--info` as `raw+08/+0C`) |
| `+0x10` | float32 | Level magnification (40, 20, 10, … down to ~0.0012) |
| `+0x14` | uint32 | Z-plane index (0 for 2-D slides) |
| `+0x20` | uint64 | Offset of C × uint64 payload offsets (one per channel) |
| `+0x28` | uint64 | Offset of C × uint32 or uint64 payload lengths (width detected automatically) |
| `+0x38` | uint32 | Tail marker `0xEEEE04FF` |
| `+0x3C` | uint32 | Delimiter `0xEEEE04F1` |

Each channel tile is a 256 × 256, 8-bit grayscale baseline JPEG. Missing tiles
in a sparse grid are treated as black. Pyramid levels are stored down to
1 × 1 px; the converter exports only levels whose larger side is at least
256 px (`--min-level-size`).

OME-XML channel colours are written as signed 32-bit **RGBA**, as the OME
schema specifies. For example, blue (DAPI) is `65535` and red is `-16776961`.

### Python API

```python
from kfbf_converter import KFBioSlide, convert_slide

with KFBioSlide("slide.kfbf") as slide:
    print(slide.summary())                      # same text as --info
    print(slide.dimensions)                     # (width, height) at level 0
    print(slide.channel_names, slide.n_levels, slide.n_z, slide.mpp)

    # One raw 256 x 256 uint8 source tile (zeros if missing)
    tile = slide.get_tile(level=0, col=10, row=4, channel=0, z=0)

    # Any region of one channel, at any pyramid level -> (height, width) uint8
    region = slide.read_region(level=2, x=0, y=0, width=1024, height=1024, channel=1)

    # All channels at once -> (C, height, width) uint8
    stack = slide.read_region_channels(level=3, x=0, y=0, width=512, height=512)

# Convert from code; returns the list of files written
written = convert_slide("slide.kfbf", "output/", mode="all", compression="deflate")
```

</details>
