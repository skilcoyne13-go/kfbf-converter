# kfbf-converter — Pure-Python KFBio Slide Converter

## Overview

**kfbf-converter** is a fast, low-memory converter for computational pathology. It
turns **KFBio whole-slide scans** into QuPath-ready open formats:

- **fluorescence** scans (`.kfbf`, proprietary) into multi-channel **OME-TIFFs**
  and standard **BigTIFFs**
- **brightfield** scans (`.svs`, as written by KFBio scanners) into an
  **RGB pyramidal BigTIFF** with standard pixel-size tags

QuPath, Bio-Formats, napari, OMERO and most other image-analysis tools can open
the results.

### Supported inputs

| Scan | Input file | Outputs |
|---|---|---|
| Fluorescence | `.kfbf` | OME-TIFF (default), composite RGB BigTIFF, split grayscale BigTIFFs |
| Brightfield | `.svs` | RGB pyramidal BigTIFF (`<name>.tif`) |

The scanner also writes QPTIFF for fluorescence itself, so that is not offered here.
Older KFBio `.kfb` files are not supported.

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
- **Brightfield RGB BigTIFF** (`.tif`) from `.svs`: the scanner's JPEG tiles are
  copied unchanged, so there is **no extra quality loss** and the file is about
  the size of the original (e.g. 94 MB `.svs` → 91 MB `.tif` in under a second).
- **Clear size errors**: before writing, every output is checked against the
  file format, the destination drive (e.g. FAT32's 4 GB file limit) and the
  free space. If it can't fit, you get
  *"Data is too big for the file format chosen."* with the reason and a fix,
  never a corrupt file. Every finished file is also read back and checked
  before it gets its final name.
- **Pure Python**: no C++/Rust compiler, vendor SDK or DLLs needed. All
  dependencies install from pip.
- **Low RAM use**: whole slides are processed tile by tile, so a full-resolution
  image never has to fit in memory. A 20 000 × 12 000 px, 3-channel slide converts
  to every format using about 300 MB of RAM, a 44 000 × 47 000 px one to
  OME-TIFF using about 370 MB.
- **Batch mode and interactive menu**: point it at a folder to convert every
  slide (fluorescence and brightfield can be mixed), or just run it and pick a
  format.

Every output is a tiled, pyramidal file. It holds the full resolution plus
downsampled levels, so it opens and zooms instantly in QuPath. The objective
magnification (e.g. 40×) is stored too, so QuPath can show "40.0x" rather
than a zoom factor.

---

## Quick Start (Zero Coding / Non-Technical Users)

1. **Install Python 3.10+** from [python.org](https://www.python.org/downloads/).
   During setup, tick **"Add Python to PATH"**.
2. **Install the dependencies.** Open a command prompt in this folder and run:
   ```bash
   pip install -r requirements.txt
   ```
3. **Copy both `kfbf_converter.py` and `run_converter.bat`** into the folder
   that contains your `.kfbf` and/or `.svs` files. The launcher looks for the
   script next to itself, so the two files must stay together.
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
   If the folder also has brightfield slides, the first line shows both counts,
   e.g. `Found 3 fluorescence (.kfbf) and 2 brightfield (.svs) slide(s).` The
   format you pick applies to the fluorescence slides; brightfield slides always
   become an RGB BigTIFF. A folder with only brightfield slides skips the menu
   and just asks you to press Enter to start.
5. **Press Enter** to choose the default QuPath OME-TIFF mode.

The converted files appear next to your slides. Existing outputs are skipped,
so running it again only converts new slides.

If a slide is too big for the drive you're saving to, you'll see
`Data is too big for the file format chosen.` followed by what to do (see
[Troubleshooting](#troubleshooting)). The other slides still convert.

---

## Command Line Usage

**Single file**
```bash
python kfbf_converter.py slide.kfbf --mode ome
python kfbf_converter.py slide.svs
```

**Entire folder (batch)**: `.kfbf` and `.svs` can be mixed; each slide gets
the output that fits its type
```bash
python kfbf_converter.py "C:\path\to\slides" --mode all -o "C:\path\to\output"
```

**Inspect a slide** (dimensions, magnification, pixel size, channels, pyramid
levels), without converting it:
```bash
python kfbf_converter.py slide.kfbf --info
python kfbf_converter.py slide.svs --info
```

### Output modes (fluorescence `.kfbf`)

| `--mode` | Creates | Best for |
|---|---|---|
| `ome` *(default)* | `slide.ome.tif` | QuPath and Bio-Formats analysis with the raw channels |
| `rgb` | `slide.rgb.tif` | Quick viewing, figures, tools that only accept RGB |
| `split` | `slide_DAPI.tif`, `slide_FITC.tif`, … | Tools that expect one channel per file |
| `all` | All of the above | |

### Brightfield output (`.svs`)

Each `.svs` gives one file, `slide.tif`: a tiled, pyramidal RGB BigTIFF
(photometric RGB, chunky samples) with the pixel size in standard TIFF
resolution tags. It uses the same page layout as an Aperio `.svs`, so viewers
that read slides with OpenSlide see the whole pyramid, the objective
magnification and the pixel size (see *Magnification* below). `--mode` does not apply; if you pass one, a single
info line says it was ignored for that slide.

By default the scanner's original JPEG tiles are **copied unchanged**: no
re-compression, no quality loss, about the same size as the `.svs`. The
thumbnail, slide label and macro photo are **not** copied, so label text
doesn't end up in the output. If you set `--compression` or `--tile-size`, the
tiles are decoded and re-encoded instead. On the 94 MB KFBio sample this gave:

| Command | Output size | Time | Pixels vs `.svs` |
|---|---|---|---|
| *(default: tile copy)* | 91 MB | < 1 s | identical |
| `--compression jpeg` (quality 90) | 111 MB | 3 s | re-compressed (lossy again) |
| `--compression deflate` | 697 MB | 8 s | identical |
| `--compression lzw` | 853 MB | 5 s | identical |

The scanner already saves brightfield tiles as JPEG (quality 30), so the default
copy is both the smallest and the most faithful option.

### Common options

| Option | Default | What it does |
|---|---|---|
| `-o, --output-dir` | next to the input | Where to save the converted files |
| `--compression` | `deflate` (`.svs`: copy tiles) | `deflate` or `lzw` (lossless), or `jpeg` (smaller, lossy) |
| `--overwrite` | off | Replace files that already exist |
| `--recursive` | off | Also convert slides in sub-folders |
| `--channel-names` | from slide | Rename the channels, e.g. `--channel-names DAPI CD3 CD8` (fluorescence only) |
| `--mpp` | `0.2528` | Pixel size (µm/px), used only if the slide doesn't contain one |
| `--magnification` | from slide | Objective magnification (e.g. `40`), if the slide's value is missing or wrong |
| `--tile-size` | `512` (`.svs` tile copy: `256`) | Output tile size in pixels |
| `--min-level-size` | `256` | Smallest pyramid level to keep, in pixels |
| `--no-preflight` | off | Skip the size check before writing (see below) |
| `-q, --quiet` | off | Hide progress bars |
| `-v, --verbose` | off | Debug output, including size estimates and error tracebacks |

Run `python kfbf_converter.py --help` for the full list.

### Magnification

The objective magnification is read from the slide (`.kfbf` header,
`.svs` `AppMag`) and written into every output. `--info` shows it, e.g.
`Magnification : 40x (header)`. If it's missing or wrong, set it with
`--magnification 40`.

| Output | How the magnification is stored | Checked with |
|---|---|---|
| `.ome.tif` | OME-XML `Instrument > Objective NominalMagnification`, linked from the image | official OME 2016-06 schema |
| `.tif` (brightfield), `.rgb.tif`, `_DAPI.tif` … | Aperio layout: levels as top-level pages, a thumbnail as page 2, and `Aperio Format … \|AppMag = 40\|MPP = 0.2528` descriptions | OpenSlide 4.0.1: all levels, `objective-power` and `mpp` read correctly |

Before 1.1, the BigTIFFs stored lower levels as SubIFDs, which OpenSlide
ignores: it saw only full resolution, with no pyramid and no magnification.
Fluorescence Z-stacks (several Z planes) keep the SubIFD layout for their
BigTIFFs, because the Aperio layout holds a single plane; their OME-TIFF still
has the magnification.

Bio-Formats (the other reader QuPath uses) has **not** been tested on the
Aperio-layout `.tif` files; it may list each level as a separate image.

### Size checks and file safety

Before each output file is written, the converter estimates its size and
checks:

- **File format**: outputs are always BigTIFF, which has no practical size
  limit (classic TIFF stops at 4 GB). Image width and height must fit TIFF's
  32-bit size fields. JPEG's 65 535 px limit applies to each tile, and tiles
  are at most a few thousand px, so JPEG compression is safe for any slide size.
- **Destination drive**: FAT32 (and older FAT) drives, which many USB sticks
  and SD cards use, can't hold files over 4 GB. Outputs expected to exceed
  that are refused. If only the worst case (no compression at all) would exceed
  it, you get a warning instead.
- **Free space**: the estimate plus a safety margin (256 MB + 5 %) must fit
  on the drive.

If an output can't fit, it is skipped, the rest of the batch continues, and
it counts as a failure in the final summary:

```text
Data is too big for the file format chosen.
  slide.ome.tif: would be about 5.20 GB, but the output drive E:\ is FAT32, which cannot store files larger than 4.00 GB.
  Fix: choose an NTFS or exFAT drive with -o, or use --compression jpeg (about 1.60 GB).
```

The same message, with the real cause, appears if the disk fills up or the
drive refuses a large file **while** writing.

Files are written to `<name>.part` first. When writing finishes, the file is
re-opened and checked: page and pyramid level counts, the size of every level,
tile data inside the file, and one decoded tile per level. Only then is it
renamed. If anything goes wrong, the partial file is deleted, so you never
get a `.tif` that can't be opened. Leftover `.part` files from a run that was
killed (e.g. by closing the window) are removed the next time you convert
that slide.

The estimate is deliberately cautious. If it ever refuses a file you know
will fit, `--no-preflight` turns the check off. Written files are still
verified.

---

## QuPath Integration Guide

**Fluorescence (`.ome.tif`)**

1. **Drag and drop** the `.ome.tif` file into the QuPath window. If QuPath asks
   for an image type, choose **Fluorescence**.
2. Open **Brightness & Contrast** (**Shift + C**) to show or hide each channel
   and set its display range.
3. Channel names, colours and the pixel size (µm) are already set, so
   measurements and scale bars are calibrated from the start.
4. The objective magnification (e.g. 40×) is stored in the OME-XML
   (`Instrument > Objective > NominalMagnification`). In the **Image** tab,
   **Magnification** should read 40, and at full resolution the toolbar zoom
   should show **40.0x** instead of "1.0x". *This was checked against the OME
   schema, but not yet in QuPath itself.*

**Brightfield (`.tif` from `.svs`)**

1. **Drag and drop** the `.tif` into QuPath.
2. Set the image type to **Brightfield (H&E)** if QuPath asks, or later under
   **Image → Set image type**. QuPath can also guess it automatically
   (*Preferences → Auto-estimate image type on opening*).
3. In the **Image** tab, check that **Pixel width** and **Pixel height** match
   the `MPP` shown by `--info` (e.g. 0.4859 µm), that **Magnification** is 20
   (KFBio brightfield scans are 20×), and that zooming out is smooth (pyramid
   levels are loaded).
4. Zoom to full resolution: the toolbar should show **20.0x**.

*Checked with OpenSlide and tifffile, but not yet in QuPath itself.* If QuPath
shows no magnification or only one resolution level, note which reader it
used (**Image** tab → **Server type**) and share it.

---

## Troubleshooting

| Problem | Fix |
|---|---|
| `Python is not installed or not in PATH` | Reinstall Python and tick **Add Python to PATH** |
| `No module named ...` | Run `pip install -r requirements.txt` |
| `No .kfbf or .svs files found` | Put the slides in the same folder you run the converter from |
| `Data is too big for the file format chosen.` | Follow the **Fix:** line under it. Usually: save to an NTFS or exFAT drive with `-o` (FAT32 USB sticks can't hold files over 4 GB), free up disk space, or use `--compression jpeg` for a smaller file. If you're sure it fits, try `--no-preflight` |
| `Output failed verification and was deleted: ...` | The written file didn't read back correctly and was removed. Check the drive (try another one with `-o`), then run again with `-v` and share the output |
| `unrecognised slide file` | Only KFBio `.kfbf` (fluorescence) and Aperio `.svs` (brightfield) slides are supported; older `.kfb` files are not |
| `already converted by kfbf_converter` | That file is an output of this tool; convert the original `.svs` / `.kfbf` instead |
| QuPath shows a zoom like "1.33x" instead of "40.0x" | Check `--info`: if it says `Magnification : unknown`, convert again with `--magnification 40` (and `--overwrite`) |
| Wrong or generic channel names | Set them with `--channel-names` |
| Image looks scrambled or wrong | Run `--info` on the slide and share the output |

## Requirements

Python 3.10+ with `tifffile`, `imagecodecs`, `numpy`, `pillow` and `tqdm`
(see `requirements.txt`).

To run the tests (they generate small synthetic slides, so no slide data is
needed):
```bash
pip install -r requirements-dev.txt
python -m pytest tests
```

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

### The `.svs` brightfield format (as written by KFBio)

KFBio scanners write brightfield scans as Aperio-style `.svs` files. SVS is a
standard tiled TIFF, so nothing needed reverse engineering. The layout below
was confirmed by inspecting three real KFBio `.svs` slides with `tifffile`:
H&E and AB-PAS stains, 11999 × 12321 to 73744 × 45156 px, all 20×. The
"sample" column refers to the 26046 × 36880 px one. Fields marked
*unverified* were not present in any of them.

**File header**

| Offset | Type | Field |
|---|---|---|
| `0x00` | 4 bytes | `49 49 2A 00` (`II*\0`): little-endian **classic** TIFF, or `49 49 2B 00` (`II+\0`): **BigTIFF**, which KFBio uses for large slides (the 73744 × 45156 px H&E) |
| `0x04` | uint32 | Offset of the first IFD |

**Pages (IFDs)**, in file order:

| Page | Content | Size (sample) | Encoding |
|---|---|---|---|
| 0 | Level 0 (full resolution) | 26046 × 36880 | tiled 256 × 256, JPEG, YCbCr |
| 1 | Thumbnail | 820 × 1164 | stripped, JPEG |
| 2 – 7 | Levels 1 – 6, each about ½ the previous (see *Level sizes*) | 13023 × 18440 … 408 × 580 | tiled 256 × 256, JPEG, YCbCr |
| 8 | Slide **label** (`subfiletype = 1`) | 1164 × 968 | LZW, RGB |
| 9 | **Macro** overview (`subfiletype = 9`) | 2940 × 968 | JPEG |

The converter reads only the pyramid (pages 0, 2–7; it uses `tifffile`'s
`series[0].levels`). The thumbnail, label and macro are never exported, so
label text (which may identify a patient) does not end up in the output.

**Tiles**: each tile is a self-contained baseline JFIF JPEG (`SOF0`, 3
components, **YCbCr 4:2:2**: luma sampling 2 × 1, chroma 1 × 1). There is
no Adobe `APP14` marker and no shared `JPEGTables` tag (abbreviated tiles with
`JPEGTables` are still supported, *unverified on KFBio files*). The TIFF
`YCbCrSubSampling` tag is absent, although the JPEG streams are 4:2:2, so the
converter reads the subsampling from each level's JPEG `SOF` header instead.
Tiles stored with photometric RGB (the "Aperio RGB JPEG" variant) are decoded
as RGB, not YCbCr (*unverified on KFBio files*). Zero-length tiles are
treated as **white** (blank glass); none of the samples had any.

**Level sizes**: each level is downsampled exactly 2× (no stretching: checked
by comparing each level with the one above it, averaged 2 × 2), but the
*stored* size isn't always an exact half, and it varies between slides:

| Slide | Example | What KFBio stored |
|---|---|---|
| 26046 × 36880 | half of 3255 × 4610 is 1628 × 2305 | **1632 × 2320**: padded with white |
| 11999 × 12321 | half of 11999 × 12321 is 5999.5 × 6160.5 | 5999 × 6160: rounded down (fine) |
| 73744 × 45156 | half of 4609 × 2822 is 2304.5 × 1411 | **2303** × 1411: 1 – 2 px short |

Readers like `tifffile` only treat a level as part of the pyramid if it is
exactly ½ (rounded either way), and the converter's own geometry is exact too.
So the converter declares the exact size: larger stored levels are cropped
(only padding is removed), and levels up to 2 % (or 2 px) short are
extended. The added pixels read as white. A tile copy keeps the scanner's own
JPEG tile, whose encoder padding repeats the edge pixel. `--info` marks
these levels, e.g. `(stored 2303x1411, extended with white)`.

**Metadata**: every page's `ImageDescription` reads like
`Aperio Format Scanned By KFBIO | 26046x36879 [0,0 26046x36879] (256x256), JPEG/RGB Q=30|AppMag = 20|MPP = 0.485896|...`.

- `MPP` is the **level-0** pixel size (µm/px) and is repeated unchanged on every
  level; `AppMag` is the objective magnification.
- There are **no** `XResolution` / `ResolutionUnit` tags, so readers that
  don't parse Aperio descriptions don't know the pixel size. The converted
  BigTIFF stores it in standard resolution tags.
- The sizes in the description don't match the actual pages (e.g. 36879 vs
  36880 rows; level 1 says 13000 × 17920 but is 13023 × 18440). The converter
  trusts only the TIFF tags.
- `Q=30` is the scanner's JPEG quality: the source is already lossy, which is why
  the converter copies the original tiles instead of re-encoding them (see
  *Brightfield output (`.svs`)* above).

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

**Brightfield and automatic reader choice**

```python
from kfbf_converter import open_slide, SVSSlide, convert_slide

# open_slide picks the reader by magic bytes first, file extension second:
# KFBioSlide for .kfbf, SVSSlide for .svs. Unknown files raise SlideFormatError.
with open_slide("slide.svs") as slide:          # -> SVSSlide
    print(slide.summary())                       # same text as --info
    print(slide.kind, slide.is_rgb, slide.n_channels, slide.mpp)  # brightfield True 3 0.485896

    # One decoded source tile -> (256, 256, 3) uint8 RGB (white if missing)
    tile = slide.get_tile(level=0, col=10, row=4)

    # Any region at any level -> (height, width, 3) uint8 RGB (white outside the image)
    rgb = slide.read_region(level=2, x=0, y=0, width=1024, height=1024)

# Brightfield always writes <stem>.tif; mode is ignored. Leave compression and
# tile_size unset to copy the original JPEG tiles, or set them to re-encode.
convert_slide("slide.svs", "output/")
convert_slide("slide.svs", "output/", compression="jpeg")
```

`convert_slide` also takes `magnification=None` (overrides the slide's
objective magnification, like `--magnification`), `preflight=True` (set it
to `False` to skip the size check) and `errors=None`. Readers and
`open_slide` accept `magnification=` too; `slide.magnification` is the
value used (0 = unknown) and `slide.magnification_source` says where it came
from (`header`, `Aperio AppMag`, `--magnification` or `missing`). By default a refused or unverifiable output raises
`DataTooBigError` / `VerificationError`. If you pass a list, those failures are
printed and appended as `(output_path, exception)`, and the remaining outputs
are still written. `DataTooBigError.lines()` returns the three user-facing
lines (headline, reason, fix).

### How the size estimate works

For every planned output:

- **uncompressed** = sum over exported levels of
  *tiles × tile_size² × samples × planes* (OME: planes = Z × channels; RGB and
  brightfield: 3 samples).
- **estimate** = uncompressed × compression factor, plus TIFF overhead
  (16 bytes per tile for the BigTIFF offset/byte-count tables, about 1 KB per
  IFD, and OME-XML). For a `.svs` tile copy, the estimate is the exact sum of
  the copied tile bytes plus that overhead.
- The factors were calibrated on real KFBio slides (worst compressed /
  uncompressed ratios seen: deflate 0.26, lzw 0.30, jpeg q90 0.073), with about
  1.5× headroom, because tissue-dense slides compress worse than those samples:
  none 1.0, deflate 0.40, lzw 0.45, jpeg 0.12 (rising towards 0.35 at quality
  100).

The free-space margin is 256 MiB + 5 % of the estimate. The FAT limit is
4 GiB − 1 per file. The file system is detected with `GetVolumeInformationW`
on Windows and `/proc/self/mounts` on Linux; on macOS it isn't detected, but
the free-space check still runs.

</details>
