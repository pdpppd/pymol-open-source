# PyMOL Open Source

PyMOL is a molecular visualization and analysis application for working with
protein structures, ligands, maps, trajectories, and publication-quality
rendering. This repository contains the open-source PyMOL codebase, including
the classic command-driven workflow, the Qt desktop interface, and newer
interactive features such as the AI Assistant and labeled grid mode.

## What Is In This Repo

- Core PyMOL rendering, selection, and molecular graphics engine
- Qt and Tk desktop interfaces
- Python command layer and scripting support
- Structure loading, styling, alignment, measurement, and export workflows
- Optional AI Assistant inside the Qt application
- Grid mode for multi-object and multi-state comparisons, now with labels

## AI Assistant

The Qt application includes an `AI -> AI Assistant...` panel that runs an
OpenAI-powered agent inside PyMOL.

### What it can do

- Execute PyMOL commands on your behalf
- Inspect the current scene and loaded objects
- Capture the current viewport and use it as visual context
- Search the Protein Data Bank before fetching structures
- Look up built-in PyMOL command documentation when it needs command syntax

### How it works

The assistant uses the OpenAI Responses API and an internal tool loop. The
model receives your chat request, can call local tools, reads the tool output,
and continues until it has a final answer.

Available tools include:

- `run_command`: runs a single PyMOL command and returns text output
- `view_screen`: captures the current viewport and sends it back as an image
- `get_scene_info`: summarizes loaded objects, chains, and residue counts
- `search_pdb`: searches the RCSB PDB and returns matching entries
- `docs_grep` and `docs_view`: search the local PyMOL command reference text

The assistant is guided to inspect the current scene first, check command
results before proceeding, and verify visual changes by capturing the screen
again after edits.

### AI requirements

The AI Assistant is optional. To use it you need:

- A Qt build of PyMOL
- The `openai` Python package
- An `OPENAI_API_KEY` environment variable

Example:

```bash
pip install openai
export OPENAI_API_KEY=your_api_key_here
```

When PyMOL starts, open the assistant from the `AI` menu in the Qt interface.
The panel lets you send prompts, cancel a running request, clear history, and
change reasoning effort between `high`, `medium`, `low`, and `none`.

## Grid Mode Labels

PyMOL grid mode is useful for comparing multiple objects or states in a tiled
layout. In this version, object grid mode includes visible labels in each grid
cell so you can immediately tell which object you are looking at.

### What changed

- Each grid cell now shows the slot number and object name
- Labels are drawn as a 2D overlay in the top-left corner of the cell
- Label size is controlled by the global `grid_label_size` setting

This is especially useful when comparing:

- multiple docking poses
- homologous structures
- conformational ensembles
- design variants
- ligand-bound versus apo forms

### Using grid mode

From the GUI, select one of the grid layouts such as:

- `By Object`
- `By State`
- `By Object-State`

From the PyMOL command line:

```pml
set grid_mode, 1
set grid_label_size, 14
```

`grid_mode = 1` enables grid mode by object, which is the mode where the new
object labels are shown.

## Installation

See also: [INSTALL](INSTALL)

### Requirements

Build/install requirements from this repository:

- C++17 compiler
- CMake 3.13+
- Python 3.9+
- `pip`
- OpenGL
- GLEW
- `libpng`
- `freetype`

Optional components:

- PyQt5, PyQt6, PySide2, or PySide6 for the Qt interface
- GLUT/freeglut for the legacy GUI
- `libxml2` for COLLADA export
- `msgpack-c` and `mmtf-cpp` for fast MMTF support
- `catch2` for tests
- `openvr` for VR support
- `libnetcdf` for VMD plugins
- `Pmw` for legacy GUI/plugins

### Install PyMOL

Standard installation:

```bash
pip install .
```

Developer-oriented installation with verbose build output and C++ tests:

```bash
pip install --verbose --no-build-isolation --config-settings testing=True .
```

### Run PyMOL

Depending on your environment, the launcher is typically installed at one of:

```bash
$PYMOL_PATH/bin/pymol
```

or

```bash
$PYTHONUSERBASE/bin/pymol
```

You can also usually launch it directly with:

```bash
pymol
```

### Install the optional AI dependencies

If you want the AI Assistant, install an OpenAI client library and set your API
key before launching the Qt app:

```bash
pip install openai
export OPENAI_API_KEY=your_api_key_here
```

If PyMOL starts without a supported Qt binding, it can fall back to the legacy
Tk interface, but the AI Assistant panel is part of the Qt application.

## Development Notes

- PyMOL can be used interactively, from Python scripts, or from the command line
- Most Python packages for the application live under `modules/`
- The Qt AI panel is implemented in `modules/pmg_qt/ai_panel.py`
- The AI tool loop is implemented in `modules/pmg_qt/ai_agent.py`

## License

Copyright (c) [Schrodinger, LLC](https://www.schrodinger.com/)

Published under a BSD-like license. See [LICENSE](LICENSE).
