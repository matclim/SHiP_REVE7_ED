# REve event display for a sampling calorimeter

A minimal, fast, web-based event display built on **REve** (ROOT 7's event
visualisation environment, the successor to TEve and the technology behind
CMS's FireworksWeb). Loads a GDML detector, keeps only the volumes you ask
for, and renders calorimeter hits as a single GPU-instanced box collection
per event.

## What it does

* Imports a GDML file via `TGeoManager::Import`.
* Walks the geometry tree and, for every node whose **volume name** matches
  any of a list of regular expressions, captures the local `TGeoShape` plus
  the cumulative global transform and wraps them in a `REveGeoShape`.
  Recursion stops at matched nodes, so you don't pay for sub-structure you
  don't want to see.
* For the example detector, the filter `^ECAL_GL\d+_` keeps the 20 ECAL
  layer envelopes (Lead absorbers + WidePVT and ThinPS scintillator layers)
  and skips the thousands of `HPL_*` hodoscope fibres entirely. That's the
  difference between a smooth display and a slideshow.
* Reads a `TTree` (default name `calo_events`) of per-event
  `std::vector<double>` branches `edep`, `x_global`, `y_global`, `z_global`
  (plus the categorical labels `type`, `section`, `layer`, `hcal`,
  `hexant`).
* Renders all hits in one event as a `REveBoxSet` with
  `kBT_AABoxFixedDim` — the most efficient mode, every box the same size,
  positions sent to the GPU as an instanced array. Colour is mapped from
  `edep` via `REveRGBAPalette`.

## Requirements

* A recent ROOT (≥ 6.30 recommended) configured with **`webgui=ON`**,
  **`root7=ON`**, **`gdml=ON`**. LCG releases ≥ LCG_103 already include
  these. Verify with `root-config --has-webgui --has-root7 --has-gdml`.
* CMake ≥ 3.16, a C++17 compiler.
* A modern browser (Chrome / Firefox / Safari) for the client side.

## Build

```bash
source /path/to/root/bin/thisroot.sh   # so CMake can find ROOTConfig.cmake
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## Run

```bash
./build/event_display detector.gdml events.root
# or with a non-default tree name:
./build/event_display detector.gdml events.root my_tree_name
```

A browser window opens at the local REve server URL
(something like `http://localhost:NNNN/?...`). You'll see the 20 ECAL
layer envelopes (Lead in grey, WidePVT in blue, ThinPS in orange) and the
hits of event 0.

The default browser UI gives you the standard REve viewer controls
(rotate, zoom, pan, transparency sliders, element tree on the left). For
event navigation you have three options:

1. The browser console: `eve_next()`, `eve_prev()`, `eve_goto(42)`.
2. From a ROOT prompt attached to the same process.
3. Wire up your own UI buttons via `REveManager`'s MIR ("Method
   Invocation Request") mechanism — see ROOT's `tutorials/eve7/`.

## Tweaking

* **Which volumes to draw.** Edit the `include` list in `main()`. Patterns
  are `std::regex` (ECMAScript flavour) and use `regex_search`, so a
  pattern matches anywhere in the volume name — the trailing hex pointer
  suffixes ROOT sometimes leaves on GDML names won't break things. To
  also draw, say, individual scintillator strips inside each layer,
  comment out the `return` after the matched-node block in
  `BuildGeoShapes` so recursion continues.
* **Colours.** The `Lead` / `WidePVT` / `ThinPS` mapping is in
  `BuildGeoShapes`. Add cases for whatever pattern you want.
* **Hit drawing.** `cell_xy` / `cell_z` in `GotoEvent` set the box size.
  For per-hit variable size use `kBT_AABox` instead of `kBT_AABoxFixedDim`
  and pass the dimensions to `AddBox`. For directional / cone-shaped
  primitives (e.g. shower axes), use `kBT_Cone`.
* **Clustering.** The hit-grouping placeholder in `GotoEvent` is just a
  comment for now. To draw real clusters, run your clustering algorithm
  on `(*fXg, *fYg, *fZg, *fEdep)`, then add either:
    * another `REveBoxSet` containing one large box per cluster centroid
      (sized e.g. by total energy), or
    * a `REvePointSet` of centroids with `SetMarkerStyle(20)` and per-point
      colour from energy, or
    * for "ellipsoidal" clusters, `REveEllipsoid` (one per cluster).
  Add them under `fEventHolder` so they're cleared on the next event.

## Why this is fast

The recipe is the same one CMS uses for FireworksWeb:

1. **Filter the geometry hard.** A few dozen `REveGeoShape` entries upload
   to the GPU once and stay there — the per-event cost is zero.
2. **Use instanced primitives for hits.** `REveBoxSet` with a fixed-dim
   box type sends positions as a packed array; rendering 10⁵ hits costs
   one draw call.
3. **Palette colouring on the GPU.** `REveRGBAPalette` does the
   value→RGBA mapping once; per-hit colour does not require per-hit
   uniforms.
4. **Clear-and-rebuild per event, leave geometry alone.** `fEventHolder`
   is wiped and refilled inside a `REveManager::ChangeGuard`; the
   `fGeoHolder` and its REveGeoShape children never change after init.

## Next steps

* Add 2D projected views: spawn `REveProjectionManager` instances with
  `kPT_RPhi` and `kPT_RhoZ` (or for a fixed-target geometry like this
  one, just spawn a second 3D viewer locked to a side-on camera).
* Cluster visualisation (see above).
* Connect a "play through events" button using REve's MIR mechanism.
* Persist user view configurations via `REveManager::SaveVizDB`.
