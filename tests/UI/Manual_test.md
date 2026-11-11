# Viewer UI checks

Build with `npm run build:dev`, serve the repository, and open
`tests/webpages/3dmol-ui.html`. The demo has a protein viewer and an independent
trajectory viewer and uses local data. Test with a keyboard, mouse, and touch
emulation. The optional UI bundle includes its own scoped styles and does not
require jQuery or an external stylesheet.

## Attaching the controls

Load `build/3Dmol.js` followed by `build/3Dmol.ui.js`:

```js
const viewer = $3Dmol.createViewer('viewer', { backgroundColor: 'white' });
const controls = new $3Dmol.StateManager(viewer, { open: true });
// Existing data-ui="true" embeddings continue to work.
// Enable just playback controls: new $3Dmol.StateManager(viewer, { edit: false }).
// Change availability without rebuilding the controls:
controls.setControls({ edit: false, trajectory: true });
controls.setControls({ edit: true });
// Minimize or restore playback without changing its settings:
controls.ui.setTrajectoryMinimized(true);
controls.ui.setTrajectoryMinimized(false);
// Programmatic model changes are detected automatically; refresh immediately:
controls.refresh();
// Repeated initialization does not create another panel.
viewer.ui.initiateUI();
// Release DOM, observers, listeners and playback timers when finished:
controls.destroy();
```

The implementation is in `src/ui/*.ts`; `ui.ts`, `state.ts`, `trajectory.ts`, and
`styles.ts` separate the controls, appearance state, playback, and styling.
`UI.Form`, `UI.Icons`, the old generated-form specification dictionaries, and
internal `ui.tools` floating panels have been replaced by native controls.
StateManager selection/style methods remain available, with typed Map records.
The surface API is now `addSurface(type, style, selection, name?, calculationSpec?)`, returning a
promise for the UI surface ID. `loadSurface` registers an existing viewer surface.
`calculationSpec` defaults to the displayed selection; `{}` includes all atoms.

## Load and models

- Start empty: selections and surface creation are disabled; the load form is open.
- Load PDB/CIF/SDF/XYZ locally; check filename, atom count, fit, and default style.
- Load a second model: the first model and its appearance remain intact.
- Try an empty file, wrong format, bad PDB ID, and failed URL: show a readable error
  without removing the current scene. A loading operation cannot be submitted twice.
- Test local `.gz` files and an explicit format for extensionless URLs.
- Test PDB IDs and PubChem CIDs when online. URL servers must allow CORS.
- Toggle Show model, focus each model, and remove the first of two models. Filters
  and playback should remain usable with the gap in model IDs.

## Selections and appearance

- Model → atom group → chain → residue name → residue number → element choices
  come from the current scene. Chains, residue names and elements accept multiple
  checked values. Search long lists; no Ctrl/Command key is needed.
  Values within one field are alternatives; different fields narrow each other.
  Existing choices remain selected if another filter makes them unavailable, with
  zero matches shown rather than silently changing the selection.
- Residue numbers use From/To bounds and endpoint sliders. Add a range, or use equal
  bounds for one residue; add more ranges for disjoint regions. Overlapping ranges
  merge. Edit/remove selected ranges, and use All to remove the residue constraint.
  Test zero, negative residue numbers, sparse numbering and invalid/reversed bounds.
- Check blank chain IDs and residue 0. No-match filters show zero and disable Add.
- Cartoon is offered for selections with protein/nucleic-acid backbone atoms.
- Choose Representation → None for a selection used only for surfaces or residue
  labels. It leaves existing atom appearances intact. Choosing None on an existing
  representation removes that representation; the empty card can add one again.
- Rename a selection directly beside its visibility checkbox; names update in the
  surface selector. Enter commits, blank names retain the previous name.
- Saved selections show their criteria and atom count. Edit opens the same filters
  in the card; Cancel discards changes, and Save selection updates the matching atoms
  while retaining the name, representations, visibility and residue-label setting.
  The New selection form follows the existing selection cards.
- Collapse individual selections with the chevron beside the name. The name,
  visibility checkbox, and summary remain visible; Edit and settings are hidden.
  Appearances do not change. Collapsed cards stay collapsed across refreshes, and
  expanding the chevron restores any unfinished edits. Keyboard Enter/Space also
  toggles the chevron.
- Add two overlapping selections, then change representation, colors, bond radius,
  and opacity. The later selection takes precedence on overlapping atoms.
- Use Add beside New representation to combine cartoons, sticks, spheres, lines or crosses in
  one selection. Each representation can be edited or removed independently;
  already-used types are omitted so two editors cannot overwrite the same style.
- Every representation has Settings containing opacity and its applicable sizes.
  Lines have opacity, sticks also have bond radius, and more detailed styles expose
  their corresponding controls. Expanding Settings persists while editing.
- Under Settings, adjust sphere and cross radii independently. Atomic radii
  use each element's van der Waals radius and a scale factor; fixed radii use Å.
  Switching to fixed radii must not retain a hidden scale multiplier. Ball & stick
  exposes separate bond and sphere sizing.
- For cartoons, try Rectangle, Oval, Parabola, Edged and Backbone trace. Check
  thickness, ribbon, helix cylinders and strand arrows. Trace should offer only
  its applicable size options; ribbon/custom width should hide arrow options.
  Editing a size keeps Settings open and preserves colors and other styles.
- Each representation has its own outlined block; Settings is visibly nested
  inside that block rather than looking like a separate selection section.
- The checkbox beside a selection name toggles all its representations and residue
  labels together. Removing/disabling a selection restores the previous appearance,
  including on every trajectory frame.
- Editing a control keeps keyboard focus on that control. Imported named colors
  retain their color when opacity changes.
- Toggle Residue labels on each selection. Repeated toggles must not accumulate
  labels or remove labels created by the host application.
- Under Color by, choose Color scheme and a built-in scheme such as greenCarbon,
  orangeCarbon or chainHetatm. Changing opacity must preserve the scheme.
- Add/remove a surface for all atoms or a selection. Choose element, chain, solid
  color, or a named color scheme before creation or on the resulting surface card.
  Recoloring must preserve opacity and reuse the existing geometry.
- Show surface for selects the displayed region. Calculate using selects the atoms
  that determine its shape: Same as shown, All atoms, or a saved selection (including
  those with representation None). Calculation atoms must include every shown atom;
  invalid choices explain why Add surface is disabled. Test a chain shown within a
  larger complex, renaming both choices, and removing a model used only for calculation.
- Opacity preview is readable;
  the Add surface button stays disabled while calculation is pending. Surfaces are
  snapshots at the frame where they were generated, not animated meshes.

## Trajectory

- Check “Load multiple models as trajectory frames” when loading multiframe PDB,
  XYZ, GRO, or LAMMPS files. Uncheck it to use the first model only.
- A model with one frame has no playback bar; multiple frames reveal it without
  starting playback. Frames added later with `setCoordinates` also reveal it.
- Minimize playback to the icon in the upper left and restore it. Playback,
  current frame, speed, Loop, and a dragged position are retained. With both
  controls enabled the icon sits beside the hamburger; with only trajectory
  controls enabled it takes the hamburger's position. Removing the trajectory
  hides both the bar and icon, including after minimizing.
- Attach with `{ edit: false, trajectory: true }`, `{ edit: true, trajectory: false }`,
  and neither option: only the requested controls should be available, with both
  enabled by default. Change availability using `controls.setControls(...)`.
  Repeat for embeddings with `data-ui="true"`, `data-ui="true" data-ui-edit="false"`,
  `data-ui="true" data-ui-trajectory="false"`, and just `data-ui-trajectory="true"`.
- Test Play/Pause, first/previous/next/last, slider arrow keys and dragging, speed,
  and loop. Displayed frame numbers are one-based. Playback affects all models;
  shorter trajectories hold their final frame.
- Play/Pause use icons with accessible names. The compact bar has a translucent
  background, while controls remain readable over light and dark structures.
- Drag anywhere on the widget background, frame counter, or gaps with mouse or
  touch to reposition it. Moving it must not rotate the molecule or change frames.
  Focus the widget itself for arrow-key movement; Shift moves farther, and Home or
  a double-click on the background restores its docked position. Resize the viewer
  and confirm the bar stays reachable. Buttons, sliders, speed and Loop controls
  (including labels and keyboard interactions) must remain independent of dragging.
  Check that the full FPS text is visible beside the speed dropdown arrow.
- Turn Loop off: stop at the final frame. Play again starts at the beginning.
- Rapidly scrub a slow remote trajectory: frame requests are serialized and the
  final requested frame wins. A failure pauses playback and shows an error; seeking
  again retries. Hiding the tab or removing the viewer stops UI playback.
- Collapse the sidebar: playback stays reachable. In viewers shorter than 320px,
  close the sidebar to reveal docked playback and recover space for settings;
  a manually positioned bar remains reachable. Restoring minimized playback in
  these short viewers closes the sidebar so the playback bar is visible.

## Responsive and accessible behavior

- Check embedded sizes 1100×760, 640×400, 360×640, 320×480, and 360×280.
- Sidebar contents scroll internally without horizontal overflow. The canvas
  remains interactive outside the controls, and controls do not rotate the molecule.
- The subtle close button overlays the corner without a header band. Check it at
  the top and bottom of the sidebar; section titles and status text remain readable.
  Touch controls remain comfortably sized despite the denser desktop spacing.
- In compact layouts, settings and playback do not overlap. Closing settings or
  pressing Escape returns focus to the hamburger button. Tab navigation has a visible focus ring.
- Two viewers have independent controls, unique panel IDs, styles and playback.
- Fit all, PNG, background and projection work. There is no VRML export button. Host labels/shapes survive
  appearance edits. Remove/recreate controls: no duplicate UI or lingering timers.

## Automated checks

```sh
npm run build:dev
npx tsc --noEmit --skipLibCheck
npx jest tests/UI --runInBand
node tests/UI/browser-smoke.cjs
```

`skipLibCheck` bypasses an existing missing `BinaryData` type in netcdfjs's
installed declaration file. The DOM suites avoid native canvas and exercise real
model selection and styling. The browser smoke uses Puppeteer (installed with
glcheck) and a real WebGL viewer; set `UI_SCREENSHOTS=1` to save screenshots in
`/tmp`. Run `npm run build:prod` to verify both production bundles.
