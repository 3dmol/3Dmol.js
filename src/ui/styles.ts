/** Scoped styles for the optional viewer UI. No host-page CSS is required. */
export const uiStyles: string = `
.molui {
  --molui-ink: #191933;
  --molui-muted: #595974;
  --molui-border: #dcdced;
  --molui-accent: #000066;
  --molui-paper: rgb(252 252 255 / 96%);
  position: absolute;
  inset: 0;
  z-index: 10;
  pointer-events: none;
  color: var(--molui-ink);
  font: 14px/1.45 system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
  text-align: left;
  color-scheme: light;
}
.molui *, .molui *::before, .molui *::after { box-sizing: border-box; font-family: inherit; }
.molui form { margin: 0; }
.molui [hidden] { display: none !important; }
.molui button, .molui input, .molui select { font: inherit; color: inherit; }
.molui button, .molui summary, .molui input, .molui select { touch-action: manipulation; }
.molui button {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  gap: 6px;
  min-height: 34px;
  padding: 6px 10px;
  margin: 0;
  border: 1px solid var(--molui-border);
  border-radius: 8px;
  background: #fff;
  color: var(--molui-ink);
  font-weight: 600;
  line-height: 1.25;
  text-transform: none;
  letter-spacing: normal;
  white-space: normal;
  cursor: pointer;
}
.molui button:hover:not(:disabled) { background: #efeffb; border-color: #b7b7df; }
.molui button:active:not(:disabled) { background: #e2e2f6; }
.molui button:disabled { opacity: .45; cursor: default; }
.molui button.molui-primary { background: var(--molui-accent); color: #fff; border-color: var(--molui-accent); }
.molui button.molui-primary:hover:not(:disabled) { background: #202080; border-color: #202080; }
.molui button.molui-danger { color: #a02d36; }
.molui button.molui-danger:hover:not(:disabled) { background: #fff0ef; border-color: #efc0c3; }
.molui button.molui-icon { width: 36px; height: 36px; padding: 0; font-size: 21px; }
.molui :is(button, input, select, summary, a):focus-visible {
  outline: 3px solid #5454c4;
  outline-offset: 3px;
}
.molui-toggle, .molui-trajectory-toggle {
  position: absolute;
  top: 12px;
  left: 12px;
  pointer-events: auto;
  box-shadow: 0 3px 14px #0000661a;
}
.molui-trajectory-toggle { left: 56px; }
.molui-no-edit .molui-trajectory-toggle { left: 12px; }
.molui-toggle svg, .molui-trajectory-toggle svg { display: block; width: 20px; height: 20px; pointer-events: none; }
.molui-panel {
  position: absolute;
  top: 56px;
  left: 12px;
  bottom: 12px;
  width: 320px;
  max-width: calc(100% - 24px);
  display: flex;
  flex-direction: column;
  min-height: 0;
  overflow: hidden;
  pointer-events: auto;
  background: var(--molui-paper);
  backdrop-filter: blur(4px);
  -webkit-backdrop-filter: blur(4px);
  border: 1px solid var(--molui-border);
  border-radius: 12px;
  box-shadow: 0 12px 36px #0000661f, 0 2px 8px #0000660d;
}
.molui-closed .molui-panel { display: none; }
/* Leave room for the scrollbar beside the floating close control. */
.molui-header { position: absolute; top: 4px; right: 12px; z-index: 2; }
.molui .molui-header > button.molui-icon {
  width: 32px;
  height: 32px;
  min-height: 32px;
  border-color: #e8e8f2;
  background: #fcfcffe8;
  color: var(--molui-muted);
  font-size: 18px;
  font-weight: 400;
}
.molui .molui-header > button:hover { background: #efeffb; color: var(--molui-accent); }
.molui-header > button:focus-visible { outline-offset: -3px; }
.molui-body {
  flex: 1 1 auto;
  min-height: 0;
  overflow: auto;
  overscroll-behavior: contain;
  scrollbar-width: thin;
  scrollbar-color: #babadb transparent;
}
.molui-section { margin: 0; padding: 10px 12px; border-bottom: 1px solid var(--molui-border); }
.molui-section:last-child { border-bottom: 0; }
.molui-section .molui-section { margin-top: 10px; padding: 8px 0 0; border-top: 1px solid var(--molui-border); border-bottom: 0; }
.molui-section > summary {
  min-height: 34px;
  display: flex;
  align-items: center;
  gap: 8px;
  margin: -4px 0;
  padding: 4px 0;
  color: var(--molui-ink);
  font-weight: 700;
  list-style: none;
  cursor: pointer;
}
.molui-section > summary::-webkit-details-marker { display: none; }
.molui-section > summary::before { content: '›'; width: 12px; font-size: 21px; line-height: 1; color: var(--molui-accent); }
.molui-section[open] > summary { margin-bottom: 8px; }
.molui-body > .molui-section:first-of-type > summary { padding-right: 40px; }
.molui-section[open] > summary::before { transform: rotate(90deg); }
.molui-section h3, .molui-section h4 { margin: 0 0 8px; font-size: 14px; }
.molui-field { display: flex; flex-direction: column; gap: 4px; min-width: 0; margin: 0 0 8px; }
.molui-filters { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); column-gap: 8px; }
.molui-selection-editor { min-width: 0; }
.molui-multiselect { min-width: 0; margin: 0 0 8px; }
.molui-multiselect[open] { grid-column: 1 / -1; }
.molui-multiselect > summary { display: block; list-style: none; cursor: pointer; }
.molui-multiselect > summary::-webkit-details-marker { display: none; }
.molui-multiselect-label { display: block; margin-bottom: 4px; color: var(--molui-muted); font-size: 12px; font-weight: 600; }
.molui-multiselect-value {
  position: relative;
  display: block;
  min-height: 34px;
  padding: 6px 26px 6px 8px;
  border: 1px solid #c3c3dd;
  border-radius: 8px;
  background: #fff;
  line-height: 1.3;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}
.molui-multiselect-value::after { content: '⌄'; position: absolute; right: 8px; top: 6px; color: var(--molui-muted); }
.molui-multiselect[open] .molui-multiselect-value { border-color: #8888bd; }
.molui-multiselect[open] .molui-multiselect-value::after { content: '⌃'; }
.molui-multiselect-body { margin-top: 6px; padding: 8px; border: 1px solid var(--molui-border); border-radius: 8px; background: #fff; }
.molui .molui-multiselect-search { margin-bottom: 6px; }
.molui-multiselect-actions { display: flex; justify-content: flex-end; margin-bottom: 4px; }
.molui-multiselect-actions > button { min-height: 30px; padding: 4px 8px; font-size: 12px; }
.molui-multiselect-options { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 0 8px; max-height: 180px; overflow: auto; overscroll-behavior: contain; scrollbar-width: thin; }
.molui-multiselect-options > .molui-check { min-width: 0; padding: 2px 4px; border-radius: 5px; }
.molui-multiselect-options > .molui-check:hover { background: #efeffb; }
.molui-multiselect-options > .molui-check > span { min-width: 0; overflow-wrap: anywhere; }
.molui-multiselect-options > .molui-unavailable { color: var(--molui-muted); font-style: italic; }
.molui-residue-available { margin: 0 0 8px; color: var(--molui-muted); font-size: 12px; }
.molui-residue-chips { display: flex; flex-wrap: wrap; gap: 6px; margin-bottom: 8px; }
.molui-residue-chips:empty { display: none; }
.molui-residue-chip { display: inline-flex; max-width: 100%; border: 1px solid #c3c3dd; border-radius: 7px; background: #ededfb; color: var(--molui-accent); }
.molui-residue-chip > button { min-height: 30px; border: 0; border-radius: 0; background: transparent; color: inherit; font-size: 12px; font-variant-numeric: tabular-nums; }
.molui-residue-chip > button:first-child { min-width: 0; padding: 4px 8px; border-radius: 6px 0 0 6px; overflow-wrap: anywhere; }
.molui-residue-chip > button:last-child { flex: 0 0 28px; min-width: 28px; padding: 0; border-left: 1px solid #c3c3dd; border-radius: 0 6px 6px 0; }
.molui-residue-endpoints { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 8px; }
.molui-residue-endpoints > .molui-field { margin: 0; }
.molui-residue-endpoints input[type=number] { font-variant-numeric: tabular-nums; }
.molui-residue-actions { display: flex; flex-wrap: wrap; align-items: center; gap: 6px; margin-top: 8px; }
.molui-residue-actions > button { flex: 0 0 auto; padding: 6px 8px; font-size: 12px; }
.molui-residue-actions > button:first-child { flex: 1 1 auto; }
.molui-residue-error { margin: 6px 0 0; color: #982e37; font-size: 12px; }
.molui-field > span, .molui-field > label { color: var(--molui-muted); font-size: 12px; font-weight: 600; }
.molui input:not([type=checkbox]):not([type=range]), .molui select {
  width: 100%;
  min-width: 0;
  min-height: 34px;
  padding: 6px 8px;
  margin: 0;
  border: 1px solid #c3c3dd;
  border-radius: 8px;
  background: #fff;
  line-height: 1.3;
}
.molui select { padding-right: 24px; cursor: pointer; }
.molui select:disabled, .molui input:disabled { background: #f2f2f8; color: #74748b; cursor: default; }
.molui input[type=color] { padding: 4px !important; cursor: pointer; }
.molui input[type=file] { overflow: hidden; font-size: 12px; }
.molui input[type=file]::file-selector-button {
  padding: 6px 9px;
  margin: 0 8px 0 0;
  border: 0;
  border-radius: 5px;
  background: #ebebfa;
  color: #000066;
  font: inherit;
  cursor: pointer;
}
.molui input::placeholder { color: #77778e; opacity: 1; }
.molui input[type=range] { width: 100%; min-width: 0; min-height: 26px; margin: 0; accent-color: var(--molui-accent); cursor: pointer; }
.molui input[type=checkbox] { width: 17px; height: 17px; flex: 0 0 auto; margin: 0; accent-color: var(--molui-accent); }
.molui-row, .molui-actions { display: flex; flex-wrap: wrap; align-items: center; gap: 8px; }
.molui-row > .molui-field { flex: 1 1 100px; }
.molui-actions { margin: 8px 0 0; }
.molui-actions > button { flex: 1 1 auto; }
.molui-check { display: inline-flex; align-items: center; gap: 6px; min-height: 34px; margin: 0; cursor: pointer; }
.molui-card {
  padding: 10px;
  margin: 8px 0 0;
  border: 1px solid var(--molui-border);
  border-radius: 10px;
  background: #fff;
  overflow-wrap: anywhere;
}
.molui-card-header { display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap; gap: 6px; margin: 0 0 6px; }
.molui-card-header > strong { flex: 1; min-width: 70px; font-size: 13px; }
.molui-card-header > button { padding: 5px 8px; font-size: 12px; }
.molui-card .molui-field:last-child { margin-bottom: 0; }
.molui-card > strong { display: block; margin-bottom: 6px; }
.molui-selection-header { display: flex; flex-wrap: nowrap; align-items: center; gap: 6px; min-width: 0; }
.molui-selection-visibility { display: inline-flex; flex: 0 0 32px; align-items: center; justify-content: center; min-height: 34px; margin: 0; border-radius: 7px; cursor: pointer; }
.molui-selection-visibility:hover { background: #efeffb; }
.molui .molui-selection-header input.molui-selection-name { flex: 1 1 0; width: 0; min-width: 0; padding: 6px; color: var(--molui-accent); font-weight: 700; }
.molui button.molui-selection-collapse { flex: 0 0 30px; min-width: 30px; min-height: 34px; padding: 0; border-color: transparent; color: var(--molui-muted); background: transparent; }
.molui-selection-collapse::before { content: ''; width: 7px; height: 7px; border: solid currentColor; border-width: 0 1.5px 1.5px 0; transform: translateY(-2px) rotate(45deg); }
.molui-selection-collapse[aria-expanded=false]::before { transform: translateX(-2px) rotate(-45deg); }
.molui-selection-summary-row { display: flex; align-items: flex-start; gap: 8px; margin-top: 6px; }
.molui-selection-summary { flex: 1 1 auto; min-width: 0; margin: 0; color: var(--molui-muted); font-size: 12px; overflow-wrap: anywhere; }
.molui-selection-summary-row > button { flex: 0 0 auto; padding: 5px 8px; font-size: 12px; }
.molui-selection-create { margin-top: 14px; padding-top: 12px; border-top: 1px solid var(--molui-border); }
.molui-selection-create > h3 { color: var(--molui-accent); }
.molui-selection-edit { margin-top: 8px; padding: 8px; border: 1px solid #c3c3dd; border-radius: 8px; background: #f5f5fd; }
.molui-selection-edit > .molui-actions { margin-top: 8px; }
.molui-style { margin-top: 8px; padding: 8px; border: 1px solid var(--molui-border); border-radius: 8px; background: #f8f8fd; }
.molui-style-header { display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap; gap: 6px; margin-bottom: 8px; }
.molui-style-header > strong { flex: 1 1 auto; min-width: 0; color: var(--molui-accent); font-size: 13px; }
.molui-style-header > button { flex: 0 1 auto; min-height: 30px; padding: 5px 8px; font-size: 12px; }
.molui-add-representation { display: flex; flex-wrap: wrap; align-items: end; gap: 8px; margin-top: 8px; padding-top: 8px; border-top: 1px dashed #c3c3dd; }
.molui-add-representation > .molui-field { flex: 1 1 120px; margin: 0; }
.molui-add-representation > select { flex: 1 1 120px; width: 0; }
.molui-add-representation > button { flex: 0 0 auto; padding: 6px 12px; }
.molui-card > .molui-check, .molui-selection-content > .molui-check { margin-top: 4px; }
.molui-style-grid { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 8px; margin-bottom: 8px; }
.molui-style-grid > .molui-field { margin: 0; }
.molui-style-options { display: flex; flex-wrap: wrap; gap: 0 12px; margin-bottom: 4px; }
.molui-style .molui-geometry { margin: 8px 0 0 3px; padding: 6px 8px; border: 0; border-left: 2px solid #b7b7df; border-radius: 0 7px 7px 0; background: #eeeefa; }
.molui-style .molui-geometry > summary { min-height: 30px; margin: 0; font-size: 12px; font-weight: 600; }
.molui-style .molui-geometry[open] > summary { margin-bottom: 8px; }

.molui-list { display: flex; flex-direction: column; gap: 0; margin: 0; padding: 0; list-style: none; }
.molui-muted { color: var(--molui-muted); font-size: 12px; margin: 6px 0; }
.molui-surface-validation[data-error=true] { color: #982e37; }
.molui-count { display: inline-block; margin: 0 0 8px; padding: 3px 8px; color: #000066; background: #ededfb; border-radius: 20px; font-size: 12px; font-weight: 600; }
.molui-empty { padding: 10px; margin: 8px 0; border-radius: 9px; background: #f2f2fa; color: var(--molui-muted); font-size: 13px; }
.molui-status { min-height: 40px; padding: 8px 52px 8px 12px; margin: 0; border-top: 1px solid var(--molui-border); background: #efeffa; color: #292970; font-size: 12px; overflow-wrap: anywhere; }
.molui-status:empty { display: none; }
.molui-status[data-error=true] { color: #982e37; background: #fff1ef; border-top-color: #efd0cf; }
.molui-trajectory {
  position: absolute;
  z-index: 3;
  right: 12px;
  bottom: 12px;
  left: 348px;
  display: grid;
  grid-template-columns: minmax(128px, 1fr) 118px;
  align-items: center;
  gap: 4px 6px;
  max-width: 350px;
  max-height: calc(100% - 24px);
  overflow: auto;
  overscroll-behavior: contain;
  margin-left: auto;
  padding: 6px;
  border: 1px solid #dcdcedbb;
  border-radius: 10px;
  background: rgb(252 252 255 / 84%);
  backdrop-filter: blur(8px);
  -webkit-backdrop-filter: blur(8px);
  box-shadow: 0 7px 25px #0000661f;
  pointer-events: auto;
  cursor: grab;
  touch-action: none;
  user-select: none;
}
.molui-trajectory:focus-visible { outline: 3px solid #5454c4; outline-offset: 2px; }
.molui-trajectory-dragging { cursor: grabbing; }
.molui-trajectory-header { grid-row: 1; grid-column: 2; min-width: 0; display: flex; align-items: center; justify-content: flex-end; gap: 4px; }
.molui button.molui-trajectory-minimize { flex: 0 0 24px; width: 24px; height: 26px; min-width: 24px; min-height: 26px; padding: 0; border-color: transparent; background: transparent; color: var(--molui-muted); }
.molui-frame-count { min-width: 0; overflow: hidden; text-overflow: ellipsis; color: var(--molui-muted); font-size: 11px; font-variant-numeric: tabular-nums; white-space: nowrap; }
.molui-trajectory svg { display: block; width: 16px; height: 16px; flex: 0 0 auto; pointer-events: none; }
.molui-trajectory-transport { grid-row: 1; grid-column: 1; min-width: 0; display: flex; justify-content: center; gap: 1px; }
.molui .molui-trajectory-transport button { flex: 1 1 24px; min-width: 24px; min-height: 30px; padding: 4px 0; line-height: 1; }
.molui .molui-trajectory-transport button.molui-primary { flex: 1.2 1 28px; min-width: 28px; }
.molui .molui-frame-slider { grid-row: 2; grid-column: 1; min-width: 0; min-height: 28px; }
.molui-trajectory-options { grid-row: 2; grid-column: 2; min-width: 0; display: flex; align-items: center; justify-content: space-between; gap: 4px; font-size: 11px; }
.molui-trajectory-options > label { display: inline-flex; align-items: center; gap: 4px; }
.molui-trajectory-options > label:first-child { position: relative; flex: 0 0 auto; min-width: 66px; font-size: 0; gap: 0; }
.molui-trajectory-options > label:first-child::after { content: ''; position: absolute; right: 8px; top: calc(50% - 4px); width: 6px; height: 6px; border: solid var(--molui-muted); border-width: 0 1.5px 1.5px 0; transform: rotate(45deg); pointer-events: none; }
.molui-trajectory-options select { appearance: none; -webkit-appearance: none; width: 66px; min-width: 66px; min-height: 28px; padding: 4px 18px 4px 6px; border-radius: 6px; font-size: 11px; cursor: pointer; }
.molui-trajectory-options .molui-check { min-height: 28px; gap: 4px; }
.molui-trajectory-options input[type=checkbox] { width: 14px; height: 14px; }
.molui-trajectory-error { grid-row: 3; grid-column: 1 / -1; margin: 2px 0 0; color: #982e37; font-size: 12px; }
.molui-trajectory-floating .molui-trajectory { right: auto; bottom: auto; width: 350px; max-width: calc(100% - 16px); margin: 0; }
.molui-closed .molui-trajectory { left: 12px; }
.molui-compact .molui-panel { right: 12px; width: auto; max-width: 340px; }
.molui-compact.molui-has-trajectory:not(.molui-trajectory-floating) .molui-panel { bottom: var(--molui-trajectory-space, 156px); }
.molui-compact .molui-trajectory { left: 12px; }
.molui-small .molui-panel { top: 56px; bottom: 8px; }
.molui-small.molui-compact.molui-has-trajectory:not(.molui-trajectory-floating) .molui-panel { bottom: var(--molui-trajectory-space, 156px); }
.molui.molui-tiny:not(.molui-closed) .molui-panel { bottom: 8px; }
.molui-tiny:not(.molui-closed):not(.molui-trajectory-floating) .molui-trajectory { visibility: hidden; pointer-events: none; }
@media (pointer: coarse) {
  .molui button, .molui select, .molui input:not([type=checkbox]):not([type=range]),
  .molui .molui-style-header > button, .molui-section > summary, .molui-style .molui-geometry > summary,
  .molui .molui-check, .molui-selection-visibility { min-height: 40px; }
  .molui-multiselect-value { min-height: 40px; padding-top: 9px; padding-bottom: 9px; }
  .molui-multiselect-value::after { top: 9px; }
  .molui-residue-chip > button:last-child { flex-basis: 40px; min-width: 40px; }
  .molui-selection-visibility { flex-basis: 40px; }
  .molui button.molui-selection-collapse { flex-basis: 40px; min-width: 40px; min-height: 40px; }
  .molui button.molui-icon { width: 40px; height: 40px; }
  .molui-trajectory-toggle { left: 60px; }
  .molui button.molui-trajectory-minimize { flex-basis: 40px; width: 40px; height: 40px; min-width: 40px; min-height: 40px; }
  .molui .molui-header > button.molui-icon { width: 40px; height: 40px; min-height: 40px; }
  .molui-body > .molui-section:first-of-type > summary { padding-right: 48px; }
  .molui-status { min-height: 48px; padding-right: 60px; }
  .molui input[type=range] { min-height: 32px; }
  .molui-trajectory { grid-template-columns: minmax(0, 1fr) 118px; }
  .molui-trajectory-header { grid-row: 1; grid-column: 1; justify-content: flex-start; min-height: 40px; }
  .molui-trajectory-options { grid-row: 1; grid-column: 2; }
  .molui-trajectory-transport { grid-row: 2; grid-column: 1 / -1; }
  .molui .molui-trajectory-transport button { min-width: 40px; min-height: 40px; }
  .molui .molui-trajectory-transport button.molui-primary { min-width: 40px; }
  .molui .molui-frame-slider { grid-row: 3; grid-column: 1 / -1; min-height: 32px; }
  .molui-trajectory-error { grid-row: 4; }
  .molui-panel, .molui-small .molui-panel { top: 60px; }
}
@media (prefers-reduced-motion: no-preference) {
  .molui button { transition: background-color .12s, border-color .12s; }
}
@media (forced-colors: active) {
  .molui button, .molui-panel, .molui-trajectory { border-color: ButtonText; }
}
`;
