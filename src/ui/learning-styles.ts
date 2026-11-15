export const learningStyles = `
.molui-learning { position:absolute; top:12px; right:12px; z-index:12; max-width:calc(100% - 24px); font:13px/1.4 system-ui,sans-serif; color:#171733; pointer-events:auto; }
.molui-learning [hidden] { display:none !important; }
.molui-learning button,.molui-learning input { font:inherit; color:inherit; border:1px solid #d9d9ed; border-radius:7px; padding:6px 9px; background:rgba(255,255,255,.85); }
.molui-learning button { cursor:pointer; }
.molui-learning button:disabled { opacity:.5; cursor:default; }
.molui-learning :is(button,input,[tabindex]):focus-visible { outline:2px solid #5454c4; outline-offset:2px; }
.molui-learning-launch { width:36px; height:36px; padding:6px !important; box-shadow:0 2px 8px #0001; }
.molui-learning-launch svg { width:22px; height:22px; }
.molui-learning-dialog { width:244px; max-width:100%; box-sizing:border-box; border:1px solid #d9d9ed; border-radius:9px; background:rgba(255,255,255,.82); backdrop-filter:blur(8px); box-shadow:0 3px 12px #0001; overflow:hidden; }
.molui-learning-header { display:flex; align-items:center; gap:8px; padding:4px 6px; background:rgba(0,0,102,.78); color:white; }
.molui-learning-drag { flex:1; min-width:0; cursor:move; touch-action:none; font-weight:600; overflow-wrap:anywhere; }
.molui-learning-drag::before { content:'⠿'; margin-right:6px; opacity:.7; }
.molui-learning-header button { border:0; background:transparent; color:white; padding:0 4px; font-size:20px; line-height:24px; }
.molui-learning-body { padding:8px; }
.molui-learning-row { display:flex; gap:6px; align-items:center; flex-wrap:wrap; }
.molui-learning-connect { align-items:flex-end; }
.molui-learning-connect > label { display:flex; flex-direction:column; gap:4px; flex:1 1 120px; min-width:0; margin:0; }
.molui-learning-connect input { box-sizing:border-box; width:100%; min-width:0; }
.molui-learning-connect > button { flex:0 0 auto; }
.molui-learning-stats { font-size:12px; display:flex; flex-direction:column; gap:2px; }
.molui-learning-error { color:#982e37; margin:6px 0 0; max-width:260px; }
`;
