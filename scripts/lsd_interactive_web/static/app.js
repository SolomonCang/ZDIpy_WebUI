/**
 * app.js — LSD Interactive Fit frontend.
 *
 * Communicates with the FastAPI backend to load observations and compute
 * Voigt / Unno-Rachkovsky model profiles in real time.
 */

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------
let _velObs = [];
let _medianI = [];
let _individualI = [];
let _velModel = [];
let _currentModel = 'voigt';
let _params = {};
let _ranges = {};
let _debounceTimer = null;

// Slider definitions:  key → { label, unit, fmt_decimals, ur_only }
const SLIDER_META = {
  line_strength:          { label: 'd  (line depth)',              unit: '',      decimals: 3, ur_only: false },
  gauss_width_kms:        { label: 'σ_G  (Gauss width)',           unit: 'km/s',  decimals: 1, ur_only: false },
  lorentz_width_fraction: { label: 'a_L  (Lorentz fraction)',      unit: '',      decimals: 3, ur_only: false },
  v_center:               { label: 'v_c  (center velocity)',       unit: 'km/s',  decimals: 1, ur_only: false },
  inst_resolution:        { label: 'R  (spectral resolution)',     unit: '',      decimals: 0, ur_only: false },
  limb_darkening:         { label: 'ε  (limb darkening)',          unit: '',      decimals: 2, ur_only: false },
  beta:                   { label: 'β  (source function slope)',   unit: '',      decimals: 3, ur_only: true  },
  fI:                     { label: 'fI  (filling factor I)',       unit: '',      decimals: 2, ur_only: true  },
};

// ---------------------------------------------------------------------------
// Plotly config
// ---------------------------------------------------------------------------
const PLOTLY_LAYOUT = {
  paper_bgcolor: '#1e1e2e',
  plot_bgcolor: '#1e1e2e',
  font: { color: '#cdd6f4', family: 'Inter, Noto Sans, sans-serif', size: 13 },
  margin: { t: 40, r: 20, b: 50, l: 60 },
  xaxis: {
    title: 'Velocity (km/s)',
    gridcolor: '#313244',
    zerolinecolor: '#585b70',
    color: '#cdd6f4',
  },
  yaxis: {
    title: 'Stokes I',
    gridcolor: '#313244',
    zerolinecolor: '#585b70',
    color: '#cdd6f4',
  },
  legend: {
    bgcolor: '#313244',
    bordercolor: '#585b70',
    font: { color: '#cdd6f4', size: 11 },
    x: 0.01, y: 0.99,
  },
  hovermode: 'x unified',
};

const PLOTLY_CONFIG = {
  responsive: true,
  displaylogo: false,
  modeBarButtonsToRemove: ['select2d', 'lasso2d'],
};

// ---------------------------------------------------------------------------
// API helpers
// ---------------------------------------------------------------------------
async function apiGet(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${r.status}: ${r.statusText}`);
  return r.json();
}

async function apiPost(url, body) {
  const r = await fetch(url, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  });
  if (!r.ok) throw new Error(`${r.status}: ${r.statusText}`);
  return r.json();
}

// ---------------------------------------------------------------------------
// Slider DOM generation
// ---------------------------------------------------------------------------
function buildSliders() {
  const container = document.getElementById('sliders');
  container.innerHTML = '';
  for (const [key, meta] of Object.entries(SLIDER_META)) {
    const range = _ranges[key] || { min: 0, max: 1, step: 0.01 };
    const val = _params[key] ?? range.min;

    const row = document.createElement('div');
    row.className = 'slider-row' + (meta.ur_only ? ' ur-only' : '');
    row.dataset.key = key;

    const header = document.createElement('div');
    header.className = 'slider-header';

    const label = document.createElement('span');
    label.className = 'slider-label';
    label.textContent = meta.label + (meta.unit ? ` (${meta.unit})` : '');

    const valSpan = document.createElement('span');
    valSpan.className = 'slider-value';
    valSpan.id = `val-${key}`;
    valSpan.textContent = formatVal(val, meta.decimals);

    header.appendChild(label);
    header.appendChild(valSpan);

    const input = document.createElement('input');
    input.type = 'range';
    input.id = `slider-${key}`;
    input.min = range.min;
    input.max = range.max;
    input.step = range.step;
    input.value = clamp(val, range.min, range.max);
    input.addEventListener('input', () => onSliderInput(key));

    row.appendChild(header);
    row.appendChild(input);
    container.appendChild(row);
  }
  updateURVisibility();
}

function formatVal(v, decimals) {
  return Number(v).toFixed(decimals);
}

function clamp(v, lo, hi) {
  return Math.min(hi, Math.max(lo, v));
}

// ---------------------------------------------------------------------------
// UR visibility toggle
// ---------------------------------------------------------------------------
function updateURVisibility() {
  const isUR = _currentModel === 'unno';
  document.querySelectorAll('.slider-row.ur-only').forEach(el => {
    el.classList.toggle('visible', isUR);
  });
  // Adjust line_strength slider range for UR vs Voigt
  const lsSlider = document.getElementById('slider-line_strength');
  if (lsSlider && _ranges.line_strength) {
    const r = _ranges.line_strength;
    const newMax = isUR ? (r.max_unno || r.max) : r.max;
    const newStep = isUR ? (r.step_unno || r.step) : r.step;
    lsSlider.max = newMax;
    lsSlider.step = newStep;
  }
}

// ---------------------------------------------------------------------------
// Model selector
// ---------------------------------------------------------------------------
function initModelSelector() {
  document.querySelectorAll('.model-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      _currentModel = btn.dataset.model;
      document.querySelectorAll('.model-btn').forEach(b =>
        b.classList.toggle('active', b === btn));
      updateURVisibility();
      requestCompute();
    });
  });
}

function setActiveModelButton(model) {
  document.querySelectorAll('.model-btn').forEach(b =>
    b.classList.toggle('active', b.dataset.model === model));
}

// ---------------------------------------------------------------------------
// Slider change handler
// ---------------------------------------------------------------------------
function onSliderInput(key) {
  const slider = document.getElementById(`slider-${key}`);
  const meta = SLIDER_META[key];
  const val = parseFloat(slider.value);
  _params[key] = val;
  document.getElementById(`val-${key}`).textContent = formatVal(val, meta.decimals);
  requestCompute();
}

// ---------------------------------------------------------------------------
// Compute model (debounced)
// ---------------------------------------------------------------------------
function requestCompute() {
  if (_debounceTimer) clearTimeout(_debounceTimer);
  _debounceTimer = setTimeout(doCompute, 30);
}

async function doCompute() {
  try {
    const body = {
      model_type: _currentModel,
      line_strength: _params.line_strength,
      gauss_width_kms: _params.gauss_width_kms,
      lorentz_width_fraction: _params.lorentz_width_fraction,
      v_center: _params.v_center,
      inst_resolution: _params.inst_resolution,
      limb_darkening: _params.limb_darkening,
      beta: _params.beta,
      fI: _params.fI,
    };
    const data = await apiPost('/api/compute', body);
    updatePlotModel(data.I_model);
    updateStatus(data.rms);
    updateParamsPanel();
  } catch (e) {
    console.error('Compute error:', e);
  }
}

// ---------------------------------------------------------------------------
// Plotly chart
// ---------------------------------------------------------------------------
function initPlot() {
  // Individual observations (gray)
  const traces = [];
  for (let i = 0; i < _individualI.length; i++) {
    traces.push({
      x: _velObs, y: _individualI[i],
      type: 'scatter', mode: 'lines',
      line: { color: '#585b70', width: 0.8 },
      opacity: 0.5,
      showlegend: i === 0,
      name: i === 0 ? 'Individual phases' : '',
      hoverinfo: 'skip',
    });
  }

  // Median Stokes I
  traces.push({
    x: _velObs, y: _medianI,
    type: 'scatter', mode: 'lines',
    line: { color: '#89b4fa', width: 2 },
    name: 'Median Stokes I',
  });

  // Continuum reference
  traces.push({
    x: [_velObs[0], _velObs[_velObs.length - 1]],
    y: [1.0, 1.0],
    type: 'scatter', mode: 'lines',
    line: { color: '#585b70', width: 1, dash: 'dash' },
    showlegend: false,
    hoverinfo: 'skip',
  });

  // Model (to be updated)
  traces.push({
    x: _velModel, y: new Array(_velModel.length).fill(1.0),
    type: 'scatter', mode: 'lines',
    line: { color: '#f38ba8', width: 2.5 },
    name: 'Model',
  });

  Plotly.newPlot('plot-container', traces, PLOTLY_LAYOUT, PLOTLY_CONFIG);
}

function updatePlotModel(I_model) {
  const modelLabel = _currentModel === 'unno' ? 'UR Model' : 'Voigt Model';
  // Model is the last trace
  const numTraces = _individualI.length + 3;  // individual + median + continuum + model
  Plotly.restyle('plot-container', {
    y: [I_model],
    name: [modelLabel],
  }, [numTraces - 1]);
}

// ---------------------------------------------------------------------------
// Status bar
// ---------------------------------------------------------------------------
function updateStatus(rms) {
  const el = document.getElementById('status-bar');
  const tag = _currentModel === 'unno' ? 'UR' : 'Voigt';
  if (rms != null) {
    el.innerHTML = `<span class="model-tag">[${tag}]</span> RMS in ±3σ = <span class="rms">${rms.toFixed(5)}</span>`;
  } else {
    el.innerHTML = `<span class="model-tag">[${tag}]</span>`;
  }
}

// ---------------------------------------------------------------------------
// Params panel
// ---------------------------------------------------------------------------
function updateParamsPanel() {
  const lines = [];
  lines.push(`model_type: ${_currentModel}`);
  lines.push(`line_strength: ${formatVal(_params.line_strength, 5)}`);
  lines.push(`gauss_width_kms: ${formatVal(_params.gauss_width_kms, 3)}`);
  lines.push(`lorentz_width_fraction: ${formatVal(_params.lorentz_width_fraction, 5)}`);
  lines.push(`v_center: ${formatVal(_params.v_center, 3)} km/s`);
  lines.push(`inst_resolution: ${formatVal(_params.inst_resolution, 0)}`);
  lines.push(`limb_darkening: ${formatVal(_params.limb_darkening, 4)}`);
  if (_currentModel === 'unno') {
    lines.push(`beta: ${formatVal(_params.beta, 5)}`);
    lines.push(`fI: ${formatVal(_params.fI, 4)}`);
  }
  document.getElementById('params-panel').textContent = lines.join('\n');
}

// ---------------------------------------------------------------------------
// Buttons
// ---------------------------------------------------------------------------
function initButtons() {
  document.getElementById('btn-reset').addEventListener('click', resetParams);
  document.getElementById('btn-copy').addEventListener('click', copyParams);
  document.getElementById('btn-snapshot').addEventListener('click', saveSnapshot);
  document.getElementById('btn-auto-estimate').addEventListener('click', autoEstimate);
  document.getElementById('btn-save-config').addEventListener('click', saveToConfig);
}

async function saveToConfig() {
  const btn = document.getElementById('btn-save-config');
  btn.disabled = true;
  btn.textContent = 'Saving…';
  try {
    const body = {
      model_type: _currentModel,
      line_strength: _params.line_strength,
      gauss_width_kms: _params.gauss_width_kms,
      lorentz_width_fraction: _params.lorentz_width_fraction,
      limb_darkening: _params.limb_darkening,
      unno_beta: _params.beta,
      unno_filling_factor_I: _params.fI,
    };
    const data = await apiPost('/api/save_config', body);
    const el = document.getElementById('status-bar');
    const prev = el.innerHTML;
    el.innerHTML = '<span style="color: var(--green);">Saved to config.json ✓</span>';
    setTimeout(() => { el.innerHTML = prev; }, 2500);
  } catch (e) {
    console.error('Save config error:', e);
    const el = document.getElementById('status-bar');
    el.innerHTML = '<span style="color: var(--red);">Save failed: ' + e.message + '</span>';
  } finally {
    btn.disabled = false;
    btn.textContent = 'Save to Config';
  }
}

let _initParams = {};

function resetParams() {
  _params = { ..._initParams };
  // Reset slider DOM values
  for (const [key, meta] of Object.entries(SLIDER_META)) {
    const slider = document.getElementById(`slider-${key}`);
    if (slider) {
      slider.value = _params[key];
      document.getElementById(`val-${key}`).textContent = formatVal(_params[key], meta.decimals);
    }
  }
  requestCompute();
}

function copyParams() {
  const lines = [];
  lines.push(`"model_type": "${_currentModel}",`);
  lines.push(`"line_strength": ${formatVal(_params.line_strength, 5)},`);
  lines.push(`"gauss_width_kms": ${formatVal(_params.gauss_width_kms, 3)},`);
  lines.push(`"lorentz_width_fraction": ${formatVal(_params.lorentz_width_fraction, 5)},`);
  lines.push(`"limb_darkening": ${formatVal(_params.limb_darkening, 4)},`);
  if (_currentModel === 'unno') {
    lines.push(`"unno_beta": ${formatVal(_params.beta, 5)},`);
    lines.push(`"unno_filling_factor_I": ${formatVal(_params.fI, 4)},`);
  }
  const text = lines.join('\n');
  navigator.clipboard.writeText(text).then(() => {
    const el = document.getElementById('status-bar');
    const prev = el.innerHTML;
    el.innerHTML = '<span style="color: var(--green);">Copied to clipboard!</span>';
    setTimeout(() => { el.innerHTML = prev; }, 1500);
  });
}

function saveSnapshot() {
  Plotly.downloadImage('plot-container', {
    format: 'png',
    width: 1200,
    height: 700,
    filename: 'lsd_interactive_snapshot',
  });
}

async function autoEstimate() {
  const btn = document.getElementById('btn-auto-estimate');
  btn.disabled = true;
  btn.textContent = 'Estimating…';
  try {
    const body = {
      model_type: _currentModel,
      gauss_width_kms: _params.gauss_width_kms,
      lorentz_width_fraction: _params.lorentz_width_fraction,
      v_center: _params.v_center,
      inst_resolution: _params.inst_resolution,
      limb_darkening: _params.limb_darkening,
      beta: _params.beta,
      fI: _params.fI,
    };
    const data = await apiPost('/api/auto_estimate', body);
    if (data.line_strength != null) {
      _params.line_strength = data.line_strength;
      // Update slider DOM
      const slider = document.getElementById('slider-line_strength');
      const meta = SLIDER_META.line_strength;
      if (slider) {
        slider.value = clamp(data.line_strength, parseFloat(slider.min), parseFloat(slider.max));
        document.getElementById('val-line_strength').textContent = formatVal(data.line_strength, meta.decimals);
      }
      // Update chart & status
      if (data.I_model) updatePlotModel(data.I_model);
      updateStatus(data.rms);
      updateParamsPanel();
      // Flash feedback
      const el = document.getElementById('status-bar');
      const prev = el.innerHTML;
      el.innerHTML = `<span style="color: var(--green);">Auto d = ${formatVal(data.line_strength, 5)}  (RMS ${data.rms != null ? data.rms.toFixed(5) : '—'})</span>`;
      setTimeout(() => { el.innerHTML = prev; }, 2500);
    }
  } catch (e) {
    console.error('Auto-estimate error:', e);
  } finally {
    btn.disabled = false;
    btn.textContent = 'Auto d';
  }
}

// ---------------------------------------------------------------------------
// Init
// ---------------------------------------------------------------------------
async function init() {
  try {
    const data = await apiGet('/api/init');
    _velObs = data.vel_obs;
    _medianI = data.median_I;
    _individualI = data.individual_I;
    _velModel = data.vel_model;
    _params = data.params;
    _ranges = data.slider_ranges;
    _initParams = { ...data.params };
    _currentModel = data.params.model_type;

    setActiveModelButton(_currentModel);
    buildSliders();
    initModelSelector();
    initButtons();
    initPlot();
    await doCompute();
  } catch (e) {
    console.error('Init failed:', e);
    document.getElementById('status-bar').textContent = 'Error: ' + e.message;
  }
}

init();
