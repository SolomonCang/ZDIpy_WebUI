/**
 * config.js — Configuration form for ZDIpy WebUI.
 *
 * Renders accordion sections from CONFIG_SCHEMA, syncs with GET/PUT /api/config.
 * The observations.files array is preserved from the server during save.
 */

// ---------------------------------------------------------------------------
// Schema definition
// ---------------------------------------------------------------------------
const CONFIG_SCHEMA = [
  {
    section: 'star', label: 'Stellar Parameters', icon: '⭐', open: true,
    fields: [
      { key: 'inclination_deg',                    label: 'Inclination (°)',                    type: 'number', min: 0,   max: 90,  step: 0.5,  tooltip: 'Angle between the rotation axis and the line of sight (0° = pole-on, 90° = equator-on).' },
      { key: 'vsini_kms',                          label: 'v sin i  (km/s)',                    type: 'number', min: 0,   max: 500, step: 0.1,  tooltip: 'Projected equatorial rotation velocity (km/s). Determines the width of the spectral line profile.' },
      { key: 'period_days',                        label: 'Period  (days)',                     type: 'number', min: 0,             step: 1e-4, tooltip: 'Stellar rotation period in days. Used to convert observation timestamps into rotation phases.' },
      { key: 'differential_rotation_rad_per_day',  label: 'Differential rotation  (rad/day)',   type: 'number',                     step: 1e-3, tooltip: 'Differential rotation coefficient dΩ (rad/day): angular velocity difference between equator and poles.' },
      { key: 'mass_msun',                          label: 'Mass  (M☉)',                         type: 'number', min: 0,             step: 0.01, tooltip: 'Stellar mass in solar masses. Used for gravity darkening calculations.' },
      { key: 'radius_rsun',                        label: 'Radius  (R☉)',                       type: 'number', min: 0,             step: 0.01, tooltip: 'Stellar radius in solar radii. Used for gravity darkening and surface grid scaling.' },
    ],
  },
  {
    section: 'grid', label: 'Grid', icon: '🌐', open: true,
    fields: [
      { key: 'nRings', label: 'Grid rings  (nRings)', type: 'number', min: 5, max: 200, step: 1, tooltip: 'Number of latitude rings in the stellar surface grid. Higher values give finer resolution but increase computation time.' },
    ],
  },
  {
    section: 'inversion', label: 'Inversion Control', icon: '🔄', open: true,
    fields: [
      { key: 'target_form',    label: 'Target form',   type: 'select', options: [{value: 'C', label: 'Chi2'}, {value: 'E', label: 'Entropy'}], tooltip: 'MEM convergence target: Chi2 drives the solution toward a target χ² value; Entropy drives it toward a target entropy value.' },
      { key: 'target_value',   label: 'Target value',  type: 'number',             step: 0.01,  tooltip: 'Numerical target for the selected form. For Chi2: reduced χ² ≈ 1.0 is ideal. For Entropy: target entropy value.' },
      { key: 'num_iterations', label: 'Max iterations', type: 'number', min: 0, max: 1000, step: 1, tooltip: 'Maximum number of MEM iterations. The solver stops early if convergence is reached before this limit.' },
      { key: 'test_aim',       label: 'Convergence test_aim', type: 'number',      step: 1e-5,  tooltip: 'Convergence threshold: iteration stops when the aim function change falls below this value. Smaller = stricter.' },
    ],
  },
  {
    section: 'magnetic', label: 'Magnetic Field', icon: '🧲', open: false,
    fields: [
      { key: 'fit_magnetic',   label: 'Fit magnetic map',                    type: 'checkbox', tooltip: 'Enable inversion of the magnetic field map from Stokes V profiles.' },
      { key: 'l_max',          label: 'ℓ_max  (spherical harmonic order)',   type: 'number', min: 1, max: 30, step: 1, tooltip: 'Maximum spherical harmonic order ℓ for the magnetic field expansion. Higher values allow finer magnetic structures but increase parameters.' },
      { key: 'default_bent',   label: 'Default B entropy slope  (G)',        type: 'number',             step: 1,    tooltip: 'Default entropy reference level for the magnetic field (Gauss). Controls the regularisation strength.' },
      { key: 'geometry_type',  label: 'Geometry type',                       type: 'select', options: ['Full', 'Poloidal', 'PotTor', 'Potential'], tooltip: 'Magnetic field topology: Full = all components, Poloidal = poloidal only, PotTor = potential + toroidal, Potential = potential field.' },
      { key: 'init_from_file', label: 'Initialise from file',                type: 'checkbox', tooltip: 'Load initial magnetic field coefficients from the file specified below instead of starting from zero.' },
      { key: 'init_file',      label: 'Init coefficients file',              type: 'text',     tooltip: 'Path to a file containing initial spherical harmonic coefficients for the magnetic field.' },
    ],
  },
  {
    section: 'brightness', label: 'Brightness Map', icon: '☀️', open: false,
    fields: [
      { key: 'fit_brightness', label: 'Fit brightness map',              type: 'checkbox', tooltip: 'Enable inversion of the surface brightness map from Stokes I profiles.' },
      { key: 'chi2_scale_I',   label: 'χ² scale for Stokes I',          type: 'number', step: 0.01,  tooltip: 'Weighting factor applied to the Stokes I χ² contribution during inversion.' },
      { key: 'chi2_scale_V',   label: 'χ² scale for Stokes V',          type: 'number', step: 0.01,  tooltip: 'Weighting factor applied to the Stokes V χ² contribution during joint inversion. Increase this when brightness fitting overwhelms the magnetic update.' },
      { key: 'entropy_scale',  label: 'Entropy scale',                   type: 'number', step: 0.01,  tooltip: 'Weighting factor for the entropy regularisation term in the brightness inversion.' },
      { key: 'entropy_form',   label: 'Entropy form  (1=image, 2=fill)', type: 'select', options: [1, 2], tooltip: 'Entropy functional form: 1 = image entropy (favours uniform map), 2 = fill entropy (favours uniform deviation from default).' },
      { key: 'default_bright', label: 'Default brightness',              type: 'number', step: 0.01,  tooltip: 'Reference (default) brightness for the entropy regularisation. Pixels are regularised toward this value.' },
      { key: 'max_bright',     label: 'Max brightness',                  type: 'number', step: 0.001, tooltip: 'Maximum allowed brightness per pixel. Acts as an upper bound constraint during the inversion.' },
      { key: 'init_from_file', label: 'Initialise from file',            type: 'checkbox', tooltip: 'Load an initial brightness map from file rather than starting from the default brightness.' },
      { key: 'init_file',      label: 'Init brightness file',            type: 'text',     tooltip: 'Path to a file containing the initial per-pixel brightness values.' },
    ],
  },
  {
    section: 'line_model', label: 'Spectral Line Model', icon: '📡', open: false,
    fields: [
      // ── Model selector (spans full width) ──
      { key: 'model_type', label: 'Line profile model', type: 'select', fullWidth: true,
        options: [
          { value: 'voigt',            label: 'Voigt  (weak-field approximation)' },
          { value: 'unno',             label: 'Unno-Rachkovsky  (Milne-Eddington, full polarised RT)' },
          { value: 'halpha_compound',  label: 'H-alpha compound double-Voigt  (weak-field)' },
        ],
        tooltip: 'Spectral line model: "Voigt" uses the weak-field approximation (fast). "Unno-Rachkovsky" solves full polarised RT. "H-alpha compound" models emission + absorption double-Voigt for H\u03b1 profiles.',
      },
      // ── Common parameters ──
      { key: 'estimate_strength',       label: 'Auto-estimate line strength',         type: 'checkbox', tooltip: 'Automatically estimate the LSD line strength from the observed profiles instead of using the fixed value below.' },
      { key: 'wavelength_nm',           label: 'Wavelength  (nm)',                    type: 'number', step: 0.001,  tooltip: 'Central wavelength of the LSD line (nm). Used for computing the Zeeman splitting in Stokes V.' },
      { key: 'line_strength',           label: 'Line strength',                       type: 'number', step: 1e-4,   tooltip: 'LSD equivalent line strength (depth weight). Ignored when auto-estimate is enabled.' },
      { key: 'gauss_width_kms',         label: 'Gaussian width  (km/s)',              type: 'number', step: 0.01,   tooltip: 'Gaussian component of the Voigt profile width (km/s), representing thermal/turbulent broadening.' },
      { key: 'lorentz_width_fraction',  label: 'Lorentz width fraction',              type: 'number', min: 0, max: 1, step: 0.001, tooltip: 'Fraction of the total Voigt width contributed by the Lorentzian component (0 = pure Gaussian, 1 = pure Lorentzian).' },
      { key: 'lande_g',                 label: 'Landé g factor',                      type: 'number', step: 0.001,  tooltip: 'Effective Landé g factor of the LSD line. Scales the magnetic splitting in the Stokes V profile.' },
      { key: 'limb_darkening',          label: 'Limb darkening coefficient',          type: 'number', min: 0, max: 1, step: 0.01, tooltip: 'Linear limb darkening coefficient ε (0–1). Higher values darken the stellar limb more strongly.' },
      { key: 'gravity_darkening',       label: 'Gravity darkening coefficient',       type: 'number', min: 0, max: 1, step: 0.01, tooltip: 'Gravity darkening exponent β (0–1). Applies von Zeipel gravity darkening to surface elements.' },
      // ── Unno-Rachkovsky specific parameters ──
      { key: '_unno_heading', type: 'subheader', label: 'Unno-Rachkovsky parameters',
        visibleWhen: { key: 'model_type', value: 'unno' } },
      { key: 'unno_beta',             label: 'Source function slope  β',       type: 'number', step: 0.01,
        visibleWhen: { key: 'model_type', value: 'unno' },
        tooltip: 'Slope of the Planck function with optical depth in the Milne-Eddington atmosphere: B(τ) = B₀(1 + β·τ). Set ≤ 0 to derive automatically from the limb darkening coefficient.' },
      { key: 'unno_filling_factor_I', label: 'Stokes I filling factor  f_I', type: 'number', min: 0, max: 1, step: 0.01,
        visibleWhen: { key: 'model_type', value: 'unno' },
        tooltip: 'Fraction of each surface element covered by the magnetised atmosphere that contributes to Stokes I. The remainder contributes a field-free (quiet) Stokes I profile.' },
      { key: 'unno_filling_factor_V', label: 'Stokes V filling factor  f_V', type: 'number', min: 0, step: 0.01,
        visibleWhen: { key: 'model_type', value: 'unno' },
        tooltip: 'Filling factor applied to the Stokes V profile. Values < 1 reduce the effective Zeeman splitting (mimicking unresolved mixed-polarity elements). Values > 1 are unusual but technically permitted.' },
      // ── H-alpha compound specific parameters ──
      { key: '_halpha_heading', type: 'subheader', label: 'H\u03b1 compound model parameters',
        visibleWhen: { key: 'model_type', value: 'halpha_compound' } },
      { key: 'emission_strength',        label: 'Emission strength',                   type: 'number', min: 0, step: 0.01,
        visibleWhen: { key: 'model_type', value: 'halpha_compound' },
        tooltip: 'Peak amplitude of the emission component in the H\u03b1 compound profile. Larger values give a stronger filled-in emission core.' },
      { key: 'emission_gauss_kms',       label: 'Emission Gaussian width  (km/s)',      type: 'number', min: 0, step: 1,
        visibleWhen: { key: 'model_type', value: 'halpha_compound' },
        tooltip: 'Gaussian half-width of the broad emission component (km/s). Controls how far the chromospheric emission extends in velocity.' },
      { key: 'emission_lorentz_ratio',   label: 'Emission Lorentz fraction',            type: 'number', min: 0, max: 1, step: 0.01,
        visibleWhen: { key: 'model_type', value: 'halpha_compound' },
        tooltip: 'Lorentzian fraction of the emission Voigt profile (0 = pure Gaussian, 1 = pure Lorentzian).' },
      { key: 'absorption_strength',      label: 'Absorption strength',                  type: 'number', min: 0, step: 0.01,
        visibleWhen: { key: 'model_type', value: 'halpha_compound' },
        tooltip: 'Peak depth of the narrow absorption component superimposed on the emission core. Set to 0 to disable.' },
      { key: 'absorption_gauss_kms',     label: 'Absorption Gaussian width  (km/s)',    type: 'number', min: 0, step: 0.5,
        visibleWhen: { key: 'model_type', value: 'halpha_compound' },
        tooltip: 'Gaussian half-width of the narrow absorption component (km/s).' },
      { key: 'absorption_lorentz_ratio', label: 'Absorption Lorentz fraction',           type: 'number', min: 0, max: 1, step: 0.01,
        visibleWhen: { key: 'model_type', value: 'halpha_compound' },
        tooltip: 'Lorentzian fraction of the absorption Voigt profile (0 = pure Gaussian, 1 = pure Lorentzian).' },
      { key: 'filling_factor_V',         label: 'Stokes V filling factor  f_V',         type: 'number', min: 0, step: 0.01,
        visibleWhen: { key: 'model_type', value: 'halpha_compound' },
        tooltip: 'Filling factor applied to the H\u03b1 Stokes V profile. Values < 1 reduce the effective Zeeman signal.' },
      // ── H-alpha pre-processing flags ──
      { key: '_halpha_preproc_heading', type: 'subheader', label: 'H\u03b1 预处理选项',
        visibleWhen: { key: 'model_type', value: 'halpha_compound' } },
      { key: 'halpha_normalize_emission', label: '发射强度归一化',       type: 'checkbox',
        visibleWhen: { key: 'model_type', value: 'halpha_compound' },
        tooltip: '对各历元 Stokes I 发射峰高度做中值归一化，同步缩放 V/N 及对应噪声数组，使不同活动水平的观测可比。' },
      { key: 'halpha_auto_init',         label: '自动估算复合模型初始参数', type: 'checkbox',
        visibleWhen: { key: 'model_type', value: 'halpha_compound' },
        tooltip: '对归一化后的 Stokes I 取中值谱，拟合双 Voigt 模型，自动设定发射/吸收成分的初始参数，并在 Run 面板中绘图显示。' },
    ],
  },
  {
    section: 'instrument', label: 'Instrument', icon: '🔬', open: false,
    fields: [
      { key: 'spectral_resolution', label: 'Spectral resolution  R', type: 'number', min: 1000, step: 1000, tooltip: 'Instrument spectral resolving power R = λ/Δλ. Used to convolve the synthetic profiles with the instrumental profile.' },
    ],
  },
  {
    section: 'velocity_grid', label: 'Velocity Grid', icon: '📈', open: false,
    fields: [
      { key: 'vel_start_kms', label: 'Grid start  (km/s)', type: 'number', step: 1, tooltip: 'Start of the velocity grid in km/s. Should extend sufficiently beyond ±v sin i to capture the full line profile.' },
      { key: 'vel_end_kms',   label: 'Grid end  (km/s)',   type: 'number', step: 1, tooltip: 'End of the velocity grid in km/s. Must be symmetric with Grid start and cover the entire spectral line.' },
    ],
  },
  {
    // jdate_ref only — the files array is managed by observations.js
    section: 'observations', label: 'Observations (reference date)', icon: '🕰️', open: false,
    fields: [
      { key: 'jdate_ref', label: 'Reference Julian Date  (jdate_ref)', type: 'number', step: 0.001, tooltip: 'Reference Julian Date used to compute rotation phases: phase = (JD − jdate_ref) / period.' },
    ],
  },
  {
    section: 'output', label: 'Output Files', icon: '📁', open: false,
    fields: [
      { key: 'mag_coeff_file',       label: 'Magnetic coefficients file',              type: 'text', tooltip: 'Output file path for the spherical harmonic coefficients of the recovered magnetic field.' },
      { key: 'bright_map_file',      label: 'Brightness map file',                     type: 'text', tooltip: 'Output file path for the recovered surface brightness map (per-pixel values).' },
      { key: 'bright_map_gdark_file',label: 'Brightness map (gravity dark.) file',     type: 'text', tooltip: 'Output file path for the brightness map corrected for gravity darkening.' },
      { key: 'line_models_file',     label: 'Line models file',                        type: 'text', tooltip: 'Output file path for the synthetic (model) line profiles at each observed phase.' },
      { key: 'observed_used_file',   label: 'Observed (used) file',                    type: 'text', tooltip: 'Output file path listing the observed profiles that were actually used in the inversion.' },
      { key: 'fit_summary_file',     label: 'Fit summary file',                        type: 'text', tooltip: 'Output file path for the fit summary: χ², entropy, convergence status per iteration.' },
    ],
  },
];

// ---------------------------------------------------------------------------
// Conditional field visibility
// ---------------------------------------------------------------------------

/**
 * Walk every element that carries [data-cond-field] / [data-cond-val] and
 * show/hide it based on the current value of its controller element.
 * Called after _buildForm and after _populate.
 */
function _applyAllConditionals() {
  document.querySelectorAll('[data-cond-field][data-cond-val]').forEach(el => {
    const ctrl = document.getElementById(el.dataset.condField);
    const show = ctrl != null && String(ctrl.value) === el.dataset.condVal;
    el.classList.toggle('field-row--hidden', !show);
  });
}

/** Attach 'change' listeners to every controller so UI updates live. */
function _setupConditionalListeners() {
  const seen = new Set();
  document.querySelectorAll('[data-cond-field]').forEach(el => {
    const id = el.dataset.condField;
    if (seen.has(id)) return;
    seen.add(id);
    const ctrl = document.getElementById(id);
    if (!ctrl) return;
    ctrl.addEventListener('change', _applyAllConditionals);
    ctrl.addEventListener('input',  _applyAllConditionals);
  });
}

// ---------------------------------------------------------------------------
// DOM rendering
// ---------------------------------------------------------------------------
function _fieldId(section, key) { return `field-${section}-${key}`; }

function _helpIcon(tooltip) {
  if (!tooltip) return '';
  const escaped = tooltip.replace(/"/g, '&quot;');
  return `<span class="param-help" data-tip="${escaped}">?</span>`;
}

function _renderField(section, f) {
  const id = _fieldId(section, f.key);

  // Build conditional data attributes (appended to wrapper)
  let condAttrs = '';
  let condClass = '';
  if (f.visibleWhen) {
    const ctrlId = _fieldId(section, f.visibleWhen.key);
    condAttrs = ` data-cond-field="${ctrlId}" data-cond-val="${f.visibleWhen.value}"`;
    condClass = ' field-row--conditional field-row--hidden';
  }
  const fullWidthClass = f.fullWidth ? ' field-row--full' : '';

  // Subheader / visual separator
  if (f.type === 'subheader') {
    return `<div class="field-subheader${condClass}${fullWidthClass}"${condAttrs}>${f.label}</div>`;
  }

  const wrapClass = `field-row${condClass}${fullWidthClass}`;

  if (f.type === 'checkbox') {
    return `
      <div class="${wrapClass}"${condAttrs}>
        <label class="checkbox-label" for="${id}">
          <input type="checkbox" id="${id}" />
          ${f.label}
          ${_helpIcon(f.tooltip)}
        </label>
      </div>`;
  }

  if (f.type === 'select') {
    const opts = f.options.map(o => {
      const val = (typeof o === 'object') ? o.value : o;
      const lbl = (typeof o === 'object') ? o.label : o;
      return `<option value="${val}">${lbl}</option>`;
    }).join('');
    return `
      <div class="${wrapClass}"${condAttrs}>
        <label class="field-label" for="${id}">${f.label}${_helpIcon(f.tooltip)}</label>
        <select id="${id}" class="field-select">${opts}</select>
      </div>`;
  }

  // number / text
  const attrs = [
    `type="${f.type}"`,
    `id="${id}"`,
    `class="field-input"`,
    f.min  !== undefined ? `min="${f.min}"`   : '',
    f.max  !== undefined ? `max="${f.max}"`   : '',
    f.step !== undefined ? `step="${f.step}"` : '',
  ].filter(Boolean).join(' ');

  return `
    <div class="${wrapClass}"${condAttrs}>
      <label class="field-label" for="${id}">${f.label}${_helpIcon(f.tooltip)}</label>
      <input ${attrs} />
    </div>`;
}

function _buildForm() {
  const container = document.getElementById('cfg-form');
  if (!container) return;

  const html = CONFIG_SCHEMA.map(sec => {
    const bodyHtml = `<div class="field-grid">${sec.fields.map(f => _renderField(sec.section, f)).join('')}</div>`;
    return `
      <div class="accordion${sec.open ? ' open' : ''}">
        <div class="accordion-header">
          <span class="accordion-icon">${sec.icon}</span>
          <span class="accordion-label">${sec.label}</span>
          <span class="accordion-arrow">▶</span>
        </div>
        <div class="accordion-body">${bodyHtml}</div>
      </div>`;
  }).join('');

  container.innerHTML = html;
  _applyAllConditionals();       // set initial visibility before any load
  _setupConditionalListeners();  // wire change events
}

// ---------------------------------------------------------------------------
// Populate / collect
// ---------------------------------------------------------------------------
function _populate(cfg) {
  for (const sec of CONFIG_SCHEMA) {
    const secData = cfg[sec.section] || {};
    for (const f of sec.fields) {
      const el = document.getElementById(_fieldId(sec.section, f.key));
      if (!el) continue;
      const val = secData[f.key];
      if (val === undefined || val === null) continue;
      if (f.type === 'checkbox') {
        el.checked = Boolean(val);
      } else {
        el.value = val;
      }
    }
  }
  // Re-apply conditional visibility after values have been set
  _applyAllConditionals();
}

function _collect() {
  const cfg = {};
  for (const sec of CONFIG_SCHEMA) {
    cfg[sec.section] = cfg[sec.section] || {};
    for (const f of sec.fields) {
      // Skip pseudo-fields (subheaders have no input element)
      if (f.type === 'subheader') continue;
      const el = document.getElementById(_fieldId(sec.section, f.key));
      if (!el) continue;
      if (f.type === 'checkbox') {
        cfg[sec.section][f.key] = el.checked ? 1 : 0;
      } else if (f.type === 'number') {
        cfg[sec.section][f.key] = parseFloat(el.value);
      } else if (f.type === 'select') {
        // preserve numeric option types
        const raw = el.value;
        cfg[sec.section][f.key] = isNaN(Number(raw)) ? raw : Number(raw);
      } else {
        cfg[sec.section][f.key] = el.value;
      }
    }
  }
  return cfg;
}

// ---------------------------------------------------------------------------
// Load / save actions
// ---------------------------------------------------------------------------
async function loadConfig() {
  const statusEl = document.getElementById('cfg-status');
  try {
    const cfg = await apiFetch('/api/config');
    _populate(cfg);
    setStatus(statusEl, 'Config loaded.', 'ok');
    // clear status after 2 s
    setTimeout(() => setStatus(statusEl, ''), 2000);
  } catch (err) {
    setStatus(statusEl, `Error loading config: ${err.message}`, 'error');
  }
}

async function saveConfig() {
  const statusEl = document.getElementById('cfg-status');
  setStatus(statusEl, 'Saving…');
  try {
    // Read-modify-write: preserve observations.files from server
    const current = await apiFetch('/api/config');
    const formCfg = _collect();

    // Merge observations: keep the server's files array
    formCfg.observations = {
      ...formCfg.observations,
      files: (current.observations || {}).files || [],
    };

    await apiFetch('/api/config', { method: 'PUT', body: formCfg });
    setStatus(statusEl, '✔ Config saved.', 'ok');
    setTimeout(() => setStatus(statusEl, ''), 3000);
  } catch (err) {
    setStatus(statusEl, `❌ ${err.message}`, 'error');
  }
}

// ---------------------------------------------------------------------------
// Initialise
// ---------------------------------------------------------------------------
document.addEventListener('DOMContentLoaded', () => {
  _buildForm();
  loadConfig();

  document.getElementById('cfg-load-btn')?.addEventListener('click', loadConfig);
  document.getElementById('cfg-save-btn')?.addEventListener('click', saveConfig);
});
