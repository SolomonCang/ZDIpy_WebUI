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
      { key: 'inclination_deg',                    label: 'Inclination (°)',                    type: 'number', min: 0,   max: 90,  step: 0.5 },
      { key: 'vsini_kms',                          label: 'v sin i  (km/s)',                    type: 'number', min: 0,   max: 500, step: 0.1 },
      { key: 'period_days',                        label: 'Period  (days)',                     type: 'number', min: 0,             step: 1e-4 },
      { key: 'differential_rotation_rad_per_day',  label: 'Differential rotation  (rad/day)',   type: 'number',                     step: 1e-3 },
      { key: 'mass_msun',                          label: 'Mass  (M☉)',                         type: 'number', min: 0,             step: 0.01 },
      { key: 'radius_rsun',                        label: 'Radius  (R☉)',                       type: 'number', min: 0,             step: 0.01 },
    ],
  },
  {
    section: 'grid', label: 'Grid', icon: '🔭', open: true,
    fields: [
      { key: 'nRings', label: 'Grid rings  (nRings)', type: 'number', min: 5, max: 200, step: 1 },
    ],
  },
  {
    section: 'inversion', label: 'Inversion Control', icon: '🔄', open: true,
    fields: [
      { key: 'target_form',    label: 'Target form  (C = χ²,  E = entropy)', type: 'select', options: ['C', 'E'] },
      { key: 'target_value',   label: 'Target value',                        type: 'number',             step: 0.01 },
      { key: 'num_iterations', label: 'Max iterations',                      type: 'number', min: 0, max: 1000, step: 1 },
      { key: 'test_aim',       label: 'Convergence test_aim',                type: 'number',             step: 1e-5 },
    ],
  },
  {
    section: 'magnetic', label: 'Magnetic Field', icon: '🧲', open: false,
    fields: [
      { key: 'fit_magnetic',   label: 'Fit magnetic map',                    type: 'checkbox' },
      { key: 'l_max',          label: 'ℓ_max  (spherical harmonic order)',   type: 'number', min: 1, max: 30, step: 1 },
      { key: 'default_bent',   label: 'Default B entropy slope  (G)',        type: 'number',             step: 1 },
      { key: 'geometry_type',  label: 'Geometry type',                       type: 'select', options: ['Full', 'Poloidal', 'PotTor', 'Potential'] },
      { key: 'init_from_file', label: 'Initialise from file',                type: 'checkbox' },
      { key: 'init_file',      label: 'Init coefficients file',              type: 'text' },
    ],
  },
  {
    section: 'brightness', label: 'Brightness Map', icon: '☀️', open: false,
    fields: [
      { key: 'fit_brightness', label: 'Fit brightness map',              type: 'checkbox' },
      { key: 'chi2_scale_I',   label: 'χ² scale for Stokes I',          type: 'number', step: 0.01 },
      { key: 'entropy_scale',  label: 'Entropy scale',                   type: 'number', step: 0.01 },
      { key: 'entropy_form',   label: 'Entropy form  (1=image, 2=fill)', type: 'select', options: [1, 2] },
      { key: 'default_bright', label: 'Default brightness',              type: 'number', step: 0.01 },
      { key: 'max_bright',     label: 'Max brightness',                  type: 'number', step: 0.001 },
      { key: 'init_from_file', label: 'Initialise from file',            type: 'checkbox' },
      { key: 'init_file',      label: 'Init brightness file',            type: 'text' },
    ],
  },
  {
    section: 'line_model', label: 'Spectral Line Model', icon: '📡', open: false,
    fields: [
      { key: 'estimate_strength',       label: 'Auto-estimate line strength',         type: 'checkbox' },
      { key: 'model_file',              label: 'Line model file',                     type: 'text' },
      { key: 'wavelength_nm',           label: 'Wavelength  (nm)',                    type: 'number', step: 0.001 },
      { key: 'line_strength',           label: 'Line strength',                       type: 'number', step: 1e-4 },
      { key: 'gauss_width_kms',         label: 'Gaussian width  (km/s)',              type: 'number', step: 0.01 },
      { key: 'lorentz_width_fraction',  label: 'Lorentz width fraction',              type: 'number', min: 0, max: 1, step: 0.001 },
      { key: 'lande_g',                 label: 'Landé g factor',                      type: 'number', step: 0.001 },
      { key: 'limb_darkening',          label: 'Limb darkening coefficient',          type: 'number', min: 0, max: 1, step: 0.01 },
      { key: 'gravity_darkening',       label: 'Gravity darkening coefficient',       type: 'number', min: 0, max: 1, step: 0.01 },
    ],
  },
  {
    section: 'instrument', label: 'Instrument', icon: '🔬', open: false,
    fields: [
      { key: 'spectral_resolution', label: 'Spectral resolution  R', type: 'number', min: 1000, step: 1000 },
    ],
  },
  {
    section: 'velocity_grid', label: 'Velocity Grid', icon: '📈', open: false,
    fields: [
      { key: 'vel_start_kms', label: 'Grid start  (km/s)', type: 'number', step: 1 },
      { key: 'vel_end_kms',   label: 'Grid end  (km/s)',   type: 'number', step: 1 },
    ],
  },
  {
    // jdate_ref only — the files array is managed by observations.js
    section: 'observations', label: 'Observations (reference date)', icon: '🕰️', open: false,
    fields: [
      { key: 'jdate_ref', label: 'Reference Julian Date  (jdate_ref)', type: 'number', step: 0.001 },
    ],
  },
  {
    section: 'output', label: 'Output Files', icon: '📁', open: false,
    fields: [
      { key: 'mag_coeff_file',       label: 'Magnetic coefficients file',              type: 'text' },
      { key: 'bright_map_file',      label: 'Brightness map file',                     type: 'text' },
      { key: 'bright_map_gdark_file',label: 'Brightness map (gravity dark.) file',     type: 'text' },
      { key: 'line_models_file',     label: 'Line models file',                        type: 'text' },
      { key: 'observed_used_file',   label: 'Observed (used) file',                    type: 'text' },
      { key: 'fit_summary_file',     label: 'Fit summary file',                        type: 'text' },
    ],
  },
];

// ---------------------------------------------------------------------------
// DOM rendering
// ---------------------------------------------------------------------------
function _fieldId(section, key) { return `field-${section}-${key}`; }

function _renderField(section, f) {
  const id = _fieldId(section, f.key);
  const wrapClass = 'field-row';

  if (f.type === 'checkbox') {
    return `
      <div class="${wrapClass}">
        <label class="checkbox-label" for="${id}">
          <input type="checkbox" id="${id}" />
          ${f.label}
        </label>
      </div>`;
  }

  if (f.type === 'select') {
    const opts = f.options.map(o => `<option value="${o}">${o}</option>`).join('');
    return `
      <div class="${wrapClass}">
        <label class="field-label" for="${id}">${f.label}</label>
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
    <div class="${wrapClass}">
      <label class="field-label" for="${id}">${f.label}</label>
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
}

function _collect() {
  const cfg = {};
  for (const sec of CONFIG_SCHEMA) {
    cfg[sec.section] = cfg[sec.section] || {};
    for (const f of sec.fields) {
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
