/**
 * app.js — Tab navigation and shared utilities for ZDIpy WebUI.
 */

// ---------------------------------------------------------------------------
// Tab navigation
// ---------------------------------------------------------------------------
document.querySelectorAll('.tab-btn').forEach(btn => {
  btn.addEventListener('click', () => {
    const target = btn.dataset.tab;

    document.querySelectorAll('.tab-btn').forEach(b => {
      b.classList.toggle('active', b.dataset.tab === target);
      b.setAttribute('aria-selected', b.dataset.tab === target);
    });
    document.querySelectorAll('.tab-panel').forEach(p => {
      p.classList.toggle('active', p.id === `tab-${target}`);
    });
  });
});

// ---------------------------------------------------------------------------
// Shared fetch wrapper
// ---------------------------------------------------------------------------
/**
 * Thin wrapper around fetch that handles JSON response/error automatically.
 * @param {string} url
 * @param {RequestInit} [opts]
 * @returns {Promise<any>}
 */
async function apiFetch(url, opts = {}) {
  const defaults = {
    headers: { 'Content-Type': 'application/json' },
  };
  const merged = { ...defaults, ...opts };
  if (merged.body && typeof merged.body === 'object' && !(merged.body instanceof FormData)) {
    merged.body = JSON.stringify(merged.body);
  }
  if (merged.body instanceof FormData) {
    delete merged.headers['Content-Type'];  // let browser set multipart boundary
  }
  const resp = await fetch(url, merged);
  if (!resp.ok) {
    let detail = resp.statusText;
    try { const j = await resp.json(); detail = j.detail || JSON.stringify(j); } catch (_) {}
    throw new Error(`${resp.status}: ${detail}`);
  }
  // Some endpoints return 204 No Content
  if (resp.status === 204) return null;
  return resp.json();
}

// ---------------------------------------------------------------------------
// Status bar helpers
// ---------------------------------------------------------------------------
/**
 * Set status bar content with optional severity class.
 * @param {HTMLElement} el
 * @param {string} msg
 * @param {'ok'|'error'|'warn'|''} [cls]
 */
function setStatus(el, msg, cls = '') {
  if (!el) return;
  el.textContent = msg;
  el.className = 'status-bar' + (cls ? ` ${cls}` : '');
}

// ---------------------------------------------------------------------------
// Accordion toggle (delegated)
// ---------------------------------------------------------------------------
document.addEventListener('click', e => {
  const hdr = e.target.closest('.accordion-header');
  if (hdr) {
    const acc = hdr.closest('.accordion');
    if (acc) acc.classList.toggle('open');
  }
});

// ---------------------------------------------------------------------------
// Plotly layout fragments (theme-aware)
// ---------------------------------------------------------------------------
const _PLOTLY_DARK = {
  paper_bgcolor: '#21262d',
  plot_bgcolor:  '#0d1117',
  font: { family: "'JetBrains Mono', 'Consolas', monospace", color: '#f0f6fc', size: 10 },
  margin: { t: 28, r: 16, b: 40, l: 50 },
  xaxis: { gridcolor: '#30363d', zerolinecolor: '#484f58', color: '#8b949e' },
  yaxis: { gridcolor: '#30363d', zerolinecolor: '#484f58', color: '#8b949e' },
};
const _PLOTLY_LIGHT = {
  paper_bgcolor: '#ffffff',
  plot_bgcolor:  '#f8f8fa',
  font: { family: "'JetBrains Mono', 'Consolas', monospace", color: '#1c1c1e', size: 10 },
  margin: { t: 28, r: 16, b: 40, l: 50 },
  xaxis: { gridcolor: '#e0e0e8', zerolinecolor: '#c0c0c8', color: '#5a5a7a' },
  yaxis: { gridcolor: '#e0e0e8', zerolinecolor: '#c0c0c8', color: '#5a5a7a' },
};

// PLOTLY_BASE_LAYOUT is a live object — update its properties on theme change
const PLOTLY_BASE_LAYOUT = { ..._PLOTLY_DARK };

const PLOTLY_CONFIG = {
  responsive: true,
  displayModeBar: true,
  modeBarButtonsToRemove: ['sendDataToCloud', 'autoScale2d', 'resetScale2d'],
  displaylogo: false,
};

// ---------------------------------------------------------------------------
// Light / Dark theme toggle
// ---------------------------------------------------------------------------
const _THEME_KEY = 'zdipy-theme';

function _applyTheme(theme) {
  const isLight = theme === 'light';
  document.documentElement.dataset.theme = isLight ? 'light' : '';
  const btn = document.getElementById('theme-toggle-btn');
  if (btn) {
    btn.textContent = isLight ? '☀️' : '🌙';
    btn.title = isLight ? 'Switch to dark mode' : 'Switch to light mode';
    btn.setAttribute('aria-label', btn.title);
  }
  // Update Plotly layout tokens so newly-rendered charts pick up current theme
  Object.assign(PLOTLY_BASE_LAYOUT, isLight ? _PLOTLY_LIGHT : _PLOTLY_DARK);
}

// Initialise from stored preference (default: dark)
_applyTheme(localStorage.getItem(_THEME_KEY) || 'dark');

document.addEventListener('DOMContentLoaded', () => {
  const btn = document.getElementById('theme-toggle-btn');
  if (!btn) return;
  btn.addEventListener('click', () => {
    const next = document.documentElement.dataset.theme === 'light' ? 'dark' : 'light';
    localStorage.setItem(_THEME_KEY, next);
    _applyTheme(next);
  });
});

// Expose globals used by other modules
window.apiFetch       = apiFetch;
window.setStatus      = setStatus;
window.PLOTLY_BASE_LAYOUT = PLOTLY_BASE_LAYOUT;
window.PLOTLY_CONFIG      = PLOTLY_CONFIG;
