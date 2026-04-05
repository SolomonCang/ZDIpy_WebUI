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
// Plotly common layout fragment (dark theme)
// ---------------------------------------------------------------------------
const PLOTLY_BASE_LAYOUT = {
  paper_bgcolor: '#21262d',
  plot_bgcolor:  '#0d1117',
  font: { family: "'JetBrains Mono', 'Consolas', monospace", color: '#f0f6fc', size: 10 },
  margin: { t: 28, r: 16, b: 40, l: 50 },
  xaxis: { gridcolor: '#30363d', zerolinecolor: '#484f58', color: '#8b949e' },
  yaxis: { gridcolor: '#30363d', zerolinecolor: '#484f58', color: '#8b949e' },
};

const PLOTLY_CONFIG = {
  responsive: true,
  displayModeBar: true,
  modeBarButtonsToRemove: ['sendDataToCloud', 'autoScale2d', 'resetScale2d'],
  displaylogo: false,
};

// Expose globals used by other modules
window.apiFetch       = apiFetch;
window.setStatus      = setStatus;
window.PLOTLY_BASE_LAYOUT = PLOTLY_BASE_LAYOUT;
window.PLOTLY_CONFIG      = PLOTLY_CONFIG;
