/**
 * run.js — Run controls and real-time SSE log streaming.
 *
 * POST /api/run           → starts background ZDI thread
 * GET  /api/run/stream    → SSE log stream (EventSource)
 * GET  /api/run/status    → poll status when SSE is not active
 * GET  /api/run/halpha_init_plot → Plotly JSON for Hα pre-processing result
 */

let _sse = null;          // active EventSource or null
let _isRunning = false;

// ---------------------------------------------------------------------------
// DOM references
// ---------------------------------------------------------------------------
function _el(id) { return document.getElementById(id); }

function _setBadge(cls, text) {
  const badge = _el('run-badge');
  if (!badge) return;
  badge.className = `badge badge-${cls}`;
  badge.textContent = text;
}

function _appendLog(text) {
  const log = _el('run-log');
  if (!log) return;
  log.textContent += text + '\n';
  log.scrollTop = log.scrollHeight;
}

// ---------------------------------------------------------------------------
// Hα pre-processing plot
// ---------------------------------------------------------------------------
async function _loadHalphaInitPlot() {
  const card    = _el('halpha-preproc-card');
  const plotDiv = _el('halpha-preproc-plot');
  const noData  = _el('halpha-preproc-no-data');
  const badge   = _el('halpha-preproc-badge');
  if (!card || !plotDiv) return;

  try {
    const resp = await apiFetch('/api/run/halpha_init_plot');
    // resp is { data: [...], layout: {...} }
    card.style.display = '';
    if (noData) noData.style.display = 'none';
    if (badge) badge.textContent = '已更新';

    const darkLayout = {
      ...resp.layout,
      paper_bgcolor: 'var(--surface, #1e1e2e)',
      plot_bgcolor:  'var(--surface, #1e1e2e)',
      font:          { color: 'var(--text-primary, #cdd6f4)' },
    };
    Plotly.react(plotDiv, resp.data, darkLayout,
      { responsive: true, displayModeBar: true });
  } catch (err) {
    if (err.message && err.message.startsWith('404')) {
      // 预处理未启用，静默忽略
      if (card) card.style.display = 'none';
    } else {
      console.warn('[halpha_init_plot]', err.message);
    }
  }
}

// ---------------------------------------------------------------------------
// SSE stream
// ---------------------------------------------------------------------------
function _startStream() {
  if (_sse) { _sse.close(); _sse = null; }

  _sse = new EventSource('/api/run/stream');

  _sse.onmessage = e => {
    let msg;
    try { msg = JSON.parse(e.data); } catch (_) { return; }

    if (msg.type === 'log') {
      // Internal signal: Hα init plot is ready in state — fetch & render it
      if (msg.text === '[HALPHA_INIT_PLOT_READY]') {
        _loadHalphaInitPlot();
        return; // do not display this internal marker in the log
      }
      _appendLog(msg.text);
    } else if (msg.type === 'complete') {
      _setBadge('done', 'Done');
      setStatus(_el('run-status-bar'), '✔ Run complete.', 'ok');
      _isRunning = false;
      _setRunBtnState(false);
      _sse.close(); _sse = null;
      window.loadAllResults?.();
      // Also try to load the plot in case the signal was missed
      _loadHalphaInitPlot();
    } else if (msg.type === 'error') {
      _setBadge('error', 'Error');
      setStatus(_el('run-status-bar'), '❌ Run failed — see log for details.', 'error');
      _isRunning = false;
      _setRunBtnState(false);
      _sse.close(); _sse = null;
    }
  };

  _sse.onerror = () => {
    // SSE connection closed by server after run finishes — this is expected
    if (_sse) { _sse.close(); _sse = null; }
    if (_isRunning) {
      _setBadge('error', 'Disconnected');
      _isRunning = false;
      _setRunBtnState(false);
    }
  };
}

// ---------------------------------------------------------------------------
// Run button
// ---------------------------------------------------------------------------
function _setRunBtnState(running) {
  const btn = _el('run-btn');
  if (!btn) return;
  btn.disabled = running;
  btn.textContent = running ? '⏳ Running…' : '▶ Run ZDI';
}

async function startRun() {
  const statusEl = _el('run-status-bar');

  // Clear log
  const log = _el('run-log');
  if (log) log.textContent = '';

  // Hide previous Hα plot
  const card = _el('halpha-preproc-card');
  if (card) card.style.display = 'none';

  const forwardOnly = _el('run-forward-only')?.checked ?? false;
  const verbose     = parseInt(_el('run-verbose')?.value ?? '1', 10);

  try {
    const resp = await apiFetch('/api/run', {
      method: 'POST',
      body: { forward_only: forwardOnly, verbose },
    });

    if (resp.status === 'busy') {
      setStatus(statusEl, '⚠️ Another run is in progress. Please wait.', 'warn');
      return;
    }

    _isRunning = true;
    _setRunBtnState(true);
    _setBadge('running', 'Running');
    setStatus(statusEl, 'ZDI pipeline started — streaming log…');
    _startStream();
  } catch (err) {
    setStatus(statusEl, `❌ Failed to start: ${err.message}`, 'error');
  }
}

// ---------------------------------------------------------------------------
// Init
// ---------------------------------------------------------------------------
document.addEventListener('DOMContentLoaded', () => {
  _el('run-btn')?.addEventListener('click', startRun);

  _el('run-clear-btn')?.addEventListener('click', () => {
    const log = _el('run-log');
    if (log) log.textContent = '';
  });

  // Restore badge state on page reload
  apiFetch('/api/run/status').then(s => {
    if (s.status === 'running') {
      _isRunning = true;
      _setRunBtnState(true);
      _setBadge('running', 'Running');
      _startStream();
    } else if (s.status === 'done') {
      _setBadge('done', 'Done');
      if (s.log_tail.length) _el('run-log').textContent = s.log_tail.join('\n');
      // Try to load the Hα init plot if present
      _loadHalphaInitPlot();
    } else if (s.status === 'error') {
      _setBadge('error', 'Error');
    }
  }).catch(() => {});   // non-fatal on first load
});
