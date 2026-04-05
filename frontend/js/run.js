/**
 * run.js — Run controls and real-time SSE log streaming.
 *
 * POST /api/run           → starts background ZDI thread
 * GET  /api/run/stream    → SSE log stream (EventSource)
 * GET  /api/run/status    → poll status when SSE is not active
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
// SSE stream
// ---------------------------------------------------------------------------
function _startStream() {
  if (_sse) { _sse.close(); _sse = null; }

  _sse = new EventSource('/api/run/stream');

  _sse.onmessage = e => {
    let msg;
    try { msg = JSON.parse(e.data); } catch (_) { return; }

    if (msg.type === 'log') {
      _appendLog(msg.text);
    } else if (msg.type === 'complete') {
      _setBadge('done', 'Done');
      setStatus(_el('run-status-bar'), '✔ Run complete.', 'ok');
      _isRunning = false;
      _setRunBtnState(false);
      _sse.close(); _sse = null;
      window.loadAllResults?.();
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
    } else if (s.status === 'error') {
      _setBadge('error', 'Error');
    }
  }).catch(() => {});   // non-fatal on first load
});
