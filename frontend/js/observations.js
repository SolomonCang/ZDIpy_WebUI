/**
 * observations.js — Manage LSD profile files and observation config entries.
 *
 * Responsibilities:
 *   - Drag-and-drop / browse upload  → POST /api/observations/upload
 *   - List existing files            → GET  /api/observations
 *   - Delete a file                  → DELETE /api/observations/{fname}
 *   - Structured per-row obs editor  → read-modify-write via GET/PUT /api/config
 *   - File-existence validation      → POST /api/observations/validate
 */

// ---------------------------------------------------------------------------
// File list (LSDprof/ directory)
// ---------------------------------------------------------------------------
async function refreshFileList() {
  const tbody = document.getElementById('obs-file-tbody');
  const noMsg = document.getElementById('obs-no-files');
  if (!tbody) return;

  try {
    const files = await apiFetch('/api/observations');
    tbody.innerHTML = '';

    if (!files.length) {
      noMsg && (noMsg.style.display = '');
      return;
    }
    noMsg && (noMsg.style.display = 'none');

    files.forEach(f => {
      const tr = document.createElement('tr');
      const kb = (f.size_bytes / 1024).toFixed(1);
      tr.innerHTML = `
        <td>${escHtml(f.name)}</td>
        <td>${kb} KB</td>
        <td>
          <button class="btn btn-secondary btn-sm obs-link-btn" data-name="${escHtml(f.name)}" title="Add to observation entries">
            ＋ Link
          </button>
        </td>
        <td>
          <button class="btn btn-danger btn-sm obs-del-btn" data-name="${escHtml(f.name)}" title="Delete file">
            ✕
          </button>
        </td>`;
      tbody.appendChild(tr);
    });

    // Populate datalist for filename autocomplete in entry editor
    _populateDatalist(files.map(f => f.name));
  } catch (err) {
    setStatus(document.getElementById('obs-status'), `Error loading files: ${err.message}`, 'error');
  }
}

function _populateDatalist(names) {
  const dl = document.getElementById('obs-file-datalist');
  if (!dl) return;
  dl.innerHTML = '';
  names.forEach(n => {
    const opt = document.createElement('option');
    opt.value = `LSDprof/${n}`;
    dl.appendChild(opt);
  });
}

// File table button handler (delete + link-to-entry, delegated)
document.getElementById('obs-file-tbody')?.addEventListener('click', async e => {
  // ── Link button: add a new observation entry row pre-filled with this file
  const linkBtn = e.target.closest('.obs-link-btn');
  if (linkBtn) {
    const name = linkBtn.dataset.name;
    const tbody = document.getElementById('obs-entries-tbody');
    if (!tbody) return;
    tbody.appendChild(buildObsRow({ filename: `LSDprof/${name}`, jdate: '', vel_center_kms: '' }));
    _renumberRows();
    // Scroll entry table into view so user notices it was added
    document.querySelector('.obs-table-wrap')?.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
    setStatus(document.getElementById('obs-status'), `Added "${name}" to observation entries.`, 'ok');
    return;
  }

  // ── Delete button
  const delBtn = e.target.closest('.obs-del-btn');
  if (!delBtn) return;
  const name = delBtn.dataset.name;
  if (!confirm(`Delete "${name}"?`)) return;
  try {
    await apiFetch(`/api/observations/${encodeURIComponent(name)}`, { method: 'DELETE' });
    setStatus(document.getElementById('obs-status'), `Deleted ${name}.`, 'ok');
    refreshFileList();
  } catch (err) {
    setStatus(document.getElementById('obs-status'), `❌ ${err.message}`, 'error');
  }
});

// ---------------------------------------------------------------------------
// Upload helpers
// ---------------------------------------------------------------------------
async function uploadFiles(fileList) {
  const statusEl = document.getElementById('obs-status');
  if (!fileList.length) return;
  const fd = new FormData();
  for (const f of fileList) fd.append('files', f);
  try {
    setStatus(statusEl, `Uploading ${fileList.length} file(s)…`);
    const res = await apiFetch('/api/observations/upload', { method: 'POST', body: fd });
    setStatus(statusEl, `✔ Saved: ${res.saved.join(', ')}`, 'ok');
    refreshFileList();
  } catch (err) {
    setStatus(statusEl, `❌ Upload failed: ${err.message}`, 'error');
  }
}

// ---------------------------------------------------------------------------
// Drag-and-drop
// ---------------------------------------------------------------------------
function initDropZone() {
  const zone  = document.getElementById('obs-drop-zone');
  const input = document.getElementById('obs-file-input');
  if (!zone || !input) return;

  zone.addEventListener('dragover', e => { e.preventDefault(); zone.classList.add('drag-over'); });
  zone.addEventListener('dragleave', () => zone.classList.remove('drag-over'));
  zone.addEventListener('drop', e => {
    e.preventDefault();
    zone.classList.remove('drag-over');
    uploadFiles(Array.from(e.dataTransfer.files));
  });
  zone.addEventListener('click', e => {
    if (e.target.tagName !== 'LABEL') input.click();
  });
  input.addEventListener('change', () => {
    uploadFiles(Array.from(input.files));
    input.value = '';
  });
}

// ---------------------------------------------------------------------------
// Observation entries table
// ---------------------------------------------------------------------------
let _obsRowCounter = 0;

/**
 * Build a single <tr> for one observation entry.
 * @param {{filename:string, jdate:number, vel_center_kms:number, renormalize_wings:boolean}|null} entry
 */
function buildObsRow(entry) {
  const idx = ++_obsRowCounter;
  const tr = document.createElement('tr');
  tr.dataset.rowId = idx;

  const fn = entry ? escHtml(String(entry.filename))  : '';
  const jd = entry != null ? entry.jdate              : '';
  const rv = entry != null ? entry.vel_center_kms      : '';
  const rw = entry != null ? entry.renormalize_wings   : false;

  tr.innerHTML = `
    <td class="obs-row-num"></td>
    <td>
      <input type="text"
             class="field-input obs-fn"
             value="${fn}"
             list="obs-file-datalist"
             placeholder="LSDprof/file.prof"
             spellcheck="false"
             autocomplete="off" />
    </td>
    <td>
      <input type="number"
             class="field-input obs-jdate"
             value="${jd}"
             step="0.00001"
             placeholder="2456886.394" />
    </td>
    <td>
      <input type="number"
             class="field-input obs-rv"
             value="${rv}"
             step="0.1"
             placeholder="-19.8" />
    </td>
    <td style="text-align: center;">
      <input type="checkbox"
             class="obs-rw"
             ${rw ? 'checked' : ''}
             title="Renormalize Stokes I wings for this entry" />
    </td>
    <td class="obs-file-status obs-file-unknown" title="Not validated">—</td>
    <td>
      <button class="btn btn-danger btn-sm obs-del-row" title="Remove this observation">✕</button>
    </td>`;

  tr.querySelector('.obs-fn').addEventListener('change', () => _clearRowStatus(tr));
  return tr;
}

function _clearRowStatus(tr) {
  const cell = tr.querySelector('.obs-file-status');
  if (!cell) return;
  cell.className = 'obs-file-status obs-file-unknown';
  cell.textContent = '—';
  cell.title = 'Not validated';
}

function _renumberRows() {
  const rows = document.querySelectorAll('#obs-entries-tbody tr');
  rows.forEach((tr, i) => {
    const numCell = tr.querySelector('.obs-row-num');
    if (numCell) numCell.textContent = i + 1;
  });
  const noMsg = document.getElementById('obs-no-entries');
  if (noMsg) noMsg.style.display = rows.length ? 'none' : '';
}

// Delegated delete-row handler
document.getElementById('obs-entries-tbody')?.addEventListener('click', e => {
  const btn = e.target.closest('.obs-del-row');
  if (!btn) return;
  btn.closest('tr').remove();
  _renumberRows();
});

// ---------------------------------------------------------------------------
// Load observations from config into table
// ---------------------------------------------------------------------------
async function loadObsEditor() {
  const tbody   = document.getElementById('obs-entries-tbody');
  const jdateEl = document.getElementById('obs-jdate-ref');
  if (!tbody) return;

  try {
    const cfg = await apiFetch('/api/config');
    const obs = cfg.observations || {};
    if (jdateEl) jdateEl.value = obs.jdate_ref ?? 2456892.015;

    tbody.innerHTML = '';
    _obsRowCounter = 0;
    (obs.files || []).forEach(entry => tbody.appendChild(buildObsRow(entry)));
    _renumberRows();

    // Auto-validate after loading
    if ((obs.files || []).length) validateObsFiles();
  } catch (err) {
    setStatus(document.getElementById('obs-status'), `Error loading config: ${err.message}`, 'error');
  }
}

// ---------------------------------------------------------------------------
// File-existence validation
// ---------------------------------------------------------------------------
async function validateObsFiles() {
  const tbody    = document.getElementById('obs-entries-tbody');
  const statusEl = document.getElementById('obs-status');
  if (!tbody) return;

  const rows  = Array.from(tbody.querySelectorAll('tr'));
  const paths = rows.map(tr => tr.querySelector('.obs-fn')?.value.trim() ?? '');

  if (!paths.length) return;

  try {
    setStatus(statusEl, 'Validating file paths…');
    const result = await apiFetch('/api/observations/validate', {
      method: 'POST',
      body: { paths }
    });

    let missing = 0;
    rows.forEach((tr, i) => {
      const cell = tr.querySelector('.obs-file-status');
      if (!cell) return;
      const path = paths[i];
      const info = result[path];
      if (!info) {
        cell.className = 'obs-file-status obs-file-unknown';
        cell.textContent = '—';
        cell.title = 'Unknown';
        return;
      }
      if (info.exists) {
        cell.className = 'obs-file-status obs-file-ok';
        cell.textContent = '✓';
        cell.title = `File found · Stokes: ${info.stokes}`;
      } else {
        cell.className = 'obs-file-status obs-file-missing';
        cell.textContent = '✗';
        cell.title = 'File not found on server';
        missing++;
      }
    });

    if (missing) {
      setStatus(statusEl, `⚠ ${missing} file(s) not found on server.`, 'warn');
    } else {
      setStatus(statusEl, `✔ All ${rows.length} file(s) found.`, 'ok');
    }
  } catch (err) {
    setStatus(statusEl, `❌ Validation error: ${err.message}`, 'error');
  }
}

// ---------------------------------------------------------------------------
// Save observations to config
// ---------------------------------------------------------------------------
async function saveObservations() {
  const statusEl = document.getElementById('obs-status');
  const tbody    = document.getElementById('obs-entries-tbody');
  const jdateEl  = document.getElementById('obs-jdate-ref');
  if (!tbody) return;

  const rows  = Array.from(tbody.querySelectorAll('tr'));
  const files = [];
  let valid   = true;

  rows.forEach((tr, i) => {
    if (!valid) return;
    const fn = tr.querySelector('.obs-fn')?.value.trim()   ?? '';
    const jd = parseFloat(tr.querySelector('.obs-jdate')?.value ?? '');
    const rv = parseFloat(tr.querySelector('.obs-rv')?.value    ?? '');
    const rw = tr.querySelector('.obs-rw')?.checked ?? false;

    if (!fn) {
      setStatus(statusEl, `❌ Row ${i + 1}: filename is required.`, 'error');
      valid = false;
      return;
    }
    if (isNaN(jd)) {
      setStatus(statusEl, `❌ Row ${i + 1}: Julian Date must be a number.`, 'error');
      valid = false;
      return;
    }
    if (isNaN(rv)) {
      setStatus(statusEl, `❌ Row ${i + 1}: RV center must be a number.`, 'error');
      valid = false;
      return;
    }
    files.push({ filename: fn, jdate: jd, vel_center_kms: rv, renormalize_wings: rw });
  });

  if (!valid) return;

  try {
    setStatus(statusEl, 'Saving…');
    const cfg = await apiFetch('/api/config');
    cfg.observations = {
      ...(cfg.observations || {}),
      jdate_ref: parseFloat(jdateEl?.value ?? cfg.observations?.jdate_ref ?? 0),
      files,
    };
    await apiFetch('/api/config', { method: 'PUT', body: cfg });
    setStatus(statusEl, `✔ Saved ${files.length} observation(s).`, 'ok');
    if (files.length) validateObsFiles();
  } catch (err) {
    setStatus(statusEl, `❌ ${err.message}`, 'error');
  }
}

// ---------------------------------------------------------------------------
// Utility
// ---------------------------------------------------------------------------
function escHtml(str) {
  return String(str)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;');
}

// ---------------------------------------------------------------------------
// Init
// ---------------------------------------------------------------------------
document.addEventListener('DOMContentLoaded', () => {
  initDropZone();
  refreshFileList();
  loadObsEditor();

  document.getElementById('obs-refresh-btn')?.addEventListener('click', () => {
    refreshFileList();
    loadObsEditor();
  });

  document.getElementById('obs-add-btn')?.addEventListener('click', () => {
    const tbody = document.getElementById('obs-entries-tbody');
    if (!tbody) return;
    tbody.appendChild(buildObsRow(null));
    _renumberRows();
  });

  document.getElementById('obs-validate-btn')?.addEventListener('click', validateObsFiles);
  document.getElementById('obs-save-btn')?.addEventListener('click', saveObservations);
});
