/**
 * results.js — Render ZDI results using Plotly.js.
 *
 * Primary endpoints (return Plotly JSON directly):
 *   GET /api/plots/profiles                    → ProfilePlotData JSON
 *   GET /api/plots/surface_map?map_type=...    → SurfaceMapData JSON
 *   GET /api/plots/light_curve                 → LightCurvePlotData JSON
 *
 * Legacy fallback endpoints (used when primary endpoints unavailable):
 *   GET /api/results/profiles    → raw phase data
 *   GET /api/results/magnetic    → Br, Bclat, Blon on stellar surface
 *   GET /api/results/brightness  → brightness map
 */

// ---------------------------------------------------------------------------
// Shared Plotly layout helpers
// ---------------------------------------------------------------------------
function _darkLayout(overrides = {}) {
  const layout = {
    ...PLOTLY_BASE_LAYOUT,
    ...overrides,
    xaxis: { ...PLOTLY_BASE_LAYOUT.xaxis, ...(overrides.xaxis || {}) },
    yaxis: { ...PLOTLY_BASE_LAYOUT.yaxis, ...(overrides.yaxis || {}) },
  };
  // Merge yaxis2 if present in overrides (PLOTLY_BASE_LAYOUT has no yaxis2)
  if (overrides.yaxis2) {
    layout.yaxis2 = { ...PLOTLY_BASE_LAYOUT.yaxis, ...overrides.yaxis2 };
  }
  return layout;
}

// Merge a Plotly-JSON layout returned by the backend into the dark theme
function _mergeLayout(backendLayout) {
  return {
    ...backendLayout,
    paper_bgcolor: PLOTLY_BASE_LAYOUT.paper_bgcolor,
    plot_bgcolor:  PLOTLY_BASE_LAYOUT.plot_bgcolor,
    font: { ...PLOTLY_BASE_LAYOUT.font, ...(backendLayout.font || {}) },
  };
}

// ---------------------------------------------------------------------------
// Line profiles  (primary: /api/plots/profiles, fallback: /api/results/profiles)
// ---------------------------------------------------------------------------
async function loadProfiles() {
  const container = document.getElementById('profiles-container');
  const noData    = document.getElementById('profiles-no-data');
  if (!container) return;

  // ── Try primary endpoint ──
  try {
    const plotJson = await apiFetch('/api/plots/profiles');
    noData && (noData.style.display = 'none');
    container.innerHTML = '';

    // Group traces by phase index (legendgroup = "phase0", "phase1", ...)
    const phaseMap = {};
    for (const trace of plotJson.data) {
      const lg  = trace.legendgroup || 'phase0';
      const idx = parseInt(lg.replace('phase', ''), 10);
      if (!phaseMap[idx]) phaseMap[idx] = [];
      phaseMap[idx].push(trace);
    }

    const phaseKeys = Object.keys(phaseMap).map(Number).sort((a, b) => a - b);
    if (!phaseKeys.length) {
      container.style.display = 'none';
      noData && (noData.style.display = '');
      return;
    }

    const nCols = Math.min(4, phaseKeys.length);
    container.style.display = 'grid';
    container.style.gridTemplateColumns = `repeat(${nCols}, 1fr)`;

    // Shared per-phase layout: dual-axis (I on top, V on bottom)
    const phaseLayout = {
      margin: { t: 24, r: 6, b: 28, l: 50 },
      legend: {
        x: 0.02, y: 0.99, xanchor: 'left', yanchor: 'top',
        bgcolor: 'rgba(0,0,0,0)', font: { size: 7 }, orientation: 'h',
      },
      xaxis: {
        title: { text: 'v (km/s)', font: { size: 9 } },
        tickfont: { size: 8 },
      },
      yaxis: {
        title: { text: 'I/Ic', font: { size: 9 } },
        tickfont: { size: 8 },
        domain: [0.46, 1.0],
      },
      yaxis2: {
        title: { text: 'V/Ic', font: { size: 9 } },
        tickfont: { size: 8 },
        domain: [0.0, 0.40],
        anchor: 'x',
        zeroline: true,
        zerolinecolor: '#484f58',
        zerolinewidth: 1,
      },
    };

    phaseKeys.forEach((phIdx, i) => {
      const traces = phaseMap[phIdx];
      const card = document.createElement('div');
      card.className = 'phase-card';

      // Extract phase label from the first I trace name (e.g. "I obs  φ=1.230")
      const firstName = traces.find(t => !t.yaxis)?.name || traces[0]?.name || '';
      const phaseLabel = firstName.match(/φ=(.+)/)?.[1] ?? String(i + 1);

      card.innerHTML = `
        <div class="phase-title">φ = ${phaseLabel}</div>
        <div id="pp-${i}" style="height:260px"></div>`;
      container.appendChild(card);

      // Force all traces visible; show legend only for I-axis traces (not V)
      const phTraces = traces.map(t => ({
        ...t,
        visible: true,
        showlegend: t.yaxis !== 'y2',
        name: t.name ? t.name.replace(/\s*\u03c6=.*$/, '').trim() : t.name,
      }));

      Plotly.newPlot(`pp-${i}`, phTraces, _darkLayout(phaseLayout), PLOTLY_CONFIG);
    });
    return;
  } catch (_) { /* fall through to legacy */ }

  // ── Legacy fallback ──
  let data;
  try {
    data = await apiFetch('/api/results/profiles');
  } catch (err) {
    setStatus(document.getElementById('results-status'), `Profiles error: ${err.message}`, 'error');
    return;
  }

  if (!data.available || !data.phases.length) {
    container.style.display = 'none';
    noData && (noData.style.display = '');
    return;
  }
  noData && (noData.style.display = 'none');
  container.style.display = 'grid';
  container.innerHTML = '';

  const nCols = Math.min(4, data.phases.length);
  container.style.gridTemplateColumns = `repeat(${nCols}, 1fr)`;

  data.phases.forEach((phase, i) => {
    const card = document.createElement('div');
    card.className = 'phase-card';
    card.innerHTML = `
      <div class="phase-title">Phase ${i + 1}</div>
      <div id="pI-${i}" style="height:140px"></div>
      <div id="pV-${i}" style="height:110px"></div>`;
    container.appendChild(card);

    const tracesI = [];
    if (phase.stokes_I_obs) {
      tracesI.push({
        x: phase.vel_obs, y: phase.stokes_I_obs,
        error_y: phase.stokes_I_err
          ? { type: 'data', array: phase.stokes_I_err, visible: true, color: '#484f58', thickness: 1 }
          : undefined,
        mode: 'markers', marker: { color: '#f0f6fc', size: 2.5 },
        name: 'Obs', showlegend: i === 0,
      });
    }
    if (phase.stokes_I_mod) {
      tracesI.push({
        x: phase.vel_mod, y: phase.stokes_I_mod,
        mode: 'lines', line: { color: '#f85149', width: 1.5 },
        name: 'Model', showlegend: i === 0,
      });
    }
    Plotly.newPlot(`pI-${i}`, tracesI, _darkLayout({
      margin: { t: 4, r: 6, b: 24, l: 44 },
      showlegend: i === 0, legend: { x: 0, y: 1, font: { size: 8 } },
      xaxis: { title: { text: 'v (km/s)', font: { size: 9 } }, tickfont: { size: 8 } },
      yaxis: { title: { text: 'I/Ic',     font: { size: 9 } }, tickfont: { size: 8 } },
    }), PLOTLY_CONFIG);

    const tracesV = [];
    // Zero reference drawn FIRST so V_obs and V_mod render on top of it
    if (phase.vel_obs?.length || phase.vel_mod?.length) {
      const xz = phase.vel_obs?.length ? phase.vel_obs : phase.vel_mod;
      tracesV.push({
        x: [xz[0], xz[xz.length - 1]], y: [0, 0],
        mode: 'lines', line: { color: '#30363d', width: 0.8, dash: 'dash' },
        showlegend: false, hoverinfo: 'skip',
      });
    }
    if (phase.stokes_V_obs) {
      tracesV.push({
        x: phase.vel_obs, y: phase.stokes_V_obs,
        error_y: phase.stokes_V_err
          ? { type: 'data', array: phase.stokes_V_err, visible: true, color: '#484f58', thickness: 1 }
          : undefined,
        mode: 'markers', marker: { color: '#8b949e', size: 2.5 },
        name: 'Obs V', showlegend: false,
      });
    }
    if (phase.stokes_V_mod) {
      tracesV.push({
        x: phase.vel_mod, y: phase.stokes_V_mod,
        mode: 'lines', line: { color: '#f85149', width: 1.5 },
        name: 'Model V', showlegend: false,
      });
    }
    Plotly.newPlot(`pV-${i}`, tracesV, _darkLayout({
      margin: { t: 4, r: 6, b: 28, l: 44 },
      xaxis: { title: { text: 'v (km/s)', font: { size: 9 } }, tickfont: { size: 8 } },
      yaxis: { title: { text: 'V/Ic', font: { size: 9 } }, tickfont: { size: 8 }, zeroline: true, zerolinecolor: '#30363d' },
    }), PLOTLY_CONFIG);
  });
}

// ---------------------------------------------------------------------------
// Magnetic map  (primary: /api/plots/surface_map?map_type=..., fallback: /api/results/magnetic)
// ---------------------------------------------------------------------------
async function loadMagnetic(mapType = 'radial_B') {
  const chartEl = document.getElementById('magnetic-chart');
  const noData  = document.getElementById('magnetic-no-data');
  if (!chartEl) return;

  // ── Try primary endpoint ──
  try {
    const plotJson = await apiFetch(`/api/plots/surface_map?map_type=${mapType}`);
    noData && (noData.style.display = 'none');
    chartEl.style.display = '';
    Plotly.react(chartEl, plotJson.data, _mergeLayout({
      ...plotJson.layout, height: 360,
      margin: { t: 30, r: 20, b: 50, l: 55 },
    }), PLOTLY_CONFIG);
    return;
  } catch (err) {
    if (!err.message.includes('404')) {
      // 503 = grid coords unavailable, 501 = not yet implemented — show hint
      noData && (noData.style.display = '');
      const hint = err.message.includes('503') ? 'Re-run the inversion to cache grid coordinates.' :
                   err.message.includes('501') ? 'Magnetic map requires updated pipeline. Re-run.' :
                   err.message;
      noData && (noData.textContent = hint);
      chartEl.style.display = 'none';
      return;
    }
    /* 404 = no result yet — fall through to legacy */
  }

  // ── Legacy fallback ──
  let data;
  try {
    data = await apiFetch('/api/results/magnetic');
  } catch (err) {
    setStatus(document.getElementById('results-status'), `Magnetic map error: ${err.message}`, 'error');
    return;
  }

  if (!data.available) {
    chartEl.style.display = 'none';
    noData && (noData.style.display = '');
    if (data.error) setStatus(document.getElementById('results-status'), `Magnetic map: ${data.error}`, 'warn');
    return;
  }
  noData && (noData.style.display = 'none');
  chartEl.style.display = '';

  const componentKey = mapType === 'radial_B' ? 'Br' : mapType === 'meridional_B' ? 'Bclat' : 'Blon';
  const vals = data[componentKey] || data.Br || [];
  const absMax = Math.max(...vals.map(Math.abs), 1);

  Plotly.react(chartEl, [{
    type: 'scatter', mode: 'markers',
    x: data.lon, y: data.lat,
    marker: {
      color: vals, colorscale: 'RdBu', reversescale: true,
      cmin: -absMax, cmax: absMax, size: 3,
      colorbar: { title: mapType, titlefont: { size: 10 }, tickfont: { size: 9 } },
    },
    showlegend: false, name: mapType,
  }], _darkLayout({
    margin: { t: 30, r: 20, b: 50, l: 55 }, height: 360,
    xaxis: { title: 'Longitude (°)' }, yaxis: { title: 'Latitude (°)' },
  }), PLOTLY_CONFIG);
}

// ---------------------------------------------------------------------------
// Brightness map  (primary: /api/plots/surface_map?map_type=brightness, fallback: /api/results/brightness)
// ---------------------------------------------------------------------------
async function loadBrightness() {
  const chartEl = document.getElementById('brightness-chart');
  const noData  = document.getElementById('brightness-no-data');
  if (!chartEl) return;

  // ── Try primary endpoint ──
  try {
    const plotJson = await apiFetch('/api/plots/surface_map?map_type=brightness');
    noData && (noData.style.display = 'none');
    chartEl.style.display = '';
    Plotly.react(chartEl, plotJson.data, _mergeLayout({
      ...plotJson.layout, height: 320,
      margin: { t: 30, r: 20, b: 50, l: 55 },
    }), PLOTLY_CONFIG);
    return;
  } catch (err) {
    if (!err.message.includes('404')) {
      noData && (noData.style.display = '');
      chartEl.style.display = 'none';
      return;
    }
  }

  // ── Legacy fallback ──
  let data;
  try {
    data = await apiFetch('/api/results/brightness');
  } catch (err) {
    setStatus(document.getElementById('results-status'), `Brightness error: ${err.message}`, 'error');
    return;
  }

  if (!data.available) {
    chartEl.style.display = 'none';
    noData && (noData.style.display = '');
    return;
  }
  noData && (noData.style.display = 'none');
  chartEl.style.display = '';

  const maxBright = Math.max(...data.brightness, 1.5);
  Plotly.react(chartEl, [{
    type: 'scatter', mode: 'markers',
    x: data.lon, y: data.lat,
    marker: {
      color: data.brightness, colorscale: 'Hot', reversescale: true,
      cmin: 0, cmax: maxBright, size: 3.5,
      colorbar: { title: 'Brightness', titlefont: { size: 10 }, tickfont: { size: 9 } },
    },
    showlegend: false,
    hovertemplate: 'lon=%{x:.1f}°  lat=%{y:.1f}°  bri=%{marker.color:.3f}<extra></extra>',
  }], _darkLayout({
    margin: { t: 30, r: 20, b: 50, l: 55 }, height: 320,
    xaxis: { title: 'Longitude (°)' }, yaxis: { title: 'Latitude (°)' },
  }), PLOTLY_CONFIG);
}

// ---------------------------------------------------------------------------
// Light curve  (primary: /api/plots/light_curve; hidden when 404)
// ---------------------------------------------------------------------------
async function loadLightCurve() {
  const card   = document.getElementById('lc-card');
  const chartEl = document.getElementById('lc-chart');
  if (!card || !chartEl) return;

  try {
    const plotJson = await apiFetch('/api/plots/light_curve');
    card.style.display = '';
    Plotly.react(chartEl, plotJson.data, _mergeLayout({
      ...plotJson.layout, height: 260,
      margin: { t: 30, r: 20, b: 50, l: 55 },
    }), PLOTLY_CONFIG);
  } catch (_) {
    // 404 = no light curve data; silently hide the card
    card.style.display = 'none';
  }
}

// ---------------------------------------------------------------------------
// Load all results (called by "↻ Load results" button and after SSE 'result' event)
// ---------------------------------------------------------------------------
async function loadAllResults() {
  const statusEl = document.getElementById('results-status');
  setStatus(statusEl, 'Loading results…');
  await Promise.allSettled([
    loadProfiles(),
    loadBrightness(),
    loadMagnetic('radial_B'),
    loadLightCurve(),
  ]);
  setStatus(statusEl, '');
}

// Expose for external callers (e.g. run.js after SSE 'result' event)
window.loadAllResults = loadAllResults;

// ---------------------------------------------------------------------------
// Magnetic sub-tab button handling
// ---------------------------------------------------------------------------
document.querySelectorAll('.mag-tab-btn').forEach(btn => {
  btn.addEventListener('click', () => {
    document.querySelectorAll('.mag-tab-btn').forEach(b => b.classList.remove('active'));
    btn.classList.add('active');
    loadMagnetic(btn.dataset.mag);
  });
});

// ---------------------------------------------------------------------------
// Matplotlib polar magnetic map
// ---------------------------------------------------------------------------
async function loadMagneticPolar() {
  const btn = document.getElementById('magnetic-mpl-btn');
  const statusEl = document.getElementById('magnetic-mpl-status');
  const container = document.getElementById('magnetic-mpl-img-container');
  const img = document.getElementById('magnetic-mpl-img');
  if (!btn) return;

  btn.disabled = true;
  statusEl.textContent = '生成中…';
  statusEl.style.color = '#8b949e';

  try {
    const resp = await fetch('/api/plots/magnetic_polar');
    if (!resp.ok) {
      let detail = resp.statusText;
      try { const j = await resp.json(); detail = j.detail || detail; } catch (_) {}
      throw new Error(detail);
    }
    const blob = await resp.blob();
    const url = URL.createObjectURL(blob);
    // Revoke previous object URL to avoid memory leak
    if (img.src && img.src.startsWith('blob:')) {
      URL.revokeObjectURL(img.src);
    }
    img.src = url;
    container.style.display = '';
    statusEl.textContent = '';
  } catch (err) {
    statusEl.textContent = `错误: ${err.message}`;
    statusEl.style.color = '#f85149';
  } finally {
    btn.disabled = false;
  }
}

// ---------------------------------------------------------------------------
// Matplotlib polar brightness map
// ---------------------------------------------------------------------------
async function loadBrightnessPolar() {
  const btn = document.getElementById('brightness-mpl-btn');
  const statusEl = document.getElementById('brightness-mpl-status');
  const container = document.getElementById('brightness-mpl-img-container');
  const img = document.getElementById('brightness-mpl-img');
  if (!btn) return;

  btn.disabled = true;
  statusEl.textContent = '生成中…';
  statusEl.style.color = '#8b949e';

  try {
    const resp = await fetch('/api/plots/brightness_polar');
    if (!resp.ok) {
      let detail = resp.statusText;
      try { const j = await resp.json(); detail = j.detail || detail; } catch (_) {}
      throw new Error(detail);
    }
    const blob = await resp.blob();
    const url = URL.createObjectURL(blob);
    if (img.src && img.src.startsWith('blob:')) {
      URL.revokeObjectURL(img.src);
    }
    img.src = url;
    container.style.display = '';
    statusEl.textContent = '';
  } catch (err) {
    statusEl.textContent = `错误: ${err.message}`;
    statusEl.style.color = '#f85149';
  } finally {
    btn.disabled = false;
  }
}

// ---------------------------------------------------------------------------
// PFSS 3D 磁力线外推
// ---------------------------------------------------------------------------
async function loadPFSS3D() {
  const btn      = document.getElementById('pfss-run-btn');
  const statusEl = document.getElementById('pfss-status');
  const chartEl  = document.getElementById('pfss-chart');
  const noData   = document.getElementById('pfss-no-data');
  if (!btn || !chartEl) return;

  const rss       = parseFloat(document.getElementById('pfss-rss')?.value     ?? '2.5');
  const nrho      = parseInt(document.getElementById('pfss-nrho')?.value      ?? '30', 10);
  const latSeeds  = parseInt(document.getElementById('pfss-lat-seeds')?.value ?? '18', 10);
  const lonSeeds  = parseInt(document.getElementById('pfss-lon-seeds')?.value ?? '8',  10);

  btn.disabled = true;
  statusEl.textContent = '计算中，请稍候…';
  statusEl.style.color = '#8b949e';

  try {
    const params = new URLSearchParams({
      nrho, rss, n_lat_seeds: latSeeds, n_lon_seeds: lonSeeds,
    });
    const plotJson = await apiFetch(`/api/plots/pfss_3d?${params}`);

    noData && (noData.style.display = 'none');
    chartEl.style.display = '';

    const layout = {
      ...plotJson.layout,
      paper_bgcolor: '#0d1117',
      height: 520,
      margin: { t: 40, r: 10, b: 10, l: 10 },
    };
    Plotly.react(chartEl, plotJson.data, layout, { ...PLOTLY_CONFIG, responsive: true });
    statusEl.textContent = '完成';
    statusEl.style.color = '#3fb950';
  } catch (err) {
    statusEl.textContent = `错误: ${err.message}`;
    statusEl.style.color = '#f85149';
  } finally {
    btn.disabled = false;
  }
}

// ---------------------------------------------------------------------------
// Init
// ---------------------------------------------------------------------------
document.addEventListener('DOMContentLoaded', () => {
  document.getElementById('results-refresh-btn')?.addEventListener('click', loadAllResults);
  document.getElementById('magnetic-mpl-btn')?.addEventListener('click', loadMagneticPolar);
  document.getElementById('brightness-mpl-btn')?.addEventListener('click', loadBrightnessPolar);
  document.getElementById('pfss-run-btn')?.addEventListener('click', loadPFSS3D);
});
