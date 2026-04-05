/**
 * results.js — Render ZDI results using Plotly.js.
 *
 * Endpoints consumed:
 *   GET /api/results/profiles    → phase data (Stokes I & V)
 *   GET /api/results/magnetic    → Br, Bclat, Blon on stellar surface
 *   GET /api/results/brightness  → brightness map
 */

// ---------------------------------------------------------------------------
// Shared Plotly layout helpers
// ---------------------------------------------------------------------------
function _darkLayout(overrides = {}) {
  return {
    ...PLOTLY_BASE_LAYOUT,
    ...overrides,
    xaxis: { ...PLOTLY_BASE_LAYOUT.xaxis, ...(overrides.xaxis || {}) },
    yaxis: { ...PLOTLY_BASE_LAYOUT.yaxis, ...(overrides.yaxis || {}) },
  };
}

// ---------------------------------------------------------------------------
// Line profiles
// ---------------------------------------------------------------------------
async function loadProfiles() {
  const container = document.getElementById('profiles-container');
  const noData    = document.getElementById('profiles-no-data');
  if (!container) return;

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

  // Clear previous charts
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

    // ── Stokes I ──
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
      showlegend: i === 0,
      legend: { x: 0, y: 1, font: { size: 8 } },
      xaxis: { title: { text: 'v (km/s)', font: { size: 9 } }, tickfont: { size: 8 } },
      yaxis: { title: { text: 'I/Ic',     font: { size: 9 } }, tickfont: { size: 8 } },
    }), PLOTLY_CONFIG);

    // ── Stokes V ──
    const tracesV = [];
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
        mode: 'lines', line: { color: '#58a6ff', width: 1.5 },
        name: 'Model V', showlegend: false,
      });
    }
    // Zero line
    if (phase.vel_obs?.length || phase.vel_mod?.length) {
      const xz = phase.vel_obs?.length ? phase.vel_obs : phase.vel_mod;
      tracesV.push({
        x: [xz[0], xz[xz.length - 1]], y: [0, 0],
        mode: 'lines', line: { color: '#30363d', width: 0.8, dash: 'dash' },
        showlegend: false, hoverinfo: 'skip',
      });
    }
    Plotly.newPlot(`pV-${i}`, tracesV, _darkLayout({
      margin: { t: 4, r: 6, b: 28, l: 44 },
      xaxis: { title: { text: 'v (km/s)', font: { size: 9 } }, tickfont: { size: 8 } },
      yaxis: { title: { text: 'V/Ic',     font: { size: 9 } }, tickfont: { size: 8 }, zeroline: true, zerolinecolor: '#30363d' },
    }), PLOTLY_CONFIG);
  });
}

// ---------------------------------------------------------------------------
// Magnetic map
// ---------------------------------------------------------------------------
async function loadMagnetic() {
  const chartEl = document.getElementById('magnetic-chart');
  const noData  = document.getElementById('magnetic-no-data');
  if (!chartEl) return;

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
    if (data.error) {
      setStatus(document.getElementById('results-status'), `Magnetic map: ${data.error}`, 'warn');
    }
    return;
  }
  noData && (noData.style.display = 'none');
  chartEl.style.display = '';

  const absMaxBr   = Math.max(...data.Br.map(Math.abs), 1);
  const absMaxBcl  = Math.max(...data.Bclat.map(Math.abs), 1);
  const absMaxBlon = Math.max(...data.Blon.map(Math.abs), 1);

  const makeScatter = (y, z, absMax, label, domain) => ({
    type: 'scatter',
    mode: 'markers',
    x: data.lon, y,
    marker: {
      color: z, colorscale: 'RdBu', reversescale: true,
      cmin: -absMax, cmax: absMax, size: 3,
      colorbar: { title: label, titlefont: { size: 10 }, tickfont: { size: 9 }, x: domain[1] - 0.01 },
    },
    xaxis: `x${domain[2]}`, yaxis: `y${domain[2]}`,
    showlegend: false, name: label,
    hovertemplate: `lon=%{x:.1f}°  lat=%{y:.1f}°  ${label}=%{marker.color:.1f} G<extra></extra>`,
  });

  const traces = [
    makeScatter(data.lat, data.Br,    absMaxBr,   'Br (G)',    [0.00, 0.30, '']),
    makeScatter(data.lat, data.Bclat, absMaxBcl,  'Bθ (G)',   [0.35, 0.65, 2]),
    makeScatter(data.lat, data.Blon,  absMaxBlon, 'Bφ (G)',   [0.70, 1.00, 3]),
  ];

  const axisFragment = {
    gridcolor: '#30363d', zerolinecolor: '#484f58', color: '#8b949e', tickfont: { size: 9 },
  };

  const layout = {
    paper_bgcolor: '#21262d', plot_bgcolor: '#0d1117',
    font: { family: "'JetBrains Mono', monospace", color: '#f0f6fc', size: 10 },
    margin: { t: 30, r: 20, b: 50, l: 55 },
    title: { text: 'Surface magnetic field  (Gauss)', font: { size: 12 }, y: 0.97 },
    xaxis:  { ...axisFragment, domain: [0.00, 0.27], title: 'Longitude (°)' },
    yaxis:  { ...axisFragment, title: 'Latitude (°)', anchor: 'x' },
    xaxis2: { ...axisFragment, domain: [0.37, 0.63], title: 'Longitude (°)', anchor: 'y2' },
    yaxis2: { ...axisFragment, title: 'Latitude (°)', anchor: 'x2' },
    xaxis3: { ...axisFragment, domain: [0.73, 0.99], title: 'Longitude (°)', anchor: 'y3' },
    yaxis3: { ...axisFragment, title: 'Latitude (°)', anchor: 'x3' },
    height: 380,
  };

  Plotly.react(chartEl, traces, layout, PLOTLY_CONFIG);
}

// ---------------------------------------------------------------------------
// Brightness map
// ---------------------------------------------------------------------------
async function loadBrightness() {
  const chartEl = document.getElementById('brightness-chart');
  const noData  = document.getElementById('brightness-no-data');
  if (!chartEl) return;

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

  const traces = [{
    type: 'scatter',
    mode: 'markers',
    x: data.lon, y: data.lat,
    marker: {
      color: data.brightness, colorscale: 'Hot', reversescale: true,
      cmin: 0, cmax: maxBright, size: 3.5,
      colorbar: { title: 'Brightness', titlefont: { size: 10 }, tickfont: { size: 9 } },
    },
    showlegend: false,
    hovertemplate: 'lon=%{x:.1f}°  lat=%{y:.1f}°  bri=%{marker.color:.3f}<extra></extra>',
  }];

  const layout = _darkLayout({
    title: { text: 'Relative brightness distribution', font: { size: 12 }, y: 0.97 },
    margin: { t: 36, r: 20, b: 50, l: 55 },
    height: 330,
    xaxis: { title: 'Longitude (°)', tickfont: { size: 9 } },
    yaxis: { title: 'Latitude (°)',  tickfont: { size: 9 } },
  });

  Plotly.react(chartEl, traces, layout, PLOTLY_CONFIG);
}

// ---------------------------------------------------------------------------
// Refresh all
// ---------------------------------------------------------------------------
async function loadAllResults() {
  const statusEl = document.getElementById('results-status');
  setStatus(statusEl, 'Loading results…');
  try {
    await Promise.all([loadProfiles(), loadMagnetic(), loadBrightness()]);
    setStatus(statusEl, '✔ Results loaded.', 'ok');
    setTimeout(() => setStatus(statusEl, ''), 3000);
  } catch (err) {
    setStatus(statusEl, `Error: ${err.message}`, 'error');
  }
}

// ---------------------------------------------------------------------------
// Init
// ---------------------------------------------------------------------------
document.addEventListener('DOMContentLoaded', () => {
  document.getElementById('results-refresh-btn')?.addEventListener('click', loadAllResults);
});
