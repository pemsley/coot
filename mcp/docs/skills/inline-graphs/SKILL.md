---
name: coot-inline-graphs
description: >
  Create interactive inline Chart.js graphs directly in the chat from live Coot data.
  Use this skill whenever the user asks to plot, graph, chart, or visualise any
  per-residue data from Coot — B-factors, density correlations, Ramachandran
  probabilities, rotamer scores, or any other per-residue metric. Also use when
  the user asks to overlay secondary structure on a graph, or to compare metrics
  across chains. Prefer this approach over any file-based graphing (e.g. Pygal)
  — it is faster, interactive, and renders inline in the conversation.
---

# Coot Inline Graphs

Inline graphs render Chart.js directly in the chat via the `visualize:show_widget`
tool. Coot supplies the data via Python; the widget renders it with no file I/O,
no external viewer, and full interactivity.

## Core workflow

1. **Get data from Coot** — fetch per-residue metrics using the Python API
2. **Get secondary structure** — call `add_header_secondary_structure_info()` then
   `get_header_secondary_structure_info()` if overlays are wanted
3. **Render the widget** — embed data as JS literals in the Chart.js HTML

Always call `visualize:read_me` (modules: `["interactive", "chart"]`) before the
first `visualize:show_widget` call in a session.

---

## Step 1 — Fetch per-residue data from Coot

### B-factors

```python
def get_bfactor_data(imol, chain_id):
    min_res = coot.min_resno_in_chain(imol, chain_id)
    max_res = coot.max_resno_in_chain(imol, chain_id)
    results = []
    for resno in range(min_res, max_res + 1):
        atoms = coot.residue_info_py(imol, chain_id, resno, "")
        if atoms:
            resname = coot.residue_name_py(imol, chain_id, resno, "")
            bfactors = [a[1][1] for a in atoms if isinstance(a[1][1], float)]
            mean_b = round(sum(bfactors) / len(bfactors), 2) if bfactors else 0
            results.append({"resno": resno, "resname": resname, "mean_b": mean_b})
    return results
```

### Density correlation

```python
def get_correlation_data(imol, chain_id, imol_map):
    stats = coot.map_to_model_correlation_stats_per_residue_range_py(
        imol, chain_id, imol_map, 1, 0)
    results = []
    for entry in stats[0]:
        residue_spec = entry[0]   # [chain_id, resno, ins_code]
        corr_data    = entry[1]   # [n_points, correlation]
        resno = residue_spec[1]
        correlation = corr_data[1]
        resname = coot.residue_name_py(imol, chain_id, resno, "")
        results.append({
            "resno": resno,
            "resname": resname,
            "correlation": round(correlation, 4) if correlation == correlation else None
        })
    return results
```

### Ramachandran probabilities

```python
def get_rama_data(imol, chain_id):
    rama = coot.all_molecule_ramachandran_score_py(imol)
    results = []
    for entry in rama[5]:
        if entry == -1:
            continue
        phi_psi, res_spec, score, res_names = entry
        if res_spec[0] != chain_id:
            continue
        results.append({
            "resno": res_spec[1],
            "resname": res_names[1],
            "phi": round(phi_psi[0], 1),
            "psi": round(phi_psi[1], 1),
            "rama_prob": round(score, 4)
        })
    return results
```

---

## Step 2 — Fetch secondary structure

Always try `get_header_secondary_structure_info()` first. If it returns `{}` or
`False`, call `add_header_secondary_structure_info()` to compute it from geometry,
then call `get_header_secondary_structure_info()` again.

```python
def get_secondary_structure(imol, chain_id):
    ss = coot.get_header_secondary_structure_info(imol)
    if not isinstance(ss, dict) or (not ss.get('helices') and not ss.get('strands')):
        coot.add_header_secondary_structure_info(imol)
        ss = coot.get_header_secondary_structure_info(imol)
    if not isinstance(ss, dict):
        return {'helices': [], 'strands': []}
    helices = [h for h in (ss.get('helices') or []) if h['initChainID'] == chain_id]
    strands = [s for s in (ss.get('strands') or []) if s['initChainID'] == chain_id]
    return {'helices': helices, 'strands': strands}
```

**Important:** `add_header_secondary_structure_info()` will crash Coot if called
on a molecule that already has secondary structure records populated and then
`get_header_secondary_structure_info()` is called — only call it when the initial
query returns empty. (Bug reported; fix applied to `c-interface-build.cc:2876`.)

---

## Step 3 — Render the widget

### Chart.js setup

Load via CDN. Always use the UMD build:
```html
<script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/4.4.1/chart.umd.js"></script>
```

For secondary structure annotation overlays, also load:
```html
<script src="https://cdnjs.cloudflare.com/ajax/libs/chartjs-plugin-annotation/3.0.1/chartjs-plugin-annotation.min.js"></script>
```

### Data embedding

Embed Coot data as a JS literal directly in the widget HTML. Do not use fetch()
or external URLs — the data comes from Coot at render time and is baked in.

```javascript
const data = [
  {"resno": 1, "resname": "ASP", "mean_b": 34.95},
  // ... all residues
];
```

### Canvas sizing

Always wrap `<canvas>` in a `<div>` with explicit height:
```html
<div style="position: relative; width: 100%; height: 300px;">
  <canvas id="chart"></canvas>
</div>
```

Set `responsive: true, maintainAspectRatio: false` in Chart.js options.
Never set height directly on the `<canvas>` element.

---

## Secondary structure overlay

Box annotations sit at the **top** of the chart as a strip. The box height is
computed dynamically so the α/β glyph sits vertically centred:

```javascript
const boxHeightUnits = Math.round(22 * yAxisMax / 280);
const boxYMax = yAxisMax;
const boxYMin = yAxisMax - boxHeightUnits;
```

Build a `resnoToIndex` lookup first (maps residue number → bar index):
```javascript
const resnoToIndex = {};
data.forEach((d, i) => { resnoToIndex[d.resno] = i; });
```

### Annotation spec

```javascript
// Helix — purple, semi-opaque, white-ish glyph text
{
  type: 'box',
  xMin: resnoToIndex[h.initSeqNum] - 0.5,
  xMax: resnoToIndex[h.endSeqNum]  + 0.5,
  yMin: boxYMin,
  yMax: boxYMax,
  backgroundColor: 'rgba(175,169,236,0.45)',
  borderColor:     'rgba(127,119,221,0.8)',
  borderWidth: 1,
  label: {
    display: true,
    content: 'α',
    position: { x: 'center', y: 'center' },
    font:  { size: 13, weight: '500' },
    color: 'rgba(255,255,255,0.85)'
  }
}

// Strand — amber, semi-opaque, white-ish glyph text
{
  type: 'box',
  xMin: resnoToIndex[s.initSeqNum] - 0.5,
  xMax: resnoToIndex[s.endSeqNum]  + 0.5,
  yMin: boxYMin,
  yMax: boxYMax,
  backgroundColor: 'rgba(239,159,39,0.35)',
  borderColor:     'rgba(186,117,23,0.7)',
  borderWidth: 1,
  label: {
    display: true,
    content: 'β',
    position: { x: 'center', y: 'center' },
    font:  { size: 13, weight: '500' },
    color: 'rgba(255,255,255,0.85)'
  }
}
```

---

## Threshold colouring

Colour bars relative to a threshold to highlight problem residues:

```javascript
// Correlation — low is bad
backgroundColor: data.map(d => d.correlation < thresh ? '#378ADD' : '#5DCAA5')

// B-factor — high is bad
backgroundColor: data.map(d => d.mean_b > thresh ? '#378ADD' : '#5DCAA5')

// Ramachandran — low probability is bad
backgroundColor: data.map(d => d.rama_prob < thresh ? '#E24B4A' : '#5DCAA5')
```

Provide a range slider to let the user adjust threshold interactively.
When switching between metrics, update the slider range accordingly:
- Correlation: min=0, max=1, step=0.01, default=0.7
- B-factor: min=0, max=`bMax`, step=1, default=20
- Ramachandran: min=0, max=1, step=0.01, default=0.02

---

## Click-to-navigate

Wire bar clicks to `sendPrompt()` so the user can jump to a residue in Coot:

```javascript
onClick: (e, els) => {
  if (els.length) {
    const d = data[els[0].index];
    sendPrompt('Navigate to residue ' + d.resno + ' ' + d.resname +
               ' in chain ' + chainId + ' of the tutorial model');
  }
}
```

---

## Axis labels and ticks

```javascript
scales: {
  x: {
    grid: { display: false },
    ticks: {
      color: '#888780',
      font: { size: 9 },
      maxRotation: 90,
      autoSkip: true,
      maxTicksLimit: 30
    }
  },
  y: {
    min: 0,
    max: yAxisMax,
    grid: { color: 'rgba(136,135,128,0.15)' },
    ticks: {
      color: '#888780',
      font: { size: 11 },
      callback: v => v + ' Å²'   // or '.toFixed(2)' for correlations
    }
  }
}
```

---

## Stat cards

Show summary metrics above the chart using the metric card pattern:

```html
<div style="background: var(--color-background-secondary);
            border-radius: var(--border-radius-md);
            padding: 10px 12px;">
  <div style="font-size: 11px; color: var(--color-text-secondary);">Mean B</div>
  <div style="font-size: 17px; font-weight: 500; color: var(--color-text-primary);"
       id="s-meanb">—</div>
</div>
```

Use a 4-column grid: residue count, mean metric, count above/below threshold,
max or min value as appropriate.

---

## Legend

Always provide a manual legend below the chart — do not use Chart.js default:

```html
<div style="display: flex; gap: 16px; margin-top: 8px;
            font-size: 12px; color: var(--color-text-secondary); flex-wrap: wrap;">
  <span style="display:flex;align-items:center;gap:4px;">
    <span style="width:10px;height:10px;border-radius:2px;background:#5DCAA5;"></span>
    Below threshold
  </span>
  <span style="display:flex;align-items:center;gap:4px;">
    <span style="width:10px;height:10px;border-radius:2px;background:#378ADD;"></span>
    Above threshold
  </span>
  <span style="display:flex;align-items:center;gap:4px;">
    <span style="width:10px;height:10px;border-radius:2px;
                 background:rgba(175,169,236,0.45);border:1px solid #7F77DD;"></span>
    Helix
  </span>
  <span style="display:flex;align-items:center;gap:4px;">
    <span style="width:10px;height:10px;border-radius:2px;
                 background:rgba(239,159,39,0.35);border:1px solid #BA7517;"></span>
    Strand
  </span>
</div>
```

---

## Tooltips

Include both the primary metric and secondary metric in tooltips:

```javascript
tooltip: {
  callbacks: {
    title: items => items[0].label,
    label: item => 'Mean B: ' + data[item.dataIndex].mean_b.toFixed(1) + ' Å²',
    afterLabel: item => {
      const r = data[item.dataIndex].resno;
      if (helices.some(h => r >= h.initSeqNum && r <= h.endSeqNum)) return 'α-helix';
      if (strands.some(s => r >= s.initSeqNum && r <= s.endSeqNum)) return 'β-strand';
      return 'loop/coil';
    }
  }
}
```

---

## Number formatting

All numbers reaching the screen must be rounded:
- B-factors: `.toFixed(1)` + `' Å²'`
- Correlations: `.toFixed(3)`
- Ramachandran probabilities: `.toFixed(4)`
- Axis tick integers: `Math.round()`

---

## Why not Pygal?

Pygal requires file I/O, a separate viewer, and a display context. It produces
black images in headless environments and is slow. Chart.js in the browser has
none of these problems and adds interactivity for free. Do not use Pygal.
