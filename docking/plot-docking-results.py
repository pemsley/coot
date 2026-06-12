#!/usr/bin/env python3
"""
Plot docking results from dock-results.csv

Usage:
   python plot-docking-results.py [dock-results.csv]

Generates:
   dock-energy-vs-distance.png  — scatter plot of total energy vs COM distance
   dock-energy-components.png   — stacked bar chart of energy components for top results
"""

import sys
import csv
import os

def read_csv(filename):
   data = {'rank': [], 'e_total': [], 'e_vdw': [], 'e_elec': [],
           'e_air': [], 'com_dist': []}
   with open(filename) as f:
      reader = csv.DictReader(f)
      for row in reader:
         for key in data:
            data[key].append(float(row[key]))
   return data


def write_energy_vs_distance_svg(data, filename):
   """Scatter plot: total energy vs distance from native COM."""

   n = len(data['e_total'])
   if n == 0:
      return

   # SVG dimensions and margins
   width = 700
   height = 500
   margin_l = 80
   margin_r = 30
   margin_t = 40
   margin_b = 60
   plot_w = width - margin_l - margin_r
   plot_h = height - margin_t - margin_b

   # Data ranges
   x_vals = data['com_dist']
   y_vals = data['e_total']
   x_min = 0
   x_max = max(x_vals) * 1.1
   y_min = min(y_vals) - abs(min(y_vals)) * 0.1
   y_max = max(y_vals) * 1.1
   if y_min > 0:
      y_min = 0

   def sx(x):
      return margin_l + (x - x_min) / (x_max - x_min) * plot_w

   def sy(y):
      return margin_t + plot_h - (y - y_min) / (y_max - y_min) * plot_h

   lines = []
   lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">')
   lines.append('<style>')
   lines.append('  text { font-family: Helvetica, Arial, sans-serif; }')
   lines.append('  .axis-label { font-size: 14px; fill: #333; }')
   lines.append('  .title { font-size: 16px; font-weight: bold; fill: #222; }')
   lines.append('  .tick-label { font-size: 11px; fill: #555; }')
   lines.append('</style>')

   # Background
   lines.append(f'<rect width="{width}" height="{height}" fill="white"/>')

   # Plot area background
   lines.append(f'<rect x="{margin_l}" y="{margin_t}" width="{plot_w}" '
                f'height="{plot_h}" fill="#f8f8f8" stroke="#ccc"/>')

   # Title
   lines.append(f'<text x="{width/2}" y="25" text-anchor="middle" class="title">'
                f'Rigid Body Docking: Energy vs Distance from Native</text>')

   # Grid lines and tick labels — Y axis
   n_y_ticks = 6
   for i in range(n_y_ticks + 1):
      y_val = y_min + (y_max - y_min) * i / n_y_ticks
      yp = sy(y_val)
      lines.append(f'<line x1="{margin_l}" y1="{yp}" x2="{margin_l + plot_w}" '
                   f'y2="{yp}" stroke="#ddd" stroke-width="0.5"/>')
      lines.append(f'<text x="{margin_l - 8}" y="{yp + 4}" text-anchor="end" '
                   f'class="tick-label">{y_val:.0f}</text>')

   # Grid lines and tick labels — X axis
   n_x_ticks = 6
   for i in range(n_x_ticks + 1):
      x_val = x_min + (x_max - x_min) * i / n_x_ticks
      xp = sx(x_val)
      lines.append(f'<line x1="{xp}" y1="{margin_t}" x2="{xp}" '
                   f'y2="{margin_t + plot_h}" stroke="#ddd" stroke-width="0.5"/>')
      lines.append(f'<text x="{xp}" y="{margin_t + plot_h + 18}" text-anchor="middle" '
                   f'class="tick-label">{x_val:.0f}</text>')

   # Zero line if visible
   if y_min < 0 < y_max:
      yp = sy(0)
      lines.append(f'<line x1="{margin_l}" y1="{yp}" x2="{margin_l + plot_w}" '
                   f'y2="{yp}" stroke="#999" stroke-width="1" stroke-dasharray="4,3"/>')

   # Data points — colour by rank (best = red, worst = blue)
   for i in range(n):
      xp = sx(x_vals[i])
      yp = sy(y_vals[i])
      # Colour gradient: rank 0 = red, higher ranks -> blue
      t = i / max(n - 1, 1)
      r = int(220 * (1 - t) + 70 * t)
      g = int(50 * (1 - t) + 130 * t)
      b = int(50 * (1 - t) + 200 * t)
      opacity = 0.8 if i < 5 else 0.5
      radius = 6 if i < 5 else 4

      lines.append(f'<circle cx="{xp:.1f}" cy="{yp:.1f}" r="{radius}" '
                   f'fill="rgb({r},{g},{b})" opacity="{opacity}" stroke="white" stroke-width="0.5"/>')

      # Label top 3
      if i < 3:
         lines.append(f'<text x="{xp + 8}" y="{yp + 4}" class="tick-label" '
                      f'font-weight="bold">#{i}</text>')

   # Axis labels
   lines.append(f'<text x="{width/2}" y="{height - 8}" text-anchor="middle" '
                f'class="axis-label">Distance from native COM (\u00c5)</text>')
   lines.append(f'<text x="18" y="{height/2}" text-anchor="middle" '
                f'class="axis-label" transform="rotate(-90,18,{height/2})">'
                f'Total energy (kcal/mol)</text>')

   # Legend
   lx = margin_l + plot_w - 160
   ly = margin_t + 15
   lines.append(f'<rect x="{lx}" y="{ly}" width="150" height="40" '
                f'fill="white" stroke="#ccc" rx="4"/>')
   lines.append(f'<circle cx="{lx+12}" cy="{ly+14}" r="5" fill="rgb(220,50,50)"/>')
   lines.append(f'<text x="{lx+22}" y="{ly+18}" class="tick-label">Best ranked</text>')
   lines.append(f'<circle cx="{lx+12}" cy="{ly+30}" r="4" fill="rgb(70,130,200)" opacity="0.5"/>')
   lines.append(f'<text x="{lx+22}" y="{ly+34}" class="tick-label">Lower ranked</text>')

   lines.append('</svg>')

   with open(filename, 'w') as f:
      f.write('\n'.join(lines))
   print(f"  Wrote {filename}")


def write_energy_components_svg(data, filename, n_show=10):
   """Horizontal bar chart of energy components for the top results."""

   n_show = min(n_show, len(data['e_total']))
   if n_show == 0:
      return

   width = 700
   bar_height = 30
   bar_gap = 8
   margin_l = 80
   margin_r = 180  # for legend
   margin_t = 50
   margin_b = 50
   plot_w = width - margin_l - margin_r
   total_bar_h = n_show * (bar_height + bar_gap) - bar_gap
   height = margin_t + total_bar_h + margin_b

   # Find the range for the x-axis
   all_vals = []
   for i in range(n_show):
      all_vals.extend([data['e_vdw'][i], data['e_elec'][i], data['e_air'][i]])
   x_min_val = min(0, min(all_vals))
   x_max_val = max(0, max(all_vals)) * 1.1

   # Ensure some range
   if abs(x_max_val - x_min_val) < 1:
      x_max_val = x_min_val + 100

   def sx(val):
      return margin_l + (val - x_min_val) / (x_max_val - x_min_val) * plot_w

   lines = []
   lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">')
   lines.append('<style>')
   lines.append('  text { font-family: Helvetica, Arial, sans-serif; }')
   lines.append('  .axis-label { font-size: 14px; fill: #333; }')
   lines.append('  .title { font-size: 16px; font-weight: bold; fill: #222; }')
   lines.append('  .tick-label { font-size: 11px; fill: #555; }')
   lines.append('  .bar-label { font-size: 11px; fill: #333; font-weight: bold; }')
   lines.append('</style>')
   lines.append(f'<rect width="{width}" height="{height}" fill="white"/>')

   lines.append(f'<text x="{(margin_l + width - margin_r)/2}" y="25" text-anchor="middle" '
                f'class="title">Energy Components by Rank</text>')

   # Zero line
   x_zero = sx(0)
   lines.append(f'<line x1="{x_zero}" y1="{margin_t - 5}" x2="{x_zero}" '
                f'y2="{margin_t + total_bar_h + 5}" stroke="#333" stroke-width="1"/>')

   # X-axis ticks
   n_ticks = 6
   for i in range(n_ticks + 1):
      val = x_min_val + (x_max_val - x_min_val) * i / n_ticks
      xp = sx(val)
      lines.append(f'<line x1="{xp}" y1="{margin_t + total_bar_h}" x2="{xp}" '
                   f'y2="{margin_t + total_bar_h + 5}" stroke="#999"/>')
      lines.append(f'<text x="{xp}" y="{margin_t + total_bar_h + 18}" '
                   f'text-anchor="middle" class="tick-label">{val:.0f}</text>')

   # Bars for each result
   colours = {'vdw': '#4477AA', 'elec': '#66CCEE', 'air': '#EE6677'}
   component_h = bar_height / 3

   for i in range(n_show):
      y_base = margin_t + i * (bar_height + bar_gap)
      dist = data['com_dist'][i]

      # Label
      lines.append(f'<text x="{margin_l - 8}" y="{y_base + bar_height/2 + 4}" '
                   f'text-anchor="end" class="bar-label">#{i} ({dist:.0f}\u00c5)</text>')

      # VDW bar
      vdw = data['e_vdw'][i]
      x1 = sx(min(0, vdw))
      x2 = sx(max(0, vdw))
      lines.append(f'<rect x="{x1}" y="{y_base}" width="{max(1, x2-x1)}" '
                   f'height="{component_h}" fill="{colours["vdw"]}" opacity="0.8"/>')

      # Elec bar
      elec = data['e_elec'][i]
      x1 = sx(min(0, elec))
      x2 = sx(max(0, elec))
      lines.append(f'<rect x="{x1}" y="{y_base + component_h}" '
                   f'width="{max(1, x2-x1)}" height="{component_h}" '
                   f'fill="{colours["elec"]}" opacity="0.8"/>')

      # AIR bar
      air = data['e_air'][i]
      x1 = sx(min(0, air))
      x2 = sx(max(0, air))
      lines.append(f'<rect x="{x1}" y="{y_base + 2*component_h}" '
                   f'width="{max(1, x2-x1)}" height="{component_h}" '
                   f'fill="{colours["air"]}" opacity="0.8"/>')

   # X axis label
   lines.append(f'<text x="{(margin_l + width - margin_r)/2}" y="{height - 12}" '
                f'text-anchor="middle" class="axis-label">Energy (kcal/mol)</text>')

   # Legend
   lx = width - margin_r + 15
   ly = margin_t + 10
   for label, colour in [('VDW', colours['vdw']),
                          ('Electrostatic', colours['elec']),
                          ('AIR', colours['air'])]:
      lines.append(f'<rect x="{lx}" y="{ly}" width="14" height="14" fill="{colour}"/>')
      lines.append(f'<text x="{lx+20}" y="{ly+12}" class="tick-label">{label}</text>')
      ly += 22

   lines.append('</svg>')

   with open(filename, 'w') as f:
      f.write('\n'.join(lines))
   print(f"  Wrote {filename}")


def main():
   csv_file = sys.argv[1] if len(sys.argv) > 1 else "dock-results.csv"

   if not os.path.exists(csv_file):
      print(f"Error: {csv_file} not found")
      sys.exit(1)

   data = read_csv(csv_file)
   n = len(data['e_total'])
   print(f"Read {n} results from {csv_file}")

   write_energy_vs_distance_svg(data, "dock-energy-vs-distance.svg")
   write_energy_components_svg(data, "dock-energy-components.svg")

   # Print summary table
   print(f"\n{'Rank':>4s}  {'E_total':>8s}  {'E_vdw':>8s}  {'E_elec':>8s}  "
         f"{'E_air':>8s}  {'Dist':>6s}")
   print("-" * 52)
   for i in range(min(10, n)):
      print(f"{i:4d}  {data['e_total'][i]:8.1f}  {data['e_vdw'][i]:8.1f}  "
            f"{data['e_elec'][i]:8.1f}  {data['e_air'][i]:8.1f}  "
            f"{data['com_dist'][i]:6.1f}")


if __name__ == "__main__":
   main()
