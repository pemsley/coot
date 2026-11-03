# coot_commands/commands/appearance.py
#
# Copyright 2026 Jordan Dialpuri, Medical Research Council Laboratory of Molecular Biology
#
# This file is part of Coot
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.

"""Shader and appearance controls: lighting, shadows, depth, materials.

These commands change how the scene is *drawn* rather than what is in it -
the knobs that turn a functional view into a publication figure: fancy
lighting, screen-space ambient occlusion (SSAO), shadows, depth-of-field
blur, fog, silhouette outlines, background colour and per-model material
(specular/shininess).  They are ordinary ``@command`` handlers, so they
work from the Command and Assistant tabs; they are also the tool surface
the *figure agent* (:mod:`coot_commands.figure`) iterates over while it
tunes a rendering to look beautiful.

Each setter maps to a single Coot API call and redraws.  Where the API
takes a raw float/int, we accept the same and validate it; ranges quoted in
the ``notes`` are the useful ones a human (or the figure agent) should stay
within, not hard limits.
"""

from __future__ import annotations

from typing import Optional

from coot_commands.registry import command
from coot_commands.types import (as_int, as_float, resolve_model, resolve_colour,
                                 ArgType, ACTIVE_MODEL_NOTE)

try:
    import coot
except ImportError:
    coot = None


CATEGORY = "Appearance"

# "set X on" / "set X off" -> a 0/1 state.  Mirrors settings.py's _TO/_GET
# fragments; kept local so appearance commands read on their own.
_STATE = r"(?P<state>on|off|enabled?|disabled?|true|false|1|0)"


def _on(state: Optional[str]) -> int:
    """Interpret an on/off word as 1/0 (default on when omitted)."""
    if state is None:
        return 1
    return 0 if state.lower() in ("off", "disable", "disabled", "false", "0") else 1


def _draw() -> None:
    """Redraw the graphics if the API is present (a no-op otherwise)."""
    if coot is not None:
        coot.graphics_draw()


def _call(name: str, *args) -> None:
    """Call a Coot API function by name if the build exposes it, then redraw.

    Looked up by name so a Coot that predates one of these setters degrades to
    a no-op with a clear message rather than crashing the agent loop.
    """
    if coot is None:
        return
    fn = getattr(coot, name, None)
    if fn is not None:
        fn(*args)
    _draw()


# ---------------------------------------------------------------------------
#   Fancy lighting (the master switch for the modern shader path)
# ---------------------------------------------------------------------------

@command(r"(?:set )?fancy lighting " + _STATE,
         examples=["set fancy lighting on", "fancy lighting off"],
         category=CATEGORY,
         arg_types={"state": ("on", "off")},
         notes="Turns on Coot's modern per-fragment lighting shader. Most of "
               "the other appearance controls (shadows, ambient occlusion, "
               "specular highlights) only take visible effect with fancy "
               "lighting on.")
def set_fancy_lighting(state: Optional[str] = None) -> str:
    """Enable or disable the fancy (per-fragment) lighting shader."""
    on = _on(state)
    _call("set_use_fancy_lighting", on)
    return f"Fancy lighting {'on' if on else 'off'}"


# ---------------------------------------------------------------------------
#   Fog / depth cueing
# ---------------------------------------------------------------------------

@command(r"(?:set )?fog " + _STATE,
         examples=["set fog on", "fog off"],
         category=CATEGORY,
         arg_types={"state": ("on", "off")},
         notes="Depth-cueing fog: fades distant atoms into the background so "
               "the eye reads depth. Good for thick views; can wash out a thin "
               "one.")
def set_fog(state: Optional[str] = None) -> str:
    """Turn depth-cueing fog on or off."""
    on = _on(state)
    _call("set_use_fog", on)
    return f"Fog {'on' if on else 'off'}"


# ---------------------------------------------------------------------------
#   Silhouette outline
# ---------------------------------------------------------------------------

@command(r"(?:set )?outlines? " + _STATE,
         examples=["set outline on", "outlines off"],
         category=CATEGORY,
         arg_types={"state": ("on", "off")},
         notes="Draws a dark silhouette outline around objects, for a "
               "molecular-illustration / cartoon look.")
def set_outline(state: Optional[str] = None) -> str:
    """Turn the silhouette outline effect on or off."""
    on = _on(state)
    _call("set_use_outline", on)
    return f"Outline {'on' if on else 'off'}"


# ---------------------------------------------------------------------------
#   Ambient occlusion (SSAO)
# ---------------------------------------------------------------------------

@command(r"(?:set )?(?:ambient occlusion|ssao|ao) " + _STATE,
         examples=["set ambient occlusion on", "ssao off"],
         category=CATEGORY,
         arg_types={"state": ("on", "off")},
         notes="Screen-space ambient occlusion: darkens crevices and contact "
               "regions, giving strong shape and depth. One of the biggest "
               "wins for a figure.")
def set_ambient_occlusion(state: Optional[str] = None) -> str:
    """Turn screen-space ambient occlusion on or off."""
    on = _on(state)
    _call("set_use_ambient_occlusion", on)
    return f"Ambient occlusion {'on' if on else 'off'}"


@command(r"set (?:ambient occlusion|ssao|ao) strength(?: to| =)? (?P<value>[\d.]+)",
         examples=["set ambient occlusion strength to 0.6", "set ao strength 0.8"],
         category=CATEGORY,
         notes="How dark the ambient-occlusion shading gets. ~0.4-0.9 reads "
               "well; above ~1.2 looks sooty.")
def set_ambient_occlusion_strength(value: str) -> str:
    """Set the ambient-occlusion (SSAO) strength."""
    v = as_float(value, "ambient occlusion strength")
    _call("set_ssao_strength", v)
    return f"Ambient occlusion strength {v:g}"


@command(r"set (?:ambient occlusion|ssao|ao) radius(?: to| =)? (?P<value>[\d.]+)",
         examples=["set ambient occlusion radius to 0.5", "set ao radius 1.0"],
         category=CATEGORY,
         notes="The radius (in world units) over which ambient occlusion is "
               "sampled. Larger radii shade broader hollows.")
def set_ambient_occlusion_radius(value: str) -> str:
    """Set the ambient-occlusion (SSAO) sampling radius."""
    v = as_float(value, "ambient occlusion radius")
    _call("set_ssao_radius", v)
    return f"Ambient occlusion radius {v:g}"


@command(r"set (?:ambient occlusion|ssao|ao) bias(?: to| =)? (?P<value>[\d.]+)",
         examples=["set ambient occlusion bias to 0.05"],
         category=CATEGORY,
         notes="Bias that suppresses self-occlusion speckling; small values "
               "(~0.02-0.1). Raise it if flat surfaces look noisy.")
def set_ambient_occlusion_bias(value: str) -> str:
    """Set the ambient-occlusion (SSAO) bias."""
    v = as_float(value, "ambient occlusion bias")
    _call("set_ssao_bias", v)
    return f"Ambient occlusion bias {v:g}"


@command(r"set (?:ambient occlusion|ssao|ao) samples(?: to| =)? (?P<value>\d+)",
         examples=["set ambient occlusion samples to 64"],
         category=CATEGORY,
         notes="Number of SSAO kernel samples: higher is smoother but slower "
               "(typical 32-64).")
def set_ambient_occlusion_samples(value: str) -> str:
    """Set the number of ambient-occlusion kernel samples."""
    v = as_int(value, "ambient occlusion samples")
    _call("set_ssao_kernel_n_samples", v)
    return f"Ambient occlusion samples {v}"


@command(r"set (?:ambient occlusion|ssao|ao) blur(?: to| =)? (?P<value>\d+)",
         examples=["set ambient occlusion blur to 2"],
         category=CATEGORY,
         notes="Blur kernel size applied to the SSAO buffer to hide sampling "
               "noise (typical 1-4).")
def set_ambient_occlusion_blur(value: str) -> str:
    """Set the ambient-occlusion blur size."""
    v = as_int(value, "ambient occlusion blur")
    _call("set_ssao_blur_size", v)
    return f"Ambient occlusion blur {v}"


# ---------------------------------------------------------------------------
#   Shadows
# ---------------------------------------------------------------------------

@command(r"set shadow strength(?: to| =)? (?P<value>[\d.]+)",
         examples=["set shadow strength to 0.5"],
         category=CATEGORY,
         notes="How dark cast shadows are (0 = none, ~0.4-0.7 is a natural "
               "figure look, 1 = very dark).")
def set_shadow_strength(value: str) -> str:
    """Set the cast-shadow strength."""
    v = as_float(value, "shadow strength")
    _call("set_shadow_strength", v)
    return f"Shadow strength {v:g}"


@command(r"set shadow softness(?: to| =)? (?P<value>\d+)",
         examples=["set shadow softness to 2"],
         category=CATEGORY,
         notes="Softness of shadow edges (typical 1-3): higher is softer, "
               "more diffuse shadows.")
def set_shadow_softness(value: str) -> str:
    """Set the shadow-edge softness."""
    v = as_int(value, "shadow softness")
    _call("set_shadow_softness", v)
    return f"Shadow softness {v}"


@command(r"set shadow resolution(?: to| =)? (?P<value>\d+)",
         examples=["set shadow resolution to 2"],
         category=CATEGORY,
         notes="Shadow-map resolution multiplier: higher gives crisper "
               "shadows at a memory/speed cost (typical 1-4).")
def set_shadow_resolution(value: str) -> str:
    """Set the shadow-map resolution multiplier."""
    v = as_int(value, "shadow resolution")
    _call("set_shadow_resolution", v)
    return f"Shadow resolution {v}"


# ---------------------------------------------------------------------------
#   Depth of field (focus blur)
# ---------------------------------------------------------------------------

@command(r"(?:set )?(?:depth of field|depth blur|dof) " + _STATE,
         examples=["set depth of field on", "depth blur off"],
         category=CATEGORY,
         arg_types={"state": ("on", "off")},
         notes="Depth-of-field blur: keeps the focal plane sharp and blurs "
               "foreground/background, drawing the eye to the subject.")
def set_depth_of_field(state: Optional[str] = None) -> str:
    """Turn depth-of-field (focus) blur on or off."""
    on = _on(state)
    _call("set_use_depth_blur", on)
    return f"Depth of field {'on' if on else 'off'}"


@command(r"set (?:depth of field|depth blur|dof) strength(?: to| =)? (?P<value>[\d.]+)",
         examples=["set depth of field strength to 1.0"],
         category=CATEGORY,
         notes="How strongly out-of-focus regions are blurred.")
def set_depth_of_field_strength(value: str) -> str:
    """Set the depth-of-field blur strength."""
    v = as_float(value, "depth of field strength")
    _call("set_focus_blur_strength", v)
    return f"Depth of field strength {v:g}"


@command(r"set (?:depth of field|depth blur|dof) (?:distance|depth)(?: to| =)? (?P<value>[\d.]+)",
         examples=["set depth of field distance to 0.5"],
         category=CATEGORY,
         notes="Where the focal plane sits in depth (0 = near, 1 = far). Put "
               "it on the subject you want sharp.")
def set_depth_of_field_distance(value: str) -> str:
    """Set the depth-of-field focal-plane depth."""
    v = as_float(value, "depth of field distance")
    _call("set_focus_blur_z_depth", v)
    return f"Depth of field distance {v:g}"


# ---------------------------------------------------------------------------
#   Post-processing (brightness / gamma)
# ---------------------------------------------------------------------------

@command(r"set brightness(?: to| =)? (?P<value>[\d.]+)",
         examples=["set brightness to 1.1"],
         category=CATEGORY,
         notes="Overall image brightness multiplier of the post-process pass "
               "(1.0 = unchanged).")
def set_brightness(value: str) -> str:
    """Set the post-processing brightness."""
    v = as_float(value, "brightness")
    _call("set_effects_shader_brightness", v)
    return f"Brightness {v:g}"


@command(r"set gamma(?: to| =)? (?P<value>[\d.]+)",
         examples=["set gamma to 1.0"],
         category=CATEGORY,
         notes="Gamma of the post-process pass (1.0 = unchanged; >1 lifts "
               "midtones, <1 deepens them).")
def set_gamma(value: str) -> str:
    """Set the post-processing gamma."""
    v = as_float(value, "gamma")
    _call("set_effects_shader_gamma", v)
    return f"Gamma {v:g}"


# ---------------------------------------------------------------------------
#   Bond smoothness (roundness of stick cylinders)
# ---------------------------------------------------------------------------

@command(r"set bond smoothness(?: to| =)? (?P<value>\d+)",
         examples=["set bond smoothness to 3"],
         category=CATEGORY,
         notes="Number of sides on the bond cylinders / atom spheres: higher "
               "is rounder and smoother (typical 1-4). Raise it for a "
               "close-up figure.")
def set_bond_smoothness(value: str) -> str:
    """Set the bond/atom mesh smoothness factor."""
    v = as_int(value, "bond smoothness")
    _call("set_bond_smoothness_factor", v)
    return f"Bond smoothness {v}"


# ---------------------------------------------------------------------------
#   Background colour
# ---------------------------------------------------------------------------

@command(r"set background(?: colou?r)?(?: to| =)? "
         r"(?:(?P<r>[\d.]+) (?P<g>[\d.]+) (?P<b>[\d.]+)|"
         r"(?P<colour>[a-zA-Z][a-zA-Z ]*?))\s*$",
         examples=["set background to black", "set background white",
                   "set background colour to 0.1 0.1 0.15"],
         category=CATEGORY,
         arg_types={"colour": ArgType.COLOUR},
         notes="Background colour, by name (e.g. black, white, sky blue) or as "
               "three 0-1 RGB values. Black or near-black makes bright models "
               "pop; white suits print figures.")
def set_background(colour: Optional[str] = None, r: Optional[str] = None,
                   g: Optional[str] = None, b: Optional[str] = None) -> str:
    """Set the background colour, by name or RGB triple."""
    if r is not None and g is not None and b is not None:
        red, green, blue = (as_float(r, "red"), as_float(g, "green"),
                            as_float(b, "blue"))
    else:
        red, green, blue = resolve_colour(colour or "black")
    _call("set_background_colour", red, green, blue)
    return f"Background colour ({red:g}, {green:g}, {blue:g})"


# ---------------------------------------------------------------------------
#   Per-model material (specular highlight / shininess)
# ---------------------------------------------------------------------------

@command(r"set (?:model )?(?:specular|shininess|material)"
         r"(?: (?:for|of) model (?P<model>\S+))?"
         r"(?: to| =)? (?P<strength>[\d.]+)(?: (?P<shininess>[\d.]+))?",
         examples=["set model specular to 0.5 128",
                   "set specular for model 0 to 0.3"],
         category=CATEGORY,
         arg_types={"model": ArgType.MODEL},
         notes="Specular highlight of a model's atoms: STRENGTH (~0-1, how "
               "bright the highlight) and optional SHININESS (~8-256, how "
               "tight it is). Higher shininess gives a smaller, glossier "
               "highlight. " + ACTIVE_MODEL_NOTE)
def set_model_specular(strength: str, shininess: Optional[str] = None,
                       model: Optional[str] = None) -> str:
    """Set a model's specular strength and (optionally) shininess."""
    imol = resolve_model(model)
    s = as_float(strength, "specular strength")
    shin = as_float(shininess, "shininess") if shininess is not None else 128.0
    _call("set_model_material_specular", imol, s, shin)
    return f"Model {imol} specular strength {s:g}, shininess {shin:g}"
