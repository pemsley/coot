# coot_commands/tts.py
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

"""Speak a line of text aloud with a local Piper (onnxruntime) neural voice.

Coot's Assistant tab spawns this as a subprocess to read the agent's final
answer out loud, in place of the macOS ``say`` command.  The text is passed as
the first command-line argument (falling back to stdin), synthesised entirely
offline with Piper, and played through the default output device.  The process
runs until playback finishes, or until Coot force-exits it to interrupt speech.

Errors are reported as a single JSON line on stdout so the GUI can surface them
(``{"type": "error", "message": "..."}``); success is silent.  Requires
``piper-tts`` and ``sounddevice`` (pip install piper-tts sounddevice).  The
voice model is chosen by ``COOT_PIPER_VOICE`` (path to a ``.onnx`` file),
defaulting to ``~/.local/share/coot/piper/en_GB-northern_english_male-medium.onnx``.
"""

from __future__ import annotations

import json
import os
import sys

DEFAULT_VOICE = os.path.expanduser(
    "~/.local/share/coot/piper/en_GB-northern_english_male-medium.onnx")


def _emit(obj: dict) -> None:
    print(json.dumps(obj), flush=True)


def _voice_path() -> str:
    return os.environ.get("COOT_PIPER_VOICE", DEFAULT_VOICE)


def speak(text: str) -> int:
    """Synthesise *text* with Piper and play it; return a process exit code."""
    text = text.strip()
    if not text:
        return 0

    try:
        import numpy as np
        import sounddevice as sd
        from piper import PiperVoice
    except ImportError as exc:
        _emit({"type": "error",
               "message": f"text-to-speech needs piper-tts + sounddevice ({exc})"})
        return 1

    voice_path = _voice_path()
    if not os.path.exists(voice_path):
        _emit({"type": "error",
               "message": f"Piper voice not found: {voice_path} "
                          "(set COOT_PIPER_VOICE or download a voice)"})
        return 1

    try:
        voice = PiperVoice.load(voice_path)
        parts = [chunk.audio_int16_array for chunk in voice.synthesize(text)]
    except Exception as exc:  # noqa: BLE001 - report any synthesis failure
        _emit({"type": "error", "message": f"speech synthesis failed: {exc}"})
        return 1

    if not parts:
        return 0
    samples = np.concatenate(parts)
    try:
        sd.play(samples, samplerate=voice.config.sample_rate)
        sd.wait()
    except Exception as exc:  # noqa: BLE001 - playback device may be unavailable
        _emit({"type": "error", "message": f"audio playback failed: {exc}"})
        return 1
    return 0


def main() -> int:
    if len(sys.argv) > 1:
        text = sys.argv[1]
    else:
        text = sys.stdin.read()
    return speak(text)


if __name__ == "__main__":
    sys.exit(main())
