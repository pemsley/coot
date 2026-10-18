# coot_commands/mic.py
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

"""Record a voice request from the microphone and output it as base64 WAV.

Protocol (one JSON object per line on stdout):
    {"type": "recording"}               — microphone opened, recording has started
    {"type": "audio", "data": "..."}   — base64-encoded 16 kHz mono WAV, ready to send
    {"type": "error", "message": "..."}

Records until a period of silence is detected or the maximum duration is
reached.  Requires ``sounddevice`` and ``numpy`` (pip install sounddevice numpy).
"""

from __future__ import annotations

import base64
import io
import json
import sys
import wave

SAMPLERATE = 16000
SILENCE_THRESHOLD = 500      # RMS on int16 scale (0–32 768); tune to mic
SILENCE_DURATION = 1.5       # seconds of quiet that ends the recording
MAX_DURATION = 30.0          # hard cap in seconds


def _emit(obj: dict) -> None:
    print(json.dumps(obj), flush=True)


def _record(samplerate: int = SAMPLERATE,
            silence_threshold: float = SILENCE_THRESHOLD,
            silence_duration: float = SILENCE_DURATION,
            max_duration: float = MAX_DURATION):
    import numpy as np
    import sounddevice as sd

    chunk_s = 0.1                              # 100 ms chunks
    chunk_frames = int(samplerate * chunk_s)
    max_chunks = int(max_duration / chunk_s)
    silence_chunks_needed = int(silence_duration / chunk_s)

    chunks = []
    silent_count = 0
    speech_started = False

    with sd.InputStream(samplerate=samplerate, channels=1, dtype="int16") as stream:
        for _ in range(max_chunks):
            chunk, _ = stream.read(chunk_frames)
            chunks.append(chunk.copy())
            rms = float(np.sqrt(np.mean(chunk.astype(np.float32) ** 2)))
            if rms >= silence_threshold:
                speech_started = True
                silent_count = 0
            elif speech_started:
                silent_count += 1
                if silent_count >= silence_chunks_needed:
                    break

    import numpy as np
    return np.concatenate(chunks, axis=0) if chunks else np.zeros((0, 1), dtype="int16")


def _to_wav_b64(samples, samplerate: int = SAMPLERATE) -> str:
    buf = io.BytesIO()
    with wave.open(buf, "wb") as wf:
        wf.setnchannels(1)
        wf.setsampwidth(2)          # int16 = 2 bytes per sample
        wf.setframerate(samplerate)
        wf.writeframes(samples.flatten().tobytes())
    return base64.b64encode(buf.getvalue()).decode("ascii")


def main() -> int:
    try:
        import numpy  # noqa: F401
    except ImportError:
        _emit({"type": "error", "message": "numpy not installed (pip install numpy)"})
        return 1
    try:
        import sounddevice  # noqa: F401
    except ImportError:
        _emit({"type": "error", "message": "sounddevice not installed (pip install sounddevice)"})
        return 1

    _emit({"type": "recording"})
    try:
        samples = _record()
    except Exception as exc:
        _emit({"type": "error", "message": f"recording failed: {exc}"})
        return 1

    import numpy as np
    if samples.size == 0 or float(np.sqrt(np.mean(samples.astype(np.float32) ** 2))) < SILENCE_THRESHOLD / 2:
        _emit({"type": "error", "message": "no speech detected"})
        return 1

    try:
        data = _to_wav_b64(samples)
    except Exception as exc:
        _emit({"type": "error", "message": f"WAV encoding failed: {exc}"})
        return 1

    _emit({"type": "audio", "data": data})
    return 0


if __name__ == "__main__":
    sys.exit(main())
