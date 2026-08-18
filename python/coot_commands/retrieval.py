# coot_commands/retrieval.py
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

"""Pick the commands relevant to a request, so the model sees few tools.

The bridge (:mod:`coot_commands.tools`) can emit all ~90 commands as tools,
but a small local model chooses far more reliably from a dozen than from
ninety.  This module ranks the commands by semantic similarity to the
user's request and returns the top few names, which
:func:`coot_commands.tools.command_tools` then narrows to via its *names*
argument.

Ranking uses text embeddings from a local, Ollama-style endpoint
(``/api/embed``, default model ``embeddinggemma`` - a ~300 MB pull;
override with ``COOT_EMBED_MODEL``).
Each command is embedded once from a short document built out of its name,
help, examples, category and notes; the query is embedded per request; we
return the commands with the highest cosine similarity.  Embeddings are
memoised for the process, so only the first request pays to embed the
command set.

As with :mod:`coot_commands.agent`, the embedding transport is injectable
(:class:`ToolRetriever` takes an *embed_fn*), so the ranking logic is
testable without a model server, and only the standard library is used.
"""

from __future__ import annotations

import json
import math
import os
import urllib.error
import urllib.request
from typing import Callable, Dict, List, Optional, Sequence

from coot_commands.registry import Command, all_commands

DEFAULT_EMBED_URL = "http://localhost:11434/api/embed"
DEFAULT_EMBED_MODEL = "embeddinggemma"

# embeddinggemma is trained with task-specific prompt prefixes and scores poorly
# without them: it embeds a *query* and a *document* asymmetrically, so the same
# text must be prefixed differently on each side. For retrieval Google specifies
# "task: search result | query: " for the query and "title: none | text: " for
# each document. Applying these turns fuzzy near-misses into confident hits (e.g.
# "Oh no! Go to A 89" then actually surfaces go_to_residue). They are specific to
# embeddinggemma, so default_retriever() only applies them for that model; a
# different embedder gets no prefix (or override via the env vars below).
EMBEDDINGGEMMA_QUERY_PREFIX = "task: search result | query: "
EMBEDDINGGEMMA_DOCUMENT_PREFIX = "title: none | text: "

# Embed a batch of texts -> one vector per text.
EmbedFn = Callable[[Sequence[str]], List[List[float]]]


def command_document(cmd: Command) -> str:
    """The text embedded to represent a command for retrieval.

    Bundles every scrap of natural language the command carries - help,
    example phrasings, category and notes - since domain terms a user might
    use ("H-bond", "rotamer", "blur") often live in the notes rather than the
    one-line help.
    """
    parts = [cmd.name.replace("_", " "), cmd.help_text or cmd.description]
    if cmd.examples:
        parts.append("Examples: " + "; ".join(cmd.examples))
    parts.append("Category: " + cmd.category)
    if cmd.notes:
        parts.append(cmd.notes)
    return ". ".join(p for p in parts if p)


def command_documents() -> Dict[str, str]:
    """Map command name -> its retrieval document, for every command."""
    docs: Dict[str, str] = {}
    for cmd in all_commands():
        docs.setdefault(cmd.name, command_document(cmd))
    return docs


def cosine(a: Sequence[float], b: Sequence[float]) -> float:
    """Cosine similarity of two vectors; 0.0 if either is degenerate."""
    dot = sum(x * y for x, y in zip(a, b))
    na = math.sqrt(sum(x * x for x in a))
    nb = math.sqrt(sum(y * y for y in b))
    if na == 0.0 or nb == 0.0:
        return 0.0
    return dot / (na * nb)


def _normalise_embed_url(url: str) -> str:
    """Accept a full endpoint or just a base, and return the embed endpoint.

    POSTing to the Ollama base URL returns 405 Method Not Allowed, so we
    tolerate a base or a ``.../api`` root and append the ``/api/embed`` path,
    and strip a trailing slash (which otherwise redirects).
    """
    url = url.rstrip("/")
    if url.endswith("/api/embed"):
        return url
    if url.endswith("/api"):
        return url + "/embed"
    return url + "/api/embed"


def ollama_embed(texts: Sequence[str], *,
                 model: Optional[str] = None,
                 url: Optional[str] = None,
                 timeout: float = 60.0) -> List[List[float]]:
    """Embed *texts* via an Ollama ``/api/embed`` endpoint (batched request)."""
    model = model or os.environ.get("COOT_EMBED_MODEL", DEFAULT_EMBED_MODEL)
    url = _normalise_embed_url(url or os.environ.get("COOT_EMBED_URL", DEFAULT_EMBED_URL))
    payload = {"model": model, "input": list(texts)}
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url, data=data, headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            body = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        # Surface the server's reason (e.g. "model 'x' not found") rather than a
        # bare status code; this message reaches the agent's fallback log.
        detail = e.read().decode("utf-8", "replace").strip()
        raise RuntimeError(
            f"embed request to {url} with model {model!r} failed "
            f"(HTTP {e.code}): {detail}") from None
    return body["embeddings"]


class ToolRetriever:
    """Rank documents against a query using an embedding transport.

    *documents* maps a name to its text; *embed_fn* embeds a batch of texts.
    Document embeddings are computed lazily on the first :meth:`select` and
    cached for the retriever's lifetime.

    *query_prefix* / *document_prefix* are prepended to the query and to each
    document before embedding, for embedders (like embeddinggemma) that expect
    task-specific prompt prefixes; both default to empty (symmetric embedding).
    """

    def __init__(self, documents: Dict[str, str], embed_fn: EmbedFn,
                 query_prefix: str = "", document_prefix: str = "") -> None:
        self.documents = documents
        self.embed_fn = embed_fn
        self.query_prefix = query_prefix
        self.document_prefix = document_prefix
        self._names: Optional[List[str]] = None
        self._vectors: Optional[List[List[float]]] = None

    def _ensure_embedded(self) -> None:
        if self._names is None:
            self._names = list(self.documents)
            self._vectors = self.embed_fn(
                [self.document_prefix + self.documents[n] for n in self._names])

    def select_scored(self, query: str, k: int) -> List[tuple[str, float]]:
        """The *k* best ``(name, cosine similarity)`` pairs for *query*.

        The scores matter where a caller needs to know not just which
        documents rank highest but whether any of them is actually relevant -
        see :func:`coot_commands.playbooks.select_playbooks`, which drops
        everything below a floor rather than always returning its top few.
        """
        self._ensure_embedded()
        query_vec = self.embed_fn([self.query_prefix + query])[0]
        scored = [(name, cosine(query_vec, vector))
                  for name, vector in zip(self._names, self._vectors)]
        scored.sort(key=lambda ns: ns[1], reverse=True)
        return scored[:k]

    def select(self, query: str, k: int) -> List[str]:
        """Return the *k* document names most similar to *query*, best first."""
        return [name for name, _ in self.select_scored(query, k)]


# Process-wide default retriever over the registry, embedded via Ollama.  Built
# lazily so importing this module never touches the network.
_default_retriever: Optional[ToolRetriever] = None


def default_prefixes() -> tuple[str, str]:
    """The (query, document) prefixes for the configured embed model.

    embeddinggemma needs its task prefixes to rank well; other models get none.
    Either can be overridden with ``COOT_EMBED_QUERY_PREFIX`` /
    ``COOT_EMBED_DOC_PREFIX`` (set them to empty strings to disable).

    Public because :mod:`coot_commands.playbooks` builds a second retriever
    over its own corpus and must embed it the same way this one does.
    """
    model = os.environ.get("COOT_EMBED_MODEL", DEFAULT_EMBED_MODEL)
    if "embeddinggemma" in model.lower():
        query, document = EMBEDDINGGEMMA_QUERY_PREFIX, EMBEDDINGGEMMA_DOCUMENT_PREFIX
    else:
        query, document = "", ""
    query = os.environ.get("COOT_EMBED_QUERY_PREFIX", query)
    document = os.environ.get("COOT_EMBED_DOC_PREFIX", document)
    return query, document


def default_retriever() -> ToolRetriever:
    """The shared retriever over all registered commands (Ollama embeddings)."""
    global _default_retriever
    if _default_retriever is None:
        query_prefix, document_prefix = default_prefixes()
        _default_retriever = ToolRetriever(
            command_documents(), ollama_embed,
            query_prefix=query_prefix, document_prefix=document_prefix)
    return _default_retriever


def select_tools(query: str, k: int = 12,
                 retriever: Optional[ToolRetriever] = None) -> List[str]:
    """Command names most relevant to *query* (convenience over the default)."""
    return (retriever or default_retriever()).select(query, k)
