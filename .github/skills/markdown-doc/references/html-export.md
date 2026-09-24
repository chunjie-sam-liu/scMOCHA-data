# Rendering a Markdown document to a single-file HTML

This is the full contract behind `scripts/render_html.sh`. Read section 12 of `SKILL.md` first; this file covers what that section only points at: the asset contract, how to refresh the vendored Mermaid version, and troubleshooting.

## A. Why a bundled script instead of a one-off pandoc command

A hand-typed `pandoc file.md -o file.html` produces a single file, but it fails the moment the document has a Mermaid diagram: pandoc emits every fenced ` ```mermaid ` block as an inert `<pre><code>` block, not a rendered figure, and a plain `--standalone` build has no JavaScript at all to draw one. Getting this right by hand every time means remembering three independent fixes -- the code-block markup, a Mermaid renderer, and making that renderer run without a network call -- so the fix lives in one script instead.

## B. The contract

| Piece | What it does | Where it lives |
| --- | --- | --- |
| `pandoc --standalone --embed-resources` | Produces one HTML file with no external stylesheet, script, or image reference | pandoc itself |
| `assets/doc.css` | Shared, GitHub-like typography and table styling, applied via `--css` and inlined by `--embed-resources` | this skill |
| `scripts/mermaid-codeblock.lua` | A pandoc Lua filter: rewrites every ` ```mermaid ` fenced block's AST node to `<pre class="mermaid">`, the markup mermaid.js expects, instead of pandoc's default code block | this skill |
| `assets/mermaid.min.js` | A vendored copy of mermaid.js (UMD build), embedded inline as a `<script>` block only when the source has at least one Mermaid fence | this skill |
| `assets/mermaid-init.html` | The two-line init/run call (`mermaid.initialize(...)`, `mermaid.run(...)`) appended after the vendored library | this skill |

**Vendored, never a CDN reference.** A `<script src="https://cdn...">` tag would leave the file broken the moment it is opened offline, on an air-gapped machine, or years later after the CDN path changes -- exactly the situations a single-file export exists for. Embedding the library inline costs about 3 MB of file size and buys permanence.

## C. Usage

```bash
pixi run bash .github/skills/markdown-doc/scripts/render_html.sh path/to/README.md
# optional second argument redirects the output path; default is the input
# with .md replaced by .html, written beside the source.
```

The script:

1. Reads the document's first `# ` heading as the HTML `<title>`.
2. Converts with pandoc, the Lua filter, and the shared stylesheet.
3. Counts ` ```mermaid ` fences in the source; if there is at least one, builds a temporary combined `<script>` block (vendored library + init snippet) and passes it through `--include-after-body`.
4. Reports the output size and the rendered Mermaid block count, and warns on stderr if that count does not match the source -- a mismatch means the Lua filter did not fire, most often because pandoc could not find it (wrong `--lua-filter` path) or the fence used a different language tag (`mermaid` must be the first class).

## D. Verification

A byte count and a matching block count only prove the file was written and the markup is present -- neither proves the JavaScript actually drew a diagram. Two ways to check further, cheapest first:

- `grep -c '<svg' output.html` -- mermaid.js replaces each `<pre class="mermaid">` with an inline `<svg>` once it runs; a count lower than the mermaid block count means at least one diagram did not render (a syntax error inside that block is the most common cause -- check the browser console).
- Open the file in an actual browser and confirm the diagrams are visible. **An editor preview pane or a mail client's HTML preview is not sufficient** -- both commonly disable script execution, which leaves the raw diagram source sitting inertly inside a `<pre>` and looks identical to the broken case from section A.

## E. Refreshing the vendored Mermaid version

```bash
curl -sL https://cdn.jsdelivr.net/npm/mermaid@<version>/dist/mermaid.min.js \
  -o .github/skills/markdown-doc/assets/mermaid.min.js
```

Re-render a document with several diagram types afterward and re-run the section D checks -- a major-version bump can change flowchart styling or the `securityLevel` default. Mermaid is MIT-licensed; the vendored file carries no separate license text because the license is unchanged from upstream.

## F. Troubleshooting

| Symptom | Cause | Fix |
| --- | --- | --- |
| `no pandoc on PATH and pixi is unavailable` | Neither a system pandoc nor a pixi project with pandoc as a dependency is reachable | Run from inside a pixi project that lists `pandoc`, or install pandoc on PATH |
| Diagrams show as plain text, no error printed | Opened in a preview pane that disables script execution | Open the file in a real browser tab (see section D) |
| Rendered mermaid block count is lower than the source count | The Lua filter did not fire on that fence | Confirm the fence opens with exactly ` ```mermaid ` (the language tag is the first class pandoc sees); check the `--lua-filter` path resolved (the script derives it from its own location, so do not `source` it from another script under a different name) |
| Output file is huge relative to the source | Expected: the vendored library is ~3.3 MB and is embedded whenever the source has any Mermaid fence, even one | Not a defect; the file is meant to be self-contained. Confirm nothing was embedded twice by checking `grep -c 'mermaid.initialize'` equals 1 |
