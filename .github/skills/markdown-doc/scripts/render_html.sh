#!/usr/bin/env bash
# @AUTHOR: Chunjie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-09-24
# @DESCRIPTION: Render one Markdown document to a single self-contained HTML
#   file with pandoc. Mermaid fences are rewritten to <pre class="mermaid">
#   and, when the source has at least one, a vendored (offline, no CDN)
#   mermaid.js is embedded so the diagrams render client-side with no
#   network access.
# @VERSION: v0.1.0

set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
assets_dir="${script_dir}/../assets"

usage() {
  cat <<'EOF'
Usage: render_html.sh INPUT.md [OUTPUT.html]

  INPUT.md     Markdown file to render.
  OUTPUT.html  Defaults to INPUT with its extension replaced by .html.
EOF
}

if [[ $# -lt 1 || "$1" == "-h" || "$1" == "--help" ]]; then
  usage
  exit "$([[ $# -lt 1 ]] && echo 1 || echo 0)"
fi

input="$1"
[[ -f "${input}" ]] || {
  echo "render_html.sh: input not found: ${input}" >&2
  exit 1
}
output="${2:-${input%.md}.html}"

pandoc_cmd=(pandoc)
if command -v pixi >/dev/null 2>&1 && ! command -v pandoc >/dev/null 2>&1; then
  pandoc_cmd=(pixi run pandoc)
fi
command -v "${pandoc_cmd[0]}" >/dev/null 2>&1 || {
  echo "render_html.sh: no pandoc on PATH and pixi is unavailable" >&2
  exit 1
}

title=$(grep -m1 '^# ' "${input}" | sed 's/^# *//')
title="${title:-$(basename "${input}" .md)}"

args=(
  --standalone --embed-resources
  --toc --toc-depth=2
  --metadata "title=${title}"
  --css "${assets_dir}/doc.css"
  --lua-filter "${script_dir}/mermaid-codeblock.lua"
)

n_mermaid_src=$(grep -cE '^```+ *mermaid' "${input}" || true)
if [[ "${n_mermaid_src}" -gt 0 ]]; then
  combined=$(mktemp)
  trap 'rm -f "${combined}"' EXIT
  {
    printf '<script>\n'
    cat "${assets_dir}/mermaid.min.js"
    printf '\n</script>\n'
    cat "${assets_dir}/mermaid-init.html"
  } >"${combined}"
  args+=(--include-after-body "${combined}")
fi

"${pandoc_cmd[@]}" "${input}" "${args[@]}" -o "${output}"

n_mermaid_out=$(grep -c '<pre class="mermaid">' "${output}" || true)
if [[ "${n_mermaid_src}" -gt 0 && "${n_mermaid_src}" -ne "${n_mermaid_out}" ]]; then
  echo "render_html.sh: warning: ${n_mermaid_src} mermaid fence(s) in ${input} but ${n_mermaid_out} rendered in ${output}" >&2
fi

bytes=$(wc -c <"${output}")
echo "render_html.sh: wrote ${output} (${bytes} bytes, ${n_mermaid_out} mermaid block(s))"
echo "render_html.sh: open it in an actual browser to confirm the diagrams draw -- an editor or mail preview may disable script execution"
