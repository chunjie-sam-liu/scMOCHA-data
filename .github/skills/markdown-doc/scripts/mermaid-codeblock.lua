-- Turn every ```mermaid fenced block into <pre class="mermaid">, the markup
-- mermaid.js expects, instead of pandoc's default inert <pre><code> block.
function CodeBlock(el)
  if el.classes[1] == "mermaid" then
    local escaped = el.text:gsub("&", "&amp;"):gsub("<", "&lt;"):gsub(">", "&gt;")
    return pandoc.RawBlock("html", '<pre class="mermaid">' .. escaped .. "</pre>")
  end
end
