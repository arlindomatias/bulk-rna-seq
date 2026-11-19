"""Self-contained HTML report builder (no R needed).

Produces a single HTML file with a floating table of contents, Markdown prose, figures embedded
as base64 and tables rendered from pandas DataFrames. The layout follows the R Markdown
`html_document` look (toc_float, kable tables, blockquote callouts) so reports from Python and
R projects look the same.

The same file is copied into every project's `scripts/` folder so that each project stays
self-contained. Keep the copies identical.

Code folding: pass `code_files` (scripts that made the figures). Blocks wrapped in
    # <figure name="volcano">
    ...
    # </figure>
are shown behind a "Code" button above the figure whose file name (without extension) matches the
block name, like `code_folding: hide` in R Markdown. Works for R and Python scripts.

Usage:
    from html_report import Report
    r = Report("Title", subtitle="Dataset", code_files=["scripts/02_de.R"])
    r.h2("Section"); r.md("Some **Markdown** text.")
    r.figure("plots/pca.png", "Caption.")
    r.table(df, "Caption.")
    r.write("reports/report.html")
"""
import base64
import datetime
import html
import math
import mimetypes
import os
import re

import markdown as _markdown

_CSS = """
:root { --accent: #2a78d6; --ink: #0b0b0b; --ink2: #52514e; --rule: #dddddd; --soft: #f5f5f4; }
* { box-sizing: border-box; }
html { scroll-behavior: smooth; }
body { margin: 0; font-family: "Helvetica Neue", Helvetica, Arial, sans-serif; font-size: 15px;
       line-height: 1.55; color: #333; background: #fff; }
.wrap { display: flex; max-width: 1400px; margin: 0 auto; padding: 0 24px; gap: 32px; }
nav.toc { flex: 0 0 270px; position: sticky; top: 20px; align-self: flex-start;
          max-height: calc(100vh - 40px); overflow-y: auto; margin-top: 20px;
          border: 1px solid var(--rule); border-radius: 4px; font-size: 14px; }
nav.toc ul { list-style: none; margin: 0; padding: 0; }
nav.toc a { display: block; padding: 7px 14px; color: var(--accent); text-decoration: none;
            border-bottom: 1px solid #eee; }
nav.toc li.l3 a { padding-left: 30px; font-size: 13px; }
nav.toc li.l3 { display: none; }
nav.toc li.l3.show { display: block; }
nav.toc a:hover { background: var(--soft); }
nav.toc a.active { background: var(--accent); color: #fff; }
main { flex: 1 1 auto; min-width: 0; max-width: 1000px; padding-bottom: 80px; }
h1.title { font-size: 34px; font-weight: 500; margin: 26px 0 4px; color: var(--ink); line-height: 1.2; }
.subtitle { font-size: 20px; color: var(--ink2); margin: 0 0 4px; }
.date { font-size: 15px; color: var(--ink2); margin: 0 0 24px; }
h2 { font-size: 27px; font-weight: 500; margin: 40px 0 12px; color: var(--ink);
     padding-top: 8px; border-top: 1px solid #eee; }
h3 { font-size: 21px; font-weight: 500; margin: 28px 0 10px; color: var(--ink); }
h4 { font-size: 17px; font-weight: 500; margin: 20px 0 8px; color: var(--ink); }
p { margin: 0 0 12px; }
a { color: var(--accent); }
code { font-family: Menlo, Consolas, monospace; font-size: 90%; background: #f6f6f6;
       padding: 1px 4px; border-radius: 3px; color: #c7254e; }
pre { background: #f6f6f6; border: 1px solid #e4e4e4; border-radius: 4px; padding: 10px 12px;
      overflow-x: auto; font-size: 13px; line-height: 1.4; }
pre code { background: none; color: #333; padding: 0; }
a code { color: var(--accent); background: #eef4fc; }
blockquote { margin: 0 0 16px; padding: 10px 20px; border-left: 5px solid #eeeeee;
             font-size: 16px; background: #fafafa; }
blockquote p:last-child { margin-bottom: 0; }
ul, ol { margin: 0 0 12px; padding-left: 26px; }
li { margin: 2px 0; }
figure { margin: 18px 0 24px; text-align: center; }
figure img { max-width: min(100%, 860px); height: auto; }
figcaption { font-size: 14px; color: var(--ink2); margin-top: 6px; text-align: left; }
.tbl { overflow-x: auto; margin: 8px 0 22px; }
table { border-collapse: collapse; font-size: 13.5px; margin: 0 auto 0 0; }
table caption { caption-side: top; text-align: left; color: var(--ink2); padding: 4px 0 6px; font-size: 14px; }
th, td { padding: 5px 10px; border-top: 1px solid var(--rule); text-align: left; vertical-align: top; }
thead th { border-bottom: 2px solid var(--rule); border-top: none; font-weight: 600; }
td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
tbody tr:nth-child(odd) { background: #fafafa; }
details { margin: 6px 0 18px; }
details > summary { cursor: pointer; color: var(--accent); font-size: 14px; margin-bottom: 6px; }
details.code { margin: 14px 0 0; }
details.code > summary { list-style: none; display: inline-block; float: right; font-size: 12px;
  color: #333; background: #f5f5f5; border: 1px solid #ccc; border-radius: 3px; padding: 2px 9px; }
details.code > summary::-webkit-details-marker { display: none; }
details.code > summary::before { content: "Show code"; }
details.code[open] > summary::before { content: "Hide code"; }
details.code > summary:hover { background: #e6e6e6; }
details.code pre { clear: both; margin-top: 30px; }
details.code .src { font-size: 12px; color: var(--ink2); margin: 4px 0 0; }
.hl-c { color: #8e908c; font-style: italic; }
.hl-s { color: #3e7a1f; }
.code-all { float: right; margin-top: 30px; font-size: 12px; color: #333; background: #f5f5f5;
  border: 1px solid #ccc; border-radius: 3px; padding: 3px 10px; cursor: pointer; }
.kpis { display: flex; flex-wrap: wrap; gap: 12px; margin: 10px 0 22px; }
.kpi { flex: 1 1 150px; border: 1px solid var(--rule); border-radius: 6px; padding: 10px 14px; }
.kpi .v { font-size: 24px; color: var(--ink); font-variant-numeric: tabular-nums; }
.kpi .l { font-size: 13px; color: var(--ink2); }
@media (max-width: 900px) { .wrap { display: block; } nav.toc { position: static; max-height: none; } }
@media print { nav.toc { display: none; } }
"""

# Highlights the TOC entry of the section in view and expands its subsections (like tocify).
_JS = """
(function () {
  const links = Array.from(document.querySelectorAll('nav.toc a'));
  const byId = new Map(links.map(a => [a.getAttribute('href').slice(1), a]));
  const heads = Array.from(document.querySelectorAll('main h2[id], main h3[id]'));
  function update() {
    let cur = heads[0];
    for (const h of heads) { if (h.getBoundingClientRect().top < 120) cur = h; else break; }
    links.forEach(a => a.classList.remove('active'));
    document.querySelectorAll('nav.toc li.l3').forEach(li => li.classList.remove('show'));
    if (!cur) return;
    const a = byId.get(cur.id); if (a) a.classList.add('active');
    const li = a && a.parentElement;
    const group = li && (li.classList.contains('l3') ? li.dataset.parent : cur.id);
    document.querySelectorAll('nav.toc li.l3[data-parent="' + group + '"]').forEach(x => x.classList.add('show'));
  }
  document.addEventListener('scroll', update, { passive: true });
  update();
  const btn = document.querySelector('.code-all');
  if (btn) btn.addEventListener('click', () => {
    const blocks = document.querySelectorAll('details.code');
    const open = !Array.from(blocks).every(d => d.open);
    blocks.forEach(d => { d.open = open; });
    btn.textContent = open ? 'Hide all code' : 'Show all code';
  });
})();
"""


def _slug(text, used):
    s = re.sub(r"<[^>]+>", "", text)
    s = re.sub(r"[^\w\s-]", "", s.lower()).strip()
    s = re.sub(r"[\s_]+", "-", s) or "section"
    base, i = s, 2
    while s in used:
        s = f"{base}-{i}"
        i += 1
    used.add(s)
    return s


def _fmt(v, digits):
    """Formats a table cell: thousands separators for integers, `digits` significant digits for
    floats, scientific notation for very small or very large values."""
    import numbers
    if isinstance(v, bool) or not isinstance(v, numbers.Number):
        return "" if v is None else str(v)
    if isinstance(v, numbers.Integral):
        return f"{int(v):,}"
    if v != v:
        return ""
    if v == 0:
        return "0"
    if abs(v) < 1e-3 or abs(v) >= 1e9:
        return f"{v:.{max(digits - 1, 1)}e}"
    if float(v).is_integer():
        return f"{int(v):,}"
    decimals = max(digits - 1 - math.floor(math.log10(abs(v))), 0)
    return f"{v:,.{decimals}f}"


_FIG_START = re.compile(r'^\s*#\s*<figure name="([^"]+)">\s*$')
_FIG_END = re.compile(r"^\s*#\s*</figure>\s*$")


def load_code_blocks(paths):
    """Returns {name: (code, source_file)} for every `# <figure name=...>` block in `paths`."""
    blocks = {}
    for path in paths or []:
        with open(path, encoding="utf-8") as f:
            lines = f.read().splitlines()
        name, buf = None, []
        for line in lines:
            m = _FIG_START.match(line)
            if m:
                name, buf = m.group(1), []
            elif name and _FIG_END.match(line):
                indent = min((len(l) - len(l.lstrip()) for l in buf if l.strip()), default=0)
                code = "\n".join(l[indent:] for l in buf).strip("\n")
                blocks.setdefault(name, (code, os.path.basename(path)))
                name = None
            elif name:
                buf.append(line)
    return blocks


_TOKEN = re.compile(r'(?P<s>"(?:\\.|[^"\\\n])*"|\'(?:\\.|[^\'\\\n])*\')|(?P<c>#[^\n]*)')


def _highlight(code):
    """Minimal highlighting shared by R and Python: strings and comments."""
    out, pos = [], 0
    for m in _TOKEN.finditer(code):
        out.append(html.escape(code[pos:m.start()]))
        cls = "hl-s" if m.group("s") else "hl-c"
        out.append(f'<span class="{cls}">{html.escape(m.group(0))}</span>')
        pos = m.end()
    out.append(html.escape(code[pos:]))
    return "".join(out)


class Report:
    def __init__(self, title, subtitle=None, date=None, code_files=None):
        self.title = title
        self.subtitle = subtitle
        self.date = date or datetime.date.today().isoformat()
        self.body = []
        self.toc = []        # (level, id, text, parent_id)
        self._ids = set()
        self._h2 = None
        self._fig = 0
        self._tab = 0
        self.code_blocks = load_code_blocks(code_files)

    # ---- structure ---------------------------------------------------------
    def _heading(self, level, text):
        hid = _slug(text, self._ids)
        if level == 2:
            self._h2 = hid
        self.toc.append((level, hid, text, self._h2))
        self.body.append(f'<h{level} id="{hid}">{html.escape(text)}</h{level}>')
        return hid

    def h2(self, text):
        return self._heading(2, text)

    def h3(self, text):
        return self._heading(3, text)

    def h4(self, text):
        self.body.append(f"<h4>{html.escape(text)}</h4>")

    # ---- content -----------------------------------------------------------
    def md(self, text):
        """Markdown paragraph(s): lists, **bold**, `code`, links, > callouts, tables."""
        text = re.sub(r"^\n+", "", text)
        indent = min((len(l) - len(l.lstrip()) for l in text.splitlines() if l.strip()), default=0)
        text = "\n".join(l[indent:] for l in text.splitlines())
        self.body.append(_markdown.markdown(text, extensions=["tables", "sane_lists", "fenced_code"]))

    def raw(self, html_text):
        self.body.append(html_text)

    def kpis(self, items):
        """Row of headline numbers: list of (value, label)."""
        cells = "".join(f'<div class="kpi"><div class="v">{html.escape(str(v))}</div>'
                        f'<div class="l">{html.escape(l)}</div></div>' for v, l in items)
        self.body.append(f'<div class="kpis">{cells}</div>')

    def code_toggle(self, names):
        """Collapsed code block(s) taken from the marked blocks of `code_files`."""
        found = [(n, *self.code_blocks[n]) for n in names if n in self.code_blocks]
        if not found:
            return ""
        srcs = ", ".join(sorted({f"<code>{html.escape(src)}</code>" for _, _, src in found}))
        lang = "python" if all(src.endswith(".py") for _, _, src in found) else "r"
        code = "\n\n".join(c for _, c, _ in found)
        return (f'<details class="code"><summary></summary><pre><code class="{lang}">{_highlight(code)}</code></pre>'
                f'<p class="src">From {srcs}</p></details>')

    def figure(self, path, caption=None, width=None, optional=False, code="auto"):
        """Embeds an image file. Missing files are reported in the page instead of failing,
        or skipped silently when `optional=True` (plots that a step only makes when it has results).
        `code`: "auto" shows the code block named after the file, a list shows those blocks,
        None shows nothing."""
        if not os.path.exists(path):
            if optional:
                return
            self.body.append(f"<blockquote><p>Figure not found: <code>{html.escape(path)}</code></p></blockquote>")
            return
        self._fig += 1
        mime = mimetypes.guess_type(path)[0] or "image/png"
        with open(path, "rb") as f:
            data = base64.b64encode(f.read()).decode()
        style = f' style="width:{width}"' if width else ""
        cap = ""
        if caption:
            cap = f"<figcaption><b>Figure {self._fig}.</b> {_markdown.markdown(caption)[3:-4]}</figcaption>"
        if code == "auto":
            code = [os.path.splitext(os.path.basename(path))[0]]
        toggle = self.code_toggle(code) if code else ""
        self.body.append(f'{toggle}<figure><img src="data:{mime};base64,{data}" alt="{html.escape(os.path.basename(path))}"{style}>{cap}</figure>')

    def table(self, df, caption=None, digits=3, max_rows=None, collapsed=False, index=False):
        """Renders a pandas DataFrame. `collapsed=True` hides it behind a "Show table" toggle."""
        if index:
            df = df.reset_index()
        shown = df if max_rows is None else df.head(max_rows)
        num = [str(c) for c in shown.columns
               if str(shown[c].dtype).startswith(("int", "float", "Int", "Float"))]
        head = "".join(f'<th class="{"num" if str(c) in num else ""}">{html.escape(str(c))}</th>' for c in shown.columns)
        rows = []
        for _, r in shown.iterrows():
            tds = "".join(
                f'<td class="{"num" if str(c) in num else ""}">{html.escape(_fmt(r[c], digits))}</td>'
                for c in shown.columns)
            rows.append(f"<tr>{tds}</tr>")
        self._tab += 1
        cap = ""
        if caption:
            cap = f"<caption><b>Table {self._tab}.</b> {_markdown.markdown(caption)[3:-4]}</caption>"
        note = ""
        if max_rows is not None and len(df) > max_rows:
            note = f'<p style="font-size:13px;color:#52514e">Showing {max_rows} of {len(df):,} rows.</p>'
        block = f'<div class="tbl"><table>{cap}<thead><tr>{head}</tr></thead><tbody>{"".join(rows)}</tbody></table>{note}</div>'
        if collapsed:
            block = f"<details><summary>Show table</summary>{block}</details>"
        self.body.append(block)

    def code(self, text, lang="", collapsed=False, summary="Show"):
        block = f'<pre><code class="{lang}">{html.escape(text.strip())}</code></pre>'
        if collapsed:
            block = f"<details><summary>{html.escape(summary)}</summary>{block}</details>"
        self.body.append(block)

    # ---- appendix ----------------------------------------------------------
    def appendix(self, outputs, env_files=(), config_files=(), session_file=None):
        """Standard appendix shared by every project report.

        outputs: list of (file, description) rows.
        env_files: conda environment YAML files (envs/*.yaml); a matching *.linux-64.pin.txt is noted.
        config_files: files with the run parameters (config.yaml, or the Snakefile when there is none).
        session_file: text file (or list of files) with the R sessionInfo() of the analysis environments.
        """
        self.h2("Appendix")
        self.h3("Outputs")
        rows = "\n".join(f"| `{f}` | {d} |" for f, d in outputs)
        self.md(f"| File | Contents |\n|---|---|\n{rows}")

        self.h3("Software")
        envs = []
        for f in env_files:
            with open(f, encoding="utf-8") as fh:
                text = fh.read()
            section = text.split("dependencies:", 1)[-1]
            deps = [l.strip()[2:] for l in section.splitlines() if l.strip().startswith("- ")]
            name = os.path.basename(f).rsplit(".", 1)[0]
            pin = os.path.join(os.path.dirname(f), f"{name}.linux-64.pin.txt")
            envs.append({"environment": name,
                         "packages": ", ".join(deps),
                         "exact versions": f"envs/{os.path.basename(pin)}" if os.path.exists(pin) else "versions in the YAML"})
        if envs:
            import pandas as pd
            self.table(pd.DataFrame(envs), "Conda environments (`envs/`). Snakemake builds one per rule with `--use-conda`.")

        self.h3("Configuration")
        for f in config_files:
            with open(f, encoding="utf-8") as fh:
                lang = "yaml" if f.endswith((".yaml", ".yml")) else "python"
                self.code(fh.read(), lang, collapsed=True, summary=f"Show {os.path.basename(f)}")

        self.h3("Session info")
        files = [session_file] if isinstance(session_file, str) else list(session_file or [])
        files = [f for f in files if f and os.path.exists(f)]
        for f in files:
            label = "" if len(files) == 1 else f" ({os.path.basename(f)})"
            with open(f, encoding="utf-8") as fh:
                self.code(fh.read(), collapsed=True, summary=f"Show R sessionInfo(){label}")
        if not files:
            self.md("No session information was recorded.")

    # ---- output ------------------------------------------------------------
    def _toc_html(self):
        items = []
        for level, hid, text, parent in self.toc:
            cls = "l2" if level == 2 else "l3"
            items.append(f'<li class="{cls}" data-parent="{parent}"><a href="#{hid}">{html.escape(text)}</a></li>')
        return f'<nav class="toc"><ul>{"".join(items)}</ul></nav>'

    def write(self, path):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        sub = f'<p class="subtitle">{html.escape(self.subtitle)}</p>' if self.subtitle else ""
        doc = (
            "<!DOCTYPE html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">"
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">"
            f"<title>{html.escape(self.title)}</title><style>{_CSS}</style></head><body>"
            f'<div class="wrap">{self._toc_html()}<main>'
            + ('<button class="code-all">Show all code</button>' if self.code_blocks else "")
            + f'<h1 class="title">{html.escape(self.title)}</h1>{sub}<p class="date">{self.date}</p>'
            + "\n".join(self.body)
            + f"</main></div><script>{_JS}</script></body></html>\n"
        )
        with open(path, "w", encoding="utf-8") as f:
            f.write(doc)
