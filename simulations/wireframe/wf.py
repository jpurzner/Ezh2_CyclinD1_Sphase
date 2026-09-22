"""Tiny SVG wireframe toolkit for the v44 model architecture figure, built module by module.

Visual grammar (shared by every module):
  activation   solid line, filled arrowhead
  inhibition   solid line, flat bar head
  drug         red, flat bar head, pill icon source
  feedback     dashed, arrowhead (slow / loop)
  equation     small grey monospace callout next to the node it governs

Pure-stdlib: writes an .svg string. Render with `qlmanage -t -s <px> -o <dir> file.svg` on macOS.
"""
from __future__ import annotations
import html

# ---------- palette ----------
PAL = dict(
    ink="#2c3e50", grey="#6b7280", light="#9ca3af", eq="#4b5563",
    act="#1f8f4a", inh="#c0392b", drug="#e74c3c", fb="#7d3c98", teal="#0e6251", blue="#2471a3",
    hh_fill="#e8f5ec", hh_line="#2e8b57",
    mycn_fill="#fff3cd", mycn_line="#e67e22",
    cd_fill="#d0ece7", cd_line="#0e6251",
    prom_fill="#fef5e7", prom_line="#b9770e",
    ezh2_fill="#ece2f2", ezh2_line="#7d3c98",
    node_fill="#ffffff",
    mem_fill="#f3f4f6", mem_line="#9ca3af",
    ext_fill="#f7fbff",
)

FONT = "Helvetica, Arial, sans-serif"
MONO = "Menlo, Consolas, monospace"
SHOW_EQ = False   # equation callouts are off by default (JP: figures without equations)


def _esc(s: str) -> str:
    return html.escape(s, quote=False)


class SVG:
    def __init__(self, w: int, h: int, bg: str = "#ffffff"):
        self.w, self.h = w, h
        self.body: list[str] = []
        self.markers: dict[tuple[str, str], str] = {}
        self.body.append(f'<rect x="0" y="0" width="{w}" height="{h}" fill="{bg}"/>')

    # ---------- markers ----------
    def _marker(self, kind: str, color: str) -> str:
        key = (kind, color)
        if key in self.markers:
            return self.markers[key]
        mid = f"m{len(self.markers)}"
        if kind == "arrow":
            d = ('<marker id="%s" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="8" markerHeight="8" '
                 'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="%s"/></marker>' % (mid, color))
        elif kind == "bar":
            d = ('<marker id="%s" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
                 'orient="auto-start-reverse"><path d="M8,0 L8,10" stroke="%s" stroke-width="2.4"/></marker>' % (mid, color))
        else:
            raise ValueError(kind)
        self.markers[key] = (mid, d)  # type: ignore[assignment]
        return self.markers[key]  # type: ignore[return-value]

    # ---------- primitives ----------
    def rect(self, x, y, w, h, fill, stroke, rx=8, sw=1.4, dashed=False, opacity=1.0):
        dash = ' stroke-dasharray="6,4"' if dashed else ""
        self.body.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" '
                         f'stroke="{stroke}" stroke-width="{sw}"{dash} opacity="{opacity}"/>')

    def text(self, x, y, s, fs=12, fill=None, anchor="start", mono=False, italic=False, bold=False,
             lh=1.25, opacity=1.0):
        fill = fill or PAL["ink"]
        fam = MONO if mono else FONT
        style = f'font-family="{fam}" font-size="{fs}" fill="{fill}" text-anchor="{anchor}"'
        if italic:
            style += ' font-style="italic"'
        if bold:
            style += ' font-weight="bold"'
        lines = s.split("\n")
        out = [f'<text x="{x}" y="{y}" {style} opacity="{opacity}">']
        for i, ln in enumerate(lines):
            dy = 0 if i == 0 else fs * lh
            out.append(f'<tspan x="{x}" dy="{dy}">{_esc(ln)}</tspan>')
        out.append("</text>")
        self.body.append("".join(out))

    def node(self, x, y, w, h, label, sub=None, fill=None, stroke=None, fs=13, rx=8, sw=1.5,
             dashed=False, bold=True, italic=False, sub_fs=10):
        """Centered box. Returns dict with edge anchor points."""
        fill = fill or PAL["node_fill"]
        stroke = stroke or PAL["ink"]
        self.rect(x - w / 2, y - h / 2, w, h, fill, stroke, rx=rx, sw=sw, dashed=dashed)
        nl = label.count("\n") + 1
        ns = (sub.count("\n") + 1) if sub else 0
        block = nl * fs * 1.2 + (ns * sub_fs * 1.25 + 3 if sub else 0)
        y0 = y - block / 2 + fs * 0.95
        self.text(x, y0, label, fs=fs, anchor="middle", bold=bold, italic=italic)
        if sub:
            self.text(x, y0 + nl * fs * 1.2 + sub_fs * 0.4, sub, fs=sub_fs, anchor="middle",
                      fill=PAL["grey"], italic=True)
        return dict(x=x, y=y, w=w, h=h,
                    l=(x - w / 2, y), r=(x + w / 2, y), t=(x, y - h / 2), b=(x, y + h / 2),
                    tl=(x - w / 2, y - h / 2), tr=(x + w / 2, y - h / 2),
                    bl=(x - w / 2, y + h / 2), br=(x + w / 2, y + h / 2))

    def eq(self, x, y, s, fs=10.5, anchor="start", fill=None):
        """Equation callout: small grey monospace (suppressed unless SHOW_EQ)."""
        if not SHOW_EQ:
            return
        self.text(x, y, s, fs=fs, fill=fill or PAL["eq"], anchor=anchor, mono=True)

    # ---------- composite glyphs ----------
    def io(self, x, y, label, sub=None, w=200, h=44, fs=12):
        """Inter-module hand-off node: dashed grey box."""
        return self.node(x, y, w, h, label, sub=sub, fill=PAL["mem_fill"], stroke=PAL["mem_line"],
                         fs=fs, dashed=True, bold=True)

    def gate(self, x, y, size=9, color=None, label=None, lx=None, ly=None):
        """Gate diamond on an edge."""
        color = color or PAL["grey"]
        self.body.append(f'<path d="M{x},{y - size} L{x + size},{y} L{x},{y + size} L{x - size},{y} z" '
                         f'fill="#ffffff" stroke="{color}" stroke-width="1.6"/>')
        if label:
            self.label(lx if lx is not None else x, ly if ly is not None else y - size - 5, label, fs=9.5, color=color)

    def badge(self, x, y, n, color=None, r=11):
        color = color or PAL["fb"]
        self.body.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{color}"/>')
        self.text(x, y + 4.5, str(n), fs=12, fill="#ffffff", anchor="middle", bold=True)

    def callout(self, x, y, w, lines, title=None, color=None, fill="#ffffff", fs=10.5, lh=15):
        """Top-left anchored box with a bold title and one line per entry. Returns node-like dict."""
        color = color or PAL["ink"]
        h = 14 + (lh if title else 0) + lh * len(lines) + 6
        self.rect(x, y, w, h, fill, color, rx=8, sw=1.4)
        yy = y + 18
        if title:
            self.text(x + 12, yy, title, fs=fs + 1, fill=color, bold=True)
            yy += lh
        for ln in lines:
            self.text(x + 12, yy, ln, fs=fs)
            yy += lh
        return dict(x=x + w / 2, y=y + h / 2, w=w, h=h, l=(x, y + h / 2), r=(x + w, y + h / 2),
                    t=(x + w / 2, y), b=(x + w / 2, y + h), tl=(x, y), tr=(x + w, y), bl=(x, y + h), br=(x + w, y + h))

    def phase_bar(self, x0, x1, y, segs, h=22, fs=11):
        """segs: list of (label, weight, fill, stroke)."""
        tot = sum(s[1] for s in segs)
        x = x0
        for lab, wgt, fill, stroke in segs:
            w = (x1 - x0) * wgt / tot
            self.rect(x, y, w, h, fill, stroke, rx=4, sw=1.2)
            self.text(x + w / 2, y + h / 2 + 4, lab, fs=fs, anchor="middle", bold=True, fill=stroke)
            x += w

    def locus(self, x, y, w=240, n_nuc=5, marked=3, pol2=True, label="Ccnd1 promoter", color=None,
              mark_color=None):
        """Chromatin fibre: DNA line, nucleosome discs, H3K27me3 flags on the first `marked` nucleosomes,
        RNA Pol II parked at the TSS with a nascent transcript. Returns anchors incl. 'pol2' and 'rna'."""
        color = color or PAL["prom_line"]
        mark_color = mark_color or PAL["ezh2_line"]
        x0, x1 = x - w / 2, x + w / 2
        self.line([(x0, y), (x1, y)], color=PAL["ink"], width=2)
        step = w / (n_nuc + 1)
        r = 13
        for i in range(n_nuc):
            cx = x0 + step * (i + 1)
            self.body.append(f'<circle cx="{cx}" cy="{y}" r="{r}" fill="{PAL["prom_fill"]}" stroke="{color}" stroke-width="1.6"/>')
            if i < marked:
                for dx in (-5, 5):
                    self.line([(cx + dx, y - r), (cx + dx, y - r - 9)], color=mark_color, width=1.4)
                    self.body.append(f'<circle cx="{cx + dx}" cy="{y - r - 12}" r="3.5" fill="{mark_color}"/>')
        out = dict(x=x, y=y, w=w, h=2 * r, l=(x0, y), r=(x1, y), t=(x, y - r - 16), b=(x, y + r),
                   tl=(x0, y - r), tr=(x1, y - r), bl=(x0, y + r), br=(x1, y + r))
        if pol2:
            px = x0 + step * (marked + 0.5)
            self.body.append(f'<ellipse cx="{px}" cy="{y - 6}" rx="19" ry="12" fill="#dbeafe" stroke="{PAL["blue"]}" stroke-width="1.4"/>')
            self.text(px, y - 2.5, "Pol II", fs=8.5, anchor="middle", bold=True, fill=PAL["blue"])
            # nascent transcript
            d = f"M{px + 10},{y - 16} q6,-8 12,0 t12,0 t12,0"
            self.body.append(f'<path d="{d}" fill="none" stroke="{PAL["teal"]}" stroke-width="1.8"/>')
            out["pol2"] = (px, y - 18)
            out["rna"] = (px + 46, y - 16)
        if label:
            self.text(x, y + r + 16, label, fs=12, anchor="middle", bold=True, italic=True, fill=color)
            out["b"] = (x, y + r + 20)
        return out

    def label(self, x, y, s, fs=10.5, color=None, anchor="middle", italic=True, bold=False):
        self.text(x, y, s, fs=fs, fill=color or PAL["grey"], anchor=anchor, italic=italic, bold=bold)

    def edge(self, pts, kind="act", color=None, width=1.8, dashed=False, curve=None, opacity=1.0):
        """pts = [(x,y),(x,y)] straight, or with curve=(cx,cy) a quadratic through control point.
        kind: act | inh | drug | fb | plain"""
        color = color or {"act": PAL["act"], "inh": PAL["inh"], "drug": PAL["drug"],
                          "fb": PAL["fb"], "plain": PAL["grey"]}[kind]
        head = "bar" if kind in ("inh", "drug") else "arrow"
        mid, mdef = self._marker(head, color)
        (x1, y1), (x2, y2) = pts[0], pts[-1]
        if curve is not None:
            cx, cy = curve
            d = f"M{x1},{y1} Q{cx},{cy} {x2},{y2}"
        elif len(pts) > 2:
            d = "M" + " L".join(f"{px},{py}" for px, py in pts)
        else:
            d = f"M{x1},{y1} L{x2},{y2}"
        dash = ' stroke-dasharray="7,5"' if (dashed or kind == "fb") else ""
        self.body.append(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}"{dash} '
                         f'marker-end="url(#{mid})" opacity="{opacity}" stroke-linecap="round"/>')

    def line(self, pts, color=None, width=1.2, dashed=False):
        color = color or PAL["light"]
        d = "M" + " L".join(f"{px},{py}" for px, py in pts)
        dash = ' stroke-dasharray="4,4"' if dashed else ""
        self.body.append(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}"{dash}/>')

    # ---------- icons ----------
    def membrane(self, y, h=26, x0=0, x1=None, label=True):
        x1 = x1 if x1 is not None else self.w
        # two leaflets of small circles
        self.body.append(f'<rect x="{x0}" y="{0}" width="{x1 - x0}" height="{y}" fill="{PAL["ext_fill"]}"/>')
        for yy in (y + 5, y + h - 5):
            cx = x0 + 6
            while cx < x1:
                self.body.append(f'<circle cx="{cx}" cy="{yy}" r="3.2" fill="#d6dde6"/>')
                cx += 9
        self.line([(x0, y), (x1, y)], color="#c3ccd6", width=1)
        self.line([(x0, y + h), (x1, y + h)], color="#c3ccd6", width=1)
        if label:
            self.text(x1 - 10, y - 10, "extracellular", fs=10.5, fill=PAL["light"], italic=True, anchor="end")
            self.text(x1 - 10, y + h + 16, "cytoplasm / nucleus", fs=10.5, fill=PAL["light"], italic=True, anchor="end")

    def dot(self, x, y, color=None, r=4):
        """Bus tap / junction dot."""
        self.body.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{color or PAL["act"]}"/>')

    def receptor(self, x, ymem, hmem, passes, label, color, faded=False, w=None):
        """Multi-pass membrane protein: `passes` vertical barrels through the membrane."""
        w = w or (passes * 7 + 10)
        op = 0.35 if faded else 1.0
        x0 = x - w / 2
        self.rect(x0, ymem - 8, w, hmem + 16, "#ffffff", color, rx=6, sw=1.4, opacity=op)
        for i in range(passes):
            bx = x0 + 5 + i * 7
            self.body.append(f'<rect x="{bx}" y="{ymem - 4}" width="4.5" height="{hmem + 8}" rx="2" '
                             f'fill="{color}" opacity="{0.55 * op}"/>')
        self.text(x, ymem + hmem + 30, label, fs=12, anchor="middle", bold=True, opacity=op)
        return dict(x=x, y=ymem + hmem / 2, w=w, h=hmem + 16,
                    l=(x0, ymem + hmem / 2), r=(x0 + w, ymem + hmem / 2),
                    t=(x, ymem - 8), b=(x, ymem + hmem + 8))

    def ligand(self, x, y, label, color, r=11):
        self.body.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{color}" opacity="0.85" stroke="{PAL["ink"]}" stroke-width="1"/>')
        self.text(x + r + 6, y + 4, label, fs=12, bold=True)
        return dict(x=x, y=y, l=(x - r, y), r=(x + r, y), t=(x, y - r), b=(x, y + r))

    def pill(self, x, y, label, w=96, h=24):
        self.rect(x - w / 2, y - h / 2, w, h, "#fdecea", PAL["drug"], rx=h / 2, sw=1.5)
        self.line([(x, y - h / 2), (x, y + h / 2)], color=PAL["drug"], width=1)
        self.text(x, y + 4, label, fs=11, anchor="middle", fill=PAL["drug"], bold=True)
        return dict(x=x, y=y, l=(x - w / 2, y), r=(x + w / 2, y), t=(x, y - h / 2), b=(x, y + h / 2))

    def tf(self, x, y, label, fill, stroke, w=84, h=40, sub=None):
        """Transcription-factor icon: rounded box with a small DNA-binding notch."""
        n = self.node(x, y, w, h, label, sub=sub, fill=fill, stroke=stroke, rx=14)
        self.body.append(f'<path d="M{x - 10},{y + h / 2} q10,8 20,0" fill="none" stroke="{stroke}" stroke-width="1.6"/>')
        return n

    def dna(self, x, y, w, label=None, color=None):
        """Short double helix segment (promoter icon)."""
        import math
        color = color or PAL["prom_line"]
        amp, per = 6, 18
        pts1 = " L".join(f"{x + i:.1f},{y + amp * math.sin(2 * math.pi * i / per):.1f}" for i in range(0, w + 1, 2))
        pts2 = " L".join(f"{x + i:.1f},{y - amp * math.sin(2 * math.pi * i / per):.1f}" for i in range(0, w + 1, 2))
        self.body.append(f'<path d="M{pts1}" fill="none" stroke="{color}" stroke-width="1.8"/>')
        self.body.append(f'<path d="M{pts2}" fill="none" stroke="{color}" stroke-width="1.8"/>')
        for i in range(per // 2, w, per // 2):
            self.line([(x + i, y - amp), (x + i, y + amp)], color=color, width=0.9)
        if label:
            self.text(x + w / 2, y + amp + 16, label, fs=11, anchor="middle", italic=True, fill=color)

    def zone(self, x, y, w, h, title, color, fill, title_fs=15):
        self.rect(x, y, w, h, fill, color, rx=14, sw=1.6, opacity=0.9)
        self.text(x + 14, y + 24, title, fs=title_fs, fill=color, bold=True)

    def legend(self, x, y, items, fs=11):
        """items: list of (kind, label)."""
        self.rect(x, y, 250, 16 + 20 * len(items), "#ffffff", "#d1d5db", rx=6, sw=1)
        for i, (kind, lab) in enumerate(items):
            yy = y + 18 + i * 20
            if kind == "eq":
                self.eq(x + 12, yy + 4, "k·x/(K+x)", fs=10)
            else:
                self.edge([(x + 12, yy), (x + 52, yy)], kind=kind, width=1.8)
            self.text(x + 62, yy + 4, lab, fs=fs)

    # ---------- output ----------
    def render(self) -> str:
        defs = "".join(d for (_, d) in self.markers.values())
        return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.w}" height="{self.h}" '
                f'viewBox="0 0 {self.w} {self.h}"><defs>{defs}</defs>' + "".join(self.body) + "</svg>")

    def save(self, path: str):
        with open(path, "w") as f:
            f.write(self.render())
        return path
