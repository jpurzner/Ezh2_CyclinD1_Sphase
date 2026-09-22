"""Build every module wireframe SVG + HTML wrappers for the browser pane, and (macOS, if Chrome is
installed) a 2x PNG of each for slides.

Run:  python3 simulations/wireframe/build_all.py            # SVG + HTML + PNG
      python3 simulations/wireframe/build_all.py --no-png   # skip the Chrome render
"""
import os, sys, importlib, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

MODULES = ["module1_mitogenic_drive", "module2_ccnd1_locus", "module3_restriction_point", "module4_s_g2_m"]
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
PNG_SCALE = 2   # 1400x800 canvas -> 2800x1600 PNG


def render_png(svg_path: str, w: int, h: int) -> str | None:
    if "--no-png" in sys.argv or not os.path.exists(CHROME):
        return None
    png = svg_path[:-4] + ".png"
    subprocess.run([CHROME, "--headless=new", "--disable-gpu", "--hide-scrollbars",
                    f"--window-size={w},{h}", f"--force-device-scale-factor={PNG_SCALE}",
                    f"--screenshot={png}", "file://" + svg_path],
                   check=True, capture_output=True)
    return png


for name in MODULES:
    mod = importlib.import_module(name)
    svg = mod.build()
    out = os.path.join(HERE, name + ".svg")
    svg.save(out)
    for suffix, shift in ((".html", 0), ("_right.html", 600)):
        with open(os.path.join(HERE, name + suffix), "w") as f:
            f.write(f'<!doctype html><html><head><meta charset="utf-8"><title>{name}</title>'
                    f'<style>html,body{{margin:0;background:#fff;overflow:hidden}}</style></head><body>'
                    f'<img src="{name}.svg?v={os.path.getmtime(out):.0f}" width="{svg.w}" height="{svg.h}" '
                    f'style="margin-left:-{shift}px"></body></html>')
    print("wrote", out)
    png = render_png(out, svg.w, svg.h)
    if png:
        print("wrote", png)
