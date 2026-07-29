"""Bundle SynEdu's small client-side enhancements into a static MyST build."""

from __future__ import annotations

import argparse
import hashlib
import re
import shutil
from pathlib import Path


BUILD_ASSET_RE = re.compile(
    r"""(?:href|src)=(?P<quote>["'])(?P<url>[^"']*/build/)[^"']+(?P=quote)"""
)
SCRIPT_RE = re.compile(
    r"""<script data-synedu-header src=["'][^"']+["'] defer></script>"""
)


def inject_site_assets(build_dir: Path, script_source: Path) -> tuple[int, Path]:
    """Copy the hashed script and add it to every rendered HTML document."""
    if not build_dir.is_dir():
        raise FileNotFoundError(f"HTML build directory does not exist: {build_dir}")
    if not script_source.is_file():
        raise FileNotFoundError(f"Header script does not exist: {script_source}")

    script_bytes = script_source.read_bytes()
    digest = hashlib.sha256(script_bytes).hexdigest()[:16]
    asset_name = f"synedu-header-{digest}.js"
    asset_path = build_dir / "build" / asset_name
    asset_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(script_source, asset_path)

    accepted = 0
    for html_path in sorted(build_dir.rglob("*.html")):
        html = html_path.read_text(encoding="utf-8")
        build_asset_match = BUILD_ASSET_RE.search(html)
        if not build_asset_match:
            continue
        accepted += 1

        asset_url = f"{build_asset_match.group('url')}{asset_name}"
        script_tag = f'<script data-synedu-header src="{asset_url}" defer></script>'

        if SCRIPT_RE.search(html):
            updated = SCRIPT_RE.sub(script_tag, html)
        elif "</head>" in html:
            updated = html.replace("</head>", f"{script_tag}</head>", 1)
        else:
            raise ValueError(f"Could not find </head> in {html_path}")

        if updated != html:
            html_path.write_text(updated, encoding="utf-8")

    if accepted == 0:
        raise RuntimeError("No rendered HTML pages accepted the SynEdu script")

    for stale_asset in asset_path.parent.glob("synedu-header-*.js"):
        if stale_asset != asset_path:
            stale_asset.unlink()

    return accepted, asset_path


def main() -> None:
    """Parse command-line arguments and inject the shared browser asset."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, default=Path("_build/html"))
    parser.add_argument(
        "--script",
        type=Path,
        default=Path("docs/_static/synedu-header.js"),
    )
    args = parser.parse_args()

    accepted, asset_path = inject_site_assets(args.build_dir, args.script)
    print(f"Ensured {asset_path.name} is present in {accepted} HTML pages.")


if __name__ == "__main__":
    main()
