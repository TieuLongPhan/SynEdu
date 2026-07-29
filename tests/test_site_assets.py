"""Regression tests for the static documentation publication path."""

from pathlib import Path

import pytest

from scripts.inject_site_assets import inject_site_assets


REPO_ROOT = Path(__file__).parents[1]


def test_site_asset_injection_is_idempotent_and_removes_stale_assets(
    tmp_path: Path,
) -> None:
    build_dir = tmp_path / "html"
    build_assets = build_dir / "build"
    build_assets.mkdir(parents=True)
    page = build_dir / "index.html"
    page.write_text(
        '<html><head><script src="/build/root.js"></script></head></html>',
        encoding="utf-8",
    )
    source = tmp_path / "header.js"
    source.write_text("console.log('header');\n", encoding="utf-8")
    stale = build_assets / "synedu-header-stale.js"
    stale.write_text("stale\n", encoding="utf-8")

    accepted_first, asset = inject_site_assets(build_dir, source)
    accepted_second, second_asset = inject_site_assets(build_dir, source)

    assert accepted_first == accepted_second == 1
    assert asset == second_asset
    assert not stale.exists()
    html = page.read_text(encoding="utf-8")
    assert html.count("data-synedu-header") == 1
    assert asset.name in html


def test_site_asset_injection_rejects_non_myst_html(tmp_path: Path) -> None:
    build_dir = tmp_path / "html"
    build_dir.mkdir()
    (build_dir / "index.html").write_text("<html></html>", encoding="utf-8")
    source = tmp_path / "header.js"
    source.write_text("/* header */\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="No rendered HTML pages"):
        inject_site_assets(build_dir, source)


def test_publication_paths_inject_assets_and_set_rtd_base_url() -> None:
    rtd = (REPO_ROOT / ".readthedocs.yml").read_text(encoding="utf-8")
    workflow = (REPO_ROOT / ".github/workflows/ci.yml").read_text(encoding="utf-8")

    assert 'BASE_URL="/$READTHEDOCS_LANGUAGE/$READTHEDOCS_VERSION"' in rtd
    assert rtd.index("scripts/inject_site_assets.py") < rtd.index(
        'cp -r _build/html/. "$READTHEDOCS_OUTPUT/html/"'
    )
    assert workflow.count("scripts/inject_site_assets.py") >= 2


def test_header_resets_before_the_missing_heading_guard() -> None:
    script = (REPO_ROOT / "docs/_static/synedu-header.js").read_text(encoding="utf-8")
    reset = script.index('nav.classList.remove("synedu-page-title-visible")')
    missing_heading_guard = script.index("if (!heading)")

    assert reset < missing_heading_guard


def test_desktop_navigation_is_not_hidden_by_the_scrolled_title() -> None:
    css = (REPO_ROOT / "docs/_static/synedu.css").read_text(encoding="utf-8")
    media = css.index("@media (max-width: 1023px)")
    hide_rule = css.index(".myst-top-nav.synedu-page-title-visible .myst-top-nav-item")

    assert hide_rule > media


def test_lesson_local_links_are_prefix_safe() -> None:
    for source in sorted((REPO_ROOT / "synedu").glob("S??/notebook.md")):
        text = source.read_text(encoding="utf-8")
        assert 'href="/docs/installing"' not in text
        assert 'href="../../docs/installing"' in text


def test_build_one_delegates_to_the_single_notebook_runner() -> None:
    makefile = (REPO_ROOT / "Makefile").read_text(encoding="utf-8")
    recipe = makefile.split("build-one:", 1)[1].split("\n\n", 1)[0]

    assert "$(MAKE) execute-notebook LESSON=$(LESSON)" in recipe
    assert "jupyter book build" not in recipe


def test_release_workflow_uses_explicit_python_setup() -> None:
    workflow = (REPO_ROOT / ".github/workflows/release.yml").read_text(encoding="utf-8")

    assert "actions/setup-python@v6" in workflow
    assert 'python-version: "3.11"' in workflow
    assert "uv python find" not in workflow
    assert workflow.count("set -o pipefail") == 3
