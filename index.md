# SynEdu

```{raw} html
<div class="synedu-home-hero synedu-home-hero--landing-clean">
  <div class="synedu-home-hero__grid synedu-home-hero__grid--clean">
    <div class="synedu-home-hero__content">
      <div class="synedu-home-kicker">
        <span class="synedu-dot"></span> SynEdu · executable talktorials
      </div>

      <div class="synedu-home-title">
        Learn reaction informatics with chemical graphs.
      </div>

      <div class="synedu-home-subtitle">
        Runnable notebooks for graph modeling, reaction rules, and reproducible workflows with RDKit and NetworkX.
      </div>

      <div class="synedu-home-actions">
        <a class="synedu-btn synedu-btn--primary" href="talktorials/">
          Browse talktorials
        </a>
        <a class="synedu-btn synedu-btn--light" href="installing/">
          Run locally
        </a>
        <a class="synedu-btn" href="api/">API</a>
        <a class="synedu-btn" href="contribute/">
          Contribute
        </a>
      </div>

      <div class="synedu-hero-inline">
        <span class="synedu-contact-pill">Graphs</span>
        <span class="synedu-contact-pill">Reaction rules</span>
        <span class="synedu-contact-pill">Executable notebooks</span>
        <span class="synedu-contact-pill">Reusable workflows</span>
      </div>
    </div>

    <div class="synedu-home-hero__logo-wrap">
      <img src="images/synedu_logo_rectangle_transparent.png"
           alt="SynEdu logo"
           class="synedu-home-hero__logo-img synedu-home-hero__logo-img--clean">
    </div>
  </div>
</div>
```

## Collections

```{raw} html
<div class="synedu-cardgrid synedu-cardgrid--home">
  <a class="synedu-card synedu-card--fund" href="talktorials/#fundamentals">
    <div class="synedu-card__top">
      <div class="synedu-card__code">S01-S02</div>
      <span class="synedu-badge synedu-badge--stable">Start here</span>
    </div>
    <div class="synedu-card__title">Fundamentals</div>
    <div class="synedu-card__desc">
      Molecules as graphs, morphisms, isomorphism, automorphisms.
    </div>
    <div class="synedu-meta">
      <span>Core concepts</span> · <span>Beginner</span>
    </div>
  </a>

  <a class="synedu-card synedu-card--lib" href="talktorials/#rule-library">
    <div class="synedu-card__top">
      <div class="synedu-card__code">S03-S07, S10</div>
      <span class="synedu-badge">Build</span>
    </div>
    <div class="synedu-card__title">Rule libraries</div>
    <div class="synedu-card__desc">
      Extract, canonicalize, verify, and curate reusable reaction rules.
    </div>
    <div class="synedu-meta">
      <span>Core workflow</span> · <span>Intermediate</span>
    </div>
  </a>

  <a class="synedu-card synedu-card--app" href="talktorials/#rule-application">
    <div class="synedu-card__top">
      <div class="synedu-card__code">S08-S09</div>
      <span class="synedu-badge">Apply</span>
    </div>
    <div class="synedu-card__title">Rule application</div>
    <div class="synedu-card__desc">
      Run transformations on datasets, validate outputs, and inspect failures.
    </div>
    <div class="synedu-meta">
      <span>Applied workflows</span> · <span>Advanced</span>
    </div>
  </a>
</div>
```

## Quickstart

```bash
conda env create -f environment.yml
conda activate synedu
jupyter lab
```

Open the notebooks in order from **S01** to **S10**.

## Documentation

- [Installation](installing.md)
- [API](api.md)
- [Talktorials](talktorials/index.md)

## License

SynEdu is released under the MIT License. See [license.md](license.md).
