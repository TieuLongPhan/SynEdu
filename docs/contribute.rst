Contribute
==========

.. raw:: html

   <div class="synedu-home-hero synedu-home-hero--compact synedu-home-hero--contribute">
     <div class="synedu-hero-main">
       <div class="synedu-home-kicker"><span class="synedu-dot"></span> Community contributions</div>
       <div class="synedu-home-title">Help improve SynEdu.</div>
       <div class="synedu-home-subtitle">
         Contribute notebooks, tooling, lightweight utilities, and feedback that improve
         clarity, reproducibility, and usability.
       </div>
     </div>
   </div>

.. raw:: html

   <div class="synedu-contrib-grid">
     <div class="synedu-contrib-card">
       <div class="synedu-contrib-card__icon">📘</div>
       <div class="synedu-contrib-card__body">
         <div class="synedu-contrib-card__title">Talktorials</div>
         <div class="synedu-contrib-card__desc">
           Add notebooks, refine explanations, and create exercises or quizzes.
         </div>
       </div>
     </div>

     <div class="synedu-contrib-card">
       <div class="synedu-contrib-card__icon">🛠️</div>
       <div class="synedu-contrib-card__body">
         <div class="synedu-contrib-card__title">Infrastructure</div>
         <div class="synedu-contrib-card__desc">
           Improve docs, CI, formatting, packaging, and reproducible workflows.
         </div>
       </div>
     </div>

     <div class="synedu-contrib-card">
       <div class="synedu-contrib-card__icon">🧩</div>
       <div class="synedu-contrib-card__body">
         <div class="synedu-contrib-card__title">Utilities</div>
         <div class="synedu-contrib-card__desc">
           Add small helpers that reduce notebook boilerplate without expanding the API unnecessarily.
         </div>
       </div>
     </div>

     <div class="synedu-contrib-card">
       <div class="synedu-contrib-card__icon">🔍</div>
       <div class="synedu-contrib-card__body">
         <div class="synedu-contrib-card__title">Feedback</div>
         <div class="synedu-contrib-card__desc">
           Report issues, suggest clarifications, and help improve the learning experience.
         </div>
       </div>
     </div>
   </div>

Workflow
--------

1. Fork the repository and create a feature branch.
2. Run tests and checks locally.
3. Open a pull request with a clear description.
4. Include screenshots when visual or documentation changes are involved.

.. code-block:: bash

   pytest
   flake8

Guidelines
----------

- Keep notebooks runnable from top to bottom.
- Prefer small, composable functions in library code.
- Move reusable logic into ``synedu`` instead of duplicating it in notebooks.
- Keep helper utilities lightweight and teaching-focused.