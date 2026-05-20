SynEco
======

SynEdu is the educational entry point into the broader ``SynEco`` software
ecosystem for reaction informatics. The talktorials introduce core concepts in a
teaching-first format, then connect those ideas to production-oriented tools once
a workflow becomes useful beyond the notebook setting.

.. raw:: html

   <!-- icon snippets reused below -->
   <svg style="display:none" xmlns="http://www.w3.org/2000/svg">
     <symbol id="ico-book" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
       <path d="M2 3h6a4 4 0 0 1 4 4v14a3 3 0 0 0-3-3H2z"/><path d="M22 3h-6a4 4 0 0 0-4 4v14a3 3 0 0 1 3-3h7z"/>
     </symbol>
     <symbol id="ico-pkg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
       <path d="M21 16V8a2 2 0 0 0-1-1.73l-7-4a2 2 0 0 0-2 0l-7 4A2 2 0 0 0 3 8v8a2 2 0 0 0 1 1.73l7 4a2 2 0 0 0 2 0l7-4A2 2 0 0 0 21 16z"/>
       <polyline points="3.27 6.96 12 12.01 20.73 6.96"/><line x1="12" y1="22.08" x2="12" y2="12"/>
     </symbol>
     <symbol id="ico-gh" viewBox="0 0 24 24" fill="currentColor">
       <path d="M12 2C6.477 2 2 6.484 2 12.017c0 4.425 2.865 8.18 6.839 9.504.5.092.682-.217.682-.483 0-.237-.008-.868-.013-1.703-2.782.605-3.369-1.343-3.369-1.343-.454-1.158-1.11-1.466-1.11-1.466-.908-.62.069-.608.069-.608 1.003.07 1.531 1.032 1.531 1.032.892 1.53 2.341 1.088 2.91.832.092-.647.35-1.088.636-1.338-2.22-.253-4.555-1.113-4.555-4.951 0-1.093.39-1.988 1.029-2.688-.103-.253-.446-1.272.098-2.65 0 0 .84-.27 2.75 1.026A9.564 9.564 0 0 1 12 6.844a9.59 9.59 0 0 1 2.504.337c1.909-1.296 2.747-1.027 2.747-1.027.546 1.379.202 2.398.1 2.651.64.7 1.028 1.595 1.028 2.688 0 3.848-2.339 4.695-4.566 4.943.359.309.678.92.678 1.855 0 1.338-.012 2.419-.012 2.747 0 .268.18.58.688.482A10.02 10.02 0 0 0 22 12.017C22 6.484 17.522 2 12 2z"/>
     </symbol>
   </svg>

   <div class="synedu-cardgrid">
     <a class="synedu-card synedu-card--fund" href="https://github.com/TieuLongPhan/SynEdu" target="_blank" rel="noopener">
       <div class="synedu-card__top">
         <span class="synedu-card__code">SynEdu</span>
         <span class="synedu-badge"><svg><use href="#ico-gh"/></svg>GitHub</span>
       </div>
       <p class="synedu-card__desc">This talktorial series — molecular graphs, reaction rules, and prediction workflows.</p>
     </a>
     <a class="synedu-card synedu-card--fund" href="https://synkit.readthedocs.io/en/latest/" target="_blank" rel="noopener">
       <div class="synedu-card__top">
         <span class="synedu-card__code">SynKit</span>
         <span class="synedu-badge"><svg><use href="#ico-book"/></svg>docs</span>
       </div>
       <p class="synedu-card__desc">Core toolkit for chemical graph handling, reaction rules, and reaction-informatics utilities.</p>
     </a>
     <a class="synedu-card synedu-card--lib" href="https://pypi.org/project/synrbl/" target="_blank" rel="noopener">
       <div class="synedu-card__top">
         <span class="synedu-card__code">SynRBL</span>
         <span class="synedu-badge"><svg><use href="#ico-pkg"/></svg>PyPI</span>
       </div>
       <p class="synedu-card__desc">Rule-based reaction balancing and preprocessing workflows.</p>
     </a>
     <a class="synedu-card synedu-card--lib" href="https://syntemp.readthedocs.io/en/latest/" target="_blank" rel="noopener">
       <div class="synedu-card__top">
         <span class="synedu-card__code">SynTemp</span>
         <span class="synedu-badge"><svg><use href="#ico-book"/></svg>docs</span>
       </div>
       <p class="synedu-card__desc">Template- and rule-oriented tooling for reaction modeling.</p>
     </a>
     <a class="synedu-card synedu-card--app" href="https://github.com/phuocchung123/SynCat" target="_blank" rel="noopener">
       <div class="synedu-card__top">
         <span class="synedu-card__code">SynCat</span>
         <span class="synedu-badge"><svg><use href="#ico-gh"/></svg>GitHub</span>
       </div>
       <p class="synedu-card__desc">Reaction categorization and related cheminformatics utilities.</p>
     </a>
     <a class="synedu-card synedu-card--app" href="https://synrxn.readthedocs.io/en/latest/" target="_blank" rel="noopener">
       <div class="synedu-card__top">
         <span class="synedu-card__code">SynRXN</span>
         <span class="synedu-badge"><svg><use href="#ico-book"/></svg>docs</span>
       </div>
       <p class="synedu-card__desc">Curated reaction datasets, benchmark tasks, standardized splits, and reproducible evaluation workflows.</p>
     </a>
   </div>

Learning relationship
---------------------

``SynEdu`` should remain readable as teaching material. The broader ``SynEco``
tools should be introduced only when they clarify a concept, reduce boilerplate,
or show how a notebook idea can become a reusable research workflow.

In practice, this means that a talktorial may begin with transparent,
step-by-step code written for learning, then point to the corresponding ecosystem
package when the same idea is needed in a scalable or reproducible setting. For
example, reaction preprocessing concepts can lead naturally to ``SynRBL``,
chemical graph and rule manipulation to ``SynKit`` or ``SynTemp``, reaction
categorization to ``SynCat``, and benchmark construction or dataset evaluation
to ``SynRXN``.