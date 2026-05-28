// SynEdu small frontend helpers (kept minimal by design).
(function () {
  var storageKey = "synedu-theme";

  function setTheme(body, theme) {
    body.setAttribute("data-synedu-theme", theme);
    try {
      window.localStorage.setItem(storageKey, theme);
    } catch (error) {
      // Storage may be unavailable in strict browser modes.
    }
  }

  function storedTheme() {
    try {
      return window.localStorage.getItem(storageKey);
    } catch (error) {
      return null;
    }
  }

  function scrollTargetTop() {
    return document.documentElement || document.body;
  }

  function pageScrollHeight() {
    return Math.max(
      document.body ? document.body.scrollHeight : 0,
      document.documentElement ? document.documentElement.scrollHeight : 0
    );
  }

  function pageIsScrollable() {
    return pageScrollHeight() > window.innerHeight + 240;
  }

  function createScrollControl() {
    var control = document.createElement("div");
    var topButton = document.createElement("button");
    var bottomButton = document.createElement("button");

    control.className = "synedu-scroll-control";
    control.setAttribute("aria-label", "Page scroll controls");

    topButton.className = "synedu-scroll-control__button";
    topButton.type = "button";
    topButton.title = "Back to top";
    topButton.setAttribute("aria-label", "Back to top");
    topButton.textContent = "↑";

    bottomButton.className = "synedu-scroll-control__button";
    bottomButton.type = "button";
    bottomButton.title = "Go to bottom";
    bottomButton.setAttribute("aria-label", "Go to bottom");
    bottomButton.textContent = "↓";

    topButton.addEventListener("click", function () {
      window.scrollTo({ top: 0, behavior: "smooth" });
    });

    bottomButton.addEventListener("click", function () {
      window.scrollTo({ top: pageScrollHeight(), behavior: "smooth" });
    });

    control.appendChild(topButton);
    control.appendChild(bottomButton);
    document.body.appendChild(control);

    function syncVisibility() {
      if (!pageIsScrollable()) {
        control.classList.remove("is-visible");
        return;
      }
      control.classList.add("is-visible");
    }

    syncVisibility();
    window.addEventListener("resize", syncVisibility);
    window.addEventListener("load", syncVisibility);
    window.setTimeout(syncVisibility, 600);
    window.setTimeout(syncVisibility, 1600);
  }

  function injectTalktorialNavigation() {
    if (!/\/talktorials\/(?:S\d+|index|all_talktorials)\.html$/.test(window.location.pathname)) {
      return;
    }

    var primary = document.querySelector(".md-sidebar--primary .md-nav--primary");
    if (!primary || primary.querySelector(".synedu-talktorial-nav")) {
      return;
    }

    var links = [
      ["S01.html", "S01", "From Molecules to Labeled Graphs"],
      ["S02.html", "S02", "Graph Morphism in Reaction Informatics"],
      ["S03.html", "S03", "Maximum Common Substructure"],
      ["S04.html", "S04", "Atom Mapping as Graph Morphism"],
      ["S05.html", "S05", "Reaction Rules as Graph Rewriting"],
      ["S06.html", "S06", "Canonicalizing Mapped Reactions and Rules"],
      ["S07.html", "S07", "From Atom-Mapped Reactions to DPO Rules"],
      ["S08.html", "S08", "One-Step Reaction Prediction"],
      ["S09.html", "S09", "Context Graph Expansion"]
    ];

    var current = window.location.pathname.split("/").pop() || "index.html";
    var nav = document.createElement("nav");
    var title = document.createElement("label");
    var list = document.createElement("ul");

    nav.className = "md-nav synedu-talktorial-nav";
    nav.setAttribute("aria-label", "Talktorial series");
    title.className = "md-nav__title synedu-talktorial-nav__title";
    title.textContent = "Talktorials";
    list.className = "md-nav__list synedu-talktorial-nav__list";

    links.forEach(function (item) {
      var li = document.createElement("li");
      var a = document.createElement("a");
      var code = document.createElement("span");
      var text = document.createElement("span");

      li.className = "md-nav__item synedu-talktorial-nav__item";
      a.className = "md-nav__link synedu-talktorial-nav__link";
      if (item[0] === current) {
        a.className += " is-active";
      }
      a.href = item[0];
      code.className = "synedu-talktorial-nav__code";
      code.textContent = item[1];
      text.className = "synedu-talktorial-nav__text";
      text.textContent = item[2];

      a.appendChild(code);
      a.appendChild(text);
      li.appendChild(a);
      list.appendChild(li);
    });

    nav.appendChild(title);
    nav.appendChild(list);
    primary.appendChild(nav);
  }

  function openExternalLinksInNewTab() {
    var links = document.querySelectorAll("a[href^='http://'], a[href^='https://']");
    for (var index = 0; index < links.length; index += 1) {
      if (links[index].hostname && links[index].hostname !== window.location.hostname) {
        links[index].setAttribute("target", "_blank");
        links[index].setAttribute("rel", "noopener noreferrer");
      }
    }
  }

  document.addEventListener("DOMContentLoaded", function () {
    var body = document.body;
    var initial = storedTheme();
    var main = document.querySelector(".md-content__inner[role='main']");
    if (main && !main.id) {
      main.id = "synedu-main";
    }
    setTheme(body, initial || "light");
    createScrollControl();
    injectTalktorialNavigation();
    openExternalLinksInNewTab();

    var toggle = document.querySelector(".synedu-theme-toggle");
    if (!toggle) {
      return;
    }
    toggle.addEventListener("click", function () {
      var next = body.getAttribute("data-synedu-theme") === "dark" ? "light" : "dark";
      setTheme(body, next);
    });
  });
})();
