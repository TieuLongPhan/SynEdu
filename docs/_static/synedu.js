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

  function pageIsScrollable() {
    return scrollTargetTop().scrollHeight > window.innerHeight + 240;
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
      window.scrollTo({ top: scrollTargetTop().scrollHeight, behavior: "smooth" });
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
