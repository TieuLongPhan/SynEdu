(() => {
  const siteTitle = "SynEdu";
  let current = {};
  let frame = 0;
  let observedPath = "";
  let initialiseTimer = 0;

  const findHeading = () =>
    document.querySelector(
      "main article h1.myst-fm-block-title, main article h1"
    );

  const pageTitle = (heading) => {
    const headingText = heading?.querySelector(".heading-text");
    return (headingText?.textContent || heading?.textContent || "")
      .replace(/\s*¶\s*$/, "")
      .trim();
  };

  const update = () => {
    frame = 0;
    const { nav, label, heading, title } = current;
    if (!nav || !label || !heading || !document.contains(heading)) return;

    const navBottom = nav.getBoundingClientRect().bottom;
    const headingBottom = heading.getBoundingClientRect().bottom;
    const showPageTitle =
      title &&
      title !== siteTitle &&
      window.scrollY > 8 &&
      headingBottom <= navBottom + 8;

    if (nav.classList.contains("synedu-page-title-visible") === showPageTitle) {
      return;
    }

    nav.classList.toggle("synedu-page-title-visible", showPageTitle);
    label.textContent = showPageTitle ? title : siteTitle;
    label.title = showPageTitle ? title : "";
    label.setAttribute("aria-label", showPageTitle ? title : siteTitle);
  };

  const scheduleUpdate = () => {
    if (!frame) frame = window.requestAnimationFrame(update);
  };

  const initialise = () => {
    initialiseTimer = 0;
    const nav = document.querySelector(".myst-top-nav");
    const home = nav?.querySelector(".myst-home-link");
    const label = home?.querySelector(":scope > span");
    const heading = findHeading();
    observedPath = window.location.pathname;
    if (!nav || !home || !label) {
      current = {};
      return;
    }

    label.classList.add("synedu-header-label");
    nav.classList.remove("synedu-page-title-visible");
    label.textContent = siteTitle;
    label.title = "";
    label.setAttribute("aria-label", siteTitle);
    home.setAttribute("aria-label", "SynEdu home");
    if (!heading) {
      current = { nav, home, label, heading: null, title: "" };
      return;
    }

    const title = pageTitle(heading);
    current = { nav, home, label, heading, title };
    scheduleUpdate();
  };

  const scheduleInitialise = () => {
    window.clearTimeout(initialiseTimer);
    initialiseTimer = window.setTimeout(initialise, 0);
  };

  window.addEventListener("scroll", scheduleUpdate, { passive: true });
  window.addEventListener("resize", scheduleUpdate, { passive: true });
  window.addEventListener("popstate", scheduleInitialise);

  const observer = new MutationObserver(() => {
    if (
      window.location.pathname !== observedPath ||
      findHeading() !== current.heading
    ) {
      scheduleInitialise();
    }
  });

  observer.observe(document.body, { childList: true, subtree: true });
  initialise();
})();
