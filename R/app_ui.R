#' Application UI
#'
#' @return A \code{bslib::page_navbar} object.
#' @export
app_ui <- function() {

  # -- Resource paths --------------------------------------------------------
  shiny::addResourcePath(
    "pga_www",
    system.file("app/www", package = "pgacmdr")
  )
  shiny::addResourcePath(
    "sdb",
    find.package("shinydashboard")
  )

  # -- Footer logos (bundled in package, served as base64 data URIs) ---------
  logos <- pga_logo_uris()

  # -- Example data for Welcome tab ------------------------------------------

  data_one_col <- data.frame(
    Population = c("Boulouparis", "Boulouparis", "Gadji", "PortLaguerre", "Sarramea"),
    B12 = c("192/194", "200/200", "0/0", "145/145", "0/0"),
    C07 = c("145/192", "179/179", "92/100", "92/92", "92/92"),
    stringsAsFactors = FALSE
  )

  data_two_col <- data.frame(
    Population = c("Boulouparis", "Gadji", "PortLaguerre", "Sarramea"),
    B12   = c(192, 200,   0, 145),
    B12_1 = c(194, 200,   0, 145),
    C07   = c(145, 179,  92,  92),
    C07_1 = c(192, 179,  92, 100),
    stringsAsFactors = FALSE
  )

  # -- Theme -----------------------------------------------------------------

  # Wes Anderson palettes — full reference:
  #  Royal1:        #FAEFD1                          (cream page bg)
  #  IsleofDogs2:   #EAD3BF #AA9486 #B6854D #39312F #1C1718
  #  IsleofDogs1:   #9986A5 #79402E #CCBA72 #D9D0D3 #8D8680
  #  AsteroidCity1: #0A9F9D #CEB175 #E54E21 #6C8645 #C18748 (navbar)
  #  FantasticFox1: #DD8D29 #E2D200 #46ACC8 #E58601 #B40F20 (box accent, download btn)
  #  Zissou1:       #3B9AB2 #78B7C5 #EBCC2A #E1AF00 #F21A00 (significance gradient)
  #  Moonrise2:     #798E87 #C27D38 #CCC591 #29211F
  #  GrandBudapest2: #C6CDF7 (dark-mode links)
  pga_theme <- bslib::bs_theme(
    version   = 5,
    bg        = "#f5f7fa",   # light grey page background
    fg        = "#333a43",   # dark grey
    primary   = "#333a43",   # dark grey primary
    secondary = "#8ea1b9",   # medium grey-blue
    "navbar-bg"            = "#333a43",  # dark grey navbar
    "navbar-color"         = "#FFFFFF",
    "navbar-hover-color"   = "#8ea1b9",
    "navbar-active-color"  = "#FFFFFF",
    "navbar-brand-color"   = "#FFFFFF",
    "nav-link-font-size"   = "0.85rem",
    "nav-link-font-weight" = "500"
  )

  # -- Custom CSS (light + dark mode) ----------------------------------------

  custom_css <- shiny::HTML("

    /* ===== BODY & LAYOUT ===== */
    body {
      font-family: 'Helvetica Neue', 'Segoe UI', Arial, sans-serif;
      font-size: 16px;
      color: #333a43;
      background-color: #f5f7fa;
    }
    p, label, .form-label, .shiny-input-container label,
    .selectize-input, .form-control, .form-select {
      font-size: 15px;
      color: #333a43;
    }
    /* Many hint/status paragraphs throughout the app were set with small
       inline font sizes (11-13px) — raise the floor so they stay readable. */
    .box p, .shiny-html-output p {
      font-size: 14px !important;
    }

    .bslib-page-navbar > .tab-content {
      padding: 20px;
      min-height: calc(100vh - 56px);
    }
    /* Welcome page is now short (module grid/help/citation removed) —
       don't force it to fill the viewport, or a big empty gap appears
       above the partner-logos footer. */
    body.pga-welcome-mode .bslib-page-navbar > .tab-content {
      min-height: auto;
    }

    /* ===== WES ANDERSON COLOR MAP =====
       IsleofDogs2: #FAEFD1(Royal1-bg) #EAD3BF #AA9486 #B6854D #39312F #1C1718
       IsleofDogs1: #9986A5  #79402E  #CCBA72  #D9D0D3  #8D8680
       Moonrise2:   #798E87  #C27D38  #CCC591  #29211F
       Royal2:      #9A8822  #F5CDB4  #F8AFA8  #FDDDA0  #74A089
       GrandBudapest2: #C6CDF7 (dark-mode link only)
    ===== */

    /* ===== BOX — white cards on light-grey page ===== */
    .box {
      border-radius: 4px !important;
      border-top:   none !important;
      border-left:  none !important;
      border-right: none !important;
      border-bottom: 2px solid #333a43 !important;
      box-shadow: 0 1px 6px rgba(57,49,47,0.10) !important;
      margin-bottom: 18px;
      background: #FFFFFF;                          /* white card */
    }
    .box-header {
      background: #FFFFFF !important;
      border-bottom: 1px solid #E8E0D6 !important;
      padding: 0 !important;                        /* title div provides its own padding */
      border-radius: 4px 4px 0 0 !important;
    }
    /* Fallback: plain title without inner div — restore padding */
    .box-header .box-title {
      padding: 10px 16px;
      display: block;
    }
    .box-title {
      font-size: 15px !important;
      font-weight: 600 !important;
      color: #333a43 !important;
    }
    .box-body { padding: 14px 16px !important; }

    /* h-tags */
    .box h2 {
      font-size: 17px; font-weight: 700;
      color: #333a43;
      border-bottom: 1px solid #8ea1b9;
      padding-bottom: 8px; margin-top: 16px;
    }
    .box h3 { font-size: 15px; font-weight: 600; color: #333a43; }
    .box h4 { font-size: 14px; font-weight: 600; color: #333a43; }

    a         { color: #8ea1b9; }
    a:hover   { color: #333a43; }

    /* ===== SECTION TITLE (module h2 with class=section-title) ===== */
    .section-title {
      color: #333a43 !important;
      border-bottom: 2px solid #8ea1b9 !important;
    }

    /* ===== NAVBAR — dark grey #333a43 ===== */
    .navbar {
      padding-top: 0; padding-bottom: 0; min-height: 52px;
      background-color: #333a43 !important;
      border-bottom: 1px solid #1e242b;
    }
    .navbar-brand { padding: 6px 12px 6px 0; }
    .navbar .nav-link {
      padding-top: 14px !important; padding-bottom: 14px !important;
      color: #FFFFFF !important;
    }
    .navbar .nav-link:hover  { color: #8ea1b9 !important; }
    .navbar .nav-link.active {
      color: #FFFFFF !important;
      font-weight: 600;
      border-bottom: 2px solid #8ea1b9;
    }
    .nav-intertryp { opacity: 0.75; font-size: 0.85rem; }
    .nav-intertryp:hover { opacity: 1; }
    .bslib-input-dark-mode {
      font-size: 1rem; padding: 14px 10px; cursor: pointer;
      color: #FFFFFF !important;
    }

    /* ===== valueBox significance gradient =====
       maroon = B40F20 dark red  (FantasticFox1 — p < 0.001)
       orange = E1AF00 amber     (Zissou1       — p < 0.05)
       aqua   = 3B9AB2 teal      (Zissou1       — not significant)
    ===== */
    .bg-maroon, .small-box.bg-maroon { background-color: #B40F20 !important; }
    .bg-maroon > .inner, .bg-maroon > .icon { background-color: #B40F20 !important; }
    .bg-orange, .small-box.bg-orange { background-color: #E1AF00 !important; }
    .bg-orange > .inner, .bg-orange > .icon { background-color: #E1AF00 !important; }
    .bg-aqua,   .small-box.bg-aqua   { background-color: #3B9AB2 !important; }
    .bg-aqua   > .inner, .bg-aqua   > .icon { background-color: #3B9AB2 !important; }

    /* ===== DARK MODE =====
       bg: IsleofDogs2 #1C1718, cards: #39312F, navbar: #29211F(Moonrise2)
       text: #EAD3BF, titles: #CCBA72(IsleofDogs1), links: #C6CDF7(GrandBudapest2)
    ===== */
    [data-bs-theme='dark'] body,
    [data-bs-theme='dark'] {
      background-color: #1C1718;                    /* IsleofDogs2 near-black */
      color: #EAD3BF;                               /* IsleofDogs2 parchment */
    }
    [data-bs-theme='dark'] .bslib-page-navbar > .tab-content {
      background-color: #1C1718;
    }
    /* Light mode page bg */
    .bslib-page-navbar > .tab-content {
      background-color: #f5f7fa;
    }
    [data-bs-theme='dark'] .navbar {
      background-color: #1e242b !important;
      border-bottom-color: #333a43 !important;
    }
    [data-bs-theme='dark'] .navbar .nav-link       { color: #8ea1b9 !important; }
    [data-bs-theme='dark'] .navbar .nav-link:hover { color: #FFFFFF !important; }
    [data-bs-theme='dark'] .navbar .nav-link.active {
      color: #FFFFFF !important;
      border-bottom-color: #8ea1b9;
    }
    [data-bs-theme='dark'] .bslib-input-dark-mode  { color: #8ea1b9 !important; }

    [data-bs-theme='dark'] .box {
      background: #39312F !important;               /* IsleofDogs2 dark brown */
      color: #EAD3BF !important;
      border-bottom-color: #8D8680 !important;      /* IsleofDogs1 taupe */
      box-shadow: 0 1px 5px rgba(0,0,0,0.50) !important;
    }
    [data-bs-theme='dark'] .box-header {
      background: #39312F !important;               /* dark mode keeps dark brown */
      border-bottom-color: #8D8680 !important;
    }
    [data-bs-theme='dark'] .box-title { color: #CCBA72 !important; }
    [data-bs-theme='dark'] .box p,
    [data-bs-theme='dark'] .box li,
    [data-bs-theme='dark'] .box label { color: #D9D0D3; }
    [data-bs-theme='dark'] .box h2 {
      color: #CCBA72; border-bottom-color: #8D8680;
    }
    [data-bs-theme='dark'] .box h3 { color: #CCBA72; }
    [data-bs-theme='dark'] .box h4 { color: #EAD3BF; }
    [data-bs-theme='dark'] a       { color: #C6CDF7; }  /* GrandBudapest2 lavender */
    [data-bs-theme='dark'] a:hover { color: #CCBA72; }  /* IsleofDogs1 gold */

    [data-bs-theme='dark'] .form-control,
    [data-bs-theme='dark'] .form-select {
      background: #29211F; color: #EAD3BF; border-color: #8D8680;
    }
    [data-bs-theme='dark'] .selectize-input,
    [data-bs-theme='dark'] .selectize-dropdown {
      background: #29211F !important; color: #EAD3BF !important;
      border-color: #8D8680 !important;
    }
    [data-bs-theme='dark'] .well {
      background: #39312F; border-color: #8D8680; color: #D9D0D3;
    }

    /* ===== BUTTONS ===== */
    /* Action/compute — dark grey, white text */
    .btn-action-primary, .btn-action-primary.btn-primary {
      background-color: #333a43 !important; color: #FFFFFF !important;
      border: none !important; font-weight: 600; border-radius: 3px;
    }
    .btn-action-primary:hover { background-color: #1e242b !important; }

    /* All download buttons — dark grey */
    .btn-download-primary, .btn-download-secondary, .btn-download-info,
    .btn-primary, .btn-block.btn-primary,
    .btn-outline-primary {
      background-color: #333a43 !important; color: #FFFFFF !important;
      border: none !important; border-radius: 3px;
    }
    .btn-download-primary:hover, .btn-download-secondary:hover,
    .btn-download-info:hover, .btn-primary:hover, .btn-outline-primary:hover {
      background-color: #1e242b !important;
    }

    /* ===== NAV TABS & PILLS (inside modules) ===== */
    /* Active tab: dark grey bg, white text */
    .nav-tabs .nav-link.active,
    .nav-pills .nav-link.active {
      background-color: #333a43 !important;
      color: #FFFFFF !important;
      border-color: #333a43 !important;
    }
    /* Inactive tab: dark grey text, transparent bg */
    .nav-tabs .nav-link,
    .nav-pills .nav-link {
      color: #333a43 !important;
    }
    .nav-tabs .nav-link:hover,
    .nav-pills .nav-link:hover {
      color: #FFFFFF !important;
      background-color: #8ea1b9 !important;
      border-color: #8ea1b9 !important;
    }

    /* ===== WELCOME PAGE ===== */
    .pga-hero {
      background: linear-gradient(145deg, #1a2035 0%, #26306B 55%, #333a43 100%);
      padding: 35px 40px 35px;
      display: flex;
      align-items: center;
      justify-content: space-between;
      margin: -20px -20px 0px -20px;
    }
    .pga-hero-text {
      flex: 1;
      text-align: left;
      padding-right: 32px;
    }
    .pga-hero-logo {
      flex-shrink: 0;
      display: flex;
      align-items: center;
    }
    .pga-hero h1 {
      color: #FFFFFF !important;
      margin: 0 0 8px 0 !important;
      border: none !important;
    }
    .pga-cta {
      display: inline-block;
      background: #A9F0D8 !important;
      color: #1a2035 !important;
      font-weight: 700 !important;
      font-size: 1rem !important;
      padding: 11px 30px !important;
      border-radius: 4px !important;
      border: none !important;
      cursor: pointer;
      text-decoration: none !important;
      transition: background 0.18s, transform 0.12s;
    }
    .pga-cta:hover {
      background: #FFFFFF !important;
      transform: translateY(-1px);
    }
    .pga-module-grid {
      display: grid;
      grid-template-columns: repeat(3, 1fr);
      gap: 14px;
      margin: 0 0 24px 0;
    }
    @media (max-width: 768px) {
      .pga-module-grid { grid-template-columns: 1fr 1fr; }
    }
    .pga-module-card {
      background: #FFFFFF;
      border-radius: 5px;
      padding: 18px 16px 14px;
      border-bottom: 3px solid #8ea1b9;
      box-shadow: 0 1px 5px rgba(51,58,67,0.08);
      transition: box-shadow 0.18s, transform 0.12s;
    }
    .pga-module-card:hover {
      box-shadow: 0 4px 14px rgba(51,58,67,0.14);
      transform: translateY(-2px);
    }
    .pga-module-card .card-icon {
      font-size: 1.5rem;
      margin-bottom: 8px;
      color: #6B64EF;
    }
    .pga-module-card h5 {
      font-size: 0.88rem !important;
      font-weight: 700 !important;
      color: #333a43 !important;
      margin: 0 0 5px 0 !important;
    }
    .pga-module-card p {
      font-size: 0.80rem !important;
      color: #666 !important;
      margin: 0 !important;
      line-height: 1.4 !important;
    }
    .pga-steps {
      counter-reset: step-counter;
      list-style: none;
      padding: 0;
      margin: 0;
    }
    .pga-steps li {
      counter-increment: step-counter;
      display: flex;
      align-items: flex-start;
      gap: 14px;
      margin-bottom: 12px;
      font-size: 14px;
    }
    .pga-steps li::before {
      content: counter(step-counter);
      min-width: 26px; height: 26px;
      background: #333a43;
      color: #FFFFFF;
      border-radius: 50%;
      display: flex; align-items: center; justify-content: center;
      font-size: 0.8rem; font-weight: 700;
    }
    /* ===== HELP PAGE ===== */
    .pga-help-section {
      margin-bottom: 28px;
    }
    .pga-help-section h3 {
      font-size: 1rem !important;
      font-weight: 700 !important;
      color: #333a43 !important;
      border-bottom: 2px solid #8ea1b9 !important;
      padding-bottom: 6px !important;
      margin-bottom: 14px !important;
    }
    .pga-format-note {
      background: #f5f7fa;
      border-left: 4px solid #6B64EF;
      border-radius: 3px;
      padding: 10px 14px;
      font-size: 13px;
      margin-bottom: 14px;
    }
    .pga-tip {
      background: #f5f7fa;
      border-left: 4px solid #2CBF9F;
      border-radius: 3px;
      padding: 10px 14px;
      font-size: 13px;
      margin-bottom: 10px;
    }

    /* Module banner (replaces h2 section-title in each tab) */
    .pga-module-banner {
      background: linear-gradient(135deg, #1a2035 0%, #26306B 60%, #333a43 100%);
      padding: 22px 40px 18px;
      display: flex;
      align-items: center;
      gap: 22px;
      margin: -20px -20px 28px -20px;
      border-bottom: 4px solid #6B64EF;
    }
    .pga-banner-icon {
      font-size: 2.8rem;
      color: rgba(255,255,255,0.18);
      flex-shrink: 0;
    }
    .pga-banner-title {
      color: #FFFFFF !important;
      font-size: 1.75rem !important;
      font-weight: 800 !important;
      margin: 0 0 5px 0 !important;
      border: none !important;
      letter-spacing: -0.2px;
    }
    .pga-banner-subtitle {
      color: #A9F0D8;
      font-size: 0.95rem;
      margin: 0 0 10px 0;
    }
    /* Method note: styled description below banner */
    .pga-method-note {
      background: #f8f9fb;
      border-left: 4px solid #6B64EF;
      border-radius: 4px;
      padding: 14px 20px;
      margin-bottom: 20px;
      font-size: 15px;
      line-height: 1.65;
      color: #2c3e50;
    }

    /* Secondary action */
    .btn-action-secondary {
      background-color: #8ea1b9 !important; color: #FFFFFF !important;
      border: none !important; border-radius: 3px;
    }
    .btn-action-secondary:hover { background-color: #333a43 !important; }

    /* ── Default: navbar hidden (welcome is entry page; sidebar shows it) ── */
    .bslib-page-navbar > nav.navbar { display: none !important; }

    /* ── Left sidebar mode: embedded flex layout ─────────────────────────── */
    body.pga-sidebar-mode.bslib-page-navbar {
      display: flex !important;
      flex-direction: row !important;
      align-items: stretch !important;
      height: 100vh !important;
      gap: 0 !important;
      padding: 0 !important;
      margin: 0 !important;
    }
    body.pga-sidebar-mode.bslib-page-navbar > nav.navbar {
      display: flex !important;
      flex: 0 0 210px !important;
      width: 210px !important;
      height: 100vh !important;
      position: sticky !important;
      top: 0 !important;
      align-self: flex-start !important;
      flex-direction: column !important;
      align-items: stretch !important;
      overflow-y: auto !important;
      overflow-x: hidden !important;
      scrollbar-width: none !important;       /* Firefox: hide sidebar scrollbar */
      z-index: 100 !important;
      background: linear-gradient(180deg, #1a2035 0%, #26306B 100%) !important;
      border-right: 1px solid rgba(107,100,239,0.25) !important;
      box-shadow: 3px 0 16px rgba(0,0,0,0.18) !important;
      padding: 0 !important;
      margin: 0 !important;
    }
    /* Chrome/Safari: hide sidebar scrollbar */
    body.pga-sidebar-mode.bslib-page-navbar > nav.navbar::-webkit-scrollbar {
      display: none !important;
    }
    /* Hide brand - PGA-cmdr label is now inside the banner */
    body.pga-sidebar-mode.bslib-page-navbar > nav.navbar .navbar-brand { display: none !important; }
    /* Collapse wrapper - full-width column */
    body.pga-sidebar-mode.bslib-page-navbar > nav.navbar .navbar-collapse {
      display: flex !important;
      flex-direction: column !important;
      width: 100% !important;
      height: 100% !important;
      padding: 6px 0 !important;
    }
    body.pga-sidebar-mode.bslib-page-navbar > nav.navbar .navbar-nav {
      flex-direction: column !important;
      width: 100% !important;
    }
    body.pga-sidebar-mode.bslib-page-navbar > nav.navbar .nav-link {
      padding: 9px 16px !important;
      text-align: left !important;
      border-radius: 0 !important;
      white-space: nowrap !important;
      font-size: 0.88rem !important;
      color: rgba(255,255,255,0.82) !important;
      border-left: 3px solid transparent !important;
    }
    body.pga-sidebar-mode.bslib-page-navbar > nav.navbar .nav-link:hover {
      background: rgba(255,255,255,0.08) !important;
      color: #FFFFFF !important;
      border-left-color: #A9F0D8 !important;
    }
    body.pga-sidebar-mode.bslib-page-navbar > nav.navbar .nav-link.active {
      background: rgba(107,100,239,0.25) !important;
      color: #FFFFFF !important;
      border-left-color: #6B64EF !important;
    }
    /* Dark mode toggle pushed to bottom of sidebar */
    body.pga-sidebar-mode.bslib-page-navbar > nav.navbar .navbar-nav.ms-auto {
      margin-top: auto !important;
      border-top: 1px solid rgba(255,255,255,0.1) !important;
      padding: 8px 0 !important;
    }
    /* Body locked to exactly viewport height; all scrolling happens inside
       the content column (.container-fluid > .tab-content). */
    body.pga-sidebar-mode {
      overflow: hidden !important;
      height: 100vh !important;
    }
    /* .container-fluid is the actual direct child of body (bslib wraps
       .tab-content inside it). Make it the flex child that takes the
       remaining width and clip its own overflow so the body never grows. */
    body.pga-sidebar-mode.bslib-page-navbar > .container-fluid {
      flex: 1 1 0 !important;
      min-width: 0 !important;
      padding: 0 !important;
      margin: 0 !important;
      height: 100vh !important;
      overflow: hidden !important;
    }
    /* tab-content is the scroll container; height is constrained so
       overflow-y:auto triggers a visible scrollbar when content overflows. */
    body.pga-sidebar-mode.bslib-page-navbar > .container-fluid > .tab-content {
      height: 100vh !important;
      overflow-y: auto !important;
      overflow-x: auto !important;
      padding-top: 0 !important;
    }
    /* Welcome page: thin overlay scrollbar — visible only while scrolling,
       never reserves a permanent gutter on Linux/Windows. */
    body.pga-welcome-mode {
      scrollbar-width: thin !important;
    }
    body.pga-welcome-mode::-webkit-scrollbar { width: 6px; }
    body.pga-welcome-mode::-webkit-scrollbar-track { background: transparent; }
    body.pga-welcome-mode::-webkit-scrollbar-thumb {
      background-color: rgba(51,58,67,0.20);
      border-radius: 3px;
    }
    body.pga-welcome-mode::-webkit-scrollbar-thumb:hover {
      background-color: rgba(51,58,67,0.50);
    }

    /* Tab pane: no top padding so banner is flush */
    body.pga-sidebar-mode.bslib-page-navbar .tab-content > .tab-pane {
      padding-top: 0 !important;
      margin-top: 0 !important;
    }
    /* Banner: sticky at top of its scroll container (.tab-content) */
    body.pga-sidebar-mode .pga-module-banner {
      position: sticky !important;
      top: 0 !important;
      z-index: 100 !important;
      margin-top: 0 !important;
      margin-left: -20px !important;
      margin-right: -20px !important;
    }
    /* Welcome footer logos */
    .pga-footer {
      background: #FFFFFF;
      margin: 0px -20px 0px -20px;
      padding: 12px 40px 10px;
      display: flex;
      flex-direction: column;
      align-items: center;
      gap: 20px;
    }
    .pga-footer-logos {
      display: flex;
      align-items: center;
      justify-content: center;
      gap: 44px;
      flex-wrap: wrap;
    }
    .pga-footer-sep {
      width: 1px;
      align-self: stretch;
      background: #d8dceb;
      margin: 0 6px;
    }
    .pga-footer img { max-height: 100px; }
    .pga-footer-caption {
      color: #6b7280;
      font-size: 17px;
      letter-spacing: 0.3px;
    }

  ")

  # -- Welcome tab content ---------------------------------------------------


  welcome_content <- shiny::tagList(

    # ── Hero ──────────────────────────────────────────────────────────────────
    shiny::div(
      class = "pga-hero",
      # Left: text + CTA
      shiny::div(
        class = "pga-hero-text",
        shiny::tags$h1(
          style = "font-size: 5.2rem; font-weight: 300; color: #F4F6FF; margin: 0 0 6px 0; letter-spacing: -1.5px; line-height: 1.05;",
          "PGA-commander V1"
        ),
        shiny::tags$div(
          style = "font-size: 5.2rem; font-weight: 600; letter-spacing: -2px; line-height: 1.05; margin: 0 0 18px 0; background: linear-gradient(90deg, #8F86FF, #5AA7FF); -webkit-background-clip: text; background-clip: text; color: transparent;",
          "PGA-cmdr V1"
        ),
        #
        shiny::div(
          class = "pga-credits",
          style = "margin: 10px 0 15px 0; padding: 10px 0; border-top: 1px solid rgba(255,255,255,0.08);",
          shiny::tags$p(
            style = "color: #A9F0D8; font-size: 25px; font-family: Inter,Segoe UI,Roboto,Helvetica,Arial,sans-serif; margin: 0; line-height: 1.5;",
            "A Versatile, user-friendly and multi-OS application", shiny::tags$br(),
          "to analyse population genetic data"
          )
        ),
        #
        shiny::div(
          class = "pga-credits",
          style = "margin: 10px 0 15px 0; padding: 10px 0; border-top: 1px solid rgba(255,255,255,0.08); border-bottom: 1px solid rgba(255,255,255,0.08);",
          shiny::tags$p(
            style = "color: #a0a8c0; font-size: 25px; font-family: Inter,Segoe UI,Roboto,Helvetica,Arial,sans-serif; margin: 0; line-height: 1.5;",
            shiny::tags$strong("Programming:", style = "color: #c8d0e8;"),
            " Naffiou Kadiri and Vincent Manzanilla",
            shiny::tags$br(),
            shiny::tags$strong("Conception and supervision:", style = "color: #c8d0e8;"),
            " Thierry de Meeûs"
          )
        ),
        shiny::tags$button(
          class   = "pga-cta",
          onclick = "document.querySelector('[data-value=\"import\"]').click();",
          shiny::icon("upload"), " Launch "
        )
      ),
      # Right: big logo — reduced ~2.7x from its previous size
      shiny::div(
        class = "pga-hero-logo",
        shiny::tags$img(
          src   = "pga_www/LogoPGAcmdr.png",
          height = "200px",
          alt   = "PGA-cmdr logo",
          style = "filter: drop-shadow(0 8px 32px rgba(0,0,0,0.55));"
        )
      )
    ),

    # ── Partner logos footer ──────────────────────────────────────────────────
    shiny::div(
      class = "pga-footer",
      shiny::div(
        class = "pga-footer-logos",
        # INTERTRYP
        shiny::tags$a(
          href = "https://umr-intertryp.cirad.fr/en", target = "_blank",
          style = "text-decoration:none;",
          if (!is.null(logos$intertryp))
            shiny::tags$img(src = logos$intertryp, height = "90px",
              alt = "INTERTRYP", title = "Hosts, Vectors and Infectious Agents",
              style = "opacity:0.9;")
          else
            shiny::tags$span(style = "font-weight:700;font-size:1.3rem;color:#26306B;letter-spacing:1px;", "INTERTRYP")
        ),
        # Université de Montpellier
        shiny::tags$a(
          href = "https://www.umontpellier.fr/en/", target = "_blank",
          style = "text-decoration:none;",
          if (!is.null(logos$montpellier))
            shiny::tags$img(src = logos$montpellier, height = "90px",
              alt = "Universite de Montpellier", title = "Universite de Montpellier",
              style = "opacity:0.9;")
          else
            shiny::tags$span(style = "font-weight:700;font-size:1.3rem;color:#26306B;letter-spacing:1px;", "UM")
        ),
        # CIRAD
        shiny::tags$a(
          href = "https://www.cirad.fr/en", target = "_blank",
          style = "text-decoration:none;",
          if (!is.null(logos$cirad))
            shiny::tags$img(src = logos$cirad, height = "90px",
              alt = "CIRAD", title = "Agricultural Research for Development",
              style = "opacity:0.9;")
          else
            shiny::tags$span(style = "font-weight:700;font-size:1.3rem;color:#26306B;letter-spacing:1px;", "CIRAD")
        ),
        # IRD
        shiny::tags$a(
          href = "https://www.ird.fr/en", target = "_blank",
          style = "text-decoration:none;",
          if (!is.null(logos$ird))
            shiny::tags$img(src = logos$ird, height = "100px",
              alt = "IRD", title = "Institut de Recherche pour le Developpement",
              style = "opacity:0.9;")
          else
            shiny::tags$span(style = "font-weight:700;font-size:1.3rem;color:#26306B;letter-spacing:1px;", "IRD")
        ),
        # UCAD
        shiny::tags$a(
          href = "https://www.ucad.sn", target = "_blank",
          style = "text-decoration:none;",
          if (!is.null(logos$ucad))
            shiny::tags$img(src = logos$ucad, height = "90px",
              alt = "UCAD", title = "Universite Cheikh Anta Diop de Dakar",
              style = "opacity:0.9;")
          else
            shiny::tags$span(style = "font-weight:700;font-size:1.3rem;color:#26306B;letter-spacing:1px;", "UCAD")
        )
      ),
      shiny::tags$div(
        class = "pga-footer-caption",
        "Intertryp, Univ Montpellier, Cirad, IRD, Montpellier, France"
      )
    )
  )

  ###
  # --- Help tab content ---
  ###
  help_content <- shiny::tagList(
    module_banner("question-circle", "Help · Documentation",""),
    shiny::fluidRow(
      shinydashboard::box(
        width = 12, solidHeader = FALSE,
        title = shiny::div(
          style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
          shiny::icon("file-alt"), " Data requirements"
        ),

        # Format overview
        shiny::div(
          class = "pga-help-section",
          shiny::tags$h3(shiny::icon("table"), " Accepted file formats"),
          shiny::div(
            class = "pga-format-note",
            shiny::icon("info-circle"),
            shiny::HTML(" <strong>PGA-cmdr accepts CSV and tabulation-delimited TXT files.</strong>
              The file must contain at least one <em>population</em> column and one or more
              <em>genetic marker</em> columns. Additional metadata columns (e.g. latitude,
              longitude, individual ID) are optional and are simply ignored — the app only
              needs Population, optionally Latitude/Longitude, and the loci range.")
          ),
          shiny::tags$ul(
            style = "font-size:14px; line-height:1.9;",
            shiny::tags$li(shiny::tags$strong("Header row:"), " required (first row = column names)."),
            shiny::tags$li(shiny::tags$strong("Separator:"), " comma (,), semicolon (;), or tab (\\t) — pick the one your file uses; the file loads automatically as soon as it is chosen, no extra click needed."),
            shiny::tags$li(shiny::tags$strong("Population column:"), " any column name; one value per individual."),
            shiny::tags$li(shiny::tags$strong("Marker columns:"), " single-column (", shiny::tags$code("192/194"), ") or two-column (", shiny::tags$code("192"), " + ", shiny::tags$code("194"), ") encoding — detected automatically (see below)."),
            shiny::tags$li(shiny::tags$strong("Missing data:"), " coded as ", shiny::tags$code("0"), " or ", shiny::tags$code("0/0"), " by default — customisable on import."),
            shiny::tags$li(shiny::tags$strong("Upload limit:"), " 500 MB.")
          )
        ),

        shiny::tags$hr(),

        # Table examples
        shiny::div(
          class = "pga-help-section",
          shiny::tags$h3(shiny::icon("columns"), " Encoding formats"),
          shiny::fluidRow(
            shiny::column(6,
              shiny::div(
                class = "pga-format-note",
                shiny::icon("compress-alt"),
                shiny::HTML(" <strong>Single-column</strong> &mdash; both alleles in one cell, separated by <code>/</code>, <code>-</code> or <code>_</code>. Each locus uses <strong>1</strong> raw column.")
              ),
              shiny::div(
                style = "overflow-x:auto;",
                shiny::HTML(kableExtra::kable_styling(
                  knitr::kable(data_one_col, "html", align = "l"),
                  full_width = FALSE, position = "left",
                  bootstrap_options = c("striped", "hover", "condensed", "bordered")
                ))
              )
            ),
            shiny::column(6,
              shiny::div(
                class = "pga-format-note",
                shiny::icon("expand-alt"),
                shiny::HTML(" <strong>Paired-column</strong> &mdash; each allele in its own column. The first allele column carries the locus name (e.g. <code>B12</code>); the second must be named <code>B12_1</code> (or <code>B12.1</code>) and sit immediately next to it. Each locus uses <strong>2</strong> raw columns.")
              ),
              shiny::div(
                style = "overflow-x:auto;",
                shiny::HTML(kableExtra::kable_styling(
                  knitr::kable(data_two_col, "html", align = "l",
                               col.names = c("Population","B12","B12_1","C07","C07_1")),
                  full_width = FALSE, position = "left",
                  bootstrap_options = c("striped", "hover", "condensed", "bordered")
                ))
              )
            )
          ),  # end fluidRow columns

          shiny::tags$hr(style = "margin: 10px 0;"),

          # --- Column assignment on Import Data
          shiny::div(
            class = "pga-format-note",
            style = "background:#f0f4ff; border-left:4px solid #6B64EF; padding:10px 14px; border-radius:4px; font-size:13px; line-height:1.7;",
            shiny::tags$p(
              style = "margin-top:0; font-weight:600;",
              shiny::icon("magic"), " Setting the loci range on the Import Data page"
            ),
            shiny::tags$p(
              "Once the file is loaded, a numbered ", shiny::tags$strong("Column reference"),
              " table appears at the bottom of the page, listing every column name next to its position \u2014",
              " use it to check exactly which number corresponds to which column before filling in the fields above."
            ),
            shiny::tags$ol(
              style = "margin-bottom:0;",
              shiny::tags$li(shiny::HTML(
                "<strong>Population / Latitude / Longitude</strong> &mdash; pick the matching column from each dropdown (Latitude/Longitude are optional, needed only for the IBD module's automatic geographic-distance calculation)."
              )),
              shiny::tags$li(shiny::HTML(
                "<strong>Column number of the first locus</strong> &mdash; the column where your markers start (check the Column reference table)."
              )),
              shiny::tags$li(shiny::HTML(
                "<strong>Number of loci</strong> &mdash; the <em>actual number of genetic markers</em> (not raw columns). ",
                "As soon as the first locus column is set, this field is <strong>auto-suggested</strong> from how many complete loci fit between that column and the end of the file, accounting for single- or paired-column format automatically. Adjust it down if some trailing columns aren't loci."
              ))
            ),
            shiny::tags$p(
              style = "margin: 8px 0 0;",
              "The app shows the resulting column range live (e.g. \u201c6 loci \u00d7 2 columns/locus = columns 9 to 20\u201d) so you can verify it against the Column reference table before clicking ",
              shiny::tags$strong("Apply"), "."
            )
          )
        ),

        shiny::tags$hr(),

        # Tips
        shiny::div(
          class = "pga-help-section",
          shiny::tags$h3(shiny::icon("lightbulb"), " Tips & common pitfalls"),
          shiny::div(class = "pga-tip", shiny::icon("check-circle"),
            shiny::HTML(" <strong>Use the default dataset</strong> to familiarise yourself with the expected format before uploading your own data.")),
          shiny::div(class = "pga-tip", shiny::icon("check-circle"),
            shiny::HTML(" <strong>Check the Column reference table</strong> on the Import Data page whenever you're unsure which column number corresponds to which field \u2014 it's always visible once a file is loaded.")),
          shiny::div(class = "pga-tip", shiny::icon("check-circle"),
            shiny::HTML(" <strong>Null alleles (FreeNA)</strong> use the standard missing-data code you set on import (e.g. <code>0/0</code>); the Null Alleles module handles the EM estimation from there.")),
          shiny::div(class = "pga-tip", shiny::icon("check-circle"),
            shiny::HTML(" <strong>Every module downloads one plain .txt file (Allele Freq, Local/Global Panmixia, Subdivision, Diversities, LD) or a zip of several files (General Stats, Null Alleles, IBD, Mantel)</strong> as soon as you click Run \u2014 there is nothing to view or export separately beforehand; the exact output file name(s) are always shown just above the Run button.")),
          shiny::div(class = "pga-tip", shiny::icon("exclamation-triangle"),
            shiny::HTML(" <strong>All locus columns must use the same encoding</strong> (single- or paired-column) within a file. Mixed encoding is not supported.")),
          shiny::div(class = "pga-tip", shiny::icon("exclamation-triangle"),
            shiny::HTML(" <strong>Population codes must be consistent</strong> across rows &mdash; trailing spaces or capitalisation differences will create duplicate populations.")),
          shiny::div(class = "pga-tip", shiny::icon("exclamation-triangle"),
            shiny::HTML(" <strong>Bootstrap over sub-samples (or over individuals for within-population FIS) needs at least 5 units</strong> to be meaningful; below that threshold the app reports <code>NA</code> instead of an unreliable confidence interval, in Null Alleles, Local/Global Panmixia, Subdivision, and Diversities."))
        ),

        shiny::tags$hr(),

        # Mac installation note
        shiny::div(
          class = "pga-help-section",
          shiny::tags$h3(shiny::icon("apple"), " macOS: prerequisites for installation"),
          shiny::div(
            class = "pga-format-note",
            shiny::icon("exclamation-triangle"),
            shiny::HTML(" <strong>macOS users only.</strong> PGA-cmdr contains C++ code compiled with OpenMP.
              Apple's default clang does <em>not</em> include OpenMP or gfortran.
              Install the following before running <code>remotes::install_github()</code>:")
          ),
          shiny::tags$ol(
            style = "font-size:14px; line-height:2.0;",
            shiny::tags$li(shiny::HTML(
              "<strong>gfortran</strong> &mdash; required for R package compilation on macOS.<br>
               Download the official R-project gfortran from
               <a href='https://mac.r-project.org/tools/' target='_blank' style='color:#6B64EF;'>mac.r-project.org/tools</a>
               and install <code>gfortran-14.2-universal.pkg</code> (or the current version)."
            )),
            shiny::tags$li(shiny::HTML(
              "<strong>libomp</strong> &mdash; OpenMP runtime library, via Homebrew:<br>
               <code>brew install libomp</code>"
            )),
            shiny::tags$li(shiny::HTML(
              "(Apple Silicon only) Add to <code>~/.R/Makevars</code> (create the file if it does not exist):<br>
               <code>LDFLAGS += -L/opt/homebrew/opt/libomp/lib -lomp</code><br>
               <code>CPPFLAGS += -I/opt/homebrew/opt/libomp/include -Xclang -fopenmp</code>"
            ))
          )
        ),

        shiny::tags$hr(),

        # Workflow
        shiny::div(
          class = "pga-help-section",
          shiny::tags$h3(shiny::icon("route"), " Recommended workflow"),
          shiny::tags$ol(
            style = "font-size:14px; line-height:2.0;",
            shiny::tags$li(shiny::HTML("<strong>Import Data</strong> &mdash; load the file, assign Population/Latitude/Longitude, set the loci range (auto-suggested), and click Apply.")),
            shiny::tags$li(shiny::HTML("<strong>Allele Freq</strong> &mdash; check allele frequencies and missing-data rates per locus and population.")),
            shiny::tags$li(shiny::HTML("<strong>General Stats</strong> &mdash; obtain Ho, Hs, Ht and FIS/FST/FIT (Weir &amp; Cockerham) per locus, plus per-population (Ho and Hs weighted by the number of individuals typed at each locus, as in Genepop) and per-allele detail.")),
            shiny::tags$li(shiny::HTML("<strong>Null Alleles</strong> &mdash; if a high Ho/Hs ratio is suspected, estimate null allele frequencies and get the ENA-corrected FST.")),
            shiny::tags$li(shiny::HTML("<strong>Local Panmixia</strong> &mdash; test for HWE (FIS) by locus and by population, with bootstrap over individuals, sub-samples and loci; one file.")),
            shiny::tags$li(shiny::HTML("<strong>Subdivision</strong> &mdash; estimate FST (bootstrap over sub-samples and over loci, permutation p-values) and the G-based test, in one file.")),
            shiny::tags$li(shiny::HTML("<strong>Diversities</strong> &mdash; obtain HS/HT confidence intervals (bootstrap over individuals, sub-samples and loci) in one file.")),
            shiny::tags$li(shiny::HTML("<strong>LD</strong> &mdash; test pairwise linkage disequilibrium across loci.")),
            shiny::tags$li(shiny::HTML("<strong>IBD</strong> &mdash; test isolation by distance: Rousset's (1997) regression and an independent Mantel test (requires GPS data for the automatic geographic distance)."))
          )
        ),

        shiny::tags$hr(),

        # Statistical methods & references
        shiny::div(
          class = "pga-help-section",
          shiny::tags$h3(shiny::icon("book-open"), " Statistical methods & key references"),
          shiny::fluidRow(
            shiny::column(6,
              shiny::div(
                style = "padding:12px 14px; background:#f8f9fc; border-radius:6px; border-left:3px solid #6B64EF; margin-bottom:12px;",
                shiny::tags$strong("F-statistics (WC84)"),
                shiny::tags$p(style = "font-size:13px; margin:6px 0 0; line-height:1.7;",
                  "FIS, FST and FIT are estimated following the unbiased moment estimators of ",
                  shiny::tags$strong("Weir & Cockerham (1984)"), ". These estimators are robust to unequal sample sizes across populations and loci. Confidence intervals are computed by non-parametric bootstrap over loci and over sub-samples; p-values by permutation."),
                shiny::tags$p(style = "font-size:12px; margin:4px 0 0; color:#777;",
                  "Weir BS, Cockerham CC. 1984. Estimating F-statistics for the analysis of population structure. Evolution 38:1358-1370.")
              ),
              shiny::div(
                style = "padding:12px 14px; background:#f8f9fc; border-radius:6px; border-left:3px solid #2CBF9F; margin-bottom:12px;",
                shiny::tags$strong("Gene diversity (Nei & Chesser 1983)"),
                shiny::tags$p(style = "font-size:13px; margin:6px 0 0; line-height:1.7;",
                  "HS (within-population gene diversity) and HT (total gene diversity) use the ",
                  shiny::tags$strong("unbiased small-sample estimator of Nei & Chesser (1983)"),
                  " (Hs = n/(n\u22121) \u00d7 [1 \u2212 \u03a3p\u00b2 \u2212 Ho/(2n)]), not the simpler biased 1\u2212\u03a3p\u00b2 form. Per-locus and multilocus estimates are reported, with bootstrap CI derived by resampling over individuals, sub-samples, or loci. Resampling individuals lowers HS by about HS/(2n), so those intervals are shifted by (observed HS \u2212 bootstrap mean); HS per population is averaged over loci weighted by the number of individuals typed (as in Genepop)."),
                shiny::tags$p(style = "font-size:12px; margin:4px 0 0; color:#777;",
                  "Nei M, Chesser RK. 1983. Estimation of fixation indices and gene diversities. Ann Hum Genet 47:253-259.")
              ),
              shiny::div(
                style = "padding:12px 14px; background:#f8f9fc; border-radius:6px; border-left:3px solid #EBCC2A; margin-bottom:12px;",
                shiny::tags$strong("Linkage disequilibrium"),
                shiny::tags$p(style = "font-size:13px; margin:6px 0 0; line-height:1.7;",
                  "Pairwise LD is tested for every locus pair using a log-likelihood ratio (G) test on the genotypic contingency table within each population; genotypes at one locus are permuted and the observed G-statistic is compared to the permutation distribution, per population and combined (\u201cAll\u201d). Significant LD between physically unlinked loci may reflect selection, admixture, or small sample size."),
                shiny::tags$p(style = "font-size:12px; margin:4px 0 0; color:#777;",
                  "Rousset F. 2008. Genepop'007: a complete re-implementation of the Genepop software. Mol Ecol Resour 8:103-106.")
              )
            ),
            shiny::column(6,
              shiny::div(
                style = "padding:12px 14px; background:#f8f9fc; border-radius:6px; border-left:3px solid #B40F20; margin-bottom:12px;",
                shiny::tags$strong("Null alleles (FreeNA / EM algorithm)"),
                shiny::tags$p(style = "font-size:13px; margin:6px 0 0; line-height:1.7;",
                  "Null allele frequencies are estimated per locus \u00d7 sub-sample by the Expectation-Maximisation (EM) algorithm implemented in ",
                  shiny::tags$strong("FreeNA"), " (Chapuis & Estoup 2007). The algorithm iterates between estimating null allele frequency from observed genotype counts and updating expected genotype frequencies until convergence. Null alleles inflate apparent FIS and bias FST downward; the module reports both raw and ENA-corrected FST/DCSE."),
                shiny::tags$p(style = "font-size:12px; margin:4px 0 0; color:#777;",
                  "Chapuis MP, Estoup A. 2007. Microsatellite null alleles and estimation of population differentiation. Mol Biol Evol 24:621-631.")
              ),
              shiny::div(
                style = "padding:12px 14px; background:#f8f9fc; border-radius:6px; border-left:3px solid #9986A5; margin-bottom:12px;",
                shiny::tags$strong("Hardy-Weinberg equilibrium (FIS)"),
                shiny::tags$p(style = "font-size:13px; margin:6px 0 0; line-height:1.7;",
                  "HWE departure within populations is quantified via FIS (WC84), by locus or by population. Permutation p-values are obtained by randomising alleles within samples. Significant positive FIS indicates excess homozygosity (inbreeding, null alleles, Wahlund effect); negative FIS indicates heterozygote excess."),
                shiny::tags$p(style = "font-size:12px; margin:4px 0 0; color:#777;",
                  "Guo SW, Thompson EA. 1992. Performing the exact test of Hardy-Weinberg proportion for multiple alleles. Biometrics 48:361-372.")
              ),
              shiny::div(
                style = "padding:12px 14px; background:#f8f9fc; border-radius:6px; border-left:3px solid #3B9AB2; margin-bottom:12px;",
                shiny::tags$strong("Interpretation thresholds"),
                shiny::tags$p(style = "font-size:13px; margin:6px 0 0; line-height:1.7;",
                  shiny::tags$strong("FST:"), " <0.05 little differentiation; 0.05-0.15 moderate; 0.15-0.25 great; >0.25 very great (Wright 1978).", shiny::tags$br(),
                  shiny::tags$strong("FIS:"), " p<0.05 significant HWE departure. Bonferroni correction recommended for multiple loci.", shiny::tags$br(),
                  shiny::tags$strong("Null alleles:"), " frequency >0.10 at a locus warrants removal or correction before FST estimation.", shiny::tags$br(),
                  shiny::tags$strong("LD:"), " p-values are not corrected by default; apply FDR or Bonferroni for multiple comparisons.")
              )
            )
          )
        )
      )
    )
  )
  ###
  # --- Assemble page ---
  ###
  bslib::page_navbar(
    id    = "main_nav",
    title = "PGA-cmdr",
    theme          = pga_theme,
    window_title   = "PGA-cmdr",
    navbar_options = bslib::navbar_options(collapsible = TRUE),
    fillable       = FALSE,

    # Inject AdminLTE CSS for shinydashboard boxes + our custom overrides
    header = shiny::tags$head(
      shiny::tags$title("PGA-cmdr"),
      shiny::tags$link(
        rel  = "stylesheet",
        href = "sdb/AdminLTE/AdminLTE.min.css"
      ),
      shiny::tags$link(
        rel  = "stylesheet",
        href = "sdb/shinydashboard.css"
      ),
      shiny::tags$style(custom_css),
      shiny::tags$script(shiny::HTML("
(function () {
  var nav, lastPane, frameId;

  function init() {
    nav = document.querySelector('nav.navbar');
    if (!nav) { setTimeout(init, 150); return; }
    document.addEventListener('shown.bs.tab', schedule, true);
    var mo = new MutationObserver(function (muts) {
      if (muts.some(function (m) {
        return m.attributeName === 'class' &&
               m.target.classList.contains('tab-pane') &&
               m.target.classList.contains('active');
      })) schedule();
    });
    mo.observe(document.body, { subtree: true, attributes: true, attributeFilter: ['class'] });
    schedule();
  }

  function schedule() {
    if (frameId) cancelAnimationFrame(frameId);
    frameId = requestAnimationFrame(doLayout);
  }

  function doLayout() {
    frameId = null;
    var pane = document.querySelector('.tab-pane.active');
    if (!pane || pane === lastPane) return;
    lastPane = pane;
    if (pane.querySelector('.pga-hero')) {
      document.body.classList.remove('pga-sidebar-mode');
      document.body.classList.add('pga-welcome-mode');
    } else {
      document.body.classList.remove('pga-welcome-mode');
      document.body.classList.add('pga-sidebar-mode');
    }
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    setTimeout(init, 0);
  }
})();
      "))
    ),
    ###
    # --- Module tabs ---
    ###
    bslib::nav_panel(
      title = "Welcome",
      icon  = shiny::icon("home"),
      value = "welcome",
      welcome_content
    ),
    bslib::nav_panel(
      title = "Import Data",
      icon  = shiny::icon("upload"),
      value = "import",
      import_data_ui("import")
    ),
    bslib::nav_panel(
      title = "Allele Freq",
      icon  = shiny::icon("chart-pie"),
      value = "allele_frequencies",
      ui_allele_frequencies("allele")
    ),
    bslib::nav_panel(
      title = "General Stats",
      icon  = shiny::icon("table"),
      value = "general_stats",
      mod_general_stats_ui("general_stats")
    ),
    bslib::nav_panel(
      title = "Local Panmixia",
      icon  = shiny::icon("flask"),
      value = "local_panmixia",
      mod_local_panmixia_ui("general_stats")
    ),
    bslib::nav_panel(
      title = "Global Panmixia",
      icon  = shiny::icon("globe"),
      value = "global_panmixia",
      mod_global_panmixia_ui("general_stats")
    ),
    bslib::nav_panel(
      title = "Subdivision",
      icon  = shiny::icon("sitemap"),
      value = "subdivision",
      mod_subdivision_ui("general_stats")
    ),
    bslib::nav_panel(
      title = "Diversities",
      icon  = shiny::icon("chart-line"),
      value = "genetic_diversities",
      mod_genetic_diversities_ui("general_stats")
    ),
    bslib::nav_panel(
      title = "LD",
      icon  = shiny::icon("link"),
      value = "linkage_desequilibrium",
      linkage_desequilibrium_UI("ld")
    ),
    bslib::nav_panel(
      title = "Null Alleles",
      icon  = shiny::icon("circle-notch"),
      value = "null_alleles",
      null_alleles_UI("null_alleles")
    ),
    bslib::nav_panel(
      title = HTML('IBD'),
      icon  = shiny::icon("atom"),
      value = "isolation_by_distance",
      isolation_by_distance_UI("ibd")
    ),
    bslib::nav_panel(
      title = "Help",
      icon  = shiny::icon("question-circle"),
      value = "help",
      help_content
    ),
    ###
    # --- Right-side controls ---
    ###
    bslib::nav_spacer(),
    bslib::nav_item(
      bslib::input_dark_mode(id = "color_mode")
    )
  )
}