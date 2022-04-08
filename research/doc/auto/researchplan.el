(TeX-add-style-hook
 "researchplan"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "11pt" "a4paper")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("doclicense" "type={CC}" "modifier={by-nc-sa}" "version={3.0}") ("geometry" "top=1cm" "bottom=1cm" "left=1cm" "right=1cm") ("pdfpages" "final") ("fullpage" "myheadings") ("babel" "english") ("inputenc" "utf8") ("fontenc" "T1") ("microtype" "protrusion=true" "expansion=true") ("caption" "font=small" "labelfont=bf") ("natbib" "round")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art11"
    "doclicense"
    "geometry"
    "pdflscape"
    "pdfpages"
    "fullpage"
    "enumerate"
    "parskip"
    "dashrule"
    "enumitem"
    "sectsty"
    "latexsym"
    "framed"
    "pbox"
    "ifdraft"
    "babel"
    "fourier"
    "shadowtext"
    "inputenc"
    "fontenc"
    "microtype"
    "textcomp"
    "soul"
    "lipsum"
    "blindtext"
    "graphicx"
    "float"
    "wrapfig"
    "rotating"
    "tikz"
    "tikzpagenodes"
    "caption"
    "pgfgantt"
    "amsmath"
    "bm"
    "amssymb"
    "marvosym"
    "mathpazo"
    "mathtools"
    "siunitx"
    "tabularx"
    "longtable"
    "array"
    "makecell"
    "totpages"
    "fancyhdr"
    "natbib"
    "url"
    "booktabs"
    "hyperref"
    "setspace"
    "titling"
    "xcolor")
   (TeX-add-symbols
    '("crule" ["argument"] 2)
    '("latlonDec" 2)
    '("HRule" 1)
    '("source" 1)
    "citationNeeded"
    "htwoo"
    "otwo"
    "cotwo"
    "chfour"
    "ntwo"
    "moc"
    "nadw"
    "sip"
    "dsw"
    "aabw"
    "asf"
    "asc"
    "acc"
    "cdp"
    "cdfi"
    "mc"
    "ais"
    "fastice"
    "Fastice"
    "wmo"
    "utas"
    "aapp"
    "imas"
    "csiro"
    "bom"
    "anu"
    "primarysupervisor"
    "primaryadvisor"
    "secondarysupervisor"
    "tertiarysupervisor"
    "quarternarysupervisor"
    "grc"
    "ORGTITLE"
    "ORGAUTHOR")
   (LaTeX-add-labels
    "fig:GEBCO_bathy"
    "sec:skills")
   (LaTeX-add-bibliographies
    "../../references/refs.bib"))
 :latex)

