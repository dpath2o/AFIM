(TeX-add-style-hook
 "ganttchart"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("standalone" "border=1pt")))
   (TeX-run-style-hooks
    "latex2e"
    "standalone"
    "standalone10"
    "pgfgantt"
    "graphicx"
    "xcolor"
    "datetime"))
 :latex)

