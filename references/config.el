;; load contrib library
(add-to-list 'load-path
             (expand-file-name "/Users/dpath2o/PHD/references/"
                               (file-name-directory
				org-find-library-dir "org")))

;; manage citations
(require 'org-bibtex)

;; export citations
(require 'ox-bibtex)
(setq org-bibtex-file "refs.org")
