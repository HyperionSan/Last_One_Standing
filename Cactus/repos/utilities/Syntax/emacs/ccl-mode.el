;;; ccl-mode.el --- Cactus ccl file editing commands for Emacs
;;; $Header$

;; Author: Tom Goodale <goodale@aei-potsdam.mpg.de>
;;   with some modifications by Ian Hawke <hawke@aei.mpg.de>
;;   and Erik Schnetter <schnetter@aei.mpg.de>
;; Maintainer: Erik Schnetter <schnetter@aei.mpg.de>
;; Keywords: languages, faces, ccl, cactus

;;; Commentary:

;; A major mode for editing Cactus ccl files 
;; (param.ccl,interface.ccl,schedule.ccl).
;; Originally derived from m4-mode.el by 
;; Andrew Csillag (drew@staff.prodigy.com)
;; Additional code derived from f90.el (indentation) by
;;  Torbj\"orn Einarsson <T.Einarsson@clab.ericsson.se>
;; and c-mode.el and perl-mode.el (electric brace stuff)
;; which are part of the standard emacs distribution.

;; At the moment this isn't very flexible. The relevant entries
;; in my .emacs file are
;;
;; (require 'ccl-mode)
;; (setq auto-mode-alist (append auto-mode-alist
;;                       (list '("\\.ccl$" . ccl-mode))))
;; (setq ccl-auto-newline t)
;;
;; and, as I feel that a colour screen should be used to the full,
;; I have font-lock mode on as well.
;;
;; ccl-auto-newline is similar but not identical to the way that
;; auto-newline works in C mode. The difference is that it doesn't
;; jump to the next line after the closing brace.
;;

;;; Code:

;;thank god for make-regexp.el!
;;(make-regexp '("cctk_real" "cctk_int" "cctk_char" "implements" "inherits" "type" "gf" "array" "scalar"))

;; Matching the following (hasn't Cactus got a lot of reserved words):
;; ("accumulator-base" "accumulator" "after" "array" "as" "at" "before"
;;  "cctk_char" "cctk_fpointer" "cctk_int" "cctk_pointer" "cctk_real"
;;  "cctk_string" "compact" "else" "friend" "function" "gf" "if" "implements" "in" 
;;  "include" "includes (source|header)" "inherits" "lang" "language" "out"
;;  "private" "protected" "provides" "public" "requires" "restricted" "scalar" 
;;  "schedule" "shares" "storage" "tags" "thorns" "timelevels" "triggers" 
;;  "type" "uses" "void" "while" "with"
(defvar ccl-font-lock-keywords
  `(("\\# \\(.*\\)"  1 font-lock-comment-face t)
    ("\\$[0-9*#@]" . font-lock-variable-name-face)
;; Type declarations ("^[ \t]*\\(boolean\\|\\(cctk_\\)\?\\(char\\|fpointer\\|int\\|pointer\\(_to_const\\)?\\|real\\|string\\)\\|keyword\\|void\\)\\b" . font-lock-constant-face)
;; interface.ccl header
    ("^[ \t]*\\(friend\\|i\\(mplements\\|nherits\\)\\)\\b" . font-lock-warning-face)
;; interface.ccl include files
    ("^[ \t]*\\(\\(uses \\)\?\\(includes\?\\( header\\| source\\)\?\\)\\)\\b" . font-lock-constant-face)
;; interface.ccl function aliasing FUNCTION statements
    ("^[ \t]*\\(\\(uses \\|provides \\|requires \\)\?function\\)\\b" . font-lock-keyword-face)
;; interface.ccl function aliasing details
    ("\\b\\(with\\|language\\)\\b" . font-lock-builtin-face)
;; interface.ccl access
    ("^[ \t]*\\(p\\(rivate\\|rotected\\|ublic\\)\\)\\b" . font-lock-warning-face)
;; interface.ccl block lines
    ("\\b\\compact|\\(di\\(m\\|strib\\)\\|ghostsize\\|s\\(ize\\|tagger\\)\\|t\\(ags\\|imelevels\\|ype\\)\\)\\b" . font-lock-keyword-face)
;; interface.ccl variable types
    ("\\b\\(array\\|constant\\|gf\\|scalar\\)\\b" . font-lock-type-face)
;; param.ccl access types. Note "private" is already defined.
    ("^[ \t]*\\(global\\|restricted\\|shares\\)\\b" . font-lock-warning-face)
;; param.ccl block lines.
    ("\\b\\(a\\(ccumulator\\(-base\\)\?\\|lways\\|s\\)\\|extends\\|never\\|recover\\|steerable\\|uses\\)\\b" . font-lock-keyword-face)
;; schedule.ccl; give the schedule statement special status
    ("^[ \t]*\\(schedule\\)\\b" . font-lock-function-name-face)
;; schedule.ccl; give conditionals special status.
    ("\\b\\(else\\|if\\)\\b" . font-lock-warning-face)
;; schedule.ccl; default schedule bins.
    ("\\b\\(cctk_\\)\?\\(startup\\|recover_parameters\\|paramcheck\\|basegrid\\|initial\\|postinitial\\|recover_variables\\|post_recover_variables\\|cpinitial\\|checkpoint\\|prestep\\|evol\\|postregrid\\|postregridinitial\\|postrestrict\\|postrestrictinitial\\|poststep\\|preregrid\\|preregridinitial\\|analysis\\|terminate\\|shutdown\\|wragh\\)\\b" . font-lock-constant-face)
;; schedule.ccl block lines. Note "as" is already defined.
    ("\\b\\(after\\|at\\|before\\|group\\|in\\|while\\)\\b" . font-lock-keyword-face)
;; schedule.ccl routine options.
    ("\\b\\(global\\|lang\\|level\\|map\\|meta\\|options\?\\|storage\\|sync\\(hronise\\)\?\\|triggers\?\\)\\b" . font-lock-keyword-face)
    "default font-lock-keywords")
;; configuration.ccl thorn options.
  ("\\b\\(\\(requires\\|uses\\)[ \t]+thorns\\)\\b" . font-lock-constant-face)
)

;;Where should the top level group go? CCL is practically perl anyway...

(defgroup ccl nil
  "CCL mode"
  :group 'perl)

(defgroup ccl-indent nil
  "CCL indentation"
  :prefix "ccl-"
  :group 'ccl)

(defcustom ccl-block-indent 2
  "*Indentation applied inside ccl blocks."
  :type 'integer
  :group 'ccl-indent)

(defcustom ccl-comment-region "# "
  "*String insterted by \\[ccl-comment-region]\
 at start of each line in region."
  :type 'string
  :group 'ccl-indent)

(defcustom ccl-auto-newline nil
  "*Non-nil means automatically newline before and after braces."
  :type 'boolean
  :group 'ccl-indent)

;; A temporary position to make region operators faster
(defvar ccl-cache-position nil)
(make-variable-buffer-local 'ccl-cache-position)

(defvar ccl-mode-syntax-table nil
  "syntax table used in ccl mode")
(setq ccl-mode-syntax-table (make-syntax-table))
(modify-syntax-entry ?\" "\""  ccl-mode-syntax-table)
(modify-syntax-entry ?\' "\""  ccl-mode-syntax-table)
(modify-syntax-entry ?#  "<\n" ccl-mode-syntax-table)
(modify-syntax-entry ?\n ">#"  ccl-mode-syntax-table)
;;(modify-syntax-entry ?\( "."   ccl-mode-syntax-table)
;;(modify-syntax-entry ?\) "."   ccl-mode-syntax-table)
(modify-syntax-entry ?\( "("   ccl-mode-syntax-table)
(modify-syntax-entry ?\) ")"   ccl-mode-syntax-table)
(modify-syntax-entry ?*  "."   ccl-mode-syntax-table)
(modify-syntax-entry ?:  "."   ccl-mode-syntax-table)
;;(modify-syntax-entry ?_  "_"   ccl-mode-syntax-table)
(modify-syntax-entry ?_  "w"   ccl-mode-syntax-table)

(defvar ccl-mode-map
  (let ((map (make-sparse-keymap)))
    (define-key map "\C-c\C-c" 'ccl-comment-region)
    (define-key map "\t"       'ccl-indent-line)
    (define-key map "\M-q"     'ccl-indent-region)
    (define-key map "{"        'electric-ccl-brace)
    (define-key map "}"        'electric-ccl-brace-noadvance)
    map))

;;;###autoload
(defun ccl-mode ()
  "A major-mode to edit ccl files like param.ccl,interface.ccl or schedule.ccl.
\\{ccl-mode-map}
"
  (interactive)
  (kill-all-local-variables)
  (use-local-map ccl-mode-map)

  (make-local-variable 'comment-start)
  (setq comment-start "#")
  (make-local-variable 'parse-sexp-ignore-comments)
  (setq parse-sexp-ignore-comments t)

  (make-local-variable 'indent-line-function)
  (setq indent-line-function 'ccl-indent-line)
  (make-local-variable 'indent-region-function)
  (setq indent-region-function 'ccl-indent-region)

  (make-local-variable	'font-lock-defaults)  
  (setq major-mode 'ccl-mode)
  (setq mode-name "CCL")
  (setq font-lock-defaults `(ccl-font-lock-keywords nil t nil nil))
  (set-syntax-table ccl-mode-syntax-table)
  (run-hooks 'ccl-mode-hook))

(defun electric-ccl-brace (arg)    
  "Insert character and correct line's indentation."
  (interactive "P")
  (let (insertpos)
    (if (and (not arg)
	     (eolp)
	     (or (save-excursion
		   (skip-chars-backward " \t")
		   (bolp))
		 (if ccl-auto-newline (progn (ccl-indent-line) (newline) t) nil)))
	(progn
	  (insert last-command-char)
	  (ccl-indent-line)
	  (if ccl-auto-newline
	      (progn
		(newline)
		;; (newline) may have done auto-fill
		(setq insertpos (- (point) 2))
		(ccl-indent-line)))
	  (save-excursion
	    (if insertpos (goto-char (1+ insertpos)))
	    (delete-char -1))))
    (if insertpos
	(save-excursion
	  (goto-char insertpos)
	  (self-insert-command (prefix-numeric-value arg)))
      (self-insert-command (prefix-numeric-value arg)))))

(defun electric-ccl-brace-noadvance (arg)    
  "Insert character and correct line's indentation."
  (interactive "P")
  (let (insertpos)
    (if (and (not arg)
	     (eolp)
	     (or (save-excursion
		   (skip-chars-backward " \t")
		   (bolp))
		 (if ccl-auto-newline (progn (ccl-indent-line) (newline) t) nil)))
	(progn
	  (insert last-command-char)
	  (ccl-indent-line))
      (if insertpos
	  (save-excursion
	    (goto-char insertpos)
	    (self-insert-command (prefix-numeric-value arg)))
	(self-insert-command (prefix-numeric-value arg))))))


(defun ccl-comment-region (beg-region end-region)
  "Comment/uncomment every line in the region.
Insert ccl-comment-region at the beginning of every line in the region
or, if already present, remove it."
  (interactive "*r")
  (let ((end (make-marker)))
    (set-marker end end-region)
    (goto-char beg-region)
    (beginning-of-line)
    (if (looking-at (regexp-quote ccl-comment-region))
	(delete-region (point) (match-end 0))
      (insert ccl-comment-region))
    (while (and  (zerop (forward-line 1))
		 (< (point) (marker-position end)))
      (if (looking-at (regexp-quote ccl-comment-region))
	  (delete-region (point) (match-end 0))
	(insert ccl-comment-region)))
    (set-marker end nil)))

(defsubst ccl-get-beg-of-line ()
  (save-excursion (beginning-of-line) (point)))

(defsubst ccl-in-comment ()
  (let ((beg-pnt
	 (if (and ccl-cache-position (> (point) ccl-cache-position))
	     ccl-cache-position
	   (point-min))))
    (nth 4 (parse-partial-sexp beg-pnt (point)))))

(defsubst ccl-line-continued ()
  (save-excursion
    (let ((bol (ccl-get-beg-of-line)))
      (end-of-line)
      (while (ccl-in-comment)
	(search-backward "#" bol)
	(skip-chars-backward "#"))
      (skip-chars-backward " \t")
      (= (preceding-char) ?&))))

(defsubst ccl-current-indentation ()
  "Return indentation of current line."
  (save-excursion
    (beginning-of-line) (skip-chars-forward " \t")
    (current-column)))

(defsubst ccl-present-statement-cont ()
  "Return continuation properties of present statement."
  (let (pcont cont)
    (save-excursion
      (setq pcont (if (ccl-previous-statement) (ccl-line-continued) nil)))
    (setq cont (ccl-line-continued))
    (cond ((and (not pcont) (not cont)) 'single)
 	  ((and (not pcont) cont)       'begin)
 	  ((and pcont       (not cont)) 'end)
 	  ((and pcont       cont)       'middle)
 	  (t (error)))))

(defsubst ccl-indent-to (col)
  "Indent current line to column COL."
  (beginning-of-line)
  (delete-horizontal-space)
  (if (zerop (current-column))
      (indent-to col)
    (indent-to col 1)))

(defun ccl-indent-line ()
  "Indent current line as CCL code."
  (interactive)
  (let (indent (pos (make-marker)) (case-fold-search t))
    (set-marker pos (point))
    (beginning-of-line)
    (if (save-excursion (and (ccl-previous-statement) (ccl-line-continued)))
	(progn (skip-chars-forward " \t")))
    (setq indent (ccl-calculate-indent))
    (if (zerop (- indent (current-column)))
	nil
      (ccl-indent-to indent))
    ;; If initial point was within line's indentation,
    ;; position after the indentation.  Else stay at same point in text.
    (if (< (point) (marker-position pos))
	(goto-char (marker-position pos)))
    (set-marker pos nil)))


(defun ccl-indent-region (beg-region end-region)
  "Indent every line in region by forward parsing."
  (interactive "*r")
  (let ((end-region-mark (make-marker)) (save-point (point-marker)))
    (set-marker end-region-mark end-region)
    (goto-char beg-region)
    ;; first find a line which is not a continuation line or comment
    (beginning-of-line)
    (while (< (point) end-region-mark) (progn (ccl-indent-line) (forward-line 1)))
    (set-marker end-region-mark nil)
    (if (string-match "XEmacs" emacs-version)
	(zmacs-deactivate-region)
      (deactivate-mark))))

(defun ccl-calculate-indent ()
  "Calculate the indent column based on previous statements."
  (interactive)
  (let (icol cont (case-fold-search t) (pnt (point)))
    (save-excursion
      (if (not (ccl-previous-statement))
	  (setq icol 0)
	(setq cont (ccl-present-statement-cont))
	(if (eq cont 'end)
	    (while (not (eq 'begin nil))
	      (ccl-previous-statement)))
	(cond ((eq cont 'begin)
	       (setq icol (+ (ccl-current-indentation)
			     ccl-continuation-indent)))
	      ((eq cont 'middle) (setq icol(current-indentation)))
	      (t (setq icol (ccl-current-indentation))
		 (skip-chars-forward " \t")
		 (if (looking-at "[ \t]*$") 
		     (setq icol (ccl-calculate-indent)))
		 (if (looking-at ".*{")
		     (setq icol (+ icol ccl-block-indent)))
		 (goto-char pnt)
		 (beginning-of-line)
		 (if (looking-at "[ \t]*}")
		     (setq icol (- icol ccl-block-indent)))
		 ))))
    icol))

(defun ccl-previous-statement ()
  "Move point to beginning of the previous CCL statement.
Return nil if no previous statement is found."
  (interactive)
  (let (not-first-statement)
    (beginning-of-line)
    (while (and (setq not-first-statement (zerop (forward-line -1)))
		(looking-at "[ \t]*\\(#\\)")))
    not-first-statement))

(provide 'ccl-mode)

;;; ccl-mode.el ends here
