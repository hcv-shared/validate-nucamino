(defun pytest (root-dir buffer-name)
    (compile "python3.6 -m unittest && date"))

(defun run-pytest ()
  (interactive)
  (pytest "/home/nathaniel/tmp/validate-nucamino" "*pytest*"))

(global-set-key (kbd "<f5>") 'run-pytest)
