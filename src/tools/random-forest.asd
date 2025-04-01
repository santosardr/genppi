(asdf:defsystem :random-forest
  :serial t
  :depends-on (:cl-random-forest :cl-store)
  :components ((:file "random-forest")))
