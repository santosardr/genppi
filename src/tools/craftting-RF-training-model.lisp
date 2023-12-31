;; The random-forest.lisp can be found at https://github.com/santosardr/non-CSPs/tree/main/src/lisp
;; The random-forest.lisp depends on cl-random-forest and cl-store installed via quicklisp.lisp
;; Start SBCL using this OS command: sbcl --dynamic-space-size 12288 --control-stack-size 20 --load random-forest.lisp

(defvar arff-data (random-forest:read-arff-file "training.arff"))

(defparameter *target* (make-array (length (getf arff-data :states_num)) :element-type 'fixnum :initial-contents (getf arff-data :states_num)))

(defparameter *datamatrix* (make-array (list (length (getf arff-data :states_num)) (length (getf arff-data :attributes)) ) :element-type 'single-float :initial-contents (getf arff-data :data) ))

(defparameter *n-class* 2)

(defparameter *forest* (CL-RANDOM-FOREST::make-forest *n-class* *datamatrix* *target* :n-tree 500 :bagging-ratio 1.0 :max-depth (truncate (/ (+ 1 (length (getf arff-data :attributes))) 2)) :min-region-samples 1 :n-trial 10))

(random-forest:serialize-forest *forest* "model.dat")
