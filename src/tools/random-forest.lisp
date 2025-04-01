(defpackage :random-forest
  (:use :cl)
  (:export :file-or-directory
           :string-empty-p
           :split-sequence
           :read-arff-file
           :predict
           :train-predict
           :serialize-forest
           :deserialize-forest))

(in-package :random-forest)

(unless (find-package :cl-random-forest)
  (require :cl-random-forest))

(unless (find-package :cl-store)
  (require :cl-store ))

(defun file-or-directory (pathname)
  (if (uiop:directory-exists-p pathname)
      :directory
      (if (uiop:file-exists-p pathname)
          :regular-file
          nil)))

(defun string-empty-p (string)
  (zerop (length string)))
(defun split-sequence (separator sequence &optional remove-empty-subseqs)
  (let* ((result)
		(start 0)
		(end (position separator sequence :start start)))
	(loop while end do
		(push (subseq sequence start end) result)
		(setf start (1+ end))
		(setf end (position separator sequence :start start)))
	(push (subseq sequence start) result)
	(if remove-empty-subseqs
		(remove-if #'string-empty-p (reverse result))
	(reverse result))))

(defun read-arff-file (file-path)
  (with-open-file (file file-path :direction :input)
		(let* ((lines (loop for line = (read-line file nil)
					while line
					collect line))
		       (data-start-index (position "@data" lines
						   :test (lambda (x y) (string=
                                                            (string-upcase x)
                                                            (string-upcase y)))))
			(relation_csv (subseq (first lines) 1))
			(relation (cadr (split-sequence #\Space relation_csv)))
			(attributes_csv (subseq lines 1 data-start-index))
			
			(classes_str (reduce #'append  (mapcar (lambda (attribute)
								 (let* ((start-index (position #\{ attribute))
									(end-index (position #\} attribute :from-end t))
									(values (subseq attribute (1+ start-index) end-index)))
								   (split-sequence #\, values)))
							       ;;locate the class line                       
							       (remove-if-not
								(lambda (item)
								  (search "CLASS" (string-upcase item)))
								attributes_csv)
							       )))
			(classes (mapcar (lambda (class) (read-from-string class)) classes_str))
			
			(attributes (remove-if #'null (mapcar (lambda (att)
								(multiple-value-bind (text success)
										     (ignore-errors (read-from-string 
												     (cadr (split-sequence #\Space att))
												     );read
												    );ignore
										     (if success text 0)
										     );multiple
								);lambda
							      (butlast attributes_csv)
							      );mapcar
					       );remove
				    );attributes
			(data_csv (subseq lines (1+ data-start-index)))
			(data (mapcar (lambda (oneline)
					(butlast (mapcar (lambda (text) 
							(multiple-value-bind (number success) 
										(ignore-errors (read-from-string text)) 
										(if success number 0)))
							(split-sequence #\, oneline)))) 
				data_csv)
			)
			
			(datafloat (mapcar (lambda (oneline)
								(mapcar (lambda (number) (float number))  
										oneline
									)
								)
								data
						)
			)

			(states_nom (reduce #'append  (mapcar (lambda (oneline)
							(last (mapcar (lambda (text) 
									(multiple-value-bind (class success) 
												(ignore-errors (read-from-string text)) 
												(if success class 0)))
									(split-sequence #\, oneline)))) 
							data_csv)
					)
				    )
			(states_num (mapcar (lambda (class)  (position class classes) ) states_nom))
			);let*

			(list :relation relation
			:attributes attributes
			:data datafloat
			:states_nom states_nom
			:states_num states_num
			:classes classes))))

(defun predict (forest test-file negative-limit &optional classes)
  (let ( (TESTMATRIX) (ACCURACY) (CONFUSION) (FN) (FP) (SENSITIVITY)
	 (SPECIFICITY) (TEST-DATA) (TN) (TP) (Y-PRED)
	 (one) (two) )

  ;; Read the testing ARFF file
  (setf test-data (read-arff-file test-file ))

  (setf testmatrix (make-array
			      (list (length (getf test-data :states_num))
				    (length (getf test-data :attributes))
				    )
			      :element-type 'single-float
			      :initial-contents (getf test-data :data) ))

  ;; Predict the labels for the testing data
  (setf y-pred (loop for i from 0 to (- (length (getf test-data :states_num)) 1)
		       collect (CL-RANDOM-FOREST::predict-forest forest testmatrix i)))


  (setf confusion (list (list 0 0) (list 0 0));; Initialize the confusion matrix
	classes (if (null classes) (getf test-data :classes)  classes);; first NEG second POS
	one  (first classes)	
	two (second classes)	
	)
  (format t "Training classes: ~a~%"  classes)
  (format t "Test     classes: ~a~%" (getf test-data :classes) )
  ;; Compare the predicted labels with the actual labels and increment the confusion matrix accordingly
  (loop for i below (length y-pred) do
	(if (< i negative-limit)
	    (if (= (nth i y-pred) (position two classes) )
		(incf (first (first confusion))) ; True Negative (NEGATIVE predicted correctly)
	      (incf (second (first confusion)))) ; False Positive (POSITIVE predicted incorrectly)
	  (if (= (nth i y-pred) (position one classes))
	      (incf (second (second confusion))) ; True Positive (POSITIVE predicted correctly)
	    (incf (first (second confusion)))))) ; False Negative (NEGATIVE predicted incorrectly)

  ;; Extract the true positives (TP), false negatives (FN), true negatives (TN), and false positives (FP) from the confusion matrix
  (setf tp (second (second confusion))
        fn (second (first confusion))
        tn (first (first confusion))
        fp (first (second confusion))
	;; Calculate the sensitivity and specificity
	sensitivity (if (= (+ tp fn) 0) 0 (/ tp (+ tp fn)))
	specificity (if (= (+ tn fp) 0) 0 (/ tn (+ tn fp)))
	;; Calculate the accuracy
	accuracy (/ (+ tp tn) (+ tp tn fp fn)))
  ;; Print the confusion matrix
  (format t "~%Confusion Matrix:~%")
  (format t "  a~Cb~Cclassified as~%"  #\tab #\tab)
  (format t " ~a~C~a~C a = NEGATIVE~%"   (first (first  confusion)) #\tab (second (first confusion )) #\tab)
  (format t " ~a~C~a~C b = POSITIVE~%~%" (first (second confusion)) #\tab (second (second confusion)) #\tab)
  ;; Print the variables involved for debugging
  (format t "True Negatives:~C~a~%" #\tab tn)  
  (format t "False Negatives:~C~a~%" #\tab fn)
  (format t "True Positives:~C~a~%" #\tab tp)
  (format t "False Positives:~C~a~%~%" #\tab fp)
  ;; Print the sensitivity, specificity, and accuracy
  (format t "Sensitivity: ~a~%" (float sensitivity))
  (format t "Specificity: ~a~%" (float specificity))
  (format t "Accuracy: ~a~%" (float accuracy))
  )
)

(defun train-predict (train-file test-file negative-limit)
  (let ( (DATAMATRIX) (FOREST) (N-CLASS) (TARGET) (train-data) )
  ;; Read the training ARFF file
  (setf train-data (read-arff-file train-file ))
  ;; Separate the features (attributes) and the labels
  (setf target (make-array (length
				      (getf train-data :states_num))
				     :element-type 'fixnum
				     :initial-contents (getf train-data :states_num))
        datamatrix (make-array (list
					  (length (getf train-data :states_num));lines
					  (length (getf train-data :attributes));columns
					  ) 
					 :element-type 'single-float
					 :initial-contents (getf train-data :data) )
        n-class (length (getf train-data :classes))
        forest (CL-RANDOM-FOREST::make-forest
			  n-class
			  datamatrix
			  target
			  :n-tree 500
			  :bagging-ratio 1.0
			  :max-depth (truncate (/ (+ 1 (length (getf train-data :attributes))) 2))
			  :min-region-samples 1
			  :n-trial 10))
  (predict forest test-file negative-limit (getf train-data :classes))
  forest
  )
)

;; Serialize the `forest` object
(defun serialize-forest (forest filename)
  (with-open-file (stream filename :direction :output 
  :if-exists :supersede
  :element-type '(unsigned-byte 8))
  (cl-store:store forest stream)))

;; Deserialize the `forest` object
(defun deserialize-forest (filename)
  (with-open-file (stream filename :direction :input 
  :element-type '(unsigned-byte 8))
  (cl-store:restore stream)))	

(in-package :cl-user)

(defun main()
  (let* ((args *posix-argv*)
	 (arg-count (length args))
	 (test-file)
	 (class-border)
	 (train-file)
	 (forest))
    (if (< arg-count 3)
	(progn
	  (format t "The first argument must be a valid test file, and the second an integer.~%")
	  (format t "The integer means a position threshold separating the negatives (head) from the positives (tail) in the test file.~%")
	  (format t "Example1:~C./random-forest ../../validation1.arff 92 ../myids-filter5-89-93-90-a.arff~%" #\tab )
	  (format t "Example2:~C./random-forest ../../validation2.arff 34 ~%~CHere, the 'model.dat' will be loaded.~%" #\tab #\tab )
	  (quit)
	  ))
    (setf test-file (probe-file (nth 1 args))
	  class-border (parse-integer (nth 2 args) :junk-allowed t)
	  train-file (probe-file (if (< arg-count 4) "model.dat" (nth 3 args))))
    (if (and test-file train-file)
	(progn
	  (if (or (not (eq (random-forest:file-or-directory test-file) :REGULAR-FILE))
		  (not (eq (random-forest:file-or-directory train-file) :REGULAR-FILE)))
	      (progn
		(format t "Error: Invalid parameters: check if the path is a valid file.~%")
		(quit)
		))
	  (if (null class-border)
	      (progn
		(format t "Error: Invalid parameters: check class-border parameter.~%") 
		(quit)
		);progn null parameters
	    )
	  (if (< arg-count 4); no training data
	      (progn
		;; Read the model file
		(setf forest (random-forest:deserialize-forest "model.dat"))
		;; Using class definitions within the arff trainig file
		(random-forest:predict forest test-file class-border (list 'SECRETED 'NON-SECRETED) )
	      );prong no training data
	      (setf forest (random-forest:train-predict train-file test-file class-border))
	    );if no training
	  );progn
      (format t "Error: Invalid files: check for valid test and model/training files.~%")
      )
    forest
    );let
  );defun
