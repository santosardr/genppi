(defpackage :features
  (:use :cl)
  (:export
   :nth-string
   :parse-float
   :read-fasta
   :read-propensity
   :all-positions
   :statistics
   :convert-to-key
   :dohistogram
   :dohistogramini
   :dohistogramend
   :dohistogrammid))

(in-package :features)

					;Given a text string, a field delimiter and a position this function returns the text at the specific location.
(defun nth-string (pos tabular delim)
  (let ( (inicio 0) (fim) (retorno) (max (length tabular)) (cont 0) )
    ;do while not
    (do ()
	((or (< pos 0) (= cont pos) (> inicio max)))
       (setq fim  (position delim tabular :start inicio))
       (if (not fim) (setq fim max))
       (setq retorno (subseq tabular inicio fim) 
	     inicio (+ fim 1)
	     cont (1+ cont)
       );setq
      );do
    (if (< cont pos) nil retorno)
    );let
  );defun

(defun parse-float (s)
  (let ((readval (handler-case
                  (read-from-string s)
                  (sb-int:simple-reader-error nil)
                  (end-of-file nil))))
    (cond ((realp readval ) (+ readval 0.0))
          (t (error (concatenate 'string "not a float: " s))))))

					;Returns a list where each component has a pair of text identifier for a fasta sequence
					;and correspondent amino acids as a vector.
(defun read-fasta (filename)
  (let ((inputfile (open filename :if-does-not-exist nil))
	(fastaheader nil)
	(filecontent)
	(sequence (make-array 10 :fill-pointer 0 :adjustable t))
       );let init block
    (when inputfile
     (loop for line = (read-line inputfile nil) until (not (zerop (length line)))
	  finally (setq fastaheader line)
     );loop to skyp blanklines
     (when (char= (elt fastaheader 0) #\>) ;check for fasta format
       (loop for line = (read-line inputfile nil) 
	  while line do 
	    (when (> (length line) 0)
	      (if (char= (elt line 0) #\>) 
		  (progn
		    (setq filecontent (append filecontent (list (list (subseq fastaheader 1) sequence ))))
		    (setq fastaheader line sequence (make-array 10 :fill-pointer 0 :adjustable t))
		    );progn
		  (loop for aa across line do (vector-push-extend (intern (string-upcase  aa)) sequence))
		);if found new fasta header
	    );when there is content
       );loop to read aa sequence
       (setq filecontent (append filecontent (list (list (subseq fastaheader 1) sequence))))
       );when fasta format
     (close inputfile)
     );when inputfile
    filecontent ;return value
  );let
  );defun

(defun read-propensity (filename)
  (let ((inputfile (open filename :if-does-not-exist nil))
	(filecontent nil)
	(filelabels nil)
	);let init block
    (when inputfile
      (loop for line = (read-line inputfile nil) 
	 while line do 
	   (when (> (length line) 0)
	     (if (char/= (elt line 0) #\#)
		 (setf
		  filelabels    (append filelabels
					(list (nth-string 1 line #\Space))
					);append
		  filecontent (append filecontent
				      (list (loop for i from 2 to 21 collect
						 (parse-float (nth-string i line #\Space))
						 );loop
					    );list
				      );append
		  );setf
		 );if
	     );when there is content
	   );loop to read number sequence
      (close inputfile)
      );when inputfile
    (values filelabels filecontent) ;return value
    );let
  );defun

					; Collect a list of needle's positions in a haystack
;(all-positions 'G *sequence*)
(defun all-positions (needle haystack)
  (loop
     for element in haystack 
     and position from 0
     when (eql element needle)
     collect position))

					;Determines the size, mean and standard deviation from a numerical sample,
					;in general, a list of positions of a pattern
(defun statistics( data )
  (if data
    (let ( (n 0) (mean 0.0) (var 0.0) (diff) )
	 (loop for x across data do
	      (incf n)
	      (setq diff (- x mean)
		    var  (+ var (/ (* diff diff (- n 1)) n))
		    mean (+ mean (/ diff n))
		    );setq
	      );loop
	 (if (> n 1)
	     (values n mean (sqrt (/ var (- n 1))) )
     	     (values n mean (sqrt (/ var    n   )) )
	     )
	 );let
    );if
  );defun

					;Considering that each search pattern is a list of lisp symbols, I need to convert the sequence to a string.
(defun convert-to-key ( nucseq )
  (let ( (var "") )
    (loop for i in nucseq do
	 (setf var (concatenate 'string var (symbol-name i )))
	 )
    var
    )  
  )


(defun histo-fasta( filename )
  (let  ( (multifastalist (features:read-fasta filename ))
	 (aatargets (list  'A 'R 'N 'D 'C 'Q 'E 'G 'H 'I 'L 'K 'M 'F 'P 'S 'T 'W 'Y 'V))
	  (histolist nil)
	  (fastaname)
	  (histogram nil)
	) 
    (dolist (fastalist multifastalist histogram)
      (setq fastaname (car fastalist)
	    histolist nil
      );setq
      (dolist (aa aatargets histolist) 
	(setq histolist (append histolist (list (count aa (cadr fastalist))))) 
      );inner do 
      (setq histogram (append histogram (list (list fastaname histolist))))
    );outer do
    histogram
  );let
  )

(defun dohistogram( fastalist aatargets)
  (let  ( (histolist nil))
      (dolist (aa aatargets histolist)
	(setq histolist (append histolist (list (count aa  fastalist)))) 
      );do 
    histolist
  );let
  )

(defun dohistogramini( fastalist aatargets)
  (let  ( (histolist nil)
	  (size (length fastalist))
	  (limit 40)
	) 
      (dolist (aa aatargets histolist) 
	(setq histolist (append histolist (list (count aa  fastalist :end (if (> size limit) limit size) )))) 
      );do 
    histolist
  );let
  )

(defun dohistogramend( fastalist aatargets)
  (let  ( (histolist nil)
	  (size (- (length fastalist) 60))
	) 
      (dolist (aa aatargets histolist) 
	(setq histolist (append histolist (list (count aa  fastalist :start (if (> size 0) size 0) )))) 
      );do 
    histolist
  );let
  )

(defun dohistogrammid( fastalist aatargets)
  (let  ( (histolist nil)
	  (size (length fastalist) )
	  (ratio)
	  (start)
	  (end)
	  )
    (setf ratio (round (/ size 3))
     start (if (> size ratio) ratio size)
	    end (- size ratio)
	    end (if (> end 0) end size))
      (dolist (aa aatargets histolist) 
	(setq histolist (append histolist (list (count aa  fastalist :start start :end end )))) 
      );do 
    histolist
  );let
  )

					;The main function settles everything for the execution of the program.
					;This function also is the entry point for the executable created from this lisp code.

(in-package :cl-user)

(defvar *printlist* nil "Vector to store the final result of the program forcluster in a format to produce the final output of the program")

(defun main ()
  (if (> (length sb-ext:*posix-argv*) 1)
      (let (
					;primeiro parâmetro é o nome do arquivo
	    (fastafile nil)
	    (propensityfile nil)
	    ;(searchspace)
	    ;(milestone)
	    ;(checkmilestone)
	    (multifasta)
	    (multifasta-size)
	    (listof-sequence-names)
	    (sequence)
	    (sequence-name)
	    ;(sequence-size)
	    (sequence-number 0)
	    (sequence-histogram)
	    (pos)
	    (histoalphabet '(  A R N D C Q E G H I L K M F P S T W Y V ))
	    (alphabet)
	    (physicochemicals)
	    ); let pars
	
	(setf fastafile (if (nth 1 sb-ext:*posix-argv*) (nth 1 sb-ext:*posix-argv*) ) );setf
	(if (not (probe-file fastafile))
	    (progn
	      (format t "~%~a~%" "ERROR: Amino acid fasta file is missing")
	      (SB-EXT:EXIT)
	      )
	    );if
	(setf propensityfile (if (nth 2 sb-ext:*posix-argv*)
				 (nth 2 sb-ext:*posix-argv*)
				 "propensity.dat")
	      );setf
	(if (not (probe-file propensityfile ))
	    (progn
	      (format  t "~%~a~%" "ERROR: Propensity data file is missing")
	      (SB-EXT:EXIT)
	      )
	    );if
	
	(multiple-value-bind ( labels props) (features:read-propensity propensityfile)
	  (setf
	   labels (append labels
			  (loop for place in (list "INI" "END" "MID") append
			       (loop for label in labels collect (concatenate 'string place label))
			       )
			  )
	   alphabet (append histoalphabet  labels)
	   physicochemicals props
	   multifasta (features:read-fasta fastafile)
	   multifasta-size (length multifasta)
	   *printlist* (make-array (list  (+ multifasta-size 1)  (+ (length alphabet) 1) ) :initial-element nil)
	   );setf
	  );multiple-value-bind
	
					;replicates the labels for the begin, middle and end of each protein sequence
					;começa o processamento de todas as sequencias
	(loop for fasta in multifasta do
	     (setf sequence (cadr fasta)
		   sequence-name (car fasta)
		   listof-sequence-names (append listof-sequence-names (list sequence-name))
		   ;sequence-size (length sequence)
		   sequence-histogram (features:dohistogram sequence histoalphabet)
		   ;searchspace sequence-size
		   ;milestone (round (*  sequence-size 0.01))
		   ;milestone (if (zerop milestone) 1 milestone)
		   ;checkmilestone milestone
		   );setf
	     
					;começa o processamento para uma sequencia
	     (loop for pos from 0 to (- (length histoalphabet) 1) do
					;medidor de progresso
					;		  (when (= (mod pos checkmilestone) 0)
					;		    (setq checkmilestone (+ checkmilestone milestone))
					;(print (round (* 100 (/ checkmilestone (- searchspace 1)))))
					;		    );when checkmilestone
		  (setf (aref *printlist* (+ sequence-number 1) (+ pos 1)) (nth pos sequence-histogram))
		  );loop
	     
					;after the  histogram calculation we start to compute the physicochemicals properties
	     (setf pos (length histoalphabet))
	     (loop for physico in physicochemicals do
		  (setf (aref *printlist* (+ sequence-number 1) (+ pos 1))
		    (loop for  weigth-list in (list physico) sum
			       (loop for product in (mapcar #'* sequence-histogram weigth-list) sum product)
			       )
			  )
		  (incf pos)
		  )
	     
					;Giving a try: lets compute the number of aminoacids just at the region known to host the signal peptide
	     (setf sequence-histogram (features:dohistogramini sequence histoalphabet))
	     (loop for physico in physicochemicals do
		  (setf (aref *printlist* (+ sequence-number 1) (+ pos 1))
		    (loop for  weigth-list in (list physico) sum
			       (loop for product in (mapcar #'* sequence-histogram weigth-list) sum product)
			       )
			  )
		  (incf pos)
		  )
	     
					;Giving another  try: lets compute the number of aminoacids just at the end of the protein
	     (setf sequence-histogram (features:dohistogramend sequence histoalphabet))
	     (loop for physico in physicochemicals do
		  (setf (aref *printlist* (+ sequence-number 1) (+ pos 1))
		    (loop for  weigth-list in (list physico) sum
			       (loop for product in (mapcar #'* sequence-histogram weigth-list) sum product)
			       )
			  )
		  (incf pos)
		  )
	     
					;Giving another  try: lets compute the number of aminoacids just at the middle
	     (setf sequence-histogram (features:dohistogrammid sequence histoalphabet))
	     (loop for physico in physicochemicals do
		  (setf (aref *printlist* (+ sequence-number 1) (+ pos 1))
		    (loop for  weigth-list in (list physico) sum
			       (loop for product in (mapcar #'* sequence-histogram weigth-list) sum product)
			       )
			  )
		  (incf pos)
		  )
     
	     (incf sequence-number)
	     );loop for multifasta
	
					;final pretty print
					;write the name of the sequences in the first line of the *printlist*
					;In the first line, the first column of features:*printlist* must be an empty string. Starting the names from the second column or x=1
	(setf (aref *printlist* 0 0 ) "");empty string
					;Inserting the attribute names in first line
	(loop for y from 1 to (length histoalphabet) do
	     (setf (aref *printlist* 0 y )
		   (symbol-name (nth (- y 1) histoalphabet))))
	(loop for y from (length histoalphabet) to (length alphabet) do
	     (setf (aref *printlist* 0 y )
		   (nth (- y 1) alphabet)))
	
					;Inserting the name of the sequences in the first column of all lines
	(loop for y from 1 to multifasta-size do (setf (aref *printlist* y 0 ) (nth (- y 1) listof-sequence-names)))
					;Once we wrote the first row, now we need to write the data concerning each sequence name
	(loop for x from 0 to multifasta-size do
	
	     (loop for y from 0 to (length alphabet) do
		  (if (aref *printlist*   x y)
		      (format t "~a~a" (aref *printlist* x y) #\Tab )
		      (format t "~a~a" 0 #\Tab )
		      ))
	          (format t "~%" )
	     )

	);let
      );if
  );defun
