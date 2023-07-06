(ql:quickload "lparallel")

(define-alien-routine "sysconf" long (name int))
(defconstant +_SC_NPROCESSORS_ONLN+ 84)
(defun workers() (sysconf +_SC_NPROCESSORS_ONLN+))

(declaim (optimize (speed 3) (debug 0) (safety 0) (space 0)))
;;;;Implementation of  GenPPI software for PPI  prediction s
;;;;Master's Project in Computer Science - UFU
;;;;Student: William Ferreira dos Anjos
;;;;----------------------------------------------------------------------------

;;Global variables  used in functions
(defparameter *genomas* (make-hash-table :test #'equalp)
  "Tabela hash para armazenar os genomas de uma análise")
(defparameter *pan-genoma* (make-hash-table :test #'equalp)
  "Tabela hash para armazenar o pan-genoma dos genomas")
(defparameter *perfil-filogenetico* (make-hash-table :test #'equalp)
  "Tabela hash para armazenar o perfil filogenético das
   proteínas que estão em mais de um genoma")
(defparameter *pesos-grupos* (make-hash-table :test #'equalp)
  "Tabela hash para fazer o relatório de perfis filogenéticos")
(defparameter *relatorio-vizinhanca-genica* (make-hash-table :test #'equalp)
  "Tabela hash para armazenar o registro do número de genes conservados
   e não conservados para cada expanção de um gene pivô do pan-genoma")
(defparameter *genomes-files* (list)
  "Lista para armazenar os nomes de arquivos de genomas")
(defparameter *fusoes* (make-hash-table :test #'equalp)
  "Tabela hash para armazenar a ppi de cada genoma")
(defparameter *ppi-identified-cn* nil
  "Variável global usada na função genppi para dizer se foi
   constatado ppi por vinsinhação gênica conservada ou não")
(defparameter *ppi-identified-pp* nil
  "Variável global usada na função genppi para dizer se foi
  constatado ppi por perfil filogenético conservada ou não")
(defparameter *ppi-identified-gf* nil
  "Variável global usada na função genppi para dizer se foi
   constatado ppi por fusão gênica conservada ou não")
(defparameter *genomas-pangenoma* (list)
  "Lista usada para armazenas os genomas do pan-genoma.
  Usada nas funões de vizinhança gênica conservada"
  )
(defstruct proteina localidade similares)
(defstruct expansao localidade conservacoes)
(defstruct fusion ppi rosetta-stone)
(defstruct ppi-struct (genea "" :type string) (geneb "" :type string  ) (weight 0.0 :type short-float ) (position 0 :type integer))
;;------------------------------------------------------------------------------

;Read file and put it in a list. In this case, a list with sub-lists,
;where, each sub-list, will contain two elements, the name of the protein
;and its  amino acid sequence respectively.
;(read-fasta "/home/william/PPIs-Predictor/histograma/sequences.fasta")---------
(defun read-fasta (filename)
  ;Read file. If the same  does  not exist,  returns  nil.
  (let ((in (open filename :if-does-not-exist nil))
       	(fastalist nil)
       	(fastaname nil)
       	(proteins nil)
       	(nonaalist (list '+ '- '* '\. '\# '\ ))
        );let init block
    (when in;When in  is  not  nil,that is, when the file was read correctly.
      ;Loops
      ;This loop, reads the first line of the  file  (headerof the first protein), and assigns to the variable fastaname.
      (loop for line = (read-line in nil) until (not (zerop (length line)))
     	  finally (setq fastaname line)
        );loop to skip blank lines
      ;(return)
      ;This when only runs when the read file is a fasta, that is, it has '> at the beginning.
      (when (char= (elt fastaname 0) #\>) ;checks the fasta format.
        (loop for line = (read-line in nil)
       	  while line do
     	    (when (> (length line) 0);when there is content in the read line.
     	      (if (char= (elt line 0) #\>);if found new protein header.
          		  (progn
         		    ;loop to remove characters that are not amino acids (header name).
         		    (dolist (nonaa nonaalist fastalist)
               		      (setq fastalist (delete nonaa fastalist))
                       );dolist

               ;Assigns the  variable  proteins,a list ofproteins, where each protein is a list of two elements,
               ;the header name of the protein and its  amino  acid sequence.
         		    (setq proteins (append proteins (list (list (subseq fastaname 1) fastalist))))

         		    (setq fastaname line fastalist nil)

           		  );progn. Here ends the  actions for a new protein header found.

          		  (setq fastalist (append fastalist
                                      (loop for aa across line collect (intern (string-upcase (string aa)))
                                        );loop
                                     	);append
                    );setq
              );if found new header
         	  );when there is content
          );loop to read a sequence

        (when fastalist;when  fastalist  is  not  nil
         	;loop to remove non-amino acid characters (header name)
        	 (dolist (nonaa nonaalist fastalist)
              	   (setq fastalist (delete nonaa fastalist))
                  );dolist

          ;Assigns the  variable  proteins,a list ofproteins, where each protein is a list of two elements,
          ;the header name of the protein and its  amino  acid sequence.
        	 (setq proteins (append proteins (list (list (subseq fastaname 1) fastalist))))

          );when fastalist
        );when fasta format
      (close in)
      );when in
    proteins ;returns a list of lists, where each element  of the root list, is another list
    ;containing the name and amino  acid  sequence of the protein.
    );let
  )
;-------------------------------------------------------------------------------

;function  that generates the histogram of proteins, i.e. the frequency distribution
;in which each amino acid appears  in  proteins.
;(histo-fasta "/home/william/PPIs-Predictor/histograma/sequences.fasta");-------
(defun histo-fasta( filename )
  (let  ((multifastalist (read-fasta filename ))
       	 (aatargets (list 'A 'R 'N 'D 'C 'E 'Q 'G 'H 'I 'L 'K 'M 'F 'P 'S 'T 'W 'Y 'V 'U 'O 'X 'B 'Z 'J));26 amino acids
       	 (histolist nil)
       	 (fastaname)
       	 (histogram nil))
    ;This dolist scans the file protein list read by the read-fasta function.
    (dolist (fastalist multifastalist histogram)
            ;the  setq  below assigns  the  variable fastaname,the header name of protein n.
            (setq fastaname (concatenate 'string (string (gensym))"."(car fastalist))
             	    histolist nil
                  );setq
            ;This  dolist  scans the list of amino acids definidain  aatargets, and counts the amount oftimes
            ;that each amino acid contained  in  aatargets, appears in proteinn.
            (dolist (aa aatargets histolist)
                   	(setq histolist (append histolist (list (count aa (cadr fastalist)))))
                    );inner do                                        ;cadr pega o segundo elemento da lista, neste caso,
            ;amino  acid  sequence of protein n.

            ;Here, with  func  append,the    blending of the analyzed proteins is  being made,
            ;by creating a single  list  with  sub-lists,whereeach  sub-list is a representation
            ;of the protein histogram, containing the name of the protein and the  frequency  of times each amino acid
            ;analyzed  appears in it.
            (setq histogram (append histogram (list (list fastaname histolist))))
            );outer do
    histogram; retorna o histograma das proteinas.
    );let
  )
;-------------------------------------------------------------------------------

;Functions to break a list in a given position
(defun split-up-to (lista pos)

  (declare (type cons lista))
  (declare (type fixnum pos))

  (cond
    ((<= pos 0) ())
    ((null lista) lista)
    (t (cons (car lista) (split-up-to (cdr lista) (1- pos))))
    )
  )

(defun split-after (lista pos)

  (declare (type cons lista))
  (declare (type fixnum pos))

  (cond
    ((<= pos 0) lista)
    ((null lista) ())
    (t (split-after (cdr lista) (1- pos)))
    )
  )

(defun split (lista pos)

  (declare (type cons lista))
  (declare (type fixnum pos))

  (list (split-up-to lista pos) (split-after lista pos))
  )
;-------------------------------------------------------------------------------

;;Function to test whether two proteins are similar via histogram
(defun similar-test( seqA seqB aadifflimit checkpointminlimit )

  (declare (type cons seqA))
  (declare (type cons seqB))
  (declare (type fixnum aadifflimit))
  (declare (type fixnum checkpointminlimit))

  (let ((checkpoint 0))

    (loop
      for nA in (cadr seqA)
      for nB in (cadr seqB)

      ;;When the difference of histogram between nA (amino n of seqA) and nB (amino n of seqB)
      ;;is <= to the tolerance edift limit and defined as  function parameter(aadifflimit),  
      ;;means that the two  proteins  are similar  to    amino n.
      ;;increases the checkpoint, which accounts for the amount of amino acids that are in.  
      ;;of the difference tolerated. Finally, in if, it checks whether the amount of amino acids
      ;;considered similar, it is >= to checkpointminlimit (amount of amino acids considered in the analysis).
      do (when (<= (abs (- nA nB)) aadifflimit ) (incf checkpoint) )
      );loop
    (if (>= checkpoint checkpointminlimit) t nil)
    );let
  )
;-------------------------------------------------------------------------------

;;Block of functions to read files from a directory
(defun component-present-p (value)
  (and value (not (eql value :unspecific))))

(defun directory-pathname-p (p)
  (and
   (not (component-present-p (pathname-name p)))
   (not (component-present-p (pathname-type p)))
   p))

(defun pathname-as-directory (name)
  (let ((pathname (pathname name)))
    (when (wild-pathname-p pathname)
      (error "Can't reliably convert wild pathnames."))
    (if (not (directory-pathname-p name))
      (make-pathname
       :directory (append (or (pathname-directory pathname) (list :relative))
                          (list (file-namestring pathname)))
       :name
       nil
       :type
       nil
       :defaults pathname)
      pathname)))

(defun directory-wildcard (dirname)
  (make-pathname
   :name :wild
   :type #-clisp :wild #+clisp nil
   :defaults (pathname-as-directory dirname)))

(defun list-directory (dirname)
  (when (wild-pathname-p dirname)
    (error "Can only list concrete directory names."))
  (let ((wildcard (directory-wildcard dirname)))
    #+(or sbcl cmu lispworks)
    (directory wildcard)
    #+openmcl
    (directory wildcard :directories t)
    #+allegro
    (directory wildcard :directories-are-files nil)
    #+clisp
    (nconc
     (directory wildcard)
     (directory (clisp-subdirectories-wildcard wildcard)))))
;;------------------------------------------------------------------------------

;;Block of functions to write and read in the data bunch
(defun save-db (filename ppi)
  (with-open-file (out filename
                       :direction :output)
    (with-standard-io-syntax
      (print ppi out))))

(defun load-db (filename)
  (with-open-file (in filename)
    (with-standard-io-syntax
      (return-from load-db (read in)))))
;-------------------------------------------------------------------------------

;;Function to load  protein files and generate your histograms
(defun histo-genomas (files);files: Lista de arquivos de proteínas
  (defparameter *ppi-identified-cn* nil)
  (defparameter *ppi-identified-pp* nil)
  (defparameter *ppi-identified-gf* nil)
  (defparameter *genomas* (make-hash-table :test #'equalp))
  (defparameter *genomes-files* (list))
  (let ((laco 0)(genomenumber 0))

    (format t "Generating amino acids histogram;~%")
    (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
    (format *query-io* "[")
    (force-output *query-io*)

    (dotimes (i (length files))
             (unless (or (eql (pathname-type (pathname (elt files i))) nil)
                         (not
                          (or (string= (pathname-type (pathname (elt files i))) "fa")
                              (string= (pathname-type (pathname (elt files i))) "faa")
                              (string= (pathname-type (pathname (elt files i))) "fasta")
                              );or
                          );not
                         );or
               ;Unless the path i is a directory, or is not a fasta file, do:
               (setf *genomes-files* (append *genomes-files* (list (pathname-name (pathname (elt files i))))))
               ;;Generates the histogram for each protein file i.
               (setf (gethash (pathname-name (pathname (elt files i))) *genomas*) (histo-fasta (elt files i)))
               );unless

             (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (length files)) 100)) laco)))
                      (format *query-io* "=")
                      (force-output *query-io*)
                      )
             (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (length files)) 100)) laco))))
             (incf genomenumber)
             );dolist
    (format t "]~%")
    );let
  );defun
;;------------------------------------------------------------------------------

(defun call-similar-test (j p proteins gp aadifflimit checkpointminlimit)
    (if (similar-test p (first gp) aadifflimit checkpointminlimit)
      (progn
        (setf (third gp) (append (third gp) (list (list j (position p proteins :test #'equalp)))))
	(if (eq (type-of (second gp)) 'CONS) (second gp) (list (second gp)))
        )
      nil
    );if similar-test                                                                                                                         
 )

;;Function block to generate pangenome  
(declaim (ftype (function (fixnum fixnum ) ) pan-genoma-1))
(defun pan-genoma-1 (aadifflimit checkpointminlimit)

  (defparameter *pan-genoma* (make-hash-table :test #'equalp))
  (defparameter *perfil-filogenetico* (make-hash-table :test #'equalp))
  (defparameter *genomas-pangenoma* (list))
  (defstruct proteina localidade similares)

  (let ((pangenoma (list (list (list (string (gensym)) (list)))))
        (similares (list))
        (perfil-filo (list))
        (laco 0)
        (genomenumber 0)
	(howmanygenomes)
        )

    (declare (type fixnum aadifflimit))
    (declare (type fixnum checkpointminlimit))
    (declare (type fixnum laco))
    (declare (type fixnum genomenumber))

    (format t "~%Generating pan-genome;~%")
    (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
    (format *query-io* "[")
    (force-output *query-io*)

    (loop for j being the hash-keys in *genomas* using (hash-value proteins)
      do (progn
	   (dolist (p proteins)
	     (setf howmanygenomes (length pangenoma))
	     (if (> howmanygenomes 0)
		 (progn
		 (setf similares (append similares (lparallel:pmap 'list #'call-similar-test
								   (make-list howmanygenomes :initial-element j)
								   (make-list howmanygenomes :initial-element p)
								   (make-list howmanygenomes :initial-element proteins)
								   pangenoma
								   (make-list howmanygenomes :initial-element aadifflimit)
								   (make-list howmanygenomes :initial-element checkpointminlimit)
								   )))
		 (setf similares (remove-if (lambda (it) (not (car it))) similares))
		 ))
	     (setf pangenoma (append pangenoma (list (list p (list j (position p proteins :test #'equalp)) similares))))
	     (setf similares nil)
	     );dolist proteins
          (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *genomas*)) 100)) laco)))
                   (format *query-io* "=")
                   (force-output *query-io*)
                   )
          (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *genomas*)) 100)) laco))))
          (incf genomenumber)
          );progn
      );loop hash table *genomas*

    (setf pangenoma (rest pangenoma))
    (dolist (p pangenoma)
            (unless (eql (third p) nil)
              (setf (gethash (first (first p)) *pan-genoma*)
                    (make-proteina :localidade (second p)
                                   :similares (third p)
                                   );make-protein
                    );setf pangenome
              (setf perfil-filo (remove-duplicates (append (list (first (second p))) (loop for i in (third p) collecting (first i))) :test #'string=))
              (if (> (length perfil-filo) 1)
                (setf (gethash (first (second p)) *perfil-filogenetico*)
                      (append (gethash (first (second p)) *perfil-filogenetico*)
                              (list (list (first (first p)) perfil-filo (second (second p))))
                              );append
                      );setf phylogenetic  profile 
                );if remove-duplicates
              );unlles
            );dolist pangenome
    (loop for k being the hash-keys in *pan-genoma* using (hash-value v)
      do (progn
          (if (not (find (first (proteina-localidade v)) *genomas-pangenoma* :test #'string=))
            (setf *genomas-pangenoma* (append *genomas-pangenoma* (list (first (proteina-localidade v)))))
            )
          )
      )
    (format t "]~%")
    );let
  );defun
;-------------------------------------------------------------------------------

(declaim (ftype (function (fixnum fixnum fixnum fixnum) ) pan-genoma-2))
(defun pan-genoma-2 (aadifflimit checkpointminlimit ppaadifflimit ppcheckpointminlimit)

  (defparameter *pan-genoma* (make-hash-table :test #'equalp))
  (defparameter *perfil-filogenetico* (make-hash-table :test #'equalp))
  (defparameter *genomas-pangenoma* (list))
  (defstruct proteina localidade similares)

  (let ((pangenoma (list (list (list (string (gensym)) (list)))))
        (similares (list))
        (perfil-filo-p (list))
        (laco 0)
        (genomenumber 0)
	(howmanygenomes)
        )

    (declare (type fixnum aadifflimit))
    (declare (type fixnum checkpointminlimit))
    (declare (type fixnum ppaadifflimit))
    (declare (type fixnum ppcheckpointminlimit))
    (declare (type fixnum laco))
    (declare (type fixnum genomenumber))

    (format t "~%Generating pan-genome;~%")
    (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
    (format *query-io* "[")
    (force-output *query-io*)

    (loop for j being the hash-keys in *genomas* using (hash-value proteins)
	  do (progn
	       (dolist (p proteins)
		 (setf howmanygenomes (length pangenoma))
		 (if (> howmanygenomes 0)
		     (progn
		       (setf similares (append similares (lparallel:pmap 'list #'call-similar-test
									 (make-list howmanygenomes :initial-element j)
									 (make-list howmanygenomes :initial-element p)
									 (make-list howmanygenomes :initial-element proteins)
									 pangenoma
									 (make-list howmanygenomes :initial-element aadifflimit)
									 (make-list howmanygenomes :initial-element checkpointminlimit)
									 )))
		       (setf similares (remove-if (lambda (it) (not (car it))) similares))
		       ))
		 (setf pangenoma (append pangenoma (list (list p (list j (position p proteins :test #'equalp)) similares))))
		 (setf similares nil)
		 );dolist proteins                                                       
	       (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *genomas*)) 100)) laco)))
		 (format *query-io* "=")
		 (force-output *query-io*)
		 )
	       (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *genomas*)) 100)) laco))))
	       (incf genomenumber)
	       );dolist proteins
	  );loop hash table *genomeas*

    (setf pangenoma (rest pangenoma))
    (dolist (p pangenoma)
            (unless (eql (third p) nil)
              ;;inserting into the hash  pangenome  table.
              (setf (gethash (first (first p)) *pan-genoma*)
                    (make-proteina :localidade (second p)
                                   :similares (third p)
                                   );make-protein
                    );setf pangenome

              ;;Sweeping p-like genes to do new similarity test with p.
              (dolist (similar (third p))
		(if (and similar (first similar) (second similar))
		    (progn
                      (if (similar-test (elt (gethash (first (second p)) *genomas*) (second (second p)))
                                        (elt (gethash (first similar) *genomas*) (second similar))
                                        ppaadifflimit ppcheckpointminlimit
                                        )
                        (setf perfil-filo-p (append perfil-filo-p (list (first similar))))
                        );if similar-test
		      );progn
		  );if
		);dolist similar

              ;;inserting in the  profile-phylogenetic hash  table
              (setf perfil-filo-p (remove-duplicates (append (list (first (second p))) perfil-filo-p) :test #'string=))
              (if (> (length perfil-filo-p) 1)
                (setf (gethash (first (second p)) *perfil-filogenetico*)
                      (append (gethash (first (second p)) *perfil-filogenetico*)
                              (list (list (first (first p)) perfil-filo-p (second (second p))))
                              );append
                      );setf phylogenetic  profile 
                );if remove-duplicates
              ;;Resetting the  profile-phylo-p variable 
              (setf perfil-filo-p nil)
              );unlles
            );dolist pangenome
    (loop for k being the hash-keys in *pan-genoma* using (hash-value v)
      do (progn
          (if (not (find (first (proteina-localidade v)) *genomas-pangenoma* :test #'string=))
            (setf *genomas-pangenoma* (append *genomas-pangenoma* (list (first (proteina-localidade v)))))
            )
          )
      )
    (format t "]~%")
    );let
  );defun
;-------------------------------------------------------------------------------

;;Function block to predict  PPI's  by  conserved  gene neighborhood
(defun conserved-neighbourhood-fixed (workingdir percentage-cn w1 cw1 w2 cw2 w3 cw3 w4 cw4 aadifflimit checkpointminlimit)

  (defparameter *relatorio-vizinhanca-genica* (make-hash-table :test #'equalp)
    "Tabela hash para armazenar o registro do número de genes conservados
     e não conservados para cada expanção de um gene pivô do pan-genoma")
  (defstruct expansao localidade conservacoes)

  (let((ppi (make-hash-table :test #'equalp));Tabela hash para armazenar as PPIs constatadas
                                             (genomas (list));List for storing genomes in which gene conservation has been contacted
                                             (pivo-um);Variable  to store each  pangenome  pivot  protein
                                             (pesos);Variable to store the forces of interactions between proteins
                                             (qtd-ppi);Variable to account for the number of conserved neighborhood identified
                                             (conservacao);Boolean variable to observe a gene conservation.
                                             (interaction)
                                             (maior-peso)
                                             (indice-genoma 0)
                                             (diretorio-banco-ppi)
                                             (diretorio-banco-expansao)
                                             (ppi-db)
                                             (laco 0)
                                             (genomenumber 0)
                                             )

    (declare (type single-float percentage-cn))
    (declare (type fixnum w1))
    (declare (type fixnum cw1))
    (declare (type fixnum w2))
    (declare (type fixnum cw2))
    (declare (type fixnum w3))
    (declare (type fixnum cw3))
    (declare (type fixnum w4))
    (declare (type fixnum cw4))
    (declare (type fixnum aadifflimit))
    (declare (type fixnum checkpointminlimit))
    (declare (type fixnum laco))
    (declare (type fixnum genomenumber))

    (format t "~%Predicting ppi by conserved gene neighborhood;~%")
    (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
    (format *query-io* "[")
    (force-output *query-io*)

    (setf diretorio-banco-ppi (concatenate 'string workingdir "/database/ppi-conserved-neighborhood/"))
    (ensure-directories-exist diretorio-banco-ppi)
    (unless (eql (list-directory diretorio-banco-ppi) nil)
      (dolist (file (list-directory diretorio-banco-ppi))
              (unless (eql (pathname-type (pathname file)) nil)
                (delete-file file)
                )
              )
      )

    (setf diretorio-banco-expansao (concatenate 'string workingdir "/database/gene-neighborhood-report/"))
    (ensure-directories-exist diretorio-banco-expansao)
    (unless (eql (list-directory diretorio-banco-expansao) nil)
      (dolist (file (list-directory diretorio-banco-expansao))
              (unless (eql (pathname-type (pathname file)) nil)
                (delete-file file)
                )
              )
      )

    ;;O loop a seguir varre a tabela hash *pan-genoma* e para cada proteína verifica-se
    ;;if there is protein similarity between your right neighborhood and the right neighborhood
    ;;of the proteins similar to it.
    (loop for k being the hash-keys in *pan-genoma* using (hash-value v)
      do (progn
          (setf pivo-um (proteina-localidade v))
          (setf pesos (loop for i from 1 to w1 collecting (list 1 1)))
          (setf genomas (list))
          (setf qtd-ppi 0)

          ;;Sweeps the list of similar proteins (pivot-2) to each  pangenome protein (pivot-one)  
          (dolist (pivo-dois (proteina-similares v))
	    (if (and pivo-um pivo-dois (second pivo-um) (second pivo-dois) )
		(progn
                  ;;Expansion  Loop
                  (loop for i from 1 to w1
                    do (progn
                        (setf conservacao nil)
                        (cond
                          ;condition and action 1
                          ((and (< (+ (second pivo-um) i) (length (gethash (first pivo-um) *genomas*)))
                                (< (+ (second pivo-dois) i) (length (gethash (first pivo-dois) *genomas*))))

                           (if (similar-test (elt (gethash (first pivo-um) *genomas*) (+ (second pivo-um) i))
                                             (elt (gethash (first pivo-dois) *genomas*) (+ (second pivo-dois) i))
                                             aadifflimit checkpointminlimit)(setf conservacao t));if similar-test
                           );condition and action 1

                          ;condition and action 2
                          ((and (>= (+ (second pivo-um) i) (length (gethash (first pivo-um) *genomas*)))
                                (< (+ (second pivo-dois) i) (length (gethash (first pivo-dois) *genomas*))))

                           (if (similar-test (elt (gethash (first pivo-um) *genomas*) (- (+ (second pivo-um) i) (length (gethash (first pivo-um) *genomas*))))
                                             (elt (gethash (first pivo-dois) *genomas*) (+ (second pivo-dois) i))
                                             aadifflimit checkpointminlimit)(setf conservacao t));if similar-test
                           );condition and action 2

                          ;condition and action 3
                          ((and (< (+ (second pivo-um) i) (length (gethash (first pivo-um) *genomas*)))
                                (>= (+ (second pivo-dois) i) (length (gethash (first pivo-dois) *genomas*))))

                           (if (similar-test (elt (gethash (first pivo-um) *genomas*) (+ (second pivo-um) i))
                                             (elt (gethash (first pivo-dois) *genomas*) (- (+ (second pivo-dois) i) (length (gethash (first pivo-dois) *genomas*))))
                                             aadifflimit checkpointminlimit)(setf conservacao t));if similar-test
                           );condition and action 3

                          ;condition and action 4
                          ((and (>= (+ (second pivo-um) i) (length (gethash (first pivo-um) *genomas*)))
                                (>= (+ (second pivo-dois) i) (length (gethash (first pivo-dois) *genomas*))))

                           (if (similar-test (elt (gethash (first pivo-um) *genomas*) (- (+ (second pivo-um) i) (length (gethash (first pivo-um) *genomas*))))
                                             (elt (gethash (first pivo-dois) *genomas*) (- (+ (second pivo-dois) i) (length (gethash (first pivo-dois) *genomas*))))
                                             aadifflimit checkpointminlimit)(setf conservacao t));if similar-test
                           );condition and action 4
                          );Cond End

                        ;Block of actions for when a gene conservation is found in the  expansion
                        (if (eql conservacao t)
                          (progn
                           ;;If there is a conservation found, increase the weights of the conservation report.
                           (incf (first (elt pesos (1- i))))

                           (unless (or (equalp (first pivo-um) (first pivo-dois))
                                       (if (find-if #'(lambda (x) (equalp x (list (1- i) (first pivo-dois)))) genomas) t nil)
                                       );or
                             ;;Unless pivot gene 1 is in the same genome as pivot 2 gene, or pivot 2 gene
                             ;;is in another genome,  but  preserved gene  neighborhood  has already been
                             ;;in this genome, make:
                             ;;The weight for each preserved gene neighborhood    between genomes is increased
                             ;;different.
                             (incf (second (elt pesos (1- i))))
                             ;;The following is stored the genomes in which   the preserved genic neighbourhood was contained so that the weight does  not
                             ;;is incremented again if another  preserved genic a neighborhood is contained in that same genome.
                             (push (list (1- i) (first pivo-dois)) genomas)
                             );unless
                           );progn
                          );if gene conservation finding
                        );prong  expansion
                    );End expansion  loop
		  )
	      )
	    );End expansion   dolist  

          (setf (gethash k *relatorio-vizinhanca-genica*)
                (make-expansao :localidade pivo-um
                               :conservacoes pesos
                               );make-protein
                );setf  report-neighborhood-genic 

          ;Checking how many stored genes exist within the expansion  window
          (cond
            ;If there is at least cw1 (required number of genes stored) in the size window w1
            ((>= (count-if #'(lambda (x) (> (first x) 1)) pesos) cw1) (setf qtd-ppi w1))
            ;If there is at least cw2 (required number of genes conserved) in the size window w2
            ((>= (count-if #'(lambda (x) (> (first x) 1)) (split-up-to pesos w2)) cw2) (setf qtd-ppi w2))
            ;If there is at least cw3 (required number of genes conserved) in the size window w3
            ((>= (count-if #'(lambda (x) (> (first x) 1)) (split-up-to pesos w3)) cw3) (setf qtd-ppi w3))
            ;If there is at least cw4 (required number of genes conserved) in the size window w4
            ((>= (count-if #'(lambda (x) (> (first x) 1)) (split-up-to pesos w4)) cw4) (setf qtd-ppi w4))
            );cond

          ;;Unless you  have no gene stored in the window of 10, make:------------------------------------------------------------
          (unless (= qtd-ppi 0)
            (setf *ppi-identified-cn* t)
            ;;Creating edges between the pivot gene and the other genes in the expansion  window:---------------------------------------------------
            (dotimes (i qtd-ppi)
                     (setf (gethash (first pivo-um) ppi)
                           (append (gethash (first pivo-um) ppi)
                                   (list
                                    (make-ppi-struct
                                     :genea k
                                     :geneb (if (< (+ (second pivo-um) (1+ i)) (length (gethash (first pivo-um) *genomas*)))
                                              (first (elt (gethash (first pivo-um) *genomas*) (+ (second pivo-um) (1+ i))))
                                              (first (elt (gethash (first pivo-um) *genomas*) (- (+ (second pivo-um) (1+ i)) (length (gethash (first pivo-um) *genomas*))))) )
                                     :weight (if (> (first (elt pesos i)) 1) (* 1.0 (first (elt pesos i))) 1.0)
                                     :position (second pivo-um)
                                     )
                                    );list
                                   );append
                           );setf
                     );dotimes

            ;;Creating edges between genes not conserved within expansion  window:----------------------------------------------------
            (dotimes (i (1- qtd-ppi))
                     (unless (> (first (elt pesos i)) 1)
                       (loop for j from (1+ i) to (1- qtd-ppi)
                         do (progn
                             (setf interaction (make-ppi-struct
                                                :genea (if (< (+ (second pivo-um) (+ i 1)) (length (gethash (first pivo-um) *genomas*)))
                                                         (first (elt (gethash (first pivo-um) *genomas*) (+ (second pivo-um) (+ i 1))))
                                                         (first (elt (gethash (first pivo-um) *genomas*) (- (+ (second pivo-um) (+ i 1)) (length (gethash (first pivo-um) *genomas*))))) )
                                                :geneb (if (< (+ (second pivo-um) (+ j 1)) (length (gethash (first pivo-um) *genomas*)))
                                                         (first (elt (gethash (first pivo-um) *genomas*) (+ (second pivo-um) (+ j 1))))
                                                         (first (elt (gethash (first pivo-um) *genomas*) (- (+ (second pivo-um) (+ j 1)) (length (gethash (first pivo-um) *genomas*))))) )
                                                :weight 1.0
                                                :position (if (< (+ (second pivo-um) (+ i 1)) (length (gethash (first pivo-um) *genomas*)))
                                                            (+ (second pivo-um) (+ i 1))
                                                            (- (+ (second pivo-um) (+ i 1)) (length (gethash (first pivo-um) *genomas*)))
                                                            )
                                                )
                                   );setf interaction

                             ;;This  if is to prevent duplication of interactions in  PPI
                             (if (and (eql nil (find-if #'(lambda (x) (and
                                                                       (equalp (ppi-struct-genea interaction) (ppi-struct-genea x))
                                                                       (equalp (ppi-struct-geneb interaction) (ppi-struct-geneb x)) ) )
                                                        (gethash (first pivo-um) ppi) ) )
                                      (eql nil (gethash (ppi-struct-genea interaction)  *pan-genoma*))
                                      );and
                               ;;If the interaction of the time is no longer contained in the PPI,and if the pi protein is notcontained in the  pangenome,make:
                               (setf (gethash (first pivo-um) ppi) (append (gethash (first pivo-um) ppi) (list interaction)))
                               );if and
                             );progn
                         );loop
                       );unless
                     );dotimes
            );unless

          ;;Normalizing weights, taking  the PPI out of memory and writing to the database.
          (unless (string= (first pivo-um) (elt *genomas-pangenoma* indice-genoma))
            (setf ppi-db (gethash (elt *genomas-pangenoma* indice-genoma) ppi))
            (unless (eql ppi-db nil)
              (setf maior-peso (loop for i in ppi-db
                                 maximizing (ppi-struct-weight i) into max
                                 finally (return max)))
              (dolist (p ppi-db)
                      (setf (ppi-struct-weight p) (* (/ (ppi-struct-weight p) maior-peso) percentage-cn))
                      )
              (save-db (concatenate 'string diretorio-banco-ppi (elt *genomas-pangenoma* indice-genoma) ".db")
                       ppi-db
                       );save-db.
              ;;Taking the data saved in the ram database.
              (setf (gethash (elt *genomas-pangenoma* indice-genoma) ppi) nil)
              );unlles
            (incf indice-genoma)
            );unless

          (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *pan-genoma*)) 100)) laco)))
                   (format *query-io* "=")
                   (force-output *query-io*)
                   )
          (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *pan-genoma*)) 100)) laco))))
          (incf genomenumber)

          );Fim do bloco de acoes para cada iteracao na tabela hash *pan-genoma*
      );Fim do loop na tabela hash *pan-genoma*
    ;---------------------------------------------------------------------------

    ;;Normalizing weights of the last genome, taking  ppi out of memory and writing to the database.
    (setf ppi-db (gethash (first (last *genomas-pangenoma*)) ppi))
    (unless (eql ppi-db nil)
      (setf maior-peso (loop for i in ppi-db
                         maximizing (ppi-struct-weight i) into max
                         finally (return max)))
      (dolist (p ppi-db)
              (setf (ppi-struct-weight p) (* (/ (ppi-struct-weight p) maior-peso) percentage-cn))
              )
      (save-db (concatenate 'string diretorio-banco-ppi (first (last *genomas-pangenoma*)) ".db")
               ppi-db
               );save-db.
      (setf (gethash (first (last *genomas-pangenoma*)) ppi) nil)
      )
    ;---------------------------------------------------------------------------

    ;;Gravando a tabela hash *relatorio-vizinhanca-genica* no banco de dados.
    (save-db (concatenate 'string diretorio-banco-expansao "gene-neighborhood-report.db")
             *relatorio-vizinhanca-genica*
             );save-db.
    ;;Taking the data out of memory.
    (setf *relatorio-vizinhanca-genica* nil)

    (format t "]~%")
    ))
;-------------------------------------------------------------------------------

;;Function two to predict  PPI's  per    conserved gene neighborhood
(defun conserved-neighbourhood-dynamic (workingdir percentage-cn janela aadifflimit checkpointminlimit)

  (defparameter *relatorio-vizinhanca-genica* (make-hash-table :test #'equalp)
    "Tabela hash para armazenar o registro do número de genes conservados
   e não conservados para cada expanção de um gene pivô do pan-genoma")
  (defstruct expansao localidade conservacoes)

  (let ((ppi (make-hash-table :test #'equalp));Lista para armazenar as PPIs constatadas
                                              (pivo-um);Variable  to store each  pangenome  pivot  protein
                                              (pos); Variável para guardar a última posição de uma vizinhança gênica conservada
                                              (pesos);Variable to store the forces of interactions between proteins.
                                              (conservacao);Boolean variable to observe a gene conservation.
                                              (genomas (list));List to assist in weight incrementation
                                              (interaction)
                                              (maior-peso)
                                              (indice-genoma 0)
                                              (diretorio-banco-ppi)
                                              (diretorio-banco-expansao)
                                              (ppi-db)
                                              (laco 0)
                                              (genomenumber 0)
                                              (total-expancoes)
                                              )

    (declare (type single-float percentage-cn))
    (declare (type fixnum janela))
    (declare (type fixnum aadifflimit))
    (declare (type fixnum checkpointminlimit))
    (declare (type fixnum laco))
    (declare (type fixnum genomenumber))

    (format t "~%Predicting ppi by conserved gene neighborhood;~%")
    (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
    (format *query-io* "[")
    (force-output *query-io*)

    (setf diretorio-banco-ppi (concatenate 'string workingdir "/database/ppi-conserved-neighborhood/"))
    (ensure-directories-exist diretorio-banco-ppi)
    (unless (eql (list-directory diretorio-banco-ppi) nil)
      (dolist (file (list-directory diretorio-banco-ppi))
              (unless (eql (pathname-type (pathname file)) nil)
                (delete-file file)
                )
              )
      )

    (setf diretorio-banco-expansao (concatenate 'string workingdir "/database/gene-neighborhood-report/"))
    (ensure-directories-exist diretorio-banco-expansao)
    (unless (eql (list-directory diretorio-banco-expansao) nil)
      (dolist (file (list-directory diretorio-banco-expansao))
              (unless (eql (pathname-type (pathname file)) nil)
                (delete-file file)
                )
              )
      )

    ;;O loop a seguir varre a tabela hash *pan-genoma* e para cada proteína verifica-se
    ;;if there is protein similarity between your right neighborhood and the right neighborhood
    ;;of proteins similar to it.
    (loop for k being the hash-keys in *pan-genoma* using (hash-value v)
      do (progn
          (setf pivo-um (proteina-localidade v))
          (setf pesos (loop for i from 1 to janela collecting (list 1 1)))
          (setf genomas (list))

          (dolist (pivo-dois (proteina-similares v))
	    (if (and pivo-um pivo-dois (second pivo-um) (second pivo-dois) )
                (progn
		  (block expancao
                         (setf total-expancoes 0)
                         (setf pos 0)
                         (loop
                           (setf conservacao nil)
                           (loop for i from 1 to janela
                             do (progn
                                 (incf total-expancoes)
                                 (if (or (= total-expancoes (1- (length (gethash (first pivo-um) *genomas*))))
                                         (= total-expancoes (1- (length (gethash (first pivo-dois) *genomas*)))))
                                   (return-from expancao))
                                 (cond
                                   ;condition and action 1
                                   ((and (< (+ (second pivo-um) (+ pos i)) (length (gethash (first pivo-um) *genomas*)))
                                         (< (+ (second pivo-dois) (+ pos i)) (length (gethash (first pivo-dois) *genomas*))))

                                    (if (similar-test (elt (gethash (first pivo-um) *genomas*) (+ (second pivo-um) (+ pos i)))
                                                      (elt (gethash (first pivo-dois) *genomas*) (+ (second pivo-dois) (+ pos i)))
                                                      aadifflimit checkpointminlimit)(setf conservacao t))
                                    );condition and action 1

                                   ;condition and action 2
                                   ((and (>= (+ (second pivo-um) (+ pos i)) (length (gethash (first pivo-um) *genomas*)))
                                         (< (+ (second pivo-dois) (+ pos i)) (length (gethash (first pivo-dois) *genomas*))))

                                    (if (similar-test (elt (gethash (first pivo-um) *genomas*) (- (+ (second pivo-um) (+ pos i)) (length (gethash (first pivo-um) *genomas*))))
                                                      (elt (gethash (first pivo-dois) *genomas*) (+ (second pivo-dois) (+ pos i)))
                                                      aadifflimit checkpointminlimit)(setf conservacao t))
                                    );condition and action 2

                                   ;condition and action 3
                                   ((and (< (+ (second pivo-um) (+ pos i)) (length (gethash (first pivo-um) *genomas*)))
                                         (>= (+ (second pivo-dois) (+ pos i)) (length (gethash (first pivo-dois) *genomas*))))

                                    (if (similar-test (elt (gethash (first pivo-um) *genomas*) (+ (second pivo-um) (+ pos i)))
                                                      (elt (gethash (first pivo-dois) *genomas*) (- (+ (second pivo-dois) (+ pos i)) (length (gethash (first pivo-dois) *genomas*))))
                                                      aadifflimit checkpointminlimit)(setf conservacao t))
                                    );condition and action 3

                                   ;condition and action 4
                                   ((and (>= (+ (second pivo-um) (+ pos i)) (length (gethash (first pivo-um) *genomas*)))
                                         (>= (+ (second pivo-dois) (+ pos i)) (length (gethash (first pivo-dois) *genomas*))))

                                    (if (similar-test (elt (gethash (first pivo-um) *genomas*) (- (+ (second pivo-um) (+ pos i)) (length (gethash (first pivo-um) *genomas*))))
                                                      (elt (gethash (first pivo-dois) *genomas*) (- (+ (second pivo-dois) (+ pos i)) (length (gethash (first pivo-dois) *genomas*))))
                                                      aadifflimit checkpointminlimit)(setf conservacao t))
                                    );condition and action 4
                                   );cond

                                 ;Block of shares for when a gene conservation is found in the  expansion
                                 (if (eql conservacao t)
                                   (progn
                                    ;;Code block for reporting purposes
                                    ;;Expanding the size of the relative-de-conservation if necessary 
                                    (unless (>= (length (split-after pesos (+ pos i))) janela)
                                      (setf pesos (append pesos
                                                          (loop for i from 1 to
                                                            (- janela (length (split-after pesos (+ pos i))))
                                                            collecting (list 1 1))
                                                          );append
                                            );setf
                                      );unless
                                    ;;If there is a conservation found, increase the weights of the conservation report.
                                    (incf (first (elt pesos (+ pos (1- i)))))
                                    ;;End of code block for reporting purposes

                                    ;;Unless the similar protein is in the same genome as the protein analyzed, or if the protein
                                    ;;similar  is in another genome,  but  preserved gene  neighborhood  has already been reported
                                    ;;in this genome, make:
                                    (unless (or (equalp (first pivo-um) (first pivo-dois))
                                                (if (find-if #'(lambda (x) (equalp x (list (+ pos (1- i)) (first pivo-dois)))) genomas) t nil)
                                                );or

                                      (incf (second (elt pesos (+ pos (1- i)))))
                                      ;;The following is the genomes in which   the preserved genic neighbourhood was contained so that the weight does  not
                                      ;;is incremented again if another  preserved genic a neighborhood is contained in that same genome.
                                      (push (list (+ pos (1- i)) (first pivo-dois)) genomas)
                                      );unless
                                    (setf pos (+ pos i))
                                    (setf i janela)
                                    );progn
                                   ;;If there is no gene conservation in the  expansion,make:
                                   (if (= i janela); Se já expandiu até o limite da janela, encerra-se a expanção.
                                     (return-from expancao))
                                   );if gene conservation finding
                                 (setf conservacao nil)
                                 );progn
                             );window  loop
                           ); loop expanção
                         ); Fim bloco expanção
		  )
	      )
	  );similar dolists

          (setf (gethash k *relatorio-vizinhanca-genica*)
                (make-expansao :localidade pivo-um
                               :conservacoes pesos
                               );make-protein
                );setf  report-neighborhood-genic 

          ;(terpri)
          ;(format t "~a => ~a~%" k pesos)

          ;Unless you  have  no gene stored in the expansion  window,make:
          (unless (= (length pesos) janela)
            (setf *ppi-identified-cn* t)
            ;;Creating edges between the pivot gene and the other genes in the expansion  window:------------------------------------------------------------------------------------------------------
            (dotimes (i (length pesos))
                     (setf (gethash (first pivo-um) ppi)
                           (append (gethash (first pivo-um) ppi)
                                   (list (make-ppi-struct
                                          :genea k
                                          :geneb (if (< (+ (second pivo-um) (1+ i)) (length (gethash (first pivo-um) *genomas*)))
                                                   (first (elt (gethash (first pivo-um) *genomas*) (+ (second pivo-um) (1+ i))))
                                                   (first (elt (gethash (first pivo-um) *genomas*) (- (+ (second pivo-um) (1+ i)) (length (gethash (first pivo-um) *genomas*))))) )
                                          :weight (if (> (first (elt pesos i)) 1) (* 1.0 (first (elt pesos i))) 1.0)
                                          :position  (second pivo-um)
                                          )
                                         );list
                                   );append
                           );setf
                     );dotimes

            ;;Creating edges between genes not conserved within expansion  window:------------------------------------------------------------------------------------------------------
            (dotimes (i (1- (length pesos)))
                     (unless (> (first (elt pesos i)) 1)
                       (loop for j from (1+ i) to (1- (length pesos))
                         do (progn
                             (setf interaction (make-ppi-struct
                                                :genea (if (< (+ (second pivo-um) (+ i 1)) (length (gethash (first pivo-um) *genomas*)))
                                                         (first (elt (gethash (first pivo-um) *genomas*) (+ (second pivo-um) (+ i 1))))
                                                         (first (elt (gethash (first pivo-um) *genomas*) (- (+ (second pivo-um) (+ i 1)) (length (gethash (first pivo-um) *genomas*))))) )
                                                :geneb (if (< (+ (second pivo-um) (+ j 1)) (length (gethash (first pivo-um) *genomas*)))
                                                         (first (elt (gethash (first pivo-um) *genomas*) (+ (second pivo-um) (+ j 1))))
                                                         (first (elt (gethash (first pivo-um) *genomas*) (- (+ (second pivo-um) (+ j 1)) (length (gethash (first pivo-um) *genomas*))))) )
                                                :weight  1.0
                                                :position  (if (< (+ (second pivo-um) (+ i 1)) (length (gethash (first pivo-um) *genomas*)))
                                                             (+ (second pivo-um) (+ i 1))
                                                             (- (+ (second pivo-um) (+ i 1)) (length (gethash (first pivo-um) *genomas*)))
                                                             )
                                                )
                                   );setf interaction

                             ;;This  if is to prevent duplication of interactions in  PPI
                             (if (and (eql nil (find-if #'(lambda (x) (and
                                                                       (equalp (ppi-struct-genea interaction) (ppi-struct-genea x))
                                                              				 (equalp (ppi-struct-geneb interaction) (ppi-struct-geneb x))
                                                                       );and
                                                                  );lambda
                                                							 (gethash (first pivo-um) ppi)))
                      				            (eql nil (gethash (ppi-struct-genea interaction) *pan-genoma*))
                                      );and
                               ;;If the interaction of the time is no longer contained in the PPI,and if the pi protein is notcontained in the  pangenome,make:
                               (setf (gethash (first pivo-um) ppi) (append (gethash (first pivo-um) ppi) (list interaction)))
                               );
                             );progn
                         );loop
                       );unless
                     );dotimes
            );unless

          ;;Normalizing the weights and taking the  PPI  out of memory and writing to the database.
          (unless (string= (first pivo-um) (elt *genomas-pangenoma* indice-genoma))
            (setf ppi-db (gethash (elt *genomas-pangenoma* indice-genoma) ppi))
            (unless (eql ppi-db nil)
              (setf maior-peso (loop for i in ppi-db
                                 maximizing (ppi-struct-weight i) into max
                                 finally (return max)))
              (dolist (p ppi-db)
                      (setf (ppi-struct-weight p) (* (/ (ppi-struct-weight p) maior-peso) percentage-cn))
                      )
              (save-db (concatenate 'string diretorio-banco-ppi (elt *genomas-pangenoma* indice-genoma) ".db")
                       ppi-db
                       );save-db.
              (setf (gethash (elt *genomas-pangenoma* indice-genoma) ppi) nil)
              )
            (incf indice-genoma)
            );unless

          (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *pan-genoma*)) 100)) laco)))
                   (format *query-io* "=")
                   (force-output *query-io*)
                   )
          (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *pan-genoma*)) 100)) laco))))
          (incf genomenumber)

          );Fim do bloco de acoes para cada iteracao na tabela hash *pan-genoma*
      );Fim do loop na tabela hash *pan-genoma*---------------------------------

    ;;Normalizing weights of the last genome, taking  ppi out of memory and writing to the database.
    (setf ppi-db (gethash (first (last *genomas-pangenoma*)) ppi))
    (unless (eql ppi-db nil)
      (setf maior-peso (loop for i in ppi-db
                         maximizing (ppi-struct-weight i) into max
                         finally (return max)))
      (dolist (p ppi-db)
              (setf (ppi-struct-weight p) (* (/ (ppi-struct-weight p) maior-peso) percentage-cn))
              )
      (save-db (concatenate 'string diretorio-banco-ppi (first (last *genomas-pangenoma*)) ".db")
               ppi-db
               );save-db.
      (setf (gethash (first (last *genomas-pangenoma*)) ppi) nil)
      )
    ;---------------------------------------------------------------------------

    ;;Gravando a tabela hash *relatorio-vizinhanca-genica* no banco de dados.
    (save-db (concatenate 'string diretorio-banco-expansao "gene-neighborhood-report.db")
             *relatorio-vizinhanca-genica*
             );save-db.
    ;;Taking the data out of memory.
    (setf *relatorio-vizinhanca-genica* nil)

    (format t "]~%")
    ))
;;;;End functions for conserved gene neighborhood

;;Block of functions to predict PPI's  by  Phylogenetic Profile
(defun comparar-perfis (lista1 lista2)
  (unless (eql nil (set-exclusive-or lista1 lista2 :test #'string=))
    (length(set-exclusive-or lista1 lista2 :test #'string=))
    )
  )

(defun phylogenetic-profiles-complete (percentage-pp pptolerance workingdir)

  (defparameter *pesos-grupos* (make-hash-table :test #'equalp)
    "Tabela hash para fazer o relatório de perfis filogenéticos")

  (let ((ppi (list));;List to   store the interactions created.
                    (grupos-identicos (list));List for storing protein groups with identical profiles.
                    (grupos-similares (list));List for storing protein groups with identical and similar profiles.
                    (pesos-grupos (list))
                    (posicao-grupo);Variable to store the position of a group g.
                    ;(number-proteínas-max-group);Variável para armezar o número de proteínas do maior grupo.
                    (soma-grupos 0);Variable to store the total number of proteins between groups.
                    (peso-grupo);Variable to store the weight that will be assigned to proteins in a group g.
                    (dividendo);Variable that stores the index of a t-size group, which is closer to the largest group.
                    (divisor);Variable that stores the index of the largest group.
                    (diferenca-perfis);Variable to store the difference of  phylogenetic profiles of two groups.
                    (diretorio-banco-ppi)
                    (diretorio-banco-grupos)
                    (laco 0)
                    (genomenumber 0)
                    )

    (declare (type single-float percentage-pp))
    (declare (type fixnum laco))
    (declare (type fixnum genomenumber))
    (declare (type fixnum pptolerance))

    (setf diretorio-banco-ppi (concatenate 'string workingdir "/database/ppi-phylogenetic-profiles/"))
    (ensure-directories-exist diretorio-banco-ppi)
    (unless (eql (list-directory diretorio-banco-ppi) nil)
      (dolist (file (list-directory diretorio-banco-ppi))
              (unless (eql (pathname-type (pathname file)) nil)
                (delete-file file)
                );unless
              );dolist
      );unless

    (setf diretorio-banco-grupos (concatenate 'string workingdir "/database/phylogenetic-profile-groups/"))
    (ensure-directories-exist diretorio-banco-grupos)
    (unless (eql (list-directory diretorio-banco-grupos) nil)
      (dolist (file (list-directory diretorio-banco-grupos))
              (unless (eql (pathname-type (pathname file)) nil)
                (delete-file file)
                );unless
              );dolist
      );unless

    (format t "~%Predicting ppi by phylogenetic profiles;~%")
    (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
    (format *query-io* "[")
    (force-output *query-io*)

    (loop for k being the hash-keys in *perfil-filogenetico* using (hash-value perfis)
      do (progn
          ;;Grouping proteins by phylogenetic profile
          (dolist (p perfis)
                  (setf posicao-grupo (position-if #'(lambda (x) (eql nil (set-exclusive-or (second p) (first x) :test #'string=))) grupos-identicos))
                  (if (eql posicao-grupo nil)
                    (push (list (second p) (list (list (first p) (third p)))) grupos-identicos)
                    (push (list (first p) (third p)) (second (elt grupos-identicos posicao-grupo)))
                    );if
                  );dolist
          ;;End of The Grouping

          ;;Ordering the proteins of each group
          (dolist (g grupos-identicos)
                  ;;Ordering the proteins of each group according to their index in the  genome.
                  (setf (second g) (sort (second g) #'< :key #'second))
                  ;;Concatenando to a group g, the size of group g.
                  (setf g (nconc g (list (length (second g)))))
                  )

          (unless (= (length grupos-identicos) 0)
            (cond
              ;;If the user only wants identical profiles
              ((= pptolerance 0)
               (progn

                ;;Ordering groups by size, from smallest to largest.
                (setf grupos-identicos (sort grupos-identicos #'< :key #'third))

                ;;Writing groupings by profiles to the database.
                (save-db (concatenate 'string diretorio-banco-grupos k ".db") grupos-identicos)

                ;;by storing the  weights of each group in a hash table that will be used in the report.
                (setf (gethash k *pesos-grupos*) (remove-duplicates grupos-identicos :key #'third))

                ;;Removing groups with the same amount of genes and assigning the result to the group weights list.
                (setf pesos-grupos (remove-duplicates grupos-identicos :key #'third))

                ;;Generating interactions for each pair of proteins in each group.
                (dolist (g grupos-identicos)

                        ;;Picking up the group index of size equal to group g in the group weights list.
                        (setf dividendo (1+ (position-if #'(lambda (x) (= (third g) (third x))) pesos-grupos)))

                        ;;Picking up the index of the most populous group, that is, the last group.
                        (setf divisor (length pesos-grupos))

                        ;;Calculating the weight of group g based on their position, and the position of the last group (most populous).
                        (setf peso-grupo (* (/ dividendo divisor) percentage-pp))

                        ;;Generating  PPIs for every pair of proteins in group g.
                        (dotimes (i (1- (third g)))
                                 (loop for j from (1+ i) to (1- (third g))
                                   do (progn
                                       (push (make-ppi-struct
                                              :genea (first (elt (second g) i))
                                              :geneb (first (elt (second g) j))
                                              :weight peso-grupo
                                              :position (second (elt (second g) i)))
                                             ppi)
                                       );progn-do
                                   );loop-dotimes
                                 );dotimes
                        (setf g (nconc g (list peso-grupo)))
                        );dolist-groups-identical
                );progn
               );condition and action 1

              ;;If the user wants to consider similar profiles
              ((> pptolerance 0)
               (progn

                ;;Ordering groups by size, from largest to smallest.
                (setf grupos-identicos (sort grupos-identicos #'> :key #'third))

                ;;Adding to a group g, all groups similar to group g.
                (dotimes (i (1- (length grupos-identicos)))
                         (setf grupos-similares  (append grupos-similares (list (list (elt grupos-identicos i) (list )))))
                         (setf soma-grupos (+ soma-grupos (third (elt grupos-identicos i))))
                         (loop for j from (1+ i) to (1- (length grupos-identicos))
                           do (progn
                               (unless (> (comparar-perfis (first (elt grupos-identicos i)) (first (elt grupos-identicos j))) pptolerance)
                                 ;;Unless the difference between the profiles of a pair of groups is greater than the tolerated difference, make:
                                 ;;Adding to the gi group, the gj group, which is similar to the gi group.
                                 (push (elt grupos-identicos j) (second (elt grupos-similares i)))
                                 ;;Summing up the total amount of proteins between similar groups.
                                 (setf soma-grupos (+ soma-grupos (third (elt grupos-identicos j))))
                                 );unlles
                               );progn-do
                           );loop-dotimes
                         (setf (elt grupos-similares i) (nconc (elt grupos-similares i) (list soma-grupos)))
                         (setf soma-grupos 0)
                         );dotimes

                ;;Adding the last group that was not considered in the above dotimes, to the list of similar groups.
                (setf grupos-similares  (append grupos-similares (list (list (first (last grupos-identicos)) (list ) (length (second (first (last grupos-identicos))))))))

                ;;Ordering similar groups by size, from smallest to largest.
                (setf grupos-similares (sort grupos-similares #'< :key #'(lambda (x) (third (first x)))))

                ;;Writing groupings by profiles to the database.
                (save-db (concatenate 'string diretorio-banco-grupos k ".db") grupos-similares)

                ;;by armezenando the weights of each group in a hash table that will be used in the report.
                (setf (gethash k *pesos-grupos*) (remove-duplicates grupos-similares :key #'(lambda (x) (third (first x)))))

                ;;Removing groups with the same amount of genes and assigning the result to the group weights list.
                (setf pesos-grupos (remove-duplicates grupos-similares :key #'(lambda (x) (third (first x)))))

                ;;Sweeping the groups formed to create edges between  identical  and  similar proteins.
                (dolist (g grupos-similares)

                        ;;Picking up the group index of size equal to group g in the group weights list.
                        (setf dividendo (1+ (position-if #'(lambda (x) (= (third (first g)) (third (first x)))) pesos-grupos)))

                        ;;Picking up the index of the most populous group, i.e. the last group.
                        (setf divisor (length pesos-grupos))

                        ;;Calculating the weight of group g based on their position, and the position of the last group (most populous).
                        (setf peso-grupo (* (/ dividendo divisor) percentage-pp))

                        ;;Creating edges between proteins with identical profiles.
                        (dotimes (i (1- (third (first g))))
                                 (loop for j from (1+ i) to (1- (third (first g)))
                                   do (progn
                                       (push (make-ppi-struct
                                              :genea (first (elt (second (first g)) i))
                                              :geneb (first (elt (second (first g)) j))
                                              :weight peso-grupo
                                              :position (second (elt (second (first g)) i)))
                                             ppi)
                                       );progn-do
                                   );loop j
                                 );dotimes i

                        ;;Sweeping protein groups with profiles similar to the profile of the protein group with identical profiles.
                        (dolist (grupo-similar (second g))
                                ;Calculating the  difference  between profiles.
                                (setf diferenca-perfis (comparar-perfis (first (first g)) (first grupo-similar)))
                                ;;Sweeping the proteins of a group with a similar profile.
                                (dolist (ps (second grupo-similar))
                                        ;;Sweeping proteins   of identical profiles of a group g.
                                        (dolist (pg (second (first g)))
                                                ;;Creating  an interaction.
                                                (push (make-ppi-struct
                                                 						:genea (if (< (second pg) (second ps)) (first pg) (first ps) )
                                                 						:geneb (if (> (second pg) (second ps)) (first pg) (first ps) )
                                                 						:weight (- peso-grupo (* peso-grupo (/ (/ (* diferenca-perfis 100) (length *genomes-files*)) 100)))
                                                 						:position (if (< (second pg) (second ps))  (second pg) (second ps) ) )
                                          					       ppi);push
                                                );dolist pg
                                        );dolist ps
                                );dolist group-similar
                        );dolist groups-similar
                );progn
               );condition and action 2
              );cond
            );unless

          ;;If any PPI has  beenpredicted, the PPI identifier is  changed  to  true.
          (if (> (length ppi) 0)
            (setf *ppi-identified-pp* t)
            )

          ;;Recording the  PPI of genome k in the database.
          (save-db (concatenate 'string diretorio-banco-ppi k ".db") ppi)

          ;;Reset the variables.
          (setf grupos-identicos nil)
          (setf grupos-similares nil)
          (setf ppi nil)

          (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *perfil-filogenetico*)) 100)) laco)))
                   (format *query-io* "=")
                   (force-output *query-io*)
                   )
          (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *perfil-filogenetico*)) 100)) laco))))
          (incf genomenumber)

          );progn-loop
      );loop hash table
    (format t "]~%")
    ))
;-------------------------------------------------------------------------------

(defun phylogenetic-profiles-ppiterlimit (percentage-pp pptolerance ppiterlimit workingdir)

  (let ((ppi (list));;List  to store the interactions created.
                    (grupos-identicos (list));List for storing protein groups with identical profiles.
                    (grupos-similares (list));List for storing protein groups with identical and similar profiles.
                    (pesos-grupos (list));List to store the weight that will be assigned to each group.
                    (posicao-grupo);Variable to store the position of a group g.
                    ;(number-proteínas-max-group);Variável para armezar o número de proteínas do maior grupo.
                    ;(max-proteínas-group);Variável para armazenar o grupo com maior número de proteínas.
                    (soma-grupos 0);Variable to store the total number of proteins between groups.
                    (peso-grupo);Variable to store the weight that will be assigned to proteins in a group g.
                    (dividendo);Variable that stores the index of a t-size group, which is closer to the largest group.
                    (divisor);Variable that stores the index of the largest group.
                    (diferenca-perfis);Variable to store the  difference  of phylogenetic profiles of two groups.
                    (total-arestas);Variable to store the sum of the total edges that each group will produce.
                    ;(subtraction);Variável para armazenar a quantidade total de arestas após exclusão do grupo mais populoso.
                    ;(rest);Variável para calcular a quantidade de genes que serão aproveitados do grupo que será excluído.
                    ;(cut-off index);Variável para armazenar o ídice que marca em qual posião o maior grupo será cortato.
                    (diretorio-banco-ppi)
                    (laco 0)
                    (genomenumber 0)
                    )

    (declare (type single-float percentage-pp))
    (declare (type fixnum laco))
    (declare (type fixnum genomenumber))
    (declare (type fixnum pptolerance))
    (declare (type fixnum ppiterlimit))

    (format t "~%Predicting ppi by phylogenetic profiles;~%")
    (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
    (format *query-io* "[")
    (force-output *query-io*)

    (setf diretorio-banco-ppi (concatenate 'string workingdir "/database/ppi-phylogenetic-profiles/"))
    (ensure-directories-exist diretorio-banco-ppi)
    (unless (eql (list-directory diretorio-banco-ppi) nil)
      (dolist (file (list-directory diretorio-banco-ppi))
              (unless (eql (pathname-type (pathname file)) nil)
                (delete-file file)
                );unless
              );dolist
      );unless

    (loop for k being the hash-keys in *perfil-filogenetico* using (hash-value perfis)
      do (progn
          ;;Grouping proteins by phylogenetic profile
          (dolist (p perfis)
                  (setf posicao-grupo (position-if #'(lambda (x) (eql nil (set-exclusive-or (second p) (first x) :test #'string=))) grupos-identicos))
                  (if (eql posicao-grupo nil)
                    (push (list (second p) (list (list (first p) (third p)))) grupos-identicos)
                    (push (list (first p) (third p)) (second (elt grupos-identicos posicao-grupo)))
                    );if
                  );dolist
          ;;End of Grouping

          ;;Concatenando the size of each group g
          (dolist (g grupos-identicos)
                  ;;Concatenando to a group g, the size of group g.
                  (setf g (nconc g (list (length (second g)))))
                  )

          (unless (= (length grupos-identicos) 0)
            (cond
              ;;If the user only wants identical profiles
              ((= pptolerance 0)
               (progn

                ;;Ordering groups by size, from smallest to largest.
                (setf grupos-identicos (sort grupos-identicos #'< :key #'third))

                ;;Block to exclude groups until the total edges produced is not higher than the threshold received
                (block remove-groups
                       (loop
                         (setf total-arestas (loop for i in grupos-identicos
                                               summing (/ (* (length (second i)) (- (length (second i)) 1)) 2) into total
                                               finally (return total)))
                         (if (<= total-arestas ppiterlimit)
                           (return-from remove-groups)
                           (setf grupos-identicos (remove (first (last grupos-identicos)) grupos-identicos :test #'equalp))
                           #|(progn
                              (setf number-proteins-max-group (third (first (last grupos-identicos))))
                              (setf subtracao (- total-arestas (/ (* number-proteins-max-group (- number-proteins-max-group 1)) 2)))
                              (if (< subtracao ppiterlimit)
                                (progn
                                 (setf resto (- ppiterlimit subtracao))
                                 (do ((x 2 (incf x)))
                                   ((> (/(* x (- x 1)) 2) resto));if condition is greater than the rest, then to the loop.
                                   (setf indice-de-corte x)
                                   );of the
                                 (setf (second (first (last grupos-identicos))) (split-up-to (second (first (last grupos-identicos))) indice-de-corte))
                                 ;(setf (third (first (last grupos-identicos))) (length (second (first (last grupos-identicos)))))
                                 );progn
                                (progn
                                 (setf grupos-identicos (remove (first (last grupos-identicos)) grupos-identicos :test #'equalp))
                                 );progn
                                );if-subtraction < threshold
                              );progn|#
                           );if-total-edges < threshold
                         );loop
                       );block remove-groups
                ;;End of group removal

                ;;Ordering the proteins of each remaining group before creating interactions.
                (dolist (g grupos-identicos)
                        ;;Ordering the proteins of each group according to their index in the  genome.
                        (setf (second g) (sort (second g) #'< :key #'second))
                        )

                ;;Ordering groups by size, from smallest to largest.
                (setf grupos-identicos (sort grupos-identicos #'< :key #'third))

                ;;Removing groups with the same amount of genes and assigning the result to the group weights list.
                (setf pesos-grupos (remove-duplicates grupos-identicos :key #'third))

                ;;Generating interactions for each pair of proteins in each group.
                (dolist (g grupos-identicos)

                        ;;Picking up the group index of size equal to group g in the group weights list.
                        (setf dividendo (1+ (position-if #'(lambda (x) (= (third g) (third x))) pesos-grupos)))

                        ;;Picking up the index of the most populous group, that is, the last group.
                        (setf divisor (length pesos-grupos))

                        ;;Calculating the weight of group g based on their position, and the position of the last group (most populous).
                        (setf peso-grupo (* (/ dividendo divisor) percentage-pp))

                        ;;Generating  PPIs for every pair of proteins in group g.
                        (dotimes (i (1- (length(second g))))
                                 (loop for j from (1+ i) to (1- (length(second g)))
                                   do (progn
                                       (push (make-ppi-struct
                                              :genea (first (elt (second g) i))
                                              :geneb (first (elt (second g) j))
                                              :weight peso-grupo
                                              :position (second (elt (second g) i)))
                                             ppi)
                                       );progn-do
                                   );loop-dotimes
                                 );dotimes
                        );dolist-id  groups
                );progn
               );condition and action 1

              ;;If the user wants to consider similar profiles
              ((> pptolerance 0)
               (progn

                ;;Ordering groups by size, from largest to smallest.
                (setf grupos-identicos (sort grupos-identicos #'> :key #'third))

                ;;Adding to a group g, all groups similar to group g.
                (dotimes (i (1- (length grupos-identicos)))
                         (setf grupos-similares  (append grupos-similares (list (list (elt grupos-identicos i) (list )))))
                         (setf soma-grupos (+ soma-grupos (third (elt grupos-identicos i))))
                         (loop for j from (1+ i) to (1- (length grupos-identicos))
                           do (progn
                               (unless (> (comparar-perfis (first (elt grupos-identicos i)) (first (elt grupos-identicos j))) pptolerance)
                                 ;;Unless the difference between the profiles of a pair of groups is greater than the tolerated difference, make:
                                 ;;Adding to the gi group, the gj group, which is similar to the gi group.
                                 (push (elt grupos-identicos j) (second (elt grupos-similares i)))
                                 ;;Summing up the total amount of proteins between similar groups.
                                 (setf soma-grupos (+ soma-grupos (third (elt grupos-identicos j))))
                                 );unlles
                               );progn-do
                           );loop-dotimes
                         (setf (elt grupos-similares i) (nconc (elt grupos-similares i) (list soma-grupos)))
                         (setf soma-grupos 0)
                         );dotimes

                ;;Adding the last group that was not considered in the above dotimes, to the list of similar groups.
                (setf grupos-similares  (append grupos-similares (list (list (first (last grupos-identicos)) (list ) (length (second (first (last grupos-identicos))))))))

                ;;Ordering groups by size, from smallest to largest.
                (setf grupos-similares (sort grupos-similares #'< :key #'(lambda (x) (third (first x)))))

                (block remove-groups
                       (loop
                         (setf total-arestas (loop for g in grupos-similares
                                               summing (+
                                                        ;;Calculating how many pairs of genes with identical profiles will be generated.
                                                        (/ (* (length (second (first g))) (- (length (second (first g))) 1)) 2)
                                                        ;;Multiplying the number of genes of identical profiles by the number of genes of similar profiles.
                                                        (* (length (second (first g)))
                                                           ;;Summing up the amount of genes of similar profiles.
                                                           (loop for gs in (second g) summing (length (second gs)) into total finally (return total))
                                                           );*
                                                        );+
                                               into total
                                               finally (return total))
                               );setf full-edges
                         (if (<= total-arestas ppiterlimit)
                           (return-from remove-groups)
                           (setf grupos-similares (remove (first (last grupos-similares)) grupos-similares :test #'equalp))
                           #|(progn
                              ;(setf number-proteins-max-group (third (first (last grupos-similares))))
                              (setf max-proteins-group (first (last grupos-similares)))
                              (setf grupos-similares (remove (first (last grupos-similares)) grupos-similares :test #'equalp))
                              (setf subtracao (- total-arestas
                                                 (+
                                                  ;;Calculating how many pairs of genes with identical profiles will be generated.
                                                  (/ (* (length (second (first max-proteins-group))) (- (length (second (first max-proteins-group))) 1)) 2)
                                                  ;;Multiplying the number of genes of identical profiles by the number of genes of similar profiles.
                                                  (* (length (second (first max-proteins-group)))
                                                     ;;Summing up the amount of genes from similar profiles.
                                                     (loop for gs in (second max-proteins-group) summing (length (second gs)) into total finally (return total))
                                                     );*
                                                  );+
                                                 );-
                                    );setf subtraction
                              (if (< subtracao ppiterlimit)
                                (progn

                                 (setf resto (- ppiterlimit subtracao))
                                 (do ((x 2 (incf x)))
                                   ((> (/(* x (- x 1)) 2) resto));if condition is greater than the rest, then to the loop.
                                   (setf indice-de-corte x)
                                   );do
                                 (setf (second (first (last grupos-identicos))) (split-up-to (second (first (last grupos-identicos))) indice-de-corte))
                                 (setf (third (first (last grupos-identicos))) (length (second (first (last grupos-identicos)))))
                                 );progn
                                (progn
                                 (setf grupos-similares (remove max-proteins-group grupos-similares :test #'equalp))
                                 );progn
                                );if-subtraction < threshold
                              );progn|#
                           );if-total-< edges
                         );loop
                       );block remove-groups
                ;;End of group removal

                ;;Ordering the proteins of each group according to their index in the  genome.
                (dolist (g grupos-similares)
                        (setf (second (first g)) (sort (second (first g)) #'< :key #'second))
                        )

                ;;Ordering similar groups by size, from smallest to largest.
                (setf grupos-similares (sort grupos-similares #'< :key #'(lambda (x) (third (first x)))))

                ;;Removing groups with the same amount of genes and assigning the result to the group weights list.
                (setf pesos-grupos (remove-duplicates grupos-similares :key #'(lambda (x) (third (first x)))))

                ;;Sweeping the groups formed to create edges between  identical  and  similar proteins.
                (dolist (g grupos-similares)

                        ;;Picking up the group index of size equal to group g in the group weights list.
                        (setf dividendo (1+ (position-if #'(lambda (x) (= (third (first g)) (third (first x)))) pesos-grupos)))

                        ;;Picking up the index of the most populous group, that is, the last group.
                        (setf divisor (length pesos-grupos))

                        ;;Calculating the weight of group g based on their position, and the position of the last group (most populous).
                        (setf peso-grupo (* (/ dividendo divisor) percentage-pp))

                        ;;Creating edges between proteins with identical profiles.
                        (dotimes (i (1- (third (first g))))
                                 (loop for j from (1+ i) to (1- (third (first g)))
                                   do (progn
                                       (push (make-ppi-struct
                                              :genea (first (elt (second (first g)) i))
                                              :geneb (first (elt (second (first g)) j))
                                              :weight peso-grupo
                                              :position (second (elt (second (first g)) i)))
                                             ppi)
                                       );progn-do
                                   );loop j
                                 );dotimes i

                        ;;Sweeping protein groups with profiles similar to the profile of the protein group with identical profiles.
                        (dolist (grupo-similar (second g))
                                ;Calculating the  difference  between profiles.
                                (setf diferenca-perfis (comparar-perfis (first (first g)) (first grupo-similar)))
                                ;;Sweeping the proteins of a group with a similar profile.
                                (dolist (ps (second grupo-similar))
                                        ;;Sweeping proteins   of identical profiles of a group g.
                                        (dolist (pg (second (first g)))
                                                ;;Creating  an interaction.
                                                (push (make-ppi-struct
                                                 						:genea (if (< (second pg) (second ps)) (first pg) (first ps) )
                                                 						:geneb (if (> (second pg) (second ps)) (first pg) (first ps) )
                                                 						:weight (- peso-grupo (* peso-grupo (/ (/ (* diferenca-perfis 100) (length *genomes-files*)) 100)))
                                                 						:position (if (< (second pg) (second ps))  (second pg) (second ps) ) )
                                          					       ppi);push
                                                );dolist pg
                                        );dolist ps
                                );similar-group dolist
                        );similar groups dolist
                );progn
               );condition and action 2
              );cond
            );unless

          ;;If any PPI has  beenpredicted, the PPI identifier is  changed  to  true.
          (if (> (length ppi) 0)
            (setf *ppi-identified-pp* t)
            )

          ;;Recording the  PPI of genome k in the database.
          (save-db (concatenate 'string diretorio-banco-ppi k ".db") ppi)

          ;;Reset the variables.
          (setf grupos-identicos nil)
          (setf grupos-similares nil)
          (setf ppi nil)

          (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *perfil-filogenetico*)) 100)) laco)))
                   (format *query-io* "=")
                   (force-output *query-io*)
                   )
          (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *perfil-filogenetico*)) 100)) laco))))
          (incf genomenumber)

          );progn-loop
      );loop hash table
    (format t "]~%")
    ))
;-------------------------------------------------------------------------------

(defun phylogenetic-profiles-trim (percentage-pp pptolerance trim workingdir)

  (let ((ppi (list));;List to   store the interactions created.
                    (ppi-by-score (list));List to store  PPI's   per score.
                    (grupos-identicos (list));List for storing protein groups with identical profiles.
                    (grupos-similares (list));List for storing protein groups with identical and similar profiles.
                    (pesos-grupos (list))
                    (posicao-grupo);Variable to store the position of a group g.
                    ;(number-proteínas-max-group);Variável para armezar o número de proteínas do maior grupo.
                    (soma-grupos 0);Variable to store the total number of proteins between groups.
                    (peso-grupo);Variable to store the weight that will be assigned to proteins in a group g.
                    (dividendo);Variable that stores the index of a t-size group, which is closer to the largest group.
                    (divisor);Variable that stores the index of the largest group.
                    (diferenca-perfis);Variable to store the difference of  phylogenetic profiles of two groups.
                    (diretorio-banco-ppi)
                    (laco 0)
                    (genomenumber 0)
                    )

    (declare (type single-float percentage-pp))
    (declare (type fixnum laco))
    (declare (type fixnum genomenumber))
    (declare (type fixnum pptolerance))
    (declare (type fixnum trim))

    (format t "~%Predicting ppi by phylogenetic profiles;~%")
    (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
    (format *query-io* "[")
    (force-output *query-io*)

    (setf diretorio-banco-ppi (concatenate 'string workingdir "/database/ppi-phylogenetic-profiles/"))
    (ensure-directories-exist diretorio-banco-ppi)
    (unless (eql (list-directory diretorio-banco-ppi) nil)
      (dolist (file (list-directory diretorio-banco-ppi))
              (unless (eql (pathname-type (pathname file)) nil)
                (delete-file file)
                );unless
              );dolist
      );unless

    (loop for k being the hash-keys in *perfil-filogenetico* using (hash-value perfis)
      do (progn
          ;;Grouping proteins by phylogenetic profile
          (dolist (p perfis)
                  (setf posicao-grupo (position-if #'(lambda (x) (eql nil (set-exclusive-or (second p) (first x) :test #'string=))) grupos-identicos))
                  (if (eql posicao-grupo nil)
                    (push (list (second p) (list (list (first p) (third p)))) grupos-identicos)
                    (push (list (first p) (third p)) (second (elt grupos-identicos posicao-grupo)))
                    );if
                  );dolist
          ;;End of the Grouping

          ;;Ordering the proteins of each group and concatenating the size of the group
          (dolist (g grupos-identicos)
                  ;;Ordering the proteins of each group according to their index in the  genome.
                  (setf (second g) (sort (second g) #'< :key #'second))
                  ;;Concatenating to a group g, the size of group g.
                  (setf g (nconc g (list (length (second g)))))
                  )

          (unless (= (length grupos-identicos) 0)
            (cond
              ;;If the user only wants identical profiles
              ((= pptolerance 0)
               (progn

                ;;Ordering groups by size, from smallest to largest.
                (setf grupos-identicos (sort grupos-identicos #'< :key #'third))

                ;;Removing groups with the same amount of genes and assigning the result to the group weights list.
                (setf pesos-grupos (remove-duplicates grupos-identicos :key #'third))

                ;;Generating interactions for each pair of proteins in each group.
                (dolist (g grupos-identicos)

                        ;;Picking up the group index of size equal to group g in the group weights list.
                        (setf dividendo (1+ (position-if #'(lambda (x) (= (third g) (third x))) pesos-grupos)))

                        ;;Picking up the index of the most populous group, that is, the last group.
                        (setf divisor (length pesos-grupos))

                        ;;Calculating the weight of group g based on their position, and the position of the last group (most populous).
                        (setf peso-grupo (* (/ dividendo divisor) percentage-pp))

                        (unless (= (third g) 1)
                          (push (list peso-grupo (list)) ppi-by-score)
                          )

                        ;;Generating  PPIs for every pair of proteins in group g.
                        (dotimes (i (1- (third g)))
                                 (loop for j from (1+ i) to (1- (third g))
                                   do (progn
                                       (push (make-ppi-struct
                                              :genea (first (elt (second g) i))
                                              :geneb (first (elt (second g) j))
                                              :weight peso-grupo
                                              :position (second (elt (second g) i)))
                                             ppi)
                                       );progn-do
                                   );loop-dotimes
                                 );dotimes
                        );  dolist-identicalgroups

                (setf ppi-by-score (remove-duplicates ppi-by-score :test #'= :key #'first))

                (dolist (p ppi)
                        (setf posicao-grupo (position-if #'(lambda (x) (= (ppi-struct-weight p) (first x))) ppi-by-score))
                        (push p (second (elt ppi-by-score posicao-grupo)))
                        );dolist

                (setf ppi nil)

                (dolist (p ppi-by-score)
                        (setf (second p) (sort (second p) #'< :key #'(lambda (x) (ppi-struct-position x))))
                        (if (<= (length (second p)) trim)
                          (setf ppi (append ppi (second p)))
                          (setf ppi (append ppi (split-up-to (second p) trim)))
                          );if
                        );dolist

                (setf ppi-by-score nil)

                );progn
               );condition and action 1

              ;;If the user wants to consider similar profiles
              ((> pptolerance 0)
               (progn

                ;;Ordering groups by size, from largest to smallest.
                (setf grupos-identicos (sort grupos-identicos #'> :key #'third))

                ;;Adding to a group g, all groups similar to group g.
                (dotimes (i (1- (length grupos-identicos)))
                         (setf grupos-similares  (append grupos-similares (list (list (elt grupos-identicos i) (list )))))
                         (setf soma-grupos (+ soma-grupos (third (elt grupos-identicos i))))
                         (loop for j from (1+ i) to (1- (length grupos-identicos))
                           do (progn
                               (unless (> (comparar-perfis (first (elt grupos-identicos i)) (first (elt grupos-identicos j))) pptolerance)
                                 ;;Unless the difference between the profiles of a pair of groups is greater than the tolerated difference, make:
                                 ;;Adding to the gi group, the gj group, which is similar to the gi group.
                                 (push (elt grupos-identicos j) (second (elt grupos-similares i)))
                                 ;;Summing up the total amount of proteins between similar groups.
                                 (setf soma-grupos (+ soma-grupos (third (elt grupos-identicos j))))
                                 );unlles
                               );progn-do
                           );loop-dotimes
                         (setf (elt grupos-similares i) (nconc (elt grupos-similares i) (list soma-grupos)))
                         (setf soma-grupos 0)
                         );dotimes

                ;;Adding the last group that was not considered in the above dotimes, to the list of similar groups.
                (setf grupos-similares  (append grupos-similares (list (list (first (last grupos-identicos)) (list ) (length (second (first (last grupos-identicos))))))))

                ;;Ordering similar groups by size, from smallest to largest.
                (setf grupos-similares (sort grupos-similares #'< :key #'(lambda (x) (third (first x)))))

                ;;Removing groups with the same amount of genes and assigning the result to the group weights list.
                (setf pesos-grupos (remove-duplicates grupos-similares :key #'(lambda (x) (third (first x)))))

                ;;Sweeping the groups formed to create edges between  identical  and  similar proteins.
                (dolist (g grupos-similares)

                        ;;Picking up the group index of size equal to group g in the group weights list.
                        (setf dividendo (1+ (position-if #'(lambda (x) (= (third (first g)) (third (first x)))) pesos-grupos)))

                        ;;Taking the index of the most populous group, that is, the last group.
                        (setf divisor (length pesos-grupos))

                        ;;Calculating the weight of group g based on their position, and the position of the last group (most populous).
                        (setf peso-grupo (* (/ dividendo divisor) percentage-pp))

                        (unless (= (third g) 1)
                          (push (list peso-grupo (list)) ppi-by-score)
                          )

                        ;;Creating edges between proteins with identical profiles.
                        (dotimes (i (1- (third (first g))))
                                 (loop for j from (1+ i) to (1- (third (first g)))
                                   do (progn
                                       (push (make-ppi-struct
                                              :genea (first (elt (second (first g)) i))
                                              :geneb (first (elt (second (first g)) j))
                                              :weight peso-grupo
                                              :position (second (elt (second (first g)) i)))
                                             ppi)
                                       );progn-do
                                   );loop j
                                 );dotimes i

                        ;;Sweeping protein groups with profiles similar to the profile of the protein group with identical profiles.
                        (dolist (grupo-similar (second g))
                                ;;Calculating the  difference  between profiles.
                                (setf diferenca-perfis (comparar-perfis (first (first g)) (first grupo-similar)))
                                ;;by inserting the new weight generated into the ppi  s list by weight.
                                (push (list (- peso-grupo (* peso-grupo (/ (/ (* diferenca-perfis 100) (length *genomes-files*)) 100))) (list)) ppi-by-score)
                                ;;Sweeping the proteins of a group with similar profile.
                                (dolist (ps (second grupo-similar))
                                        ;;Sweeping proteins   of identical profiles of a group g.
                                        (dolist (pg (second (first g)))
                                                ;;Creating  an interaction.
                                                (push (make-ppi-struct
                                                 						:genea (if (< (second pg) (second ps)) (first pg) (first ps) )
                                                 						:geneb (if (> (second pg) (second ps)) (first pg) (first ps) )
                                                 						:weight (- peso-grupo (* peso-grupo (/ (/ (* diferenca-perfis 100) (length *genomes-files*)) 100)))
                                                 						:position (if (< (second pg) (second ps))  (second pg) (second ps) ) )
                                               					  ppi);push
                                                );dolist pg
                                        );dolist ps
                                );similar group-dolist
                        );dolist similar groups

                (setf ppi-by-score (remove-duplicates ppi-by-score :test #'= :key #'first))

                (dolist (p ppi)
                        (setf posicao-grupo (position-if #'(lambda (x) (= (ppi-struct-weight p) (first x))) ppi-by-score))
                        (push p (second(elt ppi-by-score posicao-grupo)))
                        );dolist

                (setf ppi nil)

                (dolist (p ppi-by-score)
                        (setf (second p) (sort (second p) #'< :key #'(lambda (x) (ppi-struct-position x))))
                        (if (<= (length (second p)) trim)
                          (setf ppi (append ppi (second p)))
                          (setf ppi (append ppi (split-up-to (second p) trim)))
                          );if
                        );dolist

                (setf ppi-by-score nil)

                );progn
               );condition and action 2
              );cond
            );unless

          ;;If any PPI has  beenpredicted, the PPI identifier is  changed  to  true.
          (if (> (length ppi) 0)
            (setf *ppi-identified-pp* t)
            )

          ;;Recording the  PPI of genome k in the database.
          (save-db (concatenate 'string diretorio-banco-ppi k ".db") ppi)

          ;;Reset the variables.
          (setf grupos-identicos nil)
          (setf grupos-similares nil)
          (setf ppi nil)

          (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *perfil-filogenetico*)) 100)) laco)))
                   (format *query-io* "=")
                   (force-output *query-io*)
                   )
          (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *perfil-filogenetico*)) 100)) laco))))
          (incf genomenumber)

          );progn-loop
      );loop hash table
    (format t "]~%")
    ))
;-------------------------------------------------------------------------------

(defun phylogenetic-profiles-threshold (percentage-pp pptolerance ppthreshold plusminus workingdir)

  (let ((ppi (list));;List to   store the interactions created.
                    (grupos-identicos (list));List for storing protein groups with identical profiles.
                    (grupos-similares (list));List for storing protein groups with identical and similar profiles.
                    (pesos-grupos (list))
                    (posicao-grupo);Variable to store the position of a group g.
                    ;(number-proteínas-max-group);Variável para armezar o número de proteínas do maior grupo.
                    (soma-grupos 0);Variable to store the total number of proteins between groups.
                    (peso-grupo);Variable to store the weight that will be assigned to proteins in a group g.
                    (dividendo);Variable that stores the index of a t-size group, which is closer to the largest group.
                    (divisor);Variable that store the index of the largest group.
                    (diferenca-perfis);Variable to store the  difference  of phylogenetic profiles of two groups.
                    (diretorio-banco-ppi)
                    (laco 0)
                    (genomenumber 0)
                    )

    (declare (type single-float percentage-pp))
    (declare (type fixnum laco))
    (declare (type fixnum genomenumber))
    (declare (type fixnum pptolerance))
    (declare (type fixnum ppthreshold))
    (declare (type (simple-array character (1)) plusminus))

    (format t "~%Predicting ppi by phylogenetic profiles;~%")
    (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
    (format *query-io* "[")
    (force-output *query-io*)

    (setf diretorio-banco-ppi (concatenate 'string workingdir "/database/ppi-phylogenetic-profiles/"))
    (ensure-directories-exist diretorio-banco-ppi)
    (unless (eql (list-directory diretorio-banco-ppi) nil)
      (dolist (file (list-directory diretorio-banco-ppi))
              (unless (eql (pathname-type (pathname file)) nil)
                (delete-file file)
                );unless
              );dolist
      );unless

    (loop for k being the hash-keys in *perfil-filogenetico* using (hash-value perfis)
      do (progn
          ;;Grouping proteins by phylogenetic profile
          (dolist (p perfis)
                  (setf posicao-grupo (position-if #'(lambda (x) (eql nil (set-exclusive-or (second p) (first x) :test #'string=))) grupos-identicos))
                  (if (eql posicao-grupo nil)
                    (progn (push (list (second p) (list (list (first p) (third p)))) grupos-identicos))
                    (progn (push (list (first p) (third p)) (second (elt grupos-identicos posicao-grupo))))
                    );if
                  );dolist
          ;;Final of the Grouping

          ;;Discarding groups less than the threshold reported by the user.
          (if (string= plusminus "<")
            (dolist (g grupos-identicos)
                    (unless (< (length (first g)) ppthreshold)
                      (setf grupos-identicos (remove g grupos-identicos :test #'equalp))
                      );unless
                    );dolist
            (dolist (g grupos-identicos)
                    (unless (> (length (first g)) ppthreshold)
                      (setf grupos-identicos (remove g grupos-identicos :test #'equalp))
                      );unless
                    );dolist
            );if string=

          ;;Ordering the proteins of each group
          (dolist (g grupos-identicos)
                  ;;Ordering the rotein  pof each group according to its index in the genome.
                  (setf (second g) (sort (second g) #'< :key #'second))
                  ;;Concatenando to a group g, the size of group g.
                  (setf g (nconc g (list (length (second g)))))
                  );dolist

          (unless (= (length grupos-identicos) 0)
            (cond
              ;;If the user only wants identical profiles
              ((= pptolerance 0)
               (progn

                ;;Ordering groups by size, from smallest to largest.
                (setf grupos-identicos (sort grupos-identicos #'< :key #'third))

                ;;Removing groups with the same amount of genes and attributing the result to the list of group weights.
                (setf pesos-grupos (remove-duplicates grupos-identicos :key #'third))

                ;;Generating interactions for each pair of proteins of each group.
                (dolist (g grupos-identicos)

                        ;;Picking up the group index of size equal to group g in the group weights list.
                        (setf dividendo (1+ (position-if #'(lambda (x) (= (third g) (third x))) pesos-grupos)))

                        ;;Picking up the index of the most populous group, that is, the last group.
                        (setf divisor (length pesos-grupos))

                        ;;Calculating the weight of group g based on their position, and the position of the last group (most populous).
                        (setf peso-grupo (* (/ dividendo divisor) percentage-pp))

                        ;;Generating  PPIs for every pair of proteins in group g.
                        (dotimes (i (1- (third g)))
                                 (loop for j from (1+ i) to (1- (third g))
                                   do (progn
                                       (push (make-ppi-struct
                                              :genea (first (elt (second g) i))
                                              :geneb (first (elt (second g) j))
                                              :weight peso-grupo
                                              :position (second (elt (second g) i)))
                                             ppi)
                                       );progn-do
                                   );loop-dotimes
                                 );dotimes
                        );dolist-groups-identicos
                );progn
               );condition and action 1

              ;;If the user wants to consider similar profiles
              ((> pptolerance 0)
               (progn

                ;;Ordering groups by size, from largest to smallest.
                (setf grupos-identicos (sort grupos-identicos #'> :key #'third))

                ;;Adding to a group g, all groups similar to group g.
                (dotimes (i (1- (length grupos-identicos)))
                         (setf grupos-similares  (append grupos-similares (list (list (elt grupos-identicos i) (list )))))
                         (setf soma-grupos (+ soma-grupos (third (elt grupos-identicos i))))
                         (loop for j from (1+ i) to (1- (length grupos-identicos))
                           do (progn
                               (unless (> (comparar-perfis (first (elt grupos-identicos i)) (first (elt grupos-identicos j))) pptolerance)
                                 ;;Unless the difference between the profiles of a pair of groups is greater than the tolerated difference, make:
                                 ;;Adding to the gi group, the gj group, which is similar to the gi group.
                                 (push (elt grupos-identicos j) (second (elt grupos-similares i)))
                                 ;;Summing up the total amount of proteins between similar groups.
                                 (setf soma-grupos (+ soma-grupos (third (elt grupos-identicos j))))
                                 );unlles
                               );progn-do
                           );loop-dotimes
                         (setf (elt grupos-similares i) (nconc (elt grupos-similares i) (list soma-grupos)))
                         (setf soma-grupos 0)
                         );dotimes

                ;;Adding the last group that was not considered in the above dotimes, to the list of similar groups.
                (setf grupos-similares  (append grupos-similares (list (list (first (last grupos-identicos)) (list ) (length (second (first (last grupos-identicos))))))))

                ;;Ordering similar groups by size, from smallest to largest.
                (setf grupos-similares (sort grupos-similares #'< :key #'(lambda (x) (third (first x)))))

                ;;Removing groups with the same amount of genes and assigning the result to the group weights list.
                (setf pesos-grupos (remove-duplicates grupos-similares :key #'(lambda (x) (third (first x)))))

                ;;Sweeping the groups formed to create edges between  identical  and  similar proteins.
                (dolist (g grupos-similares)

                        ;;Picking up the group index of size equal to group g in the group weights list.
                        (setf dividendo (1+ (position-if #'(lambda (x) (= (third (first g)) (third (first x)))) pesos-grupos)))

                        ;;Picking up the index of the most populous group, that is, the last group.
                        (setf divisor (length pesos-grupos))

                        ;;Calculating group g weight based on their position, and position of the last group (most populous).
                        (setf peso-grupo (* (/ dividendo divisor) percentage-pp))

                        ;;Creating edges between proteins with identical profiles.
                        (dotimes (i (1- (third (first g))))
                                 (loop for j from (1+ i) to (1- (third (first g)))
                                   do (progn
                                       (push (make-ppi-struct
                                              :genea (first (elt (second (first g)) i))
                                              :geneb (first (elt (second (first g)) j))
                                              :weight peso-grupo
                                              :position (second (elt (second (first g)) i)))
                                             ppi)
                                       );progn-do
                                   );loop j
                                 );dotimes i

                        ;;Sweeping protein groups with profiles similar to the profile of the protein group with identical profiles.
                        (dolist (grupo-similar (second g))
                                ;Calculating the  difference  between profiles.
                                (setf diferenca-perfis (comparar-perfis (first (first g)) (first grupo-similar)))
                                ;;Sweeping the proteins of a group with a similar profile.
                                (dolist (ps (second grupo-similar))
                                        ;;Sweeping proteins   of identical profiles of a group g.
                                        (dolist (pg (second (first g)))
                                                ;;Creating  an interaction.
                                                (push (make-ppi-struct
                                                 						:genea (if (< (second pg) (second ps)) (first pg) (first ps) )
                                                 						:geneb (if (> (second pg) (second ps)) (first pg) (first ps) )
                                                 						:weight (- peso-grupo (* peso-grupo (/ (/ (* diferenca-perfis 100) (length *genomes-files*)) 100)))
                                                 						:position (if (< (second pg) (second ps))  (second pg) (second ps) ) )
                                          					       ppi);push
                                                );dolist pg
                                        );dolist ps
                                );similar-group dolist
                        );dolist similar groups
                );progn
               );condition and action 2
              );cond
            );unless

          ;;If any PPI has  beenpredicted, the PPI identifier is  changed  to  true.
          (if (> (length ppi) 0)
            (setf *ppi-identified-pp* t)
            )

          ;;Recording the  PPI of genome k in the database.
          (save-db (concatenate 'string diretorio-banco-ppi k ".db") ppi)

          ;;Reset the variables.
          (setf grupos-identicos nil)
          (setf grupos-similares nil)
          (setf ppi nil)

          (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *perfil-filogenetico*)) 100)) laco)))
                   (format *query-io* "=")
                   (force-output *query-io*)
                   )
          (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *perfil-filogenetico*)) 100)) laco))))
          (incf genomenumber)

          );progn-loop
      );loop hash table
    (format t "]~%")
    ))
;-------------------------------------------------------------------------------

(defun phylogenetic-profiles-delete-clusters (percentage-pp pptolerance grouplimit workingdir)

  (let ((ppi (list));;List to   store the interactions created.
                    (ppi-by-score (list));List to store  PPI's   by score.
                    (grupos-identicos (list));List for storing protein groups with identical profiles.
                    (grupos-similares (list));List for storing protein groups with identical and similar profiles.
                    (pesos-grupos (list))
                    (posicao-grupo);Variable to store the position of a group g.
                    ;(number-proteínas-max-group);Variável para armezar o número de proteínas do maior grupo.
                    (soma-grupos 0);Variable to store the total number of proteins between groups.
                    (peso-grupo);Variable to store the weight that will be assigned to proteins in a group g.
                    (dividendo);Variable that stores the index of a t-size group, which is closer to the largest group.
                    (divisor);Variable that store the index of the largest group.
                    (diferenca-perfis);Variable to store the  difference  of phylogenetic profiles of two groups.
                    (diretorio-banco-ppi)
                    (laco 0)
                    (genomenumber 0)
                    )

    (declare (type single-float percentage-pp))
    (declare (type fixnum laco))
    (declare (type fixnum genomenumber))
    (declare (type fixnum pptolerance))
    (declare (type fixnum grouplimit))

    (format t "~%Predicting ppi by phylogenetic profiles;~%")
    (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
    (format *query-io* "[")
    (force-output *query-io*)

    (setf diretorio-banco-ppi (concatenate 'string workingdir "/database/ppi-phylogenetic-profiles/"))
    (ensure-directories-exist diretorio-banco-ppi)
    (unless (eql (list-directory diretorio-banco-ppi) nil)
      (dolist (file (list-directory diretorio-banco-ppi))
              (unless (eql (pathname-type (pathname file)) nil)
                (delete-file file)
                );unless
              );dolist
      );unless

    (loop for k being the hash-keys in *perfil-filogenetico* using (hash-value perfis)
      do (progn
          ;;Grouping proteins by phylogenetic profile
          (dolist (p perfis)
                  (setf posicao-grupo (position-if #'(lambda (x) (eql nil (set-exclusive-or (second p) (first x) :test #'string=))) grupos-identicos))
                  (if (eql posicao-grupo nil)
                    (push (list (second p) (list (list (first p) (third p)))) grupos-identicos)
                    (push (list (first p) (third p)) (second (elt grupos-identicos posicao-grupo)))
                    );if
                  );dolist
          ;;End of Grouping

          ;;Removing groups that would produce an unwanted number of edges with the same score.
          (dolist (g grupos-identicos)
                  (unless (<= (/ (* (length (second g)) (- (length (second g)) 1)) 2) grouplimit)
                    (setf grupos-identicos (remove g grupos-identicos :test #'equalp))
                    );unless
                  );dolist

          ;;Ordering the proteins of each group
          (dolist (g grupos-identicos)
                  ;;Ordering the proteins of each group according to their index in the  genome.
                  (setf (second g) (sort (second g) #'< :key #'second))
                  ;;Concatenando to a group g, the size of group g.
                  (setf g (nconc g (list (length (second g)))))
                  )

          (unless (= (length grupos-identicos) 0)
            (cond
              ;;If the user only wants identical profiles
              ((= pptolerance 0)
               (progn

                ;;Ordering groups by size, from smallest to largest.
                (setf grupos-identicos (sort grupos-identicos #'< :key #'third))

                ;;Removing groups with the same amount of genes and assigning the result to the group weights list.
                (setf pesos-grupos (remove-duplicates grupos-identicos :key #'third))

                ;;Generating interactions for each pair of proteins in each group.
                (dolist (g grupos-identicos)

                        ;;Picking up the group index of size equal to group g in the group weights list.
                        (setf dividendo (1+ (position-if #'(lambda (x) (= (third g) (third x))) pesos-grupos)))

                        ;;Picking up the index of the most populous group, i.e. the last group.
                        (setf divisor (length pesos-grupos))

                        ;;Calculating the weight of group g based on their position, and the position of the last group (most populous).
                        (setf peso-grupo (* (/ dividendo divisor) percentage-pp))

                        ;;Generating  PPIs for every pair of proteins in group g.
                        (dotimes (i (1- (third g)))
                                 (loop for j from (1+ i) to (1- (third g))
                                   do (progn
                                       (push (make-ppi-struct
                                              :genea (first (elt (second g) i))
                                              :geneb (first (elt (second g) j))
                                              :weight peso-grupo
                                              :position (second (elt (second g) i)))
                                             ppi)
                                       );progn-do
                                   );loop-dotimes
                                 );dotimes
                        );dolist-groups-identical
                );progn
               );condition and action 1

              ;;If the user wants to consider similar profiles
              ((> pptolerance 0)
               (progn

                ;;Ordering groups by size, from largest to smallest.
                (setf grupos-identicos (sort grupos-identicos #'> :key #'third))

                ;;Adding to a group g, all groups similar to group g.
                (dotimes (i (1- (length grupos-identicos)))
                         (setf grupos-similares  (append grupos-similares (list (list (elt grupos-identicos i) (list )))))
                         (setf soma-grupos (+ soma-grupos (third (elt grupos-identicos i))))
                         (loop for j from (1+ i) to (1- (length grupos-identicos))
                           do (progn
                               (unless (> (comparar-perfis (first (elt grupos-identicos i)) (first (elt grupos-identicos j))) pptolerance)
                                 ;;Unless the difference between the profiles of a pair of groups is greater than the tolerated difference, make:
                                 ;;Adding to the gi group, the gj group, which is similar to the gi group.
                                 (push (elt grupos-identicos j) (second (elt grupos-similares i)))
                                 ;;Summing the total amount of proteins between similar groups.
                                 (setf soma-grupos (+ soma-grupos (third (elt grupos-identicos j))))
                                 );unlles
                               );progn-do
                           );loop-dotimes
                         (setf (elt grupos-similares i) (nconc (elt grupos-similares i) (list soma-grupos)))
                         (setf soma-grupos 0)
                         );dotimes

                ;;Adding the last group that was not considered in the above dotimes, to the list of similar groups.
                (setf grupos-similares  (append grupos-similares (list (list (first (last grupos-identicos)) (list ) (length (second (first (last grupos-identicos))))))))

                ;;Ordering similar groups by size, from smallest to largest.
                (setf grupos-similares (sort grupos-similares #'< :key #'(lambda (x) (third (first x)))))

                ;;Removing groups with the same amount of genes and assigning the result to the group weights list.
                (setf pesos-grupos (remove-duplicates grupos-similares :key #'(lambda (x) (third (first x)))))

                ;;Sweeping the groups formed to create edges between  identical  and  similar proteins.
                (dolist (g grupos-similares)

                        ;;Picking up the group index of size equal to group g in the group weights list.
                        (setf dividendo (1+ (position-if #'(lambda (x) (= (third (first g)) (third (first x)))) pesos-grupos)))

                        ;;Picking up the index of the most populous group, that is, the last group.
                        (setf divisor (length pesos-grupos))

                        ;;Calculating the weight of group g based on their position, and the position of the last group (most populous).
                        (setf peso-grupo (* (/ dividendo divisor) percentage-pp))

                        (unless (= (third g) 1)
                          (push (list peso-grupo (list)) ppi-by-score)
                          )

                        ;;Creating edges between proteins with profiles identical to each other.
                        (dotimes (i (1- (third (first g))))
                                 (loop for j from (1+ i) to (1- (third (first g)))
                                   do (progn
                                       (push (make-ppi-struct
                                              :genea (first (elt (second (first g)) i))
                                              :geneb (first (elt (second (first g)) j))
                                              :weight peso-grupo
                                              :position (second (elt (second (first g)) i)))
                                             ppi)
                                       );progn-do
                                   );loop j
                                 );dotimes i

                        ;;Sweeping protein groups with profiles similar to the profile of the protein group with identical profiles.
                        (dolist (grupo-similar (second g))
                                ;Calculating the  difference  between profiles.
                                (setf diferenca-perfis (comparar-perfis (first (first g)) (first grupo-similar)))
                                ;;by inserting the new weight generated into the  ppis list  by weight.
                                (push (list (- peso-grupo (* peso-grupo (/ (/ (* diferenca-perfis 100) (length *genomes-files*)) 100))) (list)) ppi-by-score)
                                ;;Sweeping the proteins of a group with a similar profile.
                                (dolist (ps (second grupo-similar))
                                        ;;Sweeping proteins   of identical profiles of a group g.
                                        (dolist (pg (second (first g)))
                                                ;;Creating  an interaction.
                                                (push (make-ppi-struct
                                                 						:genea (if (< (second pg) (second ps)) (first pg) (first ps) )
                                                 						:geneb (if (> (second pg) (second ps)) (first pg) (first ps) )
                                                 						:weight (- peso-grupo (* peso-grupo (/ (/ (* diferenca-perfis 100) (length *genomes-files*)) 100)))
                                                 						:position (if (< (second pg) (second ps))  (second pg) (second ps) ) )
                                          					       ppi);push
                                                );dolist pg
                                        );dolist ps
                                );similar-group dolist

                        );similar group dolist

                (setf ppi-by-score (remove-duplicates ppi-by-score :test #'= :key #'first))

                (dolist (p ppi)
                        (setf posicao-grupo (position-if #'(lambda (x) (= (ppi-struct-weight p) (first x))) ppi-by-score))
                        (push p (second(elt ppi-by-score posicao-grupo)))
                        );dolist

                (setf ppi nil)

                (dolist (p ppi-by-score)
                        (setf (second p) (sort (second p) #'< :key #'(lambda (x) (ppi-struct-position x))))
                        (if (<= (length (second p)) grouplimit)
                          (setf ppi (append ppi (second p)))
                          (setf ppi (append ppi (split-up-to (second p) grouplimit)))
                          );if
                        );dolist

                (setf ppi-by-score nil)

                );progn
               );condition and action 2
              );cond
            );unless

          ;;If any PPI has  beenpredicted, the PPI identifier is  changed  to  true.
          (if (> (length ppi) 0)
            (setf *ppi-identified-pp* t)
            )

          ;;Recording the  PPI of genome k in the database.
          (save-db (concatenate 'string diretorio-banco-ppi k ".db") ppi)

          ;;Reset the variables.
          (setf grupos-identicos nil)
          (setf grupos-similares nil)
          (setf ppi nil)

          (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *perfil-filogenetico*)) 100)) laco)))
                   (format *query-io* "=")
                   (force-output *query-io*)
                   )
          (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *perfil-filogenetico*)) 100)) laco))))
          (incf genomenumber)

          );progn-loop
      );loop hash table
    (format t "]~%")
    ))

(defun phylogenetic-profiles-delete-perfil (percentage-pp pptolerance string-perfis workingdir)

  (let ((ppi (list));;List  to store the interactions created.
                    (grupos-identicos (list));List for storing protein groups with identical profiles.
                    (grupos-similares (list));List for storing protein groups with identical and similar profiles.
                    (pesos-grupos (list))
                    (perfis-indesejados (list))
                    (perfil)
                    (posicao-grupo);Variable to store the position of a group g.
                    ;(number-proteínas-max-group);Variável para armezar o número de proteínas do maior grupo.
                    (soma-grupos 0);Variable to store the total number of proteins between groups.
                    (peso-grupo);Variable to store the weight that will be assigned to proteins in a group g.
                    (dividendo);Variable that stores the index of a t-size group, which is closer to the largest group.
                    (divisor);Variable that store the index of the largest group.
                    (diferenca-perfis);Variable to store the  difference  of phylogenetic profiles of two groups.
                    (diretorio-banco-ppi)
                    (laco 0)
                    (genomenumber 0)
                    )

    (declare (type single-float percentage-pp))
    (declare (type fixnum laco))
    (declare (type fixnum genomenumber))
    (declare (type fixnum pptolerance))

    (format t "~%Predicting ppi by phylogenetic profiles;~%")
    (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
    (format *query-io* "[")
    (force-output *query-io*)

    (setf diretorio-banco-ppi (concatenate 'string workingdir "/database/ppi-phylogenetic-profiles/"))
    (ensure-directories-exist diretorio-banco-ppi)
    (unless (eql (list-directory diretorio-banco-ppi) nil)
      (dolist (file (list-directory diretorio-banco-ppi))
              (unless (eql (pathname-type (pathname file)) nil)
                (delete-file file)
                );unless
              );dolist
      );unless

    (setf string-perfis (concatenate 'string string-perfis ";"))

    ;;extracting unwanted string-profiles
    (loop for c across string-perfis
      do (progn
          (if (char= c #\;)
            (progn
             (unless (eql perfil nil)
               (setf perfis-indesejados (append perfis-indesejados (list (parse-integer perfil))))
               (setf perfil nil)
               )
             )
            (setf perfil (concatenate 'string perfil (list c)))
            );if char=
          );progn
      );loop

    (loop for k being the hash-keys in *perfil-filogenetico* using (hash-value perfis)
      do (progn
          ;;Grouping proteins by phylogenetic profile
          (dolist (p perfis)
                  (setf posicao-grupo (position-if #'(lambda (x) (eql nil (set-exclusive-or (second p) (first x) :test #'string=))) grupos-identicos))
                  (if (eql posicao-grupo nil)
                    (push (list (second p) (list (list (first p) (third p)))) grupos-identicos)
                    (push (list (first p) (third p)) (second (elt grupos-identicos posicao-grupo)))
                    );if
                  );dolist
          ;;End of Grouping

          ;;Discarding unwanted profiles.
          (dolist (g grupos-identicos)
                  (dolist (p perfis-indesejados)
                          (if (= p (length (first g)))
                            (setf grupos-identicos (remove g grupos-identicos :test #'equalp))
                            );if
                          );dolist-unwanted profiles
                  );dolist    identical groups

          ;;Ordering the proteins of each group
          (dolist (g grupos-identicos)
                  ;;Ordering the proteins of each group according to their index in the  genome.
                  (setf (second g) (sort (second g) #'< :key #'second))
                  ;;Concatenating to a group g, the size of group g.
                  (setf g (nconc g (list (length (second g)))))
                  );dolist

          (unless (= (length grupos-identicos) 0)
            (cond
              ;;If the user only wants identical profiles
              ((= pptolerance 0)
               (progn

                ;;Ordering groups by size, from smallest to largest.
                (setf grupos-identicos (sort grupos-identicos #'< :key #'third))

                ;;Removing groups with the same amount of genes and assigning the result to the group weights list.
                (setf pesos-grupos (remove-duplicates grupos-identicos :key #'third))

                ;;Generating interactions for each pair of proteins in each group.
                (dolist (g grupos-identicos)

                        ;;Picking up the group index of size equal to group g in the group weights list.
                        (setf dividendo (1+ (position-if #'(lambda (x) (= (third g) (third x))) pesos-grupos)))

                        ;;Picking up the index of the most populous group, i.e. the last group.
                        (setf divisor (length pesos-grupos))

                        ;;Calculating the weight of group g based on their position, and the position of the last group (most populous).
                        (setf peso-grupo (* (/ dividendo divisor) percentage-pp))

                        ;;Generating  PPIs  for every pair of proteins in group g.
                        (dotimes (i (1- (third g)))
                                 (loop for j from (1+ i) to (1- (third g))
                                   do (progn
                                       (push (make-ppi-struct
                                              :genea (first (elt (second g) i))
                                              :geneb (first (elt (second g) j))
                                              :weight peso-grupo
                                              :position (second (elt (second g) i)))
                                             ppi)
                                       );progn-do
                                   );loop-dotimes
                                 );dotimes
                        );dolist-groups-identicos
                );progn
               );condition and action 1

              ;;If the user wants to consider similar profiles
              ((> pptolerance 0)
               (progn

                ;;Ordering groups by size, from largest to smallest.
                (setf grupos-identicos (sort grupos-identicos #'> :key #'third))

                ;;Adding to a group g, all groups similar to group g.
                (dotimes (i (1- (length grupos-identicos)))
                         (setf grupos-similares  (append grupos-similares (list (list (elt grupos-identicos i) (list )))))
                         (setf soma-grupos (+ soma-grupos (third (elt grupos-identicos i))))
                         (loop for j from (1+ i) to (1- (length grupos-identicos))
                           do (progn
                               (unless (> (comparar-perfis (first (elt grupos-identicos i)) (first (elt grupos-identicos j))) pptolerance)
                                 ;;Unless the difference between the profiles of a pair of groups is greater than the tolerated difference, make:
                                 ;;Adding to the gi group, the gj group, which is similar to the gi group.
                                 (push (elt grupos-identicos j) (second (elt grupos-similares i)))
                                 ;;Summing up the total amount of proteins between similar groups.
                                 (setf soma-grupos (+ soma-grupos (third (elt grupos-identicos j))))
                                 );unlles
                               );progn-do
                           );loop-dotimes
                         (setf (elt grupos-similares i) (nconc (elt grupos-similares i) (list soma-grupos)))
                         (setf soma-grupos 0)
                         );dotimes

                ;;By adding the last group that was not considered in the above dotimes, to the list of similar groups.
                (setf grupos-similares  (append grupos-similares (list (list (first (last grupos-identicos)) (list ) (length (second (first (last grupos-identicos))))))))

                ;;Ordering similar groups by size, from smallest to largest.
                (setf grupos-similares (sort grupos-similares #'< :key #'(lambda (x) (third (first x)))))

                ;;Removing groups with the same amount of genes and assigning the result to the group weights list.
                (setf pesos-grupos (remove-duplicates grupos-similares :key #'(lambda (x) (third (first x)))))

                ;;Sweeping the groups formed to create edges between  identical  and  similar proteins.
                (dolist (g grupos-similares)

                        ;;Picking up the group index of size equal to group g in the group weights list.
                        (setf dividendo (1+ (position-if #'(lambda (x) (= (third (first g)) (third (first x)))) pesos-grupos)))

                        ;;Picking up the index of the most populous group, that is, the last group.
                        (setf divisor (length pesos-grupos))

                        ;;Calculating the weight of group g based on their position, and the position of the last group (most populous).
                        (setf peso-grupo (* (/ dividendo divisor) percentage-pp))

                        ;;Creating edges between proteins with identical profiles.
                        (dotimes (i (1- (third (first g))))
                                 (loop for j from (1+ i) to (1- (third (first g)))
                                   do (progn
                                       (push (make-ppi-struct
                                              :genea (first (elt (second (first g)) i))
                                              :geneb (first (elt (second (first g)) j))
                                              :weight peso-grupo
                                              :position (second (elt (second (first g)) i)))
                                             ppi)
                                       );progn-do
                                   );loop j
                                 );dotimes i

                        ;;Sweeping protein groups with profiles similar to the profile of the protein group with identical profiles.
                        (dolist (grupo-similar (second g))
                                ;Calculating the  difference  between profiles.
                                (setf diferenca-perfis (comparar-perfis (first (first g)) (first grupo-similar)))
                                ;;Sweeping the proteins of a group with a similar profile.
                                (dolist (ps (second grupo-similar))
                                        ;;Sweeping proteins   of identical profiles of a group g.
                                        (dolist (pg (second (first g)))
                                                ;;Creating  an interaction.
                                                (push (make-ppi-struct
                                                 						:genea (if (< (second pg) (second ps)) (first pg) (first ps) )
                                                 						:geneb (if (> (second pg) (second ps)) (first pg) (first ps) )
                                                 						:weight (- peso-grupo (* peso-grupo (/ (/ (* diferenca-perfis 100) (length *genomes-files*)) 100)))
                                                 						:position (if (< (second pg) (second ps))  (second pg) (second ps) ) )
                                          					       ppi);push
                                                );dolist pg
                                        );dolist ps
                                );similar-group dolist
                        );similar groups dolist
                );progn
               );condition and action 2
              );cond
            );unless

          ;;If any PPI has  beenpredicted, the PPI identifier is  changed  to  true.
          (if (> (length ppi) 0)
            (setf *ppi-identified-pp* t)
            )

          ;;Recording the  PPI of genome k in the database.
          (save-db (concatenate 'string diretorio-banco-ppi k ".db") ppi)

          ;;Reset the variables.
          (setf grupos-identicos nil)
          (setf grupos-similares nil)
          (setf ppi nil)

          (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *perfil-filogenetico*)) 100)) laco)))
                   (format *query-io* "=")
                   (force-output *query-io*)
                   )
          (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *perfil-filogenetico*)) 100)) laco))))
          (incf genomenumber)

          );progn-loop
      );loop hash table
    (format t "]~%")
    ))
;;------------------------------------------------------------------------------

;;Function to predict  PPI  by gene fusion
(defun rosetta-stone(aadifflimit checkpointminlimit percentage-rs)

  (defparameter *fusoes* (make-hash-table :test #'equalp))
  (defstruct fusion ppi rosetta-stone)

  (let ((fusoes (list))
        (gene-a)
        (gene-b)
        (fusao)
        (posicao-p)
        (qtd-genes)
        (ppi-fusoes (make-hash-table :test #'equalp))
        (laco 0)
        (genomenumber 0)
        )

    (declare (type fixnum aadifflimit))
    (declare (type fixnum checkpointminlimit))
    (declare (type single-float percentage-rs))
    (declare (type fixnum laco))
    (declare (type fixnum genomenumber))

    (format t "~%Predicting ppi by gene-fusion;~%")
    (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
    (format *query-io* "[")
    (force-output *query-io*)

    (loop for g being the hash-keys in *genomas* using (hash-value proteins)
      do (progn
          (setf qtd-genes (length proteins))
          (dolist (p proteins)
                  (setf posicao-p (position p proteins :test #'equalp))
                  (cond
                    ((< posicao-p (- qtd-genes 2))
                     (progn
                      (setf gene-a p)
                      (setf gene-b (elt proteins (+ posicao-p 1)))
                      (setf fusao (list "fusao" (mapcar #'+ (second gene-a) (second gene-b))))
                      (setf fusoes (append fusoes (list (list fusao (list g (first gene-a) (first gene-b) posicao-p) (list)))))
                      #|(setf gene-b (elt proteins (+ posicao-p 2)))
                      (setf fusao (list "fusao" (mapcar #'+ (second gene-a) (second gene-b))))
                      (setf fusoes (append fusoes (list (list fusao (list g (first gene-a) (first gene-b) posicao-p) (list)))))|#
                      )
                     );end condition 1 and action 1.
                    ((= posicao-p (- qtd-genes 2))
                     (progn
                      (setf gene-a p)
                      (setf gene-b (elt proteins (1+ posicao-p)))
                      (setf fusao (list "fusao" (mapcar #'+ (second gene-a) (second gene-b))))
                      (setf fusoes (append fusoes (list (list fusao (list g (first gene-a) (first gene-b) posicao-p) (list)))))
                      #|(setf gene-b (elt proteins 0))
                      (setf fusao (list "fusao" (mapcar #'+ (second gene-a) (second gene-b))))
                      (setf fusoes (append fusoes (list (list fusao (list g (first gene-a) (first gene-b) posicao-p) (list)))))|#
                      )
                     );end condition 2 and action 2.
                    #|((= posicao-p (- qtd-genes 1))
                       (progn
                        (setf gene-a p)
                        (setf gene-b (elt proteins 0))
                        (setf fusao (list "fusao" (mapcar #'+ (second gene-a) (second gene-b))))
                        (setf fusoes (append fusoes (list (list fusao (list g (first gene-a) (first gene-b) posicao-p) (list)))))
                        (setf gene-b (elt proteins 1))
                        (setf fusao (list "fusao" (mapcar #'+ (second gene-a) (second gene-b))))
                        (setf fusoes (append fusoes (list (list fusao (list g (first gene-a) (first gene-b) posicao-p) (list)))))
                        )
                       );condition end 3 and action 3.|#
                    );cond
                  (dolist (f fusoes)
                          (unless (or (equalp (first p) (second (second f))) (equalp (first p) (third (second f))))
                            (if (similar-test p (first f) aadifflimit checkpointminlimit)
                              (progn
                               (setf (third f) (append (third f) (list (list (first p) g))))
                               );
                              );if similar-test
                            );unless
                          );dolist mergers
                  );dolist  proteins
          (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *genomas*)) 100)) laco)))
                   (format *query-io* "=")
                   (force-output *query-io*)
                   )
          (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *genomas*)) 100)) laco))))
          (incf genomenumber)
          );progn
      );loop hash table *genomas*

    (dolist (f fusoes)
            (unless (eql (third f) nil)
              (setf *ppi-identified-gf* t)
              (setf (gethash (first (second f)) ppi-fusoes)
                    (append (gethash (first (second f)) ppi-fusoes)
                            (list
                             (make-ppi-struct
                              :genea (second (second f))
                              :geneb (third (second f))
                              :weight (* 1 percentage-rs)
                              :position (fourth (second f))
                              )
                             );list
                            );append
                    );setf
              (setf (gethash (first (second f)) *fusoes*)
                    (append
                     (gethash (first (second f)) *fusoes*)
                     (list (make-fusion :ppi (list (second (second f)) '-- (third (second f)))
                                        :rosetta-stone (third f)
                                        );make-fusion
                           );list
                     );append
                    );setf mergers
              );unlles
            );dolist pangenome
    (format t "]~%")
    (return-from rosetta-stone ppi-fusoes)
    );let
  );defun
;-------------------------------------------------------------------------------

;;Reports
(defun relatorio-perfil-filogenetico (workingdir caminho percentage-pp ppdifftolerated pphistofilter
                                                 ppaadifflimit ppaacheckminlimit aadifflimit
                                                 aacheckminlimit
                                                 )

  (let ((peso-grupo)(dividendo)(divisor)(total-arestas)(total-ppi-profile)
                    (diferenca-perfis)(perfil-cont 0)(similar-cont 0)
                    (laco 0)(genomenumber 0)(grupos))

    (ensure-directories-exist caminho)
    (with-open-file (str caminho
                         :direction :output
                         :if-exists :supersede
                         :if-does-not-exist :create)

      (format str "~%")

      (if (= ppdifftolerated 0)
        (progn
         (format str "---------------------Report of complete ppi prediction by phylogenetic profile---------------------~%")
         (if (string= pphistofilter "Y")
           (cond
             ((and (= ppaadifflimit 0) (= ppaacheckminlimit 26)) (format str "Minimum identity percentage: 100%~%"))
             ((and (= ppaadifflimit 0) (= ppaacheckminlimit 25)) (format str "Minimum identity percentage: 100%~%"))
             ((and (= ppaadifflimit 0) (= ppaacheckminlimit 24)) (format str "Minimum identity percentage: 97,96%~%"))
             ((and (= ppaadifflimit 0) (= ppaacheckminlimit 23)) (format str "Minimum identity percentage: 96,94%~%"))
             ((and (= ppaadifflimit 0) (= ppaacheckminlimit 22)) (format str "Minimum identity percentage: 96,94%~%"))
             ((and (= ppaadifflimit 0) (= ppaacheckminlimit 21)) (format str "Minimum identity percentage: 94,68%~%"))
             ((and (= ppaadifflimit 0) (= ppaacheckminlimit 20)) (format str "Minimum identity percentage: 91,75%~%"))
             ((and (= ppaadifflimit 0) (= ppaacheckminlimit 19)) (format str "Minimum identity percentage: 85,57%~%"))
             ((and (= ppaadifflimit 0) (= ppaacheckminlimit 18)) (format str "Minimum identity percentage: 50%~%"))
             ((and (= ppaadifflimit 1) (= ppaacheckminlimit 26)) (format str "Minimum identity percentage: 97,87%~%"))
             ((and (= ppaadifflimit 1) (= ppaacheckminlimit 25)) (format str "Minimum identity percentage: 92,55%~%"))
             )
           (cond
             ((and (= aadifflimit 0) (= aacheckminlimit 26)) (format str "Minimum identity percentage: 100%~%"))
             ((and (= aadifflimit 0) (= aacheckminlimit 25)) (format str "Minimum identity percentage: 100%~%"))
             ((and (= aadifflimit 0) (= aacheckminlimit 24)) (format str "Minimum identity percentage: 97,96%~%"))
             ((and (= aadifflimit 0) (= aacheckminlimit 23)) (format str "Minimum identity percentage: 96,94%~%"))
             ((and (= aadifflimit 0) (= aacheckminlimit 22)) (format str "Minimum identity percentage: 96,94%~%"))
             ((and (= aadifflimit 0) (= aacheckminlimit 21)) (format str "Minimum identity percentage: 94,68%~%"))
             ((and (= aadifflimit 0) (= aacheckminlimit 20)) (format str "Minimum identity percentage: 91,75%~%"))
             ((and (= aadifflimit 0) (= aacheckminlimit 19)) (format str "Minimum identity percentage: 85,57%~%"))
             ((and (= aadifflimit 0) (= aacheckminlimit 18)) (format str "Minimum identity percentage: 50%~%"))
             ((and (= aadifflimit 1) (= aacheckminlimit 26)) (format str "Minimum identity percentage: 97,87%~%"))
             ((and (= aadifflimit 1) (= aacheckminlimit 25)) (format str "Minimum identity percentage: 92,55%~%"))
             )
           );if pphistofilter
         (format str "---------------------------------------------------------------------------------------------------~%")

         (format t "~%Writing phylogenetic profiles report;~%")
         (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
         (format *query-io* "[")
         (force-output *query-io*)

         (dolist (genoma *genomes-files*)
                 (if (probe-file (concatenate 'string workingdir "/database/phylogenetic-profile-groups/" genoma ".db"))
                   (progn
                    (setf grupos (load-db (concatenate 'string workingdir "/database/phylogenetic-profile-groups/" genoma ".db")))
                    (setf total-arestas (loop for i in grupos
                                          summing (/ (* (length (second i)) (- (length (second i)) 1)) 2) into total
                                          finally (return total)))

                    (format str "~%Genome: ~a ~%Total profiles: ~a ~%Total ppi: ~a~%" genoma (length grupos) total-arestas)
                    (format str "---------------------------------------------------------------------------------------------------~%")
                    (dolist (g grupos)

                            ;;Incrementing or profile counter
                            (setf perfil-cont (incf perfil-cont))
                            ;;Picking up the group index of size equal to group g in the group weights list.
                            (setf dividendo (1+ (position-if #'(lambda (x) (= (third g) (third x))) (gethash genoma *pesos-grupos*))))
                            ;;Picking up the index of the most populous group, that is, the last group.
                            (setf divisor (length (gethash genoma *pesos-grupos*)))
                            ;;Calculating the weight of group g based on their position, and the position of the last group (most populous).
                            (setf peso-grupo (* (/ dividendo divisor) percentage-pp))
                            ;;Calculating the total edges created for profile g
                            (setf total-ppi-profile (/ (* (length (second g)) (- (length (second g)) 1)) 2))

                            (format str "Profile ~6a => Total genes: ~7a| Total ppi: ~9a| Weight: ~,3f | Number of genomes: ~a;~%"
                                    perfil-cont
                                    (length (second g))
                                    total-ppi-profile
                                    (if (= total-ppi-profile 0) (- 0 0) (- peso-grupo 0))
                                    (length(first g))
                                    )
                            );dotimes groups
                    (setf perfil-cont 0)
                    (format str "---------------------------------------------------------------------------------------------------~%")
                    )
                   );if
                 (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (length *genomes-files*)) 100)) laco)))
                          (format *query-io* "=")
                          (force-output *query-io*)
                          )
                 (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (length *genomes-files*)) 100)) laco))))
                 (incf genomenumber)
                 );dolist
         (format t "]~%")
         );progn

        (progn
         (format str "-----------------------Report of complete ppi prediction by phylogenetic profile-----------------------~%")
         (if (string= pphistofilter "Y")
           (cond
             ((and (= ppaadifflimit 0) (= ppaacheckminlimit 26)) (format str "Minimum identity percentage: 100%~%"))
             ((and (= ppaadifflimit 0) (= ppaacheckminlimit 25)) (format str "Minimum identity percentage: 100%~%"))
             ((and (= ppaadifflimit 0) (= ppaacheckminlimit 24)) (format str "Minimum identity percentage: 97,96%~%"))
             ((and (= ppaadifflimit 0) (= ppaacheckminlimit 23)) (format str "Minimum identity percentage: 96,94%~%"))
             ((and (= ppaadifflimit 0) (= ppaacheckminlimit 22)) (format str "Minimum identity percentage: 96,94%~%"))
             ((and (= ppaadifflimit 0) (= ppaacheckminlimit 21)) (format str "Minimum identity percentage: 94,68%~%"))
             ((and (= ppaadifflimit 0) (= ppaacheckminlimit 20)) (format str "Minimum identity percentage: 91,75%~%"))
             ((and (= ppaadifflimit 0) (= ppaacheckminlimit 19)) (format str "Minimum identity percentage: 85,57%~%"))
             ((and (= ppaadifflimit 0) (= ppaacheckminlimit 18)) (format str "Minimum identity percentage: 50%~%"))
             ((and (= ppaadifflimit 1) (= ppaacheckminlimit 26)) (format str "Minimum identity percentage: 97,87%~%"))
             ((and (= ppaadifflimit 1) (= ppaacheckminlimit 25)) (format str "Minimum identity percentage: 92,55%~%"))
             )
           (cond
             ((and (= aadifflimit 0) (= aacheckminlimit 26)) (format str "Minimum identity percentage: 100%~%"))
             ((and (= aadifflimit 0) (= aacheckminlimit 25)) (format str "Minimum identity percentage: 100%~%"))
             ((and (= aadifflimit 0) (= aacheckminlimit 24)) (format str "Minimum identity percentage: 97,96%~%"))
             ((and (= aadifflimit 0) (= aacheckminlimit 23)) (format str "Minimum identity percentage: 96,94%~%"))
             ((and (= aadifflimit 0) (= aacheckminlimit 22)) (format str "Minimum identity percentage: 96,94%~%"))
             ((and (= aadifflimit 0) (= aacheckminlimit 21)) (format str "Minimum identity percentage: 94,68%~%"))
             ((and (= aadifflimit 0) (= aacheckminlimit 20)) (format str "Minimum identity percentage: 91,75%~%"))
             ((and (= aadifflimit 0) (= aacheckminlimit 19)) (format str "Minimum identity percentage: 85,57%~%"))
             ((and (= aadifflimit 0) (= aacheckminlimit 18)) (format str "Minimum identity percentage: 50%~%"))
             ((and (= aadifflimit 1) (= aacheckminlimit 26)) (format str "Minimum identity percentage: 97,87%~%"))
             ((and (= aadifflimit 1) (= aacheckminlimit 25)) (format str "Minimum identity percentage: 92,55%~%"))
             )
           );if pphistofilter
         (format str "-------------------------------------------------------------------------------------------------------~%")

         (format t "~%Writing phylogenetic profiles report;~%")
         (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
         (format *query-io* "[")
         (force-output *query-io*)

         (dolist (genoma *genomes-files*)
                 (if (probe-file (concatenate 'string workingdir "/database/phylogenetic-profile-groups/" genoma ".db"))
                   (progn
                    (setf grupos (load-db (concatenate 'string workingdir "/database/phylogenetic-profile-groups/" genoma ".db")))
                    (setf total-arestas (loop for g in grupos
                                          summing (+
                                                   ;;Calculating how many pairs of genes with identical profiles will be generated.
                                                   (/ (* (length (second (first g))) (- (length (second (first g))) 1)) 2)
                                                   ;;Multiplying the number of genes of identical profiles by the number of genes of similar profiles.
                                                   (* (length (second (first g)))
                                                      ;;Adding up the amount of genes from similar profiles.
                                                      (loop for gs in (second g) summing (length (second gs)) into total finally (return total))
                                                      );*
                                                   );+
                                          into total
                                          finally (return total))
                          );setf full-edges

                    (format str "~%Genome: ~a ~%Total profiles: ~a ~%Total ppi: ~a~%" genoma (length grupos) total-arestas)
                    (format str "-------------------------------------------------------------------------------------------------------~%")
                    (dolist (g grupos)

                            ;;Incrementing or profile counter
                            (setf perfil-cont (incf perfil-cont))
                            ;;Picking up the group index of size equal to group g in the group weights list.
                            (setf dividendo (1+ (position-if #'(lambda (x) (= (third (first g)) (third (first x)))) (gethash genoma *pesos-grupos*))))
                            ;;Picking up the index of the most populous group, that is, the last group.
                            (setf divisor (length (gethash genoma *pesos-grupos*)))
                            ;;Calculating the weight of group g based on their position, and the position of the last group (most populous).
                            (setf peso-grupo (* (/ dividendo divisor) percentage-pp))
                            ;;Calculating the total edges created for profile g
                            (setf total-ppi-profile (/ (* (length (second (first g))) (- (length (second (first g))) 1)) 2))

                            (format str "Profile ~9a => Total genes: ~7a| Total ppi: ~9a| Weight: ~,3f | Number of genomes: ~a;~%"
                                    perfil-cont
                                    (length (second (first g)))
                                    total-ppi-profile
                                    (if (= total-ppi-profile 0) (- 0 0) (- peso-grupo 0))
                                    (length(first (first g)))
                                    );format

                            (dolist (sg (second g))
                                    ;;Incrementing or profile counter
                                    (setf similar-cont (incf similar-cont))
                                    ;Calculating the  difference  between profiles.
                                    (setf diferenca-perfis (comparar-perfis (first (first g)) (first sg)))
                                    (format str "Similar ~9a => Total genes: ~7a| Total ppi: ~9a| Weight: ~,3f | Number of genomes: ~a;~%"
                                            ;profile-cont
                                            (concatenate 'string (write-to-string perfil-cont) "." (write-to-string similar-cont))
                                            ;similar-cont
                                            (length (second sg))
                                            ;;Multiplying the number of g-profile genes   by the total number of genes in the  sg profile  (profile i similar to g)
                                            ;;It is necessarytodo this multiplication, because each gene of profile g interacts with each gene of the similar profile  sg
                                            (* (length (second (first g))) (length (second sg)))
                                            (- peso-grupo (* peso-grupo (/ (/ (* diferenca-perfis 100) (length *genomes-files*)) 100)))
                                            (length (first sg))
                                            );format
                                    );dolist sg
                            (setf similar-cont 0)
                            );dotimes groups
                    (setf perfil-cont 0)
                    (format str "------------------------------------------------------------------------------------------------------~%")
                    );progn
                   );if
                 (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (length *genomes-files*)) 100)) laco)))
                          (format *query-io* "=")
                          (force-output *query-io*)
                          )
                 (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (length *genomes-files*)) 100)) laco))))
                 (incf genomenumber)
                 );dolist
         (format t "]~%")
         );progn
        );if
      (format str "~%")
      );with-open-file
    );let
  );defun

;;Function to create a  pangenome report
(defun relatorio-pangenoma-1(caminho aadifflimit aacheckminlimit)
  (let ((protein)
        (laco 0)
        (genomenumber 0)
        (keys-hash-table (list))
        (pan-genoma (make-hash-table :test #'equalp))
        )
    (ensure-directories-exist caminho)
    (with-open-file (str caminho
                         :direction :output
                         :if-exists :supersede
                         :if-does-not-exist :create)

      (format str "~%")

      (format str "-----------------------------------Pan genoma------------------------------------~%")
      (cond
        ((and (= aadifflimit 0) (= aacheckminlimit 26)) (format str "Minimum identity percentage: 100%~%"))
        ((and (= aadifflimit 0) (= aacheckminlimit 25)) (format str "Minimum identity percentage: 100%~%"))
        ((and (= aadifflimit 0) (= aacheckminlimit 24)) (format str "Minimum identity percentage: 97,96%~%"))
        ((and (= aadifflimit 0) (= aacheckminlimit 23)) (format str "Minimum identity percentage: 96,94%~%"))
        ((and (= aadifflimit 0) (= aacheckminlimit 22)) (format str "Minimum identity percentage: 96,94%~%"))
        ((and (= aadifflimit 0) (= aacheckminlimit 21)) (format str "Minimum identity percentage: 94,68%~%"))
        ((and (= aadifflimit 0) (= aacheckminlimit 20)) (format str "Minimum identity percentage: 91,75%~%"))
        ((and (= aadifflimit 0) (= aacheckminlimit 19)) (format str "Minimum identity percentage: 85,57%~%"))
        ((and (= aadifflimit 0) (= aacheckminlimit 18)) (format str "Minimum identity percentage: 50%~%"))
        ((and (= aadifflimit 1) (= aacheckminlimit 26)) (format str "Minimum identity percentage: 97,87%~%"))
        ((and (= aadifflimit 1) (= aacheckminlimit 25)) (format str "Minimum identity percentage: 92,55%~%"))
        )
      (format str "---------------------------------------------------------------------------------")

      (format t "~%Writing pan-genome format 1;~%")
      (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
      (format *query-io* "[")
      (force-output *query-io*)

      ;;  lisp  writes  hash tables  to the database in reverse  order, sothis block of code
      ;;sort the hash  table  retrieved from the database
      (loop for k being each hash-key of *pan-genoma*
        do (push k keys-hash-table))
      (dolist (k keys-hash-table)
              (setf (gethash k pan-genoma) (gethash k *pan-genoma*))
              )
      ;;Taking unnecessary data from memory
      (setf keys-hash-table nil)
      ;;End  ordering

      (loop for k being the hash-keys in pan-genoma using (hash-value v)
        do (progn
            (format str "~%~%Protein: ~a | Genome: ~a~%" (subseq k (1+ (position #\. k))) (first (proteina-localidade v)))
            (format str "~%Similar proteins:~%~%")
            (dolist (i (proteina-similares v))
	      (if (and i (first i) (second i))
		  (progn	      
                    (setf protein (first (elt (gethash (first i) *genomas*)(second i))))
                    (format str "Protein: ~a | Genome: ~a;~%"
                            (subseq protein (1+ (position #\. protein)))
                            (first i))
                    )
		)
	      )
            (format str "---------------------------------------------------------------------------------")

            (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count pan-genoma)) 100)) laco)))
                     (format *query-io* "=")
                     (force-output *query-io*)
                     )
            (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count pan-genoma)) 100)) laco))))
            (incf genomenumber)

            );progn
        );loop hash table
      (format t "]~%")
      (format str "~%")
      );with-open-file
    )
  );----------------------------------------------------------------------------

(defun relatorio-pangenoma-2 (caminho)
  (let ((similar)
        (laco 0)
        (genomenumber 0)
        (keys-hash-table (list))
        (pan-genoma (make-hash-table :test #'equalp))
        )
    (ensure-directories-exist caminho)
    (with-open-file (str caminho
                         :direction :output
                         :if-exists :supersede
                         :if-does-not-exist :create)

      (format str "#Columns per row: Genome, Gene, Similar list~%~%")

      (format t "~%Writing pan-genome format 2;~%")
      (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
      (format *query-io* "[")
      (force-output *query-io*)

      ;;  Lisp  writes  hash tables  to the database in reverse  order, sothis block of code
      ;;sort the hash  table  retrieved from the database
      (loop for k being each hash-key of *pan-genoma*
        do (push k keys-hash-table))
      (dolist (k keys-hash-table)
              (setf (gethash k pan-genoma) (gethash k *pan-genoma*))
              )
      ;;Taking unnecessary data from memory
      (setf keys-hash-table nil)
      ;;End  ordering

      (loop for k being the hash-keys in pan-genoma using (hash-value v)
        do (progn
            (format str "~a, ~a," (first (proteina-localidade v)) (subseq k (1+ (position #\. k))))
            (dolist (i (proteina-similares v))
	      (if (and i (first i) (second i))
                  (progn	      
                    (setf similar (first (elt (gethash (first i) *genomas*)(second i))))
                    (format str " ~a" (subseq similar (1+ (position #\. similar))))
                    )
		)
	      )
            (format str ";~%")

            (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count pan-genoma)) 100)) laco)))
                     (format *query-io* "=")
                     (force-output *query-io*)
                     )
            (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count pan-genoma)) 100)) laco))))
            (incf genomenumber)

            );progn
        );loop hash table
      (format t "]~%")
      (format str "~%")
      );with-open-file
    )
  )
;-------------------------------------------------------------------------------

;;Function to create a report that shows the amount of genes stored and
;;not kept in the vicinity of pangenome  genes
(defun relatorio-conservao-genica(workingdir caminho aadifflimit aacheckminlimit expansion-type w1 window-size)

  (let ((laco 0)
        (genomenumber 0)
        (keys-hash-table (list))
        (tabela-hash-vizinhanca-genica (make-hash-table :test #'equalp))
        )

    (ensure-directories-exist caminho)
    (with-open-file (str caminho
                         :direction :output
                         :if-exists :supersede
                         :if-does-not-exist :create)

      (format str "~%")

      (format str "#-----------------------------------------------Gene neighborhood conservation report----------------------------------------------~%")
      (format str "~%")
      (cond
        ((and (= aadifflimit 0) (= aacheckminlimit 26)) (format str "#Minimum identity percentage: 100%,"))
        ((and (= aadifflimit 0) (= aacheckminlimit 25)) (format str "#Minimum identity percentage: 100%,"))
        ((and (= aadifflimit 0) (= aacheckminlimit 24)) (format str "#Minimum identity percentage: 97,96%,"))
        ((and (= aadifflimit 0) (= aacheckminlimit 23)) (format str "#Minimum identity percentage: 96,94%,"))
        ((and (= aadifflimit 0) (= aacheckminlimit 22)) (format str "#Minimum identity percentage: 96,94%,"))
        ((and (= aadifflimit 0) (= aacheckminlimit 21)) (format str "#Minimum identity percentage: 94,68%,"))
        ((and (= aadifflimit 0) (= aacheckminlimit 20)) (format str "#Minimum identity percentage: 91,75%,"))
        ((and (= aadifflimit 0) (= aacheckminlimit 19)) (format str "#Minimum identity percentage: 85,57%,"))
        ((and (= aadifflimit 0) (= aacheckminlimit 18)) (format str "#Minimum identity percentage: 50%,"))
        ((and (= aadifflimit 1) (= aacheckminlimit 26)) (format str "#Minimum identity percentage: 97,87%,"))
        ((and (= aadifflimit 1) (= aacheckminlimit 25)) (format str "#Minimum identity percentage: 92,55%,"))
        )
      (if (string= expansion-type "fixed")
        (format str " Expansion type: fixed, Window size: ~a~%" w1)
        (format str " Expansion type: dynamic, Window size: ~a~%" window-size)
        )
      (format str "~%")
      (format str "#----------------------------------------------------------------------------------------------------------------------------------~%~%")

      (format str "#Columns per row: Genome, Gene, Number of genes conserved, Number of genes not conserved, Number of genomes of each conserved gene;~%~%")

      (format t "~%Writing Gene neighborhood conservation report;~%")
      (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
      (format *query-io* "[")
      (force-output *query-io*)

      ;;  Lisp  writes  hash tables  to the database in reverse  order, sothis block of code
      ;;sort the hash  table  retrieved from the database
      (setf *relatorio-vizinhanca-genica* (load-db (concatenate 'string workingdir "/database/gene-neighborhood-report/gene-neighborhood-report.db")))
      (loop for k being each hash-key of *relatorio-vizinhanca-genica*
        do (push k keys-hash-table))
      (dolist (k keys-hash-table)
              (setf (gethash k tabela-hash-vizinhanca-genica) (gethash k *relatorio-vizinhanca-genica*))
              )
      ;;Taking data from memory
      (setf *relatorio-vizinhanca-genica* nil)
      (setf keys-hash-table nil)
      ;;End  ordering

      (loop for k being the hash-keys in tabela-hash-vizinhanca-genica using (hash-value v)
        do (progn
            (format str "~a, ~a, ~a, ~a,"
                    (first (expansao-localidade v))
                    (subseq k (1+ (position #\. k)))
                    (count-if #'(lambda (x) (> (first x) 1)) (expansao-conservacoes v))
                    (count-if #'(lambda (x) (= (first x) 1)) (expansao-conservacoes v))
                    )
            (dolist (i (expansao-conservacoes v))
                    (unless (= (first i) 1)
                      (format str " ~a" (second i))
                      )
                    )
            (unless (> (count-if #'(lambda (x) (> (first x) 1)) (expansao-conservacoes v)) 0)
              (format str " 0")
              )
            (format str ";~%")

            (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count tabela-hash-vizinhanca-genica)) 100)) laco)))
                     (format *query-io* "=")
                     (force-output *query-io*)
                     )
            (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count tabela-hash-vizinhanca-genica)) 100)) laco))))
            (incf genomenumber)

            );progn
        );loop hash table
      ;;Taking unnecessary data from memory
      (setf tabela-hash-vizinhanca-genica nil)
      (format t "]~%")
      (format str "~%")
      );with-open-file
    );let
  );----------------------------------------------------------------------------

;;Function to create the report of the amount of interactions per feature in each genome.
(defun ppi-number-report (workingdir ppi-gene-fusion gene-fusion ppcn ppcomplete pplimit pptrim ppthreshold ppdeletegroup ppdeleteprofile)

  (let ((qtd-ppi-cn)
        (qtd-ppi-pp)
        (qtd-ppi-gf)
        (laco 0)
        (genomenumber 0)
        )

    (if (= (length *genomes-files*) 1)
      (setf ppcn "N"
            ppcomplete "N"
            pplimit "N"
            pptrim "N"
            ppthreshold "N"
            ppdeletegroup "N"
            ppdeleteprofile "N"
            )
      )

    (ensure-directories-exist (concatenate 'string workingdir "/ppi-report/report.txt"))
    (with-open-file (str (concatenate 'string workingdir "/ppi-report/report.txt")
                         :direction :output
                         :if-exists :supersede
                         :if-does-not-exist :create)

      (format str "~%")

      (format str "---------------------Amount of predicted interactions---------------------~%~%")
      (format str "--------------------------------------------------------------------------~%")

      (cond
        ((and (or (string= ppcn "Y") (string= ppcomplete "Y") (string= pplimit "Y") (string= pptrim "Y")
                  (string= ppthreshold "Y") (string= ppdeletegroup "Y") (string= ppdeleteprofile "Y"))
              (string= gene-fusion "N"))

         (format t "~%Writing ppi report;~%")
         (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
         (format *query-io* "[")
         (force-output *query-io*)

         (dolist (g *genomes-files*)
                 (if (probe-file (concatenate 'string workingdir "/database/ppi-conserved-neighborhood/" g ".db"))
                   (setf qtd-ppi-cn (length (load-db (concatenate 'string workingdir "/database/ppi-conserved-neighborhood/" g ".db"))))
                   (setf qtd-ppi-cn 0)
                   )
                 (if (probe-file (concatenate 'string workingdir "/database/ppi-phylogenetic-profiles/" g ".db"))
                   (setf qtd-ppi-pp (length (load-db (concatenate 'string workingdir "/database/ppi-phylogenetic-profiles/" g ".db"))))
                   (setf qtd-ppi-pp 0)
                   )
                 (format str "Genome name: ~22t~a~%" g)
                 (format str "Number of ppi by cn: ~22t~a~%" qtd-ppi-cn)
                 (format str "Number of ppi by pp: ~22t~a~%" qtd-ppi-pp)
                 (format str "--------------------------------------------------------------------------~%")

                 (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (length *genomes-files*)) 100)) laco)))
                          (format *query-io* "=")
                          (force-output *query-io*)
                          )
                 (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (length *genomes-files*)) 100)) laco))))
                 (incf genomenumber)

                 );dolist
         (format t "]~%")
         )
        ((and (or (string= ppcn "Y") (string= ppcomplete "Y") (string= pplimit "Y") (string= pptrim "Y")
                  (string= ppthreshold "Y") (string= ppdeletegroup "Y") (string= ppdeleteprofile "Y"))
              (string= gene-fusion "Y"))

         (format t "~%Writing ppi report;~%")
         (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
         (format *query-io* "[")
         (force-output *query-io*)

         (dolist (g *genomes-files*)
                 (if (probe-file (concatenate 'string workingdir "/database/ppi-conserved-neighborhood/" g ".db"))
                   (setf qtd-ppi-cn (length (load-db (concatenate 'string workingdir "/database/ppi-conserved-neighborhood/" g ".db"))))
                   (setf qtd-ppi-cn 0)
                   )
                 (if (probe-file (concatenate 'string workingdir "/database/ppi-phylogenetic-profiles/" g ".db"))
                   (setf qtd-ppi-pp (length (load-db (concatenate 'string workingdir "/database/ppi-phylogenetic-profiles/" g ".db"))))
                   (setf qtd-ppi-pp 0)
                   )
                 (setf qtd-ppi-gf (length (gethash g ppi-gene-fusion)))
                 (format str "Genome name: ~22t~a~%" g)
                 (format str "Number of ppi by cn: ~22t~a~%" qtd-ppi-cn)
                 (format str "Number of ppi by pp: ~22t~a~%" qtd-ppi-pp)
                 (format str "Number of ppi by gf: ~22t~a~%" qtd-ppi-gf)
                 (format str "--------------------------------------------------------------------------~%")

                 (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (length *genomes-files*)) 100)) laco)))
                          (format *query-io* "=")
                          (force-output *query-io*)
                          )
                 (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (length *genomes-files*)) 100)) laco))))
                 (incf genomenumber)

                 );dolist
         (format t "]~%")
         )
        ((and (and (string= ppcn "N") (string= ppcomplete "N") (string= pplimit "N") (string= pptrim "N")
                   (string= ppthreshold "N") (string= ppdeletegroup "N") (string= ppdeleteprofile "N"))
              (string= gene-fusion "Y"))

         (format t "~%Writing ppi report;~%")
         (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
         (format *query-io* "[")
         (force-output *query-io*)

         (dolist (g *genomes-files*)
                 (if (probe-file (concatenate 'string workingdir "/database/ppi-conserved-neighborhood/" g ".db"))
                   (setf qtd-ppi-cn (length (load-db (concatenate 'string workingdir "/database/ppi-conserved-neighborhood/" g ".db"))))
                   (setf qtd-ppi-cn 0)
                   )
                 (setf qtd-ppi-gf (length (gethash g ppi-gene-fusion)))
                 (format str "Genome name: ~22t~a~%" g)
                 (format str "Number of ppi by cn: ~22t~a~%" qtd-ppi-cn)
                 (format str "Number of ppi by gf: ~22t~a~%" qtd-ppi-gf)
                 (format str "--------------------------------------------------------------------------~%")

                 (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (length *genomes-files*)) 100)) laco)))
                          (format *query-io* "=")
                          (force-output *query-io*)
                          )
                 (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (length *genomes-files*)) 100)) laco))))
                 (incf genomenumber)

                 );dolist
         (format t "]~%")
         )
        ((and (and (string= ppcn "N") (string= ppcomplete "N") (string= pplimit "N") (string= pptrim "N")
                   (string= ppthreshold "N") (string= ppdeletegroup "N") (string= ppdeleteprofile "N"))
              (string= gene-fusion "N"))

         (format t "~%Writing ppi report;~%")
         (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
         (format *query-io* "[")
         (force-output *query-io*)

         (dolist (g *genomes-files*)
                 (if (probe-file (concatenate 'string workingdir "/database/ppi-conserved-neighborhood/" g ".db"))
                   (setf qtd-ppi-cn (length (load-db (concatenate 'string workingdir "/database/ppi-conserved-neighborhood/" g ".db"))))
                   (setf qtd-ppi-cn 0)
                   )
                 (format str "Genome name: ~22t~a~%" g)
                 (format str "Number of ppi by cn: ~22t~a~%" qtd-ppi-cn)
                 (format str "--------------------------------------------------------------------------~%")

                 (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (length *genomes-files*)) 100)) laco)))
                          (format *query-io* "=")
                          (force-output *query-io*)
                          )
                 (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (length *genomes-files*)) 100)) laco))))
                 (incf genomenumber)

                 );dolist
         (format t "]~%")
         )
        );cond
      );with-open-file
    );let
  );defun
;-------------------------------------------------------------------------------

(defun gene-fusion-report (caminho aadifflimit aacheckminlimit)

  (let ((laco 0)(genomenumber 0))
    (ensure-directories-exist caminho)
    (with-open-file (str caminho
                         :direction :output
                         :if-exists :supersede
                         :if-does-not-exist :create)

      (format str "~%")

      (format str "--------------------------------------Gene fusion report---------------------------------------~%")
      (cond
        ((and (= aadifflimit 0) (= aacheckminlimit 26)) (format str "Minimum identity percentage: 100%~%"))
        ((and (= aadifflimit 0) (= aacheckminlimit 25)) (format str "Minimum identity percentage: 100%~%"))
        ((and (= aadifflimit 0) (= aacheckminlimit 24)) (format str "Minimum identity percentage: 97,96%~%"))
        ((and (= aadifflimit 0) (= aacheckminlimit 23)) (format str "Minimum identity percentage: 96,94%~%"))
        ((and (= aadifflimit 0) (= aacheckminlimit 22)) (format str "Minimum identity percentage: 96,94%~%"))
        ((and (= aadifflimit 0) (= aacheckminlimit 21)) (format str "Minimum identity percentage: 94,68%~%"))
        ((and (= aadifflimit 0) (= aacheckminlimit 20)) (format str "Minimum identity percentage: 91,75%~%"))
        ((and (= aadifflimit 0) (= aacheckminlimit 19)) (format str "Minimum identity percentage: 85,57%~%"))
        ((and (= aadifflimit 0) (= aacheckminlimit 18)) (format str "Minimum identity percentage: 50%~%"))
        ((and (= aadifflimit 1) (= aacheckminlimit 26)) (format str "Minimum identity percentage: 97,87%~%"))
        ((and (= aadifflimit 1) (= aacheckminlimit 25)) (format str "Minimum identity percentage: 92,55%~%"))
        )
      (format str "-----------------------------------------------------------------------------------------------")

      (format t "~%Writing gene fusion report;~%")
      (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
      (format *query-io* "[")
      (force-output *query-io*)

      (loop for g being the hash-keys in *fusoes* using (hash-value fusoes)
        do (progn
            (format str "~%~%Genome: ~15t~a~%~%" g)
            (dotimes (i (length fusoes))
                     (format str "Fusion ~a => ~15tPPi: ~a -- ~a~%" (1+ i)
                             (subseq (first(fusion-ppi (elt fusoes i))) (1+ (position #\. (first(fusion-ppi (elt fusoes i))))))
                             (subseq (third(fusion-ppi (elt fusoes i))) (1+ (position #\. (third(fusion-ppi (elt fusoes i))))))
                             );format
                     (dolist (f (fusion-rosetta-stone (elt fusoes i)))
                             (format str "~15tRosetta-stone: ~a in ~a~%"
                                     (subseq (first f) (1+ (position #\. (first f))))
                                     (second f)
                                     );format
                             );dolist
                     (format str "~%")
                     )
            (format str "-----------------------------------------------------------------------------------------------")

            (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *fusoes*)) 100)) laco)))
                     (format *query-io* "=")
                     (force-output *query-io*)
                     )
            (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (hash-table-count *fusoes*)) 100)) laco))))
            (incf genomenumber)

            );progn
        );loop hash table
      (format t "]~%")
      (format str "~%")
      );with-open-file
    );let
  )
;-------------------------------------------------------------------------------

;;Function that merges the  PPI networksof  each feature generating  a single PPI  
(defun genppi (
               workingdir
               gene-fusion
               percentage-gf
               percentage-pp
               ppdifftolerated
               pphistofilter
               ppaadifflimit
               ppaacheckminlimit
               ppcn
               ppcomplete
               pplimit
               ppiterlimit
               pptrim
               trim
               ppthreshold
               threshold
               plusminus
               ppdeletegroup
               grouplimit
               ppdeleteprofile
               profiles
               expansion-type
               window-size
               percentage-cn
               w1
               cw1
               w2
               cw2
               w3
               cw3
               w4
               cw4
               aadifflimit
               aacheckminlimit
               );parameters

  (let (
        (ppi-gene-fusion)
        (perfil-filo-g)
        (perfil-filo-A)
        (perfil-filo-B)
        (lista-ppi-pp (list))
        (cont 0)
        (caminho)
        (ppi-final (list))
        (laco 0)
        (genomenumber 0)
        (diretorio-banco)
        )

    ;;Loading the files and generating the  pangenome
    (histo-genomas (list-directory workingdir))
    (setf lparallel:*kernel* (lparallel:make-kernel (if (> (- (workers) 2) 0 ) (- (workers) 2) 1)))
    (if (string= pphistofilter "N")
      (pan-genoma-1 aadifflimit aacheckminlimit)
      (if (> (length *genomes-files*) 1)
        (pan-genoma-2 aadifflimit aacheckminlimit ppaadifflimit ppaacheckminlimit)
        (pan-genoma-1 aadifflimit aacheckminlimit)
        )
      )
    (lparallel:end-kernel)
    (unless (> (hash-table-count *pan-genoma*) 0)
      (format t "~%No protein in the pan-genome~%")
      )

    (if (or
         (and (= (length *genomes-files*) 1)
              (string= gene-fusion "N")
              )
         (and (and (string= ppcn "N") (string= ppcomplete "N") (string= pplimit "N") (string= pptrim "N")
                   (string= ppthreshold "N") (string= ppdeletegroup "N") (string= ppdeleteprofile "N")
                   )
              (string= gene-fusion "N")
              );and
         );or
      (setf percentage-cn (- 1.0 0))
      )

    ;;Generating  PPI   by  conserved gene neighborhood
    (unless (= (hash-table-count *pan-genoma*) 0)
      (if (string= expansion-type "fixed")
        (conserved-neighbourhood-fixed workingdir percentage-cn w1 cw1 w2 cw2 w3 cw3 w4 cw4 aadifflimit aacheckminlimit)
        (conserved-neighbourhood-dynamic workingdir percentage-cn window-size aadifflimit aacheckminlimit)
        );if
      );unless

    (if (string= gene-fusion "Y")
      (setf ppi-gene-fusion (rosetta-stone aadifflimit aacheckminlimit percentage-gf))
      )

    ;;After turning the stored gene neighborhood function, the  pangenome is
    ;;from memory and writes it to the database.
    (setf diretorio-banco (concatenate 'string workingdir "/database/pan-genome/"))
    (ensure-directories-exist diretorio-banco)
    (unless (eql (list-directory diretorio-banco) nil)
      (dolist (file (list-directory diretorio-banco))
              (unless (eql (pathname-type (pathname file)) nil)
                (delete-file file)
                )
              )
      )
    (save-db (concatenate 'string diretorio-banco "pan-genoma.db")
             *pan-genoma*
             );save-db.
    (setf *pan-genoma* nil)

    ;;Mesmo caso para a tabela hash *genomas*
    (setf diretorio-banco (concatenate 'string workingdir "/database/amino-acid-histogram/"))
    (ensure-directories-exist diretorio-banco)
    (unless (eql (list-directory diretorio-banco) nil)
      (dolist (file (list-directory diretorio-banco))
              (unless (eql (pathname-type (pathname file)) nil)
                (delete-file file)
                )
              )
      )
    (save-db (concatenate 'string diretorio-banco "amino-acid-histogram.db")
             *genomas*
             );save-db.
    (setf *genomas* nil)
    ;---------------------------------------------------------------------------

    ;;Generating  PPI by phylogenetic profile
    (cond
      ;;Start  of Condition and Action 1
      ((and (string= ppcn "Y") (> (length *genomes-files*) 1))

       (setf diretorio-banco (concatenate 'string workingdir "/database/ppi-phylogenetic-profiles/"))
       (ensure-directories-exist diretorio-banco)
       (unless (eql (list-directory diretorio-banco) nil)
         (dolist (file (list-directory diretorio-banco))
                 (unless (eql (pathname-type (pathname file)) nil)
                   (delete-file file)
                   );unless
                 );dolist
         );unless

       (format t "~%Predicting ppi by phylogenetic profiles;~%")
       (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
       (format *query-io* "[")
       (force-output *query-io*)

       ;;Generating  PPI by phylogenetic profile for interactions predicted by  conserved  gene neighborhood.
       (dolist (g *genomes-files*)
               (if (probe-file (concatenate 'string workingdir "/database/ppi-conserved-neighborhood/" g ".db"))
                 (progn
                  (setf perfil-filo-g (gethash g *perfil-filogenetico*))
                  (unless (eql perfil-filo-g nil)
                    (dolist (i (load-db (concatenate 'string workingdir "/database/ppi-conserved-neighborhood/" g ".db")))
                            (setf perfil-filo-A (find-if #'(lambda (x) (equalp (ppi-struct-genea i) (first x))) perfil-filo-g))
                            (setf perfil-filo-B (find-if #'(lambda (x) (equalp (ppi-struct-geneb i) (first x))) perfil-filo-g))
                            (unless (or (eql perfil-filo-A nil)
                                        (eql perfil-filo-B nil)
                                        );or
                              (if (eql nil (set-exclusive-or (second perfil-filo-A) (second perfil-filo-B) :test #'string=))
                                (progn
                                 (push (list (list (first perfil-filo-A) '-- (first perfil-filo-B))
                                             (list 'Weight '= percentage-pp) (third perfil-filo-A)) lista-ppi-pp)
                                 (setf cont (incf cont))
                                 );progn
                                (if (<= (length (set-exclusive-or (second perfil-filo-A) (second perfil-filo-B) :test #'string=)) ppdifftolerated)
                                  (progn
                                   (push (list (list (first perfil-filo-A) '-- (first perfil-filo-B))
                                               (list 'Weight '= percentage-pp) (third perfil-filo-A)) lista-ppi-pp)
                                   (setf cont (incf cont))
                                   );progn
                                  )
                                ;Mesclando  PPI:
                                ;(setf (third (second i)) (+ (third (second i)) (* percentage-pp 1.0)))
                                );if
                              );unless
                            );dolist
                    (setf cont 0)
                    ;recording  PPI  by phylogenetic profile in the database.
                    (save-db (concatenate 'string diretorio-banco g ".db") lista-ppi-pp)
                    (setf lista-ppi-pp nil)
                    );unlles
                  );progn
                 );if probe-file
               (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (length *genomes-files*)) 100)) laco)))
                        (format *query-io* "=")
                        (force-output *query-io*)
                        )
               (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (length *genomes-files*)) 100)) laco))))
               (incf genomenumber)
               );dolist genomas-files
       (format t "]~%")
       );End of Condition and Action 1

      ;Early  Condition and Action 2
      ((and (string= ppcomplete "Y") (> (length *genomes-files*) 1))
       ;;Generating  PPI by phylogenetic profile
       (phylogenetic-profiles-complete percentage-pp ppdifftolerated workingdir)
       );End Condition and Action 2

      ;Start  of Condition and Action 3
      ((and (string= pplimit "Y") (> (length *genomes-files*) 1))
       ;;Generating  PPI by phylogenetic profile
       (phylogenetic-profiles-ppiterlimit percentage-pp ppdifftolerated ppiterlimit workingdir)
       );End Condition and Action 3

      ;Start  of Condition and Action 4
      ((and (string= pptrim "Y") (> (length *genomes-files*) 1))
       ;;Generating  PPI by phylogenetic profile
       (phylogenetic-profiles-trim percentage-pp ppdifftolerated trim workingdir)
       );End Condition and Action 4

      ;Start  of Condition and Action 5
      ((and (string= ppthreshold "Y") (> (length *genomes-files*) 1))
       ;;Generating  PPI by phylogenetic profile
       (phylogenetic-profiles-threshold percentage-pp ppdifftolerated threshold plusminus workingdir)
       );End Condition and Action 5

      ;Early  Condition and Action 6
      ((and (string= ppdeletegroup "Y") (> (length *genomes-files*) 1))
       ;;Generating  PPI by phylogenetic profile
       (phylogenetic-profiles-delete-clusters percentage-pp ppdifftolerated grouplimit workingdir)
       );End Condition and Action 6

      ;Early  Condition and Action 7
      ((and (string= ppdeleteprofile "Y") (> (length *genomes-files*) 1))
       ;;Generating  PPI by phylogenetic profile
       (phylogenetic-profiles-delete-perfil percentage-pp ppdifftolerated profiles workingdir)
       );End Condition and Action 7
      );cond

    ;;Taking unnecessary data out of memory from now on
    (setf *perfil-filogenetico* nil)

    ;If PPI was found for one of the  features
    (if (or *ppi-identified-cn* *ppi-identified-pp* *ppi-identified-gf*)
      (progn

        ;;  PPI Report
        ;;Printing of the amounts of interactions per feature.
        (ppi-number-report workingdir ppi-gene-fusion gene-fusion ppcn ppcomplete pplimit
                           pptrim ppthreshold ppdeletegroup ppdeleteprofile)

       ;reading  database pangenome
       (setf *genomas* (load-db (concatenate 'string workingdir "/database/amino-acid-histogram/amino-acid-histogram.db")))
       (setf *pan-genoma* (load-db (concatenate 'string workingdir "/database/pan-genome/pan-genoma.db")))
       (if (> (hash-table-count *pan-genoma*) 0)
         (progn
          ;;Pangenome  report
          (setf caminho (concatenate 'string workingdir "/pan-genome/" "pan-genome-format-1.txt"))
          (relatorio-pangenoma-1 caminho aadifflimit aacheckminlimit)
          (setf caminho (concatenate 'string workingdir "/pan-genome/" "pan-genome-format-2.txt"))
          (relatorio-pangenoma-2 caminho)

          ;;Preserved gene neighborhood report
          (setf caminho (concatenate 'string workingdir "/gene-neighborhood-conservation-report/" "report.txt"))
          (relatorio-conservao-genica workingdir caminho aadifflimit aacheckminlimit expansion-type w1 window-size)
          )
         )
       ;;Taking unnecessary hash tables out of memory from now on.
       (setf *pan-genoma* nil)
       (setf *genomas* nil)
       ;------------------------------------------------------------------------

       ;;Relatório de phylogenetic profiles
       (if (and (> (length *genomes-files*) 1)
                (string= ppcomplete "Y")
                )
         (progn
          (setf caminho (concatenate 'string workingdir "/phylogenetic-profiles-report/report.txt"))
          (relatorio-perfil-filogenetico workingdir caminho percentage-pp ppdifftolerated pphistofilter
                                         ppaadifflimit ppaacheckminlimit aadifflimit aacheckminlimit)
          );progn
         );if pp report

       ;;Report on gene mergers
       (if *ppi-identified-gf*
         (progn
          (setf caminho (concatenate 'string workingdir "/gene-fusion-report/report.txt"))
          (gene-fusion-report caminho aadifflimit aacheckminlimit)
          )
         )

       ;;Resetting progress bar variables in case -ppcn is chosen.
       (setf genomenumber 0)
       (setf laco 0)

       (format t "~%Writing ppi files;~%")
       (format t "[05 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]%~%")
       (format *query-io* "[")
       (force-output *query-io*)
       ;;This  dolist generates  the  ppi-end of each genome
        (dolist (g *genomes-files*)

                (if (and *ppi-identified-cn* (probe-file (concatenate 'string workingdir "/database/ppi-conserved-neighborhood/" g ".db")))
                  (setf ppi-final (append ppi-final (load-db (concatenate 'string workingdir "/database/ppi-conserved-neighborhood/" g ".db"))))
                  )
                (if (and *ppi-identified-pp* (probe-file (concatenate 'string workingdir "/database/ppi-phylogenetic-profiles/" g ".db")))
                  (setf ppi-final (append ppi-final (load-db (concatenate 'string workingdir "/database/ppi-phylogenetic-profiles/" g ".db"))))
                  )
                (if *ppi-identified-gf*
                  (setf ppi-final (append ppi-final (gethash g ppi-gene-fusion)))
                  )

                ;This line was necessary to resolve a problem arising from the following  ordering.
                (setf ppi-final (append ppi-final (list)))

                ;;Ordering the final PPI of the genome g
                (setf ppi-final (sort ppi-final #'< :key #'(lambda (x) (ppi-struct-position x))))

                ;;Formatting the weight value of the final PPI  of the genome g
                (dolist (i ppi-final)
                        (setf (ppi-struct-weight i)
                              (with-input-from-string
                                (in (format nil "~,3f" (ppi-struct-weight i)))
                                (read in))
                              );setf
                        );dolist

                ;;Recording PPI .dot file for genome g 
                (setf caminho (concatenate 'string workingdir "/ppi-files/" g ".dot"))
                (ensure-directories-exist caminho)
                (with-open-file (str caminho
                                     :direction :output
                                     :if-exists :supersede
                                     :if-does-not-exist :create)

                  (format str (concatenate 'string "graph " g " {~%"))

                  (dolist (i ppi-final)
                          (format str "~S -- ~S [WEIGHT = ~S];~%"
                                  ;(ppi-struct-genea i)
                                  ;(ppi-struct-geneb i)
                                  (subseq (ppi-struct-genea i) (1+ (position #\. (ppi-struct-genea i))))
                                  (subseq (ppi-struct-geneb i) (1+ (position #\. (ppi-struct-geneb i))))
                                  (ppi-struct-weight i)
                                  );format
                          );dolist
                  (format str "}~%")
                  );with-open-file

                ;;Resetting the  PPIvariable -end to the next genome g
                (setf ppi-final nil)

                ;;Start of generation of  PPI files  . SIF
                ;;Recording PPI  .sif file for genome g
                (setf caminho (concatenate 'string workingdir "/ppi-files/" g ".sif"))
                (ensure-directories-exist caminho)
                (with-open-file (str caminho
                                     :direction :output
                                     :if-exists :supersede
                                     :if-does-not-exist :create)

                  (if (and *ppi-identified-cn* (probe-file (concatenate 'string workingdir "/database/ppi-conserved-neighborhood/" g ".db")))
                    (progn

                     (setf ppi-final (append ppi-final (load-db (concatenate 'string workingdir "/database/ppi-conserved-neighborhood/" g ".db"))))
                     (setf ppi-final (sort ppi-final #'< :key #'(lambda (x) (ppi-struct-position x))))

                     (dolist (i ppi-final)
                             (format str "~S~Ccn~C~S~%"
                                     ;(ppi-struct-genea i)
                                     ;(ppi-struct-geneb i)
                                     (subseq (ppi-struct-genea i) (1+ (position #\. (ppi-struct-genea i))))
                                     #\tab
                                     #\tab
                                     (subseq (ppi-struct-geneb i) (1+ (position #\. (ppi-struct-geneb i))))
                                     );format
                             );dolist
                     (setf ppi-final nil)
                     );progn
                    );if

                  (if (and *ppi-identified-pp* (probe-file (concatenate 'string workingdir "/database/ppi-phylogenetic-profiles/" g ".db")))
                    (progn
                     ;;Ordering the  PPI'sof  the genome g
                     (setf ppi-final (append ppi-final (load-db (concatenate 'string workingdir "/database/ppi-phylogenetic-profiles/" g ".db"))))
                     (setf ppi-final (sort ppi-final #'< :key #'(lambda (x) (ppi-struct-position x))))

                     (dolist (i ppi-final)
                             (format str "~S~Cpp~C~S~%"
                                     ;(ppi-struct-genea i)
                                     ;(ppi-struct-geneb i)
                                     (subseq (ppi-struct-genea i) (1+ (position #\. (ppi-struct-genea i))))
                                     #\tab
                                     #\tab
                                     (subseq (ppi-struct-geneb i) (1+ (position #\. (ppi-struct-geneb i))))
                                     );format
                             );dolist
                     (setf ppi-final nil)
                     );progn
                    );if

                  (if *ppi-identified-gf*
                    (progn
                     ;;Ordering the  PPI'sof  the genome g
                     (setf (gethash g ppi-gene-fusion) (sort (gethash g ppi-gene-fusion) #'< :key #'(lambda (x) (ppi-struct-position x))))

                     (dolist (i (gethash g ppi-gene-fusion))
                             (format str "~S~Cgf~C~S~%"
                                     ;(ppi-struct-genea i)
                                     ;(ppi-struct-geneb i)
                                     (subseq (ppi-struct-genea i) (1+ (position #\. (ppi-struct-genea i))))
                                     #\tab
                                     #\tab
                                     (subseq (ppi-struct-geneb i) (1+ (position #\. (ppi-struct-geneb i))))
                                     );format
                             );dolist
                     );progn
                    );if
                  );with-open-file
                ;;End of generation of  PPI files  . SIF

                (dotimes (j (ceiling (- (* 60 (/ (/ (* (1+ genomenumber) 100) (length *genomes-files*)) 100)) laco)))
                         (format *query-io* "=")
                         (force-output *query-io*)
                         )
                (setf laco (ceiling (+ laco (- (* 60 (/ (/ (* (1+ genomenumber) 100) (length *genomes-files*)) 100)) laco))))
                (incf genomenumber)

                );dolist
       (format t "]~%")

       (format t "~%~%PPi prediction completed successfully.~%")

       ;reading  database pangenome
       (setf *pan-genoma* (load-db (concatenate 'string workingdir "/database/pan-genome/pan-genoma.db")))
       ;;Displaying the path of the  pangenome  archives  and  conserved gene neighborhood.
       (if (> (hash-table-count *pan-genoma*) 0)
         (progn
          (if (char= #\/ (elt workingdir (1- (length workingdir))))
            (progn
             (format t "~%Pan-genome path: ~a" (concatenate 'string workingdir "pan-genome/"))
             (format t "~%Gene-neighborhood-conservation-report path: ~a" (concatenate 'string workingdir "gene-neighborhood-conservation-report/" "report.txt"))
             )
            (progn
             (format t "~%Pan-genome path: ~a" (concatenate 'string workingdir "/pan-genome/"))
             (format t "~%Gene-neighborhood-conservation-report path: ~a" (concatenate 'string workingdir "/gene-neighborhood-conservation-report/" "report.txt"))
             )
            );if char=
          );progn
         );if hash-table.
       ;;Taking the  pangenome out of ram  .
       (setf *pan-genoma* nil)
       ;------------------------------------------------------------------------

       ;;Exhibiting the path of the phylogenetic profile report.
       (if (and (> (length *genomes-files*) 1)
                (string= ppcomplete "Y")
                )
         (if (char= #\/ (elt workingdir (1- (length workingdir))))
           (format t (concatenate 'string "~%Phylogenetic-profiles-report path: "
                                  (concatenate 'string workingdir "phylogenetic-profiles-report/report.txt")))
           (format t (concatenate 'string "~%Phylogenetic-profiles-report path: "
                                  (concatenate 'string workingdir "/phylogenetic-profiles-report/report.txt")))
           )
         )

       ;;Displaying the path of the gene fusion report.
       (if *ppi-identified-gf*
         (if (char= #\/ (elt workingdir (1- (length workingdir))))
           (format t (concatenate 'string "~%Gene-fusion-report path: "
                                  (concatenate 'string workingdir "gene-fusion-report/report.txt")))
           (format t (concatenate 'string "~%Gene-fusion-report path: "
                                  (concatenate 'string workingdir "/gene-fusion-report/report.txt")))
           )
         )

       ;;Displaying the path of the .     PPI dot
       (if (char= #\/ (elt workingdir (1- (length workingdir))))
         (format t (concatenate 'string "~%PPi-files path: " (concatenate 'string workingdir "ppi-files")))
         (format t (concatenate 'string "~%PPi-files path: " (concatenate 'string workingdir "/ppi-files")))
         )

       ;;Viewing  the paths of the  PPI report
       (if (char= #\/ (elt workingdir (1- (length workingdir))))
         (format t (concatenate 'string "~%PPi-report path: " (concatenate 'string workingdir "ppi-report/report.txt")))
         (format t (concatenate 'string "~%PPi-report path: " (concatenate 'string workingdir "/ppi-report/report.txt")))
         )

       );progn  End of the stock block in case there is  PPI  of one or more  features
      ;-------------------------------------------------------------------------

      ;;If I have not been included  PPIs,make:--------------------------------
      (progn
       ;reading  database pangenome
       (setf *pan-genoma* (load-db (concatenate 'string workingdir "/database/pan-genome/pan-genoma.db")))
       (if (= (hash-table-count *pan-genoma*) 0)
         (progn
          (terpri)
          (format t "No ppi identified~%")
          )
         )
       );progn
      );if *ppi-identified-pp* *ppi-identified-cn* *ppi-identified-gf*

    ;reading  database pangenome
    (setf *pan-genoma* (load-db (concatenate 'string workingdir "/database/pan-genome/pan-genoma.db")))
    ;;Unless there were no genes in the pangenome  or no  PPI was found, make:
    (unless (or (= (hash-table-count *pan-genoma*) 0) (or *ppi-identified-pp* *ppi-identified-cn* *ppi-identified-gf*))
      (setf *genomas* (load-db (concatenate 'string workingdir "/database/amino-acid-histogram/amino-acid-histogram.db")))
      ;;Creating a pangenome  report
      (setf caminho (concatenate 'string workingdir "/pan-genome/" "pan-genome-format-1.txt"))
      (relatorio-pangenoma-1 caminho aadifflimit aacheckminlimit)
      (setf caminho (concatenate 'string workingdir "/pan-genome/" "pan-genome-format-2.txt"))
      (relatorio-pangenoma-2 caminho)
      (setf caminho (concatenate 'string workingdir "/gene-neighborhood-conservation-report/" "report.txt"))
      (relatorio-conservao-genica workingdir caminho aadifflimit aacheckminlimit expansion-type w1 window-size)
      (terpri)
      (terpri)
      (format t "No ppi identified~%")
      (if (char= #\/ (elt workingdir (1- (length workingdir))))
        (format t "~%Pan-genome path: ~a" (concatenate 'string workingdir "pan-genome/"))
        (format t "~%Pan-genome path: ~a" (concatenate 'string workingdir "/pan-genome/"))
        )
      (if (char= #\/ (elt workingdir (1- (length workingdir))))
        (format t "~%Gene-neighborhood-conservation-report path: ~a" (concatenate 'string workingdir "gene-neighborhood-conservation-report/" "report.txt"))
        (format t "~%Gene-neighborhood-conservation-report path: ~a" (concatenate 'string workingdir "/gene-neighborhood-conservation-report/" "report.txt"))
        )
      );unless
    );let genPPI
  );defun genPPI
;;------------------------------------------------------------------------------

;;Main function to receive  program parameters and run it
(defun main ()
  (let (
        (parameter 1)
        (parametertag)
        (parametervalue)
        (help T)
        (workingdir)
        (gene-fusion "N")
        (percentage-gf 0.05)
        (percentage-pp 0.30)
        (aux-percentage-pp nil)
        (ppdifftolerated 0)
        (pphistofilter "N")
        (ppaadifflimit 0)
        (ppaacheckminlimit 26)
        (ppcn "N")
        (ppcomplete "N")
        (pplimit "N")
        (ppiterlimit 500000)
        (pptrim "N")
        (trim 45000)
        (ppthreshold "N")
        (threshold 0)
        (plusminus "nil")
        (ppdeletegroup "N")
        (grouplimit 45000)
        (ppdeleteprofile "N")
        (profiles "0")
        (expansion-type "fixed")
        (window-size 1)
        (percentage-cn 0.65)
        (aux-percentage-cn nil)
        (w1 10)
        (cw1 4)
        (w2 7)
        (cw2 3)
        (w3 5)
        (cw3 2)
        (w4 3)
        (cw4 1)
        (aadifflimit 1)
        (aacheckminlimit 25)
        (directory-num 0)
        (qtd-arquivos-fasta 0)
    	   ); let pars

   	;The number 28 fixed in this code means that the program can receive up to 28 parameters on the command line.
   	;I do not worry about how many parameters have been passed, I just check if each of the six  possible
   	;was passed and if so I process the parameter.
   	(dotimes (set 31)
           		;Parameters must always be in tag pairs and tag value. Therefore, for each parameter found,
           		;the first and odd is the tag and the second and even is the value of the parameter.

             (setf parametertag (nth parameter sb-ext:*posix-argv*))

             (unless (or (string= parametertag "-pphistofilter") (string= parametertag "-ppcn")
                         (string= parametertag "-ppcomplete") (string= parametertag "-help")
                         (string= parametertag "-genefusion") (string= parametertag "-version")
                         )
               (setf parametervalue (nth (+ parameter 1) sb-ext:*posix-argv*))
               )

          	  ;(if (and parametertag parametervalue)
             ;Unless the tag does not exist, it turns out who it is.
             (unless (not parametertag)
               (setf help nil)
          	    (cond
                 ((string= parametertag "-version")
                  (progn
                   (terpri)
                   (format t "GENPPI VERSION: 1.1~%")
                   (format t "RELEASE NUMBER: 55e436b354ba1f3575529aaccc51002e7bf68ee6~%")
                   (format t "REPOSITORY: https://github.com/santosardr/genppi ~%")
                   (terpri)
                   (SB-EXT:EXIT)
                   )
                  )
                 ((string= parametertag "-help")
                  (progn
                   (terpri)
                   (format t "Mandatory parameter~%~%")
                   (format t "Directory parameter:~%")
                   (format t "-dir <workingdir> directory~%~%")
                   (format t "Optional parameters~%~%")
                   (format t "Parameters of ppi by conserved neighbourhood:~%")
                   (format t "-cnp <conserved neighbourhood score percentage> (65-default)~%")
                   (format t "-expt <type of expansion of the conserved gene neighborhood> 'fixed' or 'dynamic' (fixed-default)~%~%")
                   (format t "Parameters of conserved neighbourhood by fixed expansion:~%")
                   (format t "-w1 <window size 1> (10-default)~%")
                   (format t "-cw1 <gene conservation required for window 1> (4-default)~%")
                   (format t "-w2 <window size 2> (7-default)~%")
                   (format t "-cw2 <gene conservation required for window 2> (3-default)~%")
                   (format t "-w3 <window size 3> (5-default)~%")
                   (format t "-cw3 <gene conservation required for window 3> (2-default)~%")
                   (format t "-w4 <window size 4> (3-default)~%")
                   (format t "-cw4 <gene conservation required for window 4> (1-default)~%~%")
                   (format t "Parameters of conserved neighbourhood by dynamic expansion:~%")
                   (format t "-ws <dynamic expansion window size> (1-default)~%~%")
                   (format t "Parameters of ppi by phylogenetic profile:~%")
                   (format t "-ppp <phylogenetic profiles score percentage> (30-default)~%")
                   (format t "-ppdifftolerated <difference in phylogenetic profiles tolerated to infer ppi> (0-default)~%~%")
                   (format t "Amino acid histogram parameter settings for the phylogenetic profile:~%")
                   (format t "-pphistofilter <build the phylogenetic profile of genes with a higher percentage of identity>~%")
                   (format t "-ppaadifflimit <amino acid difference limit> (0-default)~%")
                   (format t "-ppaacheckminlimit <minimum amount of amino acids to check> (26-default)~%~%")
                   (format t "Methods of ppi prediction by phylogenetic profile and its parameters:~%~%")
                   (format t "Method 1 - Predict ppi by phylogenetic profile only for interactions predicted by conserved neighborhood~%")
                   (format t "-ppcn~%~%")
                   (format t "Method 2 - PPi prediction by phylogenetic profile without filters 2~%")
                   (format t "-ppcomplete~%~%")
                   (format t "Method 3 - Prediction of ppi by phylogenetic profile with a limit of interactions~%")
                   (format t "-ppiterlimit <maximum number of interactions desired> (500000-default)~%~%")
                   (format t "Method 4 - Prediction of ppi by phylogenetic profile with interactions limit by weight~%")
                   (format t "-trim <maximum number of interactions by weight> (45000-default)~%~%")
                   (format t "Method 5 - Prediction of ppi by the phylogenetic profile only for genes with profiles that cover a greater or lesser number of genomes than an informed threshold~%")
                   (format t "-threshold <phylogenetic profiles threshold>~%")
                   (format t "-plusminus <parameter that receives the greater than or less than sign to apply the ppthreshod filter> '<' or '>' (signs greater than and less than, must be enclosed in single or double quotes)~%~%")
                   (format t "Method 6 - Delete groups of ppi predicted by phylogenetic profile that exceed a limit of interactions by weight~%")
                   (format t "-grouplimit <limit of tolerated interactions to maintain a group of ppi with the same weight> (45000-default)~%~%")
                   (format t "Method 7 - To exclude genes with unwanted profiles in predicting ppi by phylogenetic profile~%")
                   (format t (concatenate 'string "-profiles <number of genomes in the unwanted profiles> Entry example: 7 (genes that co-occur in a total of 7 genomes will be excluded). To insert more than one profile, the entry must be enclosed in single or double quotes, and the values separated by semicolons. Example: ""\"7; 15; 21\"~%~%"))
                   (format t "Parameters of ppi by gene fusion:~%")
                   (format t "-genefusion <make ppi predictions by gene fusion>~%")
                   (format t "-gfp <gene fusion score percentage> (5-default)~%")
                   (format t "Note: if gene fusion is included, there will be an increase of more than 100% in the execution time~%~%")
                   (format t "Amino acid histogram parameters:~%")
                   (format t "-aadifflimit <amino acid difference limit> (1-default)~%")
                   (format t "-aacheckminlimit <minimum amount of amino acids to check> (25-default)~%~%")
                   (format t "Possible configurations of the histogram parameters:~%")
                   (format t "-------------------------------------------------------------~%")
                   (format t "-aadifflimit   -aacheckminlimit   minimal protein identity %~%")
                   (format t "-------------------------------------------------------------~%")
                   (format t "      0               26                     100%~%")
                   (format t "-------------------------------------------------------------~%")
                   (format t "      0               25                     100%~%")
                   (format t "-------------------------------------------------------------~%")
                   (format t "      0               24                    97,96%~%")
                   (format t "-------------------------------------------------------------~%")
                   (format t "      0               23                    96,94%~%")
                   (format t "-------------------------------------------------------------~%")
                   (format t "      0               22                    96,94%~%")
                   (format t "-------------------------------------------------------------~%")
                   (format t "      0               21                    94,68%~%")
                   (format t "-------------------------------------------------------------~%")
                   (format t "      0               20                    91,75%~%")
                   (format t "-------------------------------------------------------------~%")
                   (format t "      0               19                    85,57%~%")
                   (format t "-------------------------------------------------------------~%")
                   (format t "      0               18                    50,00%~%")
                   (format t "-------------------------------------------------------------~%")
                   (format t "      1               26                    97,87%~%")
                   (format t "-------------------------------------------------------------~%")
                   (format t "      1               25                    92,55%~%")
                   (format t "-------------------------------------------------------------~%")
                   (terpri)
                   (SB-EXT:EXIT)
                   )
                  )
                	((string= parametertag "-dir")(setf workingdir (if parametervalue
                                                                  (if (char= (elt parametervalue 0) #\-) "./" parametervalue)
                                                                  "./"
                                                                  );if parametervalue
                                                     ));execution directory
                 ((string= parametertag "-genefusion") (setf gene-fusion "Y"))
                 ((string= parametertag "-gfp")
                  (setf percentage-gf  (if parametervalue
                                         (if (numberp (read-from-string parametervalue))
                                           (if (or (< (parse-integer parametervalue) 0) (> (parse-integer parametervalue) 100))
                                             (progn
                                              (terpri)
                                              (format t "Warning:~%")
                                              (format t "Invalid parameter: -gfp <gene fusion score percentage> ~%")
                                              (format t "The prediction percentage by gene fusion can only take values from 0 to 100.")
                                              (terpri)
                                              (SB-EXT:EXIT)
                                              )
                                             (/ (parse-integer parametervalue) 100)
                                             )
                                           (progn
                                            (terpri)
                                            (format t "Warning:~%")
                                            (format t"Invalid parameter: -gfp <gene fusion score percentage>~%")
                                            (format t "The prediction percentage by gene fusion can only take values from 0 to 100.")
                                            (terpri)
                                     	      (SB-EXT:EXIT)
                                            )
                                           )
                                         percentage-gf
                                         )
                        ))
                 ((string= parametertag "-ppp")
                  (setf percentage-pp  (if parametervalue
                                         (if (numberp (read-from-string parametervalue))
                                           (if (or (< (parse-integer parametervalue) 0) (> (parse-integer parametervalue) 100))
                                             (progn
                                              (terpri)
                                              (format t "Warning:~%")
                                              (format t "Invalid parameter: -ppp <phylogenetic profiles score percentage>~%")
                                              (format t "The prediction percentage by phylogenetic profile can only take values from 0 to 100.")
                                              (terpri)
                                              (SB-EXT:EXIT)
                                              )
                                             (progn
                                              (setf aux-percentage-pp t)
                                              (/ (parse-integer parametervalue) 100)
                                              )
                                             )
                                           (progn
                                            (terpri)
                                            (format t "Warning:~%")
                                            (format t"Invalid parameter: -ppp <phylogenetic profiles score percentage>~%")
                                            (format t "The prediction percentage by phylogenetic profile can only take values from 0 to 100.")
                                            (terpri)
                                     	      (SB-EXT:EXIT)
                                            )
                                           )
                                         percentage-pp
                                         )
                        ))
                 ((string= parametertag "-ppdifftolerated")
                  (setf ppdifftolerated  (if parametervalue
                                           (if (numberp (read-from-string parametervalue))
                                             (if (< (parse-integer parametervalue) 0)
                                               (progn
                                                (terpri)
                                                (format t "Warning:~%")
                                                (format t "Invalid parameter: -ppdifftolerated <difference in phylogenetic profiles tolerated to infer ppi>~%")
                                                (format t "The difference in phylogenetic profiles tolerated to infer ppi cannot be less than zero.")
                                                (terpri)
                                                (SB-EXT:EXIT)
                                                )
                                               (- (parse-integer parametervalue) 0)
                                               )
                                             (progn
                                              (terpri)
                                              (format t "Warning:~%")
                                              (format t "Invalid parameter: -ppdifftolerated <difference in phylogenetic profiles tolerated to infer ppi>~%")
                                              (format t "The difference in phylogenetic profiles tolerated to infer ppi must be an integer.")
                                              (terpri)
                                       	      (SB-EXT:EXIT)
                                              )
                                             )
                                           ppdifftolerated
                                           )
                        ))
                 ((string= parametertag "-pphistofilter") (setf pphistofilter  "Y"))
                 ((string= parametertag "-ppaadifflimit")
                  (setf ppaadifflimit  (if parametervalue
                                         (if (numberp (read-from-string parametervalue))
                                           (if (or (= (parse-integer parametervalue) 0) (= (parse-integer parametervalue) 1))
                                             (parse-integer parametervalue)
                                             (progn
                                              (terpri)
                                              (format t "Warning:~%")
                                              (format t "Invalid parameter: -ppaadifflimit <amino acid difference limit>~%~%")
                                              (terpri)
                                              (SB-EXT:EXIT)
                                              )
                                             )
                                           (progn
                                            (terpri)
                                            (format t "Warning:~%")
                                            (format t "Invalid parameter: -ppaadifflimit <amino acid difference limit>~%~%")
                                            (terpri)
                                            (SB-EXT:EXIT)
                                            )
                                           )
                                         ppaadifflimit
                                         )
                        ))
                 ((string= parametertag "-ppaacheckminlimit")
                  (setf ppaacheckminlimit  (if parametervalue
                                             (if (numberp (read-from-string parametervalue))
                                               (if (or (= (parse-integer parametervalue) 26) (= (parse-integer parametervalue) 25)
                                                       (= (parse-integer parametervalue) 24) (= (parse-integer parametervalue) 23)
                                                       (= (parse-integer parametervalue) 22) (= (parse-integer parametervalue) 21)
                                                       (= (parse-integer parametervalue) 20) (= (parse-integer parametervalue) 19)
                                                       (= (parse-integer parametervalue) 18)
                                                       );or
                                                 (parse-integer parametervalue)
                                                 (progn
                                                  (terpri)
                                                  (format t "Warning:~%")
                                                  (format t "Invalid parameter: -ppaacheckminlimit <minimum amount of amino acids to check>~%~%")
                                                  (terpri)
                                                  (SB-EXT:EXIT)
                                                  )
                                                 )
                                               (progn
                                                (terpri)
                                                (format t "Warning:~%")
                                                (format t "Invalid parameter: -ppaacheckminlimit <minimum amount of amino acids to check>~%~%")
                                                (terpri)
                                                (SB-EXT:EXIT)
                                                )
                                               )
                                             ppaacheckminlimit
                                             )
                        ))
                 ((string= parametertag "-ppcn")
                  (setf ppcn "Y")
                  (unless (and (string= ppcomplete "N") (string= pplimit "N") (string= pptrim "N")
                               (string= ppthreshold "N") (string= ppdeletegroup "N") (string= ppdeleteprofile "N"))
                    (terpri)
                    (format t "Warning:~%")
                    (format t "The use of methods of prediction of ppi by phylogenetic profile is limited to only 1 among the 7 available methods")
                    (terpri)
                    (SB-EXT:EXIT)
                    ); unless or
                  )
                 ((string= parametertag "-ppcomplete")
                  (setf ppcomplete "Y")
                  (unless (and (string= ppcn "N") (string= pplimit "N") (string= pptrim "N")
                               (string= ppthreshold "N") (string= ppdeletegroup "N") (string= ppdeleteprofile "N"))
                    (terpri)
                    (format t "Warning:~%")
                    (format t "The use of methods of prediction of ppi by phylogenetic profile is limited to only 1 among the 7 available methods")
                    (terpri)
                    (SB-EXT:EXIT)
                    ); unless or
                  )
                 ((string= parametertag "-ppiterlimit")
                  (setf pplimit "Y")
                  (unless (and (string= ppcomplete "N") (string= ppcn "N") (string= pptrim "N")
                               (string= ppthreshold "N") (string= ppdeletegroup "N") (string= ppdeleteprofile "N"))
                    (terpri)
                    (format t "Warning:~%")
                    (format t "The use of methods of prediction of ppi by phylogenetic profile is limited to only 1 among the 7 available methods")
                    (terpri)
                    (SB-EXT:EXIT)
                    ); unless or
                  (setf ppiterlimit  (if parametervalue
                                       (if (char= (elt parametervalue 0) #\-)
                                         ppiterlimit
                                         (if (numberp (read-from-string parametervalue))
                                           (if (<= (parse-integer parametervalue) 0)
                                             (progn
                                              (terpri)
                                              (format t "Warning:~%")
                                              (format t "Invalid parameter: -ppiterlimit <maximum number of interactions desired>~%")
                                              (format t "The -ppiterlimit parameter must be greater than zero.")
                                              (terpri)
                                              (SB-EXT:EXIT)
                                              )
                                             (- (parse-integer parametervalue) 0)
                                             )
                                           (progn
                                            (terpri)
                                            (format t "Warning:~%")
                                            (format t "Invalid parameter: -ppiterlimit <maximum number of interactions desired>~%")
                                            (format t "Invalid character in the -ppiterlimit parameter.")
                                            (terpri)
                                            (SB-EXT:EXIT)
                                            )
                                           )
                                         );if
                                       ppiterlimit
                                       )
                        ))
                 ((string= parametertag "-trim")
                  (setf pptrim "Y")
                  (unless (and (string= ppcomplete "N") (string= pplimit "N") (string= ppcn "N")
                               (string= ppthreshold "N") (string= ppdeletegroup "N") (string= ppdeleteprofile "N"))
                    (terpri)
                    (format t "Warning:~%")
                    (format t "The use of methods of prediction of ppi by phylogenetic profile is limited to only 1 among the 7 available methods")
                    (terpri)
                    (SB-EXT:EXIT)
                    ); unless or

                  (setf trim  (if parametervalue
                                (if (char= (elt parametervalue 0) #\-)
                                  trim
                                  (if (numberp (read-from-string parametervalue))
                                    (if (<= (parse-integer parametervalue) 0)
                                      (progn
                                       (terpri)
                                       (format t "Warning:~%")
                                       (format t "Invalid parameter: -trim <maximum number of interactions by weight>~%")
                                       (format t "The -trim parameter must be greater than zero.")
                                       (terpri)
                                       (SB-EXT:EXIT)
                                       )
                                      (- (parse-integer parametervalue) 0)
                                      )
                                    (progn
                                     (terpri)
                                     (format t "Warning:~%")
                                     (format t "Invalid parameter: -trim <maximum number of interactions by weight>~%")
                                     (format t "Invalid character in the -trim parameter.")
                                     (terpri)
                              	      (SB-EXT:EXIT)
                                     )
                                    )
                                  );if
                                trim
                                )
                        ))
                 ((string= parametertag "-threshold")
                  (setf ppthreshold "Y")
                  (unless (and (string= ppcomplete "N") (string= pplimit "N") (string= pptrim "N")
                               (string= ppcn "N") (string= ppdeletegroup "N") (string= ppdeleteprofile "N"))
                    (terpri)
                    (format t "Warning:~%")
                    (format t "The use of methods of prediction of ppi by phylogenetic profile is limited to only 1 among the 7 available methods")
                    (terpri)
                    (SB-EXT:EXIT)
                    ); unless or
                  (setf threshold  (if parametervalue
                                     (if (numberp (read-from-string parametervalue))
                                       (if (<= (parse-integer parametervalue) 0)
                                         (progn
                                          (terpri)
                                          (format t "Warning:~%")
                                          (format t "Invalid parameter: -threshold <phylogenetic profiles threshold>~%")
                                          (format t "The -threshold parameter must be greater than zero.")
                                          (terpri)
                                          (SB-EXT:EXIT)
                                          )
                                         (- (parse-integer parametervalue) 0)
                                         )
                                       (progn
                                        (terpri)
                                        (format t "Warning:~%")
                                        (format t "Invalid parameter: -threshold <phylogenetic profiles threshold>~%")
                                        (format t "Invalid character in the -threshold parameter.")
                                        (terpri)
                                 	      (SB-EXT:EXIT)
                                        )
                                       )
                                     (progn
                                      (terpri)
                                      (format t "Warning:~%")
                                      (format t"When using the method 5, it is necessary to enter a value greater than zero through the -threshold parameter")
                                      (terpri)
                                      (SB-EXT:EXIT)
                                      )
                                     )
                        ))
                 ((string= parametertag "-plusminus")
                  (setf ppthreshold "Y")
                  (unless (and (string= ppcomplete "N") (string= pplimit "N") (string= pptrim "N")
                               (string= ppcn "N") (string= ppdeletegroup "N") (string= ppdeleteprofile "N"))
                    (terpri)
                    (format t "Warning:~%")
                    (format t "The use of methods of prediction of ppi by phylogenetic profile is limited to only 1 among the 7 available methods")
                    (terpri)
                    (SB-EXT:EXIT)
                    ); unless or
                  (setf plusminus  (if parametervalue
                                     (progn
                                      (unless (or (string= (string parametervalue) "<") (string= (string parametervalue) ">"))
                                        (progn
                                         (terpri)
                                         (format t "Warning:~%")
                                         (format t "Invalid parameter: -plusminus <parameter that receives the greater than or less than sign to apply the ppthreshod filter>~%")
                                         (format t "The -plusminus parameter can take '<' or '>'")
                                         (terpri)
                                         (SB-EXT:EXIT)
                                         );progn
                                        );unless
                                      (string parametervalue)
                                      );progn
                                     (progn
                                      (terpri)
                                      (format t "Warning:~%")
                                      (format t"When using the method 5, it is necessary to insert the less than (<) or greater than (>) sign through the -plusminus parameter (signs greater than and less than, must be enclosed in single or double quotes)")
                                      (terpri)
                                      (SB-EXT:EXIT)
                                      )
                                     )))
                 ((string= parametertag "-grouplimit")
                  (setf ppdeletegroup "Y")
                  (unless (and (string= ppcomplete "N") (string= pplimit "N") (string= pptrim "N")
                               (string= ppthreshold "N") (string= ppcn "N") (string= ppdeleteprofile "N"))
                    (terpri)
                    (format t "Warning:~%")
                    (format t "The use of methods of prediction of ppi by phylogenetic profile is limited to only 1 among the 7 available methods")
                    (terpri)
                    (SB-EXT:EXIT)
                    ); unless or
                  (setf grouplimit  (if parametervalue
                                      (if (char= (elt parametervalue 0) #\-)
                                        grouplimit
                                        (if (numberp (read-from-string parametervalue))
                                          (if (<= (parse-integer parametervalue) 0)
                                            (progn
                                             (terpri)
                                             (format t "Warning:~%")
                                             (format t "Invalid parameter: -grouplimit <limit of tolerated interactions to maintain a group of ppi with the same weight>~%")
                                             (format t "The -grouplimit parameter must be greater than zero.")
                                             (terpri)
                                             (SB-EXT:EXIT)
                                             )
                                            (- (parse-integer parametervalue) 0)
                                            )
                                          (progn
                                           (terpri)
                                           (format t "Warning:~%")
                                           (format t "Invalid parameter: -grouplimit <limit of tolerated interactions to maintain a group of ppi with the same weight>~%")
                                           (format t "Invalid character in the -grouplimit parameter.")
                                           (terpri)
                                    	      (SB-EXT:EXIT)
                                           )
                                          )
                                        );unless
                                      grouplimit
                                      )
                        ))
                 ((string= parametertag "-profiles")
                  (setf ppdeleteprofile "Y")
                  (unless (and (string= ppcomplete "N") (string= pplimit "N") (string= pptrim "N")
                               (string= ppthreshold "N") (string= ppdeletegroup "N") (string= ppcn "N"))
                    (terpri)
                    (format t "Warning:~%")
                    (format t "The use of methods of prediction of ppi by phylogenetic profile is limited to only 1 among the 7 available methods")
                    (terpri)
                    (SB-EXT:EXIT)
                    ); unless or
                  (setf profiles  (if parametervalue
                                    (progn
                                     (loop for c across parametervalue
                                       do (progn
                                           (unless (or (char= c #\;)
                                                       (char= c #\0) (char= c #\1) (char= c #\2) (char= c #\3) (char= c #\4)
                                                       (char= c #\5) (char= c #\6) (char= c #\7) (char= c #\8) (char= c #\9)
                                                       );or
                                             (progn
                                              (terpri)
                                              (format t "Warning:~%")
                                              (format t "Invalid parameter: -profiles <number of genomes in the unwanted profiles. Entry example: 7 (genes that co-occur in a total of 7 genomes will be excluded). To enter more than one profile, simply separate the input values with a semicolon. Example: 7; 15; 21>~%")
                                              (format t "Invalid character in the -profiles parameter.")
                                              (terpri)
                                              (SB-EXT:EXIT)
                                              )
                                             );unless or
                                           );progn
                                       );loop
                                     (string parametervalue)
                                     );progn
                                    (progn
                                     (terpri)
                                     (format t "Warning:~%")
                                     (format t (concatenate 'string "When using the method 7, it is necessary to enter a value greater than zero through the -profiles parameter. Entry example: 7 (genes that co-occur in a total of 7 genomes will be excluded). ""To insert more than one profile, the entry must be enclosed in single or double quotes, and the values separated by semicolons. Example: ""\"7; 15; 21\"~%~%"))
                                     (terpri)
                                     (SB-EXT:EXIT)
                                     )
                                    )))
                 ((string= parametertag "-cnp")
                  (setf percentage-cn  (if parametervalue
                                         (if (numberp (read-from-string parametervalue))
                                           (if (or (< (parse-integer parametervalue) 0) (> (parse-integer parametervalue) 100))
                                             (progn
                                              (terpri)
                                              (format t "Warning:~%")
                                              (format t "Invalid parameter: -cpp <conserved neighbourhood score percentage>~%")
                                              (format t "The prediction percentage by conserved neighbourhood can only take values from 0 to 100.")
                                              (terpri)
                                              (SB-EXT:EXIT)
                                              )
                                             (progn
                                              (setf aux-percentage-cn T)
                                              (/ (parse-integer parametervalue) 100)
                                              )
                                             )
                                           (progn
                                            (terpri)
                                            (format t "Warning:~%")
                                            (format t"Invalid parameter: -cnp <conserved neighbourhood score percentage>~%")
                                            (format t "The prediction percentage by conserved neighbourhood can only take values from 0 to 100.")
                                            (terpri)
                                     	      (SB-EXT:EXIT)
                                            )
                                           )
                                         percentage-cn)
                        ))
                 ((string= parametertag "-expt")
                  (if parametervalue
                    (if (numberp (read-from-string parametervalue))
                      (progn
                       (terpri)
                       (format t "Warning:~%")
                       (format t"Invalid parameter: -expt <type of expansion of the conserved gene neighborhood>~%")
                       (format t "The expansion of the conserved gene neighborhood can be 'fixed' or 'dynamic'")
                       (terpri)
                       (SB-EXT:EXIT)
                       )
                      (if (or (string= (string-downcase parametervalue) "fixed") (string= (string-downcase parametervalue) "dynamic"))
                        (setf expansion-type (string-downcase parametervalue))
                        (progn
                         (terpri)
                         (format t "Warning:~%")
                         (format t"Invalid parameter: -expt <type of expansion of the conserved gene neighborhood>~%")
                         (format t "The expansion of the conserved gene neighborhood can be 'fixed' or 'dynamic'")
                         (terpri)
                         (SB-EXT:EXIT)
                         )
                        )
                      )
                    (progn
                     (terpri)
                     (format t "Warning:~%")
                     (format t"Invalid parameter: -expt <type of expansion of the conserved gene neighborhood>~%")
                     (format t "The expansion of the conserved gene neighborhood can be 'fixed' or 'dynamic'")
                     (terpri)
                     (SB-EXT:EXIT)
                     )))
                 ((string= parametertag "-ws")
                  (setf window-size  (if parametervalue
                                       (if (numberp (read-from-string parametervalue))
                                         (if (<= (parse-integer parametervalue) 0)
                                           (progn
                                            (terpri)
                                            (format t "Warning:~%")
                                            (format t "The -ws parameter cannot be less than or equal to zero")
                                            (terpri)
                                            (SB-EXT:EXIT)
                                            )
                                           (parse-integer parametervalue)
                                           )
                                         (progn
                                          (terpri)
                                          (format t "Warning:~%")
                                          (format t "Invalid parameter: -ws <dynamic expansion window size>")
                                          (terpri)
                                   	      (SB-EXT:EXIT)
                                          )
                                         );if numberp
                                       (progn
                                        (terpri)
                                        (format t "Warning:~%")
                                        (format t"Invalid parameter: -ws <dynamic expansion window size>")
                                        (terpri)
                                        (SB-EXT:EXIT)
                                        )
                                       )
                        ))
                 ((string= parametertag "-w1")
                  (setf w1  (if parametervalue
                              (if (numberp (read-from-string parametervalue))
                                (if (<= (parse-integer parametervalue) 0)
                                  (progn
                                   (terpri)
                                   (format t "Warning:~%")
                                   (format t "The -w1 parameter cannot be less than or equal to zero")
                                   (terpri)
                                   (SB-EXT:EXIT)
                                   )
                                  (parse-integer parametervalue)
                                  )
                                (progn
                                 (terpri)
                                 (format t "Warning:~%")
                                 (format t"Invalid parameter: -w1 <window size 1>")
                                 (terpri)
                          	      (SB-EXT:EXIT)
                                 )
                                )
                              (progn
                               (terpri)
                               (format t "Warning:~%")
                               (format t"Invalid parameter: -w1 <window size 1>")
                               (terpri)
                               (SB-EXT:EXIT)
                               )
                              )
                        ))
                 ((string= parametertag "-cw1")
                  (setf cw1  (if parametervalue
                               (if (numberp (read-from-string parametervalue))
                                 (if (<= (parse-integer parametervalue) 0)
                                   (progn
                                    (terpri)
                                    (format t "Warning:~%")
                                    (format t "The -cw1 parameter cannot be less than or equal to zero")
                                    (terpri)
                                    (SB-EXT:EXIT)
                                    )
                                   (parse-integer parametervalue)
                                   )
                                 (progn
                                  (terpri)
                                  (format t "Warning:~%")
                                  (format t"Invalid parameter: -cw1 <gene conservation required for window 1>")
                                  (terpri)
                           	      (SB-EXT:EXIT)
                                  )
                                 )
                               (progn
                                (terpri)
                                (format t "Warning:~%")
                                (format t"Invalid parameter: -cw1 <gene conservation required for window 1>")
                                (terpri)
                                (SB-EXT:EXIT)
                                )
                               )
                        ))
                 ((string= parametertag "-w2")
                  (setf w2  (if parametervalue
                              (if (numberp (read-from-string parametervalue))
                                (if (<= (parse-integer parametervalue) 0)
                                  (progn
                                   (terpri)
                                   (format t "Warning:~%")
                                   (format t "The -w2 parameter cannot be less than or equal to zero")
                                   (terpri)
                                   (SB-EXT:EXIT)
                                   )
                                  (parse-integer parametervalue)
                                  )
                                (progn
                                 (terpri)
                                 (format t "Warning:~%")
                                 (format t"Invalid parameter: -w2 <window size 2>")
                                 (terpri)
                          	      (SB-EXT:EXIT)
                                 )
                                )
                              (progn
                               (terpri)
                               (format t "Warning:~%")
                               (format t"Invalid parameter: -w2 <window size 2>")
                               (terpri)
                               (SB-EXT:EXIT)
                               )
                              )
                        ))
                 ((string= parametertag "-cw2")
                  (setf cw2  (if parametervalue
                               (if (numberp (read-from-string parametervalue))
                                 (if (<= (parse-integer parametervalue) 0)
                                   (progn
                                    (terpri)
                                    (format t "Warning:~%")
                                    (format t "The -cw2 parameter cannot be less than or equal to zero")
                                    (terpri)
                                    (SB-EXT:EXIT)
                                    )
                                   (parse-integer parametervalue)
                                   )
                                 (progn
                                  (terpri)
                                  (format t "Warning:~%")
                                  (format t"Invalid parameter: -cw2 <gene conservation required for window 2>")
                                  (terpri)
                           	      (SB-EXT:EXIT)
                                  )
                                 )
                               (progn
                                (terpri)
                                (format t "Warning:~%")
                                (format t"Invalid parameter: -cw2 <gene conservation required for window 2>")
                                (terpri)
                                (SB-EXT:EXIT)
                                )
                               )
                        ))
                 ((string= parametertag "-w3")
                  (setf w3  (if parametervalue
                              (if (numberp (read-from-string parametervalue))
                                (if (<= (parse-integer parametervalue) 0)
                                  (progn
                                   (terpri)
                                   (format t "Warning:~%")
                                   (format t "The -w3 parameter cannot be less than or equal to zero")
                                   (terpri)
                                   (SB-EXT:EXIT)
                                   )
                                  (parse-integer parametervalue)
                                  )
                                (progn
                                 (terpri)
                                 (format t "Warning:~%")
                                 (format t"Invalid parameter: -w3 <window size 3>")
                                 (terpri)
                          	      (SB-EXT:EXIT)
                                 )
                                )
                              (progn
                               (terpri)
                               (format t "Warning:~%")
                               (format t"Invalid parameter: -w3 <window size 3>")
                               (terpri)
                               (SB-EXT:EXIT)
                               )
                              )
                        ))
                 ((string= parametertag "-cw3")
                  (setf cw3  (if parametervalue
                               (if (numberp (read-from-string parametervalue))
                                 (if (<= (parse-integer parametervalue) 0)
                                   (progn
                                    (terpri)
                                    (format t "Warning:~%")
                                    (format t "The -cw3 parameter cannot be less than or equal to zero")
                                    (terpri)
                                    (SB-EXT:EXIT)
                                    )
                                   (parse-integer parametervalue)
                                   )
                                 (progn
                                  (terpri)
                                  (format t "Warning:~%")
                                  (format t"Invalid parameter: -cw3 <gene conservation required for window 3>")
                                  (terpri)
                           	      (SB-EXT:EXIT)
                                  )
                                 )
                               (progn
                                (terpri)
                                (format t "Warning:~%")
                                (format t"Invalid parameter: -cw3 <gene conservation required for window 3>")
                                (terpri)
                                (SB-EXT:EXIT)
                                )
                               )
                        ))
                 ((string= parametertag "-w4")
                  (setf w4  (if parametervalue
                              (if (numberp (read-from-string parametervalue))
                                (if (<= (parse-integer parametervalue) 0)
                                  (progn
                                   (terpri)
                                   (format t "Warning:~%")
                                   (format t "The -w4 parameter cannot be less than or equal to zero")
                                   (terpri)
                                   (SB-EXT:EXIT)
                                   )
                                  (parse-integer parametervalue)
                                  )
                                (progn
                                 (terpri)
                                 (format t "Warning:~%")
                                 (format t"Invalid parameter: -w4 <window size 4>")
                                 (terpri)
                          	      (SB-EXT:EXIT)
                                 )
                                )
                              (progn
                               (terpri)
                               (format t "Warning:~%")
                               (format t"Invalid parameter: -w4 <window size 4>")
                               (terpri)
                               (SB-EXT:EXIT)
                               )
                              )
                        ))
                 ((string= parametertag "-cw4")
                  (setf cw4  (if parametervalue
                               (if (numberp (read-from-string parametervalue))
                                 (if (<= (parse-integer parametervalue) 0)
                                   (progn
                                    (terpri)
                                    (format t "Warning:~%")
                                    (format t "The -cw4 parameter cannot be less than or equal to zero")
                                    (terpri)
                                    (SB-EXT:EXIT)
                                    )
                                   (parse-integer parametervalue)
                                   )
                                 (progn
                                  (terpri)
                                  (format t "Warning:~%")
                                  (format t"Invalid parameter: -cw4 <gene conservation required for window 4>")
                                  (terpri)
                           	      (SB-EXT:EXIT)
                                  )
                                 )
                               (progn
                                (terpri)
                                (format t "Warning:~%")
                                (format t"Invalid parameter: -cw4 <gene conservation required for window 4>")
                                (terpri)
                                (SB-EXT:EXIT)
                                )
                               )
                        ))
                 ((string= parametertag "-aadifflimit")
                  (setf aadifflimit  (if parametervalue
                                       (if (numberp (read-from-string parametervalue))
                                         (if (or (= (parse-integer parametervalue) 0) (= (parse-integer parametervalue) 1))
                                           (parse-integer parametervalue)
                                           (progn
                                            (terpri)
                                            (format t "Warning:~%")
                                            (format t "Invalid parameter: -aadifflimit <amino acid difference limit>~%~%")
                                            (format t "Possible configurations of the histogram parameters:~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "-aadifflimit   -aacheckminlimit   minimal protein identity %~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      0               26                     100%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      0               25                     100%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      0               24                    97,96%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      0               23                    96,94%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      0               22                    96,94%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      0               21                    94,68%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      0               20                    91,75%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      0               19                    85,57%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      0               18                    50,00%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      1               26                    97,87%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      1               25                    92,55%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (terpri)
                                            (SB-EXT:EXIT)
                                            )
                                           )
                                         (progn
                                          (terpri)
                                          (format t "Warning:~%")
                                          (format t "Invalid parameter: -aadifflimit <amino acid difference limit>~%~%")
                                          (format t "Possible configurations of the histogram parameters:~%")
                                          (format t "-------------------------------------------------------------~%")
                                          (format t "-aadifflimit   -aacheckminlimit   minimal protein identity %~%")
                                          (format t "-------------------------------------------------------------~%")
                                          (format t "      0               26                     100%~%")
                                          (format t "-------------------------------------------------------------~%")
                                          (format t "      0               25                     100%~%")
                                          (format t "-------------------------------------------------------------~%")
                                          (format t "      0               24                    97,96%~%")
                                          (format t "-------------------------------------------------------------~%")
                                          (format t "      0               23                    96,94%~%")
                                          (format t "-------------------------------------------------------------~%")
                                          (format t "      0               22                    96,94%~%")
                                          (format t "-------------------------------------------------------------~%")
                                          (format t "      0               21                    94,68%~%")
                                          (format t "-------------------------------------------------------------~%")
                                          (format t "      0               20                    91,75%~%")
                                          (format t "-------------------------------------------------------------~%")
                                          (format t "      0               19                    85,57%~%")
                                          (format t "-------------------------------------------------------------~%")
                                          (format t "      0               18                    50,00%~%")
                                          (format t "-------------------------------------------------------------~%")
                                          (format t "      1               26                    97,87%~%")
                                          (format t "-------------------------------------------------------------~%")
                                          (format t "      1               25                    92,55%~%")
                                          (format t "-------------------------------------------------------------~%")
                                          (terpri)
                                          (SB-EXT:EXIT)
                                          )
                                         )
                                       (progn
                                        (terpri)
                                        (format t "Warning:~%")
                                        (format t "Invalid parameter: -aadifflimit <amino acid difference limit>~%~%")
                                        (format t "Possible configurations of the histogram parameters:~%")
                                        (format t "-------------------------------------------------------------~%")
                                        (format t "-aadifflimit   -aacheckminlimit   minimal protein identity %~%")
                                        (format t "-------------------------------------------------------------~%")
                                        (format t "      0               26                     100%~%")
                                        (format t "-------------------------------------------------------------~%")
                                        (format t "      0               25                     100%~%")
                                        (format t "-------------------------------------------------------------~%")
                                        (format t "      0               24                    97,96%~%")
                                        (format t "-------------------------------------------------------------~%")
                                        (format t "      0               23                    96,94%~%")
                                        (format t "-------------------------------------------------------------~%")
                                        (format t "      0               22                    96,94%~%")
                                        (format t "-------------------------------------------------------------~%")
                                        (format t "      0               21                    94,68%~%")
                                        (format t "-------------------------------------------------------------~%")
                                        (format t "      0               20                    91,75%~%")
                                        (format t "-------------------------------------------------------------~%")
                                        (format t "      0               19                    85,57%~%")
                                        (format t "-------------------------------------------------------------~%")
                                        (format t "      0               18                    50,00%~%")
                                        (format t "-------------------------------------------------------------~%")
                                        (format t "      1               26                    97,87%~%")
                                        (format t "-------------------------------------------------------------~%")
                                        (format t "      1               25                    92,55%~%")
                                        (format t "-------------------------------------------------------------~%")
                                        (terpri)
                                        (SB-EXT:EXIT)
                                        )
                                       )
                        ))
                 ((string= parametertag "-aacheckminlimit")
                  (setf aacheckminlimit  (if parametervalue
                                           (if (numberp (read-from-string parametervalue))
                                             (if (or (= (parse-integer parametervalue) 26) (= (parse-integer parametervalue) 25)
                                                     (= (parse-integer parametervalue) 24) (= (parse-integer parametervalue) 23)
                                                     (= (parse-integer parametervalue) 22) (= (parse-integer parametervalue) 21)
                                                     (= (parse-integer parametervalue) 20) (= (parse-integer parametervalue) 19)
                                                     (= (parse-integer parametervalue) 18)
                                                     );or
                                               (parse-integer parametervalue)
                                               (progn
                                                (terpri)
                                                (format t "Warning:~%")
                                                (format t "Invalid parameter: -aacheckminlimit <minimum amount of amino acids to check>~%~%")
                                                (format t "Possible configurations of the histogram parameters:~%")
                                                (format t "-------------------------------------------------------------~%")
                                                (format t "-aadifflimit   -aacheckminlimit   minimal protein identity %~%")
                                                (format t "-------------------------------------------------------------~%")
                                                (format t "      0               26                     100%~%")
                                                (format t "-------------------------------------------------------------~%")
                                                (format t "      0               25                     100%~%")
                                                (format t "-------------------------------------------------------------~%")
                                                (format t "      0               24                    97,96%~%")
                                                (format t "-------------------------------------------------------------~%")
                                                (format t "      0               23                    96,94%~%")
                                                (format t "-------------------------------------------------------------~%")
                                                (format t "      0               22                    96,94%~%")
                                                (format t "-------------------------------------------------------------~%")
                                                (format t "      0               21                    94,68%~%")
                                                (format t "-------------------------------------------------------------~%")
                                                (format t "      0               20                    91,75%~%")
                                                (format t "-------------------------------------------------------------~%")
                                                (format t "      0               19                    85,57%~%")
                                                (format t "-------------------------------------------------------------~%")
                                                (format t "      0               18                    50,00%~%")
                                                (format t "-------------------------------------------------------------~%")
                                                (format t "      1               26                    97,87%~%")
                                                (format t "-------------------------------------------------------------~%")
                                                (format t "      1               25                    92,55%~%")
                                                (format t "-------------------------------------------------------------~%")
                                                (terpri)
                                                (SB-EXT:EXIT)
                                                )
                                               )
                                             (progn
                                              (terpri)
                                              (format t "Warning:~%")
                                              (format t "Invalid parameter: -aacheckminlimit <minimum amount of amino acids to check>~%~%")
                                              (format t "Possible configurations of the histogram parameters:~%")
                                              (format t "-------------------------------------------------------------~%")
                                              (format t "-aadifflimit   -aacheckminlimit   minimal protein identity %~%")
                                              (format t "-------------------------------------------------------------~%")
                                              (format t "      0               26                     100%~%")
                                              (format t "-------------------------------------------------------------~%")
                                              (format t "      0               25                     100%~%")
                                              (format t "-------------------------------------------------------------~%")
                                              (format t "      0               24                    97,96%~%")
                                              (format t "-------------------------------------------------------------~%")
                                              (format t "      0               23                    96,94%~%")
                                              (format t "-------------------------------------------------------------~%")
                                              (format t "      0               22                    96,94%~%")
                                              (format t "-------------------------------------------------------------~%")
                                              (format t "      0               21                    94,68%~%")
                                              (format t "-------------------------------------------------------------~%")
                                              (format t "      0               20                    91,75%~%")
                                              (format t "-------------------------------------------------------------~%")
                                              (format t "      0               19                    85,57%~%")
                                              (format t "-------------------------------------------------------------~%")
                                              (format t "      0               18                    50,00%~%")
                                              (format t "-------------------------------------------------------------~%")
                                              (format t "      1               26                    97,87%~%")
                                              (format t "-------------------------------------------------------------~%")
                                              (format t "      1               25                    92,55%~%")
                                              (format t "-------------------------------------------------------------~%")
                                              (terpri)
                                              (SB-EXT:EXIT)
                                              )
                                             )
                                           (progn
                                            (terpri)
                                            (format t "Warning:~%")
                                            (format t "Invalid parameter: -aadifflimit <amino acid difference limit>~%~%")
                                            (format t "Possible configurations of the histogram parameters:~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "-aadifflimit   -aacheckminlimit   minimal protein identity %~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      0               26                     100%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      0               25                     100%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      0               24                    97,96%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      0               23                    96,94%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      0               22                    96,94%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      0               21                    94,68%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      0               20                    91,75%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      0               19                    85,57%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      0               18                    50,00%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      1               26                    97,87%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (format t "      1               25                    92,55%~%")
                                            (format t "-------------------------------------------------------------~%")
                                            (terpri)
                                            (SB-EXT:EXIT)
                                            )
                                           )
                        ))
                 );cond
        	      );if

             (if (or (string= parametertag "-pphistofilter") (string= parametertag "-ppcn")
                     (string= parametertag "-ppcomplete") (string= parametertag "-help")
                     (string= parametertag "-genefusion")
                     )
               (setf parameter (+ parameter 1))
               #|(if (char= (elt parametervalue 0) #\-)
                   (setf parameter (+ parameter 1))
                   (setf parameter (+ parameter 2))
                   )|#
               (if (and (string= parametertag "-dir") parametervalue)
                 (if (eql (pathname-directory (pathname parametervalue)) nil)
                   (if (char= (elt parametervalue 0) #\-)
                     (setf parameter (+ parameter 1))
                     (setf parameter (+ parameter 2))
                     )
                   (setf parameter (+ parameter 2))
                   )
                 (if (and (or (string= parametertag "-ppiterlimit") (string= parametertag "-trim")
                              (string= parametertag "-grouplimit")
                              )
                          parametervalue)
                   (if (char= (elt parametervalue 0) #\-)
                     (setf parameter (+ parameter 1))
                     (setf parameter (+ parameter 2))
                     )
                   (setf parameter (+ parameter 2))
                   );if or
                 ); if and
               );if or
          	  );dotimes

    (unless (eql help nil)
      (progn
       (terpri)
       (format t "Mandatory parameter~%~%")
       (format t "Directory parameter:~%")
       (format t "-dir <workingdir> directory~%~%")
       (format t "Optional parameters~%~%")
       (format t "Parameters of ppi by conserved neighbourhood:~%")
       (format t "-cnp <conserved neighbourhood score percentage> (65-default)~%")
       (format t "-expt <type of expansion of the conserved gene neighborhood> 'fixed' or 'dynamic' (fixed-default)~%~%")
       (format t "Parameters of conserved neighbourhood by fixed expansion:~%")
       (format t "-w1 <window size 1> (10-default)~%")
       (format t "-cw1 <gene conservation required for window 1> (4-default)~%")
       (format t "-w2 <window size 2> (7-default)~%")
       (format t "-cw2 <gene conservation required for window 2> (3-default)~%")
       (format t "-w3 <window size 3> (5-default)~%")
       (format t "-cw3 <gene conservation required for window 3> (2-default)~%")
       (format t "-w4 <window size 4> (3-default)~%")
       (format t "-cw4 <gene conservation required for window 4> (1-default)~%~%")
       (format t "Parameters of conserved neighbourhood by dynamic expansion:~%")
       (format t "-ws <dynamic expansion window size> (1-default)~%~%")
       (format t "Parameters of ppi by phylogenetic profile:~%")
       (format t "-ppp <phylogenetic profiles score percentage> (30-default)~%")
       (format t "-ppdifftolerated <difference in phylogenetic profiles tolerated to infer ppi> (0-default)~%~%")
       (format t "Amino acid histogram parameter settings for the phylogenetic profile:~%")
       (format t "-pphistofilter <build the phylogenetic profile of genes with a higher percentage of identity>~%")
       (format t "-ppaadifflimit <amino acid difference limit> (0-default)~%")
       (format t "-ppaacheckminlimit <minimum amount of amino acids to check> (26-default)~%~%")
       (format t "Methods of ppi prediction by phylogenetic profile and its parameters:~%~%")
       (format t "Method 1 - Predict ppi by phylogenetic profile only for interactions predicted by conserved neighborhood~%")
       (format t "-ppcn~%~%")
       (format t "Method 2 - PPi prediction by phylogenetic profile without filters 2~%")
       (format t "-ppcomplete~%~%")
       (format t "Method 3 - Prediction of ppi by phylogenetic profile with a limit of interactions~%")
       (format t "-ppiterlimit <maximum number of interactions desired> (500000-default)~%~%")
       (format t "Method 4 - Prediction of ppi by phylogenetic profile with interactions limit by weight~%")
       (format t "-trim <maximum number of interactions by weight> (45000-default)~%~%")
       (format t "Method 5 - Prediction of ppi by the phylogenetic profile only for genes with profiles that cover a greater or lesser number of genomes than an informed threshold~%")
       (format t "-threshold <phylogenetic profiles threshold>~%")
       (format t "-plusminus <parameter that receives the greater than or less than sign to apply the ppthreshod filter> '<' or '>' (signs greater than and less than, must be enclosed in single or double quotes)~%~%")
       (format t "Method 6 - Delete groups of ppi predicted by phylogenetic profile that exceed a limit of interactions by weight~%")
       (format t "-grouplimit <limit of tolerated interactions to maintain a group of ppi with the same weight> (45000-default)~%~%")
       (format t "Method 7 - To exclude genes with unwanted profiles in predicting ppi by phylogenetic profile~%")
       (format t (concatenate 'string "-profiles <number of genomes in the unwanted profiles> Entry example: 7 (genes that co-occur in a total of 7 genomes will be excluded). To insert more than one profile, the entry must be enclosed in single or double quotes, and the values separated by semicolons. Example: ""\"7; 15; 21\"~%~%"))
       (format t "Parameters of ppi by gene fusion:~%")
       (format t "-genefusion <make ppi predictions by gene fusion>~%")
       (format t "-gfp <gene fusion score percentage> (5-default)~%")
       (format t "Note: if gene fusion is included, there will be an increase of more than 100% in the execution time~%~%")
       (format t "Amino acid histogram parameters:~%")
       (format t "-aadifflimit <amino acid difference limit> (1-default)~%")
       (format t "-aacheckminlimit <minimum amount of amino acids to check> (25-default)~%~%")
       (format t "Possible configurations of the histogram parameters:~%")
       (format t "-------------------------------------------------------------~%")
       (format t "-aadifflimit   -aacheckminlimit   minimal protein identity %~%")
       (format t "-------------------------------------------------------------~%")
       (format t "      0               26                     100%~%")
       (format t "-------------------------------------------------------------~%")
       (format t "      0               25                     100%~%")
       (format t "-------------------------------------------------------------~%")
       (format t "      0               24                    97,96%~%")
       (format t "-------------------------------------------------------------~%")
       (format t "      0               23                    96,94%~%")
       (format t "-------------------------------------------------------------~%")
       (format t "      0               22                    96,94%~%")
       (format t "-------------------------------------------------------------~%")
       (format t "      0               21                    94,68%~%")
       (format t "-------------------------------------------------------------~%")
       (format t "      0               20                    91,75%~%")
       (format t "-------------------------------------------------------------~%")
       (format t "      0               19                    85,57%~%")
       (format t "-------------------------------------------------------------~%")
       (format t "      0               18                    50,00%~%")
       (format t "-------------------------------------------------------------~%")
       (format t "      1               26                    97,87%~%")
       (format t "-------------------------------------------------------------~%")
       (format t "      1               25                    92,55%~%")
       (format t "-------------------------------------------------------------~%")
       (terpri)
       (SB-EXT:EXIT)
       )
      )

    ;If you have missed the main data file parameter or the working directory then generate an error and abort the execution of the program
   	(if (not workingdir)
 	    (progn
       (terpri)
       (format t "Warning:~%")
 	     (format t"There are missing parameters: -dir <workingdir>")
       (terpri)
	      (SB-EXT:EXIT)
	      )
 	    )

   	;If there is a problem accessing the main data directory then abort and generate an error message
   	(if (not (probe-file workingdir))
 	    (progn
       (terpri)
       (format t "Warning:~%")
 	     (format t "~%~a~%" "Working directory is missing. Assuming current path")
	      ;(setf workingdir "./")
	      )
      (progn
       (dolist (path (list-directory workingdir))

               (if (eql (pathname-type (pathname path)) nil)
                 ;;If the path represents a directory, increment  the number of folders in thedirectory. 
                 (setf directory-num (incf directory-num))
                 ;;If the path is not from a folder, make:
                 (progn
                  (if (or  (equalp (pathname-type (pathname path)) "fasta")
                           (equalp (pathname-type (pathname path)) "fa")
                           (equalp (pathname-type (pathname path)) "faa")
                           );or
                    (incf qtd-arquivos-fasta)
                    );if or
                  );progn
                 );if
               );dolist
       (if (or (= directory-num (length (list-directory workingdir))) (= qtd-arquivos-fasta 0))
         ;;If there are only folders in the specified directory,  do:
         (progn
          (terpri)
          (format t "Warning:~%")
          (format t "There are no valid files in this directory.")
          (terpri)
          (format t "Protein sequence files must be in fasta format (.fasta, .faa or .fa).")
          (terpri)
          (SB-EXT:EXIT)
          )
         (progn
          (if (or (string= ppcn "Y") (string= ppcomplete "Y") (string= pplimit "Y") (string= pptrim "Y")
                  (string= ppthreshold "Y") (string= ppdeletegroup "Y") (string= ppdeleteprofile "Y")
                  );or
            (progn
             (cond
               ((string= ppthreshold "Y")
                (unless (> threshold 0)
                  (terpri)
                  (format t "Warning:~%")
                  (format t"When using the method 5, it is necessary to enter a value greater than zero through the -threshold parameter")
                  (terpri)
                  (SB-EXT:EXIT)
                  );unless threshold
                (unless (not (string= plusminus "nil"))
                  (terpri)
                  (format t "Warning:~%")
                  (format t"When using the method 5, it is necessary to insert the less than (<) or greater than (>) sign through the -plusminus parameter (signs greater than and less than, must be enclosed in single or double quotes)")
                  (terpri)
                  (SB-EXT:EXIT)
                  );unless plusminus
                )
               ((string= ppdeleteprofile "Y")
                (unless (not (string= profiles "0"))
                  (terpri)
                  (format t "Warning:~%")
                  (format t (concatenate 'string "When using the method 7, it is necessary to enter a value greater than zero through the -profiles parameter. Entry example: 7 (genes that co-occur in a total of 7 genomes will be excluded). ""To insert more than one profile, the entry must be enclosed in single or double quotes, and the values separated by semicolons. Example: ""\"7; 15; 21\"~%~%"))
                  (terpri)
                  (SB-EXT:EXIT)
                  );unless
                )
               );cond
             );progn
            );if or

          (terpri)
          (format t "GENPPI VERSION: 1.1~%")
          (format t "RELEASE NUMBER: 55e436b354ba1f3575529aaccc51002e7bf68ee6~%")
          (format t "REPOSITORY: https://github.com/santosardr/genppi ~%~%")
          (format t "Directory parameter:~%")
          (format t "-dir <workingdir> = ~a~%~%" workingdir)
          (format t "Parameters of ppi by conserved neighbourhood:~%")

          (cond
            ((and (string= gene-fusion "Y") (eql aux-percentage-cn nil))
             (format t "-cnp <conserved neighbourhood score percentage> = 65% (default)~%")
             )
            ((and (string= gene-fusion "Y") (not (eql aux-percentage-cn nil)))
             (format t "-cnp <conserved neighbourhood score percentage> = ~a%~%" (* percentage-cn 100))
             )
            ((and (string= gene-fusion "N")
                  (eql aux-percentage-cn nil)
                  (or (string= ppcn "Y") (string= ppcomplete "Y") (string= pplimit "Y") (string= pptrim "Y")
                      (string= ppthreshold "Y") (string= ppdeletegroup "Y") (string= ppdeleteprofile "Y")
                      );or
                  (> qtd-arquivos-fasta 1)
                  );and
                   (format t "-cnp <conserved neighbourhood score percentage> = 70%~%")
                   (setf percentage-cn (/ 70 100))
                   )
            ((and (string= gene-fusion "N")
                  (eql aux-percentage-cn nil)
                  (or (string= ppcn "Y") (string= ppcomplete "Y") (string= pplimit "Y") (string= pptrim "Y")
                      (string= ppthreshold "Y") (string= ppdeletegroup "Y") (string= ppdeleteprofile "Y")
                      );or
                  (= qtd-arquivos-fasta 1)
                  );and
                   (format t "-cnp <conserved neighbourhood score percentage> = 100%~%")
                   )
            ((and (string= gene-fusion "N") (not (eql aux-percentage-cn nil)))
             (format t "-cnp <conserved neighbourhood score percentage> = ~a%~%" (* percentage-cn 100))
             )
            ((and (string= gene-fusion "N")
                  (and (string= ppcn "N") (string= ppcomplete "N") (string= pplimit "N") (string= pptrim "N")
                       (string= ppthreshold "N") (string= ppdeletegroup "N") (string= ppdeleteprofile "N")
                       );and
                  );and
                   (format t "-cnp <conserved neighbourhood score percentage> = 100%~%")
                   )
            );cond

          (if (string= expansion-type "fixed")
            (progn
             (format t "-expt <type of expansion of the conserved gene neighborhood> = ~a~%" "fixed (default)")
             (format t "-w1 <window size 1> = ~a~%" (if (= w1 10) (string "10 (default)") (- w1 0)))
             (format t "-cw1 <gene conservation required for window 1> = ~a~%" (if (= cw1 4) (string "4 (default)") (- cw1 0)))
             (format t "-w2 <window size 2> = ~a~%" (if (= w2 7) (string "7 (default)") (- w2 0)))
             (format t "-cw2 <gene conservation required for window 2> = ~a~%" (if (= cw2 3) (string "3 (default)") (- cw2 0)))
             (format t "-w3 <window size 3> = ~a~%" (if (= w3 5) (string "5 (default)") (- w3 0)))
             (format t "-cw3 <gene conservation required for window 3> = ~a~%" (if (= cw3 2) (string "2 (default)") (- cw3 0)))
             (format t "-w4 <window size 4> = ~a~%" (if (= w4 3) (string "3 (default)") (- w4 0)))
             (format t "-cw4 <gene conservation required for window 4> = ~a~%~%" (if (= cw4 1) (string "1 (default)") (- cw4 0)))

             (if (> w2 w1)
               (progn
                (terpri)
                (format t "Warning:~%")
                (format t"The expansion window size w2 cannot be larger than the maximum expansion (w1)")
                (terpri)
                (SB-EXT:EXIT)
                );progn
               );if
             (if (> w3 w1)
               (progn
                (terpri)
                (format t "Warning:~%")
                (format t"The expansion window size w3 cannot be larger than the maximum expansion (w1)")
                (terpri)
                (SB-EXT:EXIT)
                );progn
               );if
             (if (> w4 w1)
               (progn
                (terpri)
                (format t "Warning:~%")
                (format t"The expansion window size w4 cannot be larger than the maximum expansion (w1)")
                (terpri)
                (SB-EXT:EXIT)
                );progn
               );if
             (if (> cw1 w1)
               (progn
                (terpri)
                (format t "Warning:~%")
                (format t"The cw1 parameter cannot be larger than the window size w1")
                (terpri)
                (SB-EXT:EXIT)
                );progn
               );if
             (if (> cw2 w2)
               (progn
                (terpri)
                (format t "Warning:~%")
                (format t"The cw2 parameter cannot be larger than the window size w2")
                (terpri)
                (SB-EXT:EXIT)
                );progn
               );if
             (if (> cw3 w3)
               (progn
                (terpri)
                (format t "Warning:~%")
                (format t"The cw3 parameter cannot be larger than the window size w3")
                (terpri)
                (SB-EXT:EXIT)
                );progn
               );if
             (if (> cw4 w4)
               (progn
                (terpri)
                (format t "Warning:~%")
                (format t"The cw4 parameter cannot be larger than the window size w4")
                (terpri)
                (SB-EXT:EXIT)
                );progn
               );if

             );Prong
            (progn
             (format t "-expt <type of expansion of the conserved gene neighborhood> = ~a~%" expansion-type)
             (format t "-ws <dynamic expansion window size> = ~a~%~%" (if (= window-size 1)(string "1 (default)")(- window-size 0)))
             )
            );if

          (if (or (string= ppcn "Y") (string= ppcomplete "Y") (string= pplimit "Y") (string= pptrim "Y")
                  (string= ppthreshold "Y") (string= ppdeletegroup "Y") (string= ppdeleteprofile "Y")
                  );or
            (progn
             (cond
               ((string= ppcn "Y")
                (format t "Parameters of ppi by phylogenetic profile:~%")
                (format t "Method 1 - Predict ppi by phylogenetic profile only for interactions predicted by conserved neighborhood~%")
                )
               ((string= ppcomplete "Y")
                (format t "Parameters of ppi by phylogenetic profile:~%")
                (format t "Method 2 - PPi prediction by phylogenetic profile without filters~%")
                )
               ((string= pplimit "Y")
                (format t "Parameters of ppi by phylogenetic profile:~%")
                (format t "Method 3 - Prediction of ppi by phylogenetic profile with a limit of interactions~%")
                (if (= ppiterlimit 500000)
                  (format t "-ppiterlimit <maximum number of interactions desired> = 500000 (default)~%")
                  (format t "-ppiterlimit <maximum number of interactions desired> = ~a~%"ppiterlimit)
                  )
                )
               ((string= pptrim "Y")
                (format t "Parameters of ppi by phylogenetic profile:~%")
                (format t "Method 4 - Prediction of ppi by phylogenetic profile with interactions limit by weight~%")
                (if (= trim 45000)
                  (format t "-trim <maximum number of interactions by weight> = 45000 (default)~%")
                  (format t "-trim <maximum number of interactions by weight> = ~a~%"trim)
                  )
                )
               ((string= ppthreshold "Y")
                (format t "Parameters of ppi by phylogenetic profile:~%")
                (format t "Method 5 - Prediction of ppi by the phylogenetic profile only for genes with profiles that cover a greater or lesser number of genomes than an informed threshold~%")
                (format t "-threshold <phylogenetic profiles threshold> = ~a~%"threshold)
                (format t (concatenate 'string "-plusminus <only genes that occurred in " (if (string= plusminus "<") (string "less ") (string "more ")) "than ~a genomes will be considered in the analysis>~%") threshold)
                )
               ((string= ppdeletegroup "Y")
                (format t "Parameters of ppi by phylogenetic profile:~%")
                (format t "Method 6 - Delete groups of ppi predicted by phylogenetic profile that exceed a limit of interactions by weight~%")
                (if (= grouplimit 45000)
                  (format t "-grouplimit <limit of tolerated interactions to maintain a group of ppi with the same weight> = 45000 (default)~%")
                  (format t "-grouplimit <limit of tolerated interactions to maintain a group of ppi with the same weight> = ~a~%"grouplimit)
                  )
                )
               ((string= ppdeleteprofile "Y")
                (format t "Parameters of ppi by phylogenetic profile:~%")
                (format t "Method 7 - To exclude genes with unwanted profiles in predicting ppi by phylogenetic profile~%")
                (format t "-profiles <number of genomes in the unwanted profiles> Entry = ~a (genes that co-occur in a total of ~a genomes will be excluded)~%" profiles profiles)
                )
               );cond

             (cond
               ((eql aux-percentage-pp nil)
                (format t "-ppp <phylogenetic profiles score percentage> = 30% (default)~%")
                )
               ((not (eql aux-percentage-pp nil))
                (format t "-ppp <phylogenetic profiles score percentage> = ~a%~%" (* percentage-pp 100))
                )
               );cond

             (if (= ppdifftolerated 0)
               (format t "-ppdifftolerated <difference in phylogenetic profiles tolerated to infer ppi> = 0 (default)~%")
               (format t "-ppdifftolerated <difference in phylogenetic profiles tolerated to infer ppi> = ~a~%" ppdifftolerated)
               )
             (if (string= pphistofilter "Y")
               (progn
                (format t "-pphistofilter <build the phylogenetic profile of genes with a higher percentage of identity>~%")
                (if (= ppaadifflimit 0)
                  (format t "-ppaadifflimit <amino acid difference limit> = 0 (default)~%")
                  (format t "-ppaadifflimit <amino acid difference limit> = ~a~%" ppaadifflimit)
                  )
                (if (= ppaacheckminlimit 26)
                  (format t "-ppaacheckminlimit <minimum amount of amino acids to check> = 26 (default)~%~%")
                  (format t "-ppaacheckminlimit <minimum amount of amino acids to check> = ~a~%~%" ppaacheckminlimit)
                  )
                );progn
               (format t "~%")
               );if pphistofilter
             );progn
            );if

          (cond
            ((and (string= gene-fusion "Y") (= (float percentage-gf) 0.05))
             (format t "Parameters of ppi by gene fusion:~%")
             (format t "-gfp <gene fusion score percentage> = 5% (default)~%~%")
             )
            ((and (string= gene-fusion "Y") (not (= (float percentage-gf) 0.05)))
             (format t "Parameters of ppi by gene fusion:~%")
             (format t "-gfp <gene fusion score percentage> = ~a~%~%" (* percentage-gf 100))
             )
            );cond

          (format t "Amino acid histogram parameters:~%")
          (if (= aadifflimit 1)
            (format t "-aadifflimit <amino acid difference limit> 1 (default)~%")
            (format t "-aadifflimit <amino acid difference limit> ~a~%" aadifflimit)
            )
          (if (= aacheckminlimit 25)
            (format t "-aacheckminlimit <minimum amount of amino acids to check> 25 (default)~%")
            (format t "-aacheckminlimit <minimum amount of amino acids to check> ~a~%" aacheckminlimit)
            )
          (terpri)
          (terpri)
          (format t "Making ppi prediction, please wait.~%")
          (terpri)
          (genppi
           workingdir
           gene-fusion
           (float percentage-gf)
           (float percentage-pp)
           ppdifftolerated
           pphistofilter
           ppaadifflimit
           ppaacheckminlimit
           ppcn
           ppcomplete
           pplimit
           ppiterlimit
           pptrim
           trim
           ppthreshold
           threshold
           plusminus
           ppdeletegroup
           grouplimit
           ppdeleteprofile
           profiles
           expansion-type
           window-size
           (float percentage-cn)
           w1
           cw1
           w2
           cw2
           w3
           cw3
           w4
           cw4
           aadifflimit
           aacheckminlimit
           );parameters genPPI
          (terpri)
          );progn
         );if
       );progn
 	    );if
   	);let
  );defun
;;------------------------------------------------------------------------------
