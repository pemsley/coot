;;;; Copyright 2008 by The University of Oxford

;;;; This program is free software; you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation; either version 3 of the License, or (at
;;;; your option) any later version.
 
;;;; This program is distributed in the hope that it will be useful, but
;;;; WITHOUT ANY WARRANTY; without even the implied warranty of
;;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;;;; General Public License for more details.
 
;;;; You should have received a copy of the GNU General Public License
;;;; along with this program; if not, write to the Free Software
;;;; Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
;;;; 02110-1301, USA

(use-modules (sxml simple))

;; (message-file-name "loopy.msg")
(define (write-arp-warp-xml-file imol-map XMLOutputFilename MessageFilename InputPDB ccp4i-project-dir)

  (define (cell-string imol-map)
    (let ((cell (map-cell imol-map)))
      (if (not cell)
	  #f
	  (string-append-with-spaces (map number->string cell)))))

  (let ((arp-warp-bin (getenv "xx"))
	(arp-warp-prefix arp-warp-bin))
    (if (not arp-warp-bin)
	(begin
	  (format #t "arp/warp bin not set aborting now~%")
	  (info-dialog "arp/warp bin not set - aborting now"))

	(let* ((map-params (map-parameters imol)))
	  (if map-parameters
	      (let ((InputMTZ (car map-parameters))
		    (FWTLabel (cadr map-parameters)
			      (PHIWTLabel (caddr map-parameters))
			      (pdb-file-prefix (strip-path (file-name-sans-extension InputPDB)))
			      ;; (XtalCell "64.897 78.323 38.792 90.00 90.00 90.00")
			      (XtalCell (cell-string imol-map))
			      (DictionaryFilename (string-append arp/warp-prefix "/AAbis_extension.XYZ"))
			      (SymmetryFilename (string-append arp/warp-prefix "/syminfo.lib"))
			      (StructureFilenameToC (string-append arp/warp-prefix "/set_1200_2C_my_0.3b_cos_log.llh"))
			      (StructureFilenameToN (string-append "/set_1200_2N_my_0.3b_cos_log.llh"))
			      (OutputMapFromMTZ "/tmp/paule/rnasa-1.8-all_loopy1.map")
			      (SaveLoopsDir ccp4i-project-dir)
			      (SaveLoopsBasename (string-append pdb-file-prefix "-loop"))
			      (SaveProposedPDBDir ccp4i-project-dir)
			      (SaveProposedPDBBasename (string-append "-proposed")))
		    
		    (call-with-output-file "loopy.xml"
		      (lambda (port)
			(let ((expr 
			       `(arp-warp
				 (loop
				  (MessageFilename ,message-file-name) "\n"
				  (XMLOutputFilename ,XMLOutputFilename) "\n"
				  (AbortLevel 8) "\n"
				  (MessageLevel 4) "\n"
				  (ProgramName loopy) "\n"
				  (DictionaryFilename ,DictionaryFilename) "\n"
				  (UniformAtomBFactor 20) "\n"
				  (UniformAtomRadius 0.74) "\n"
				  (AntiBumpFactor 1) "\n"
				  (SymmetryFilename ,SymmetryFilename) "\n"
				  (SpaceGroup 19) "\n"
				  (XtalCell ,XtalCell) "\n"
				  (KeepSideChain 1) "\n"
				  (KeepFragmentLongerThan 0) "\n"
				  (UseCDummyResidue 1) "\n"
				  (StructureFilenameToC ,StructureFilenameToC)  "\n"
				  (StructureFilenameToN ,StructureFilenameToN) "\n"
				  (CADistance 3.8) "\n"
				  (RandomGrid 0) "\n"
				  (GridNumber 378) "\n"
				  (DisregardNegativeDensity 1) "\n"
				  (KeepNegDensityHalfway 0) "\n"
				  (NormalizeDensity 1) "\n"
				  (BFactorForBuildSideChain 26) "\n"
				  (LikelihoodThreshold -4.5) "\n"
				  (WeightDensityLLH 1.) "\n"
				  (WeightDistanceLLH 1.) "\n"
				  (WeightStructuralLLH 1.) "\n"
				  (CheckDistanceDuringBuild yes) "\n"
				  (MaxDistanceBetweenAnchors 20.) "\n"
				  (OldComputation 0) "\n"
				  (CheckFirstAngle 1) "\n"
				  (UseCubicFitTarget 1) "\n"
				  (OverlapRejectionThreshold 0.1) "\n"
				  (DummyRejectionThreshold 0.1) "\n"
				  (LoopBuildAll 0) "\n"
				  (ExtendGapSmallerThan 0) "\n"
				  (LoopBothWays 1) "\n"
				  (LoopToC 1) "\n"
				  ;; LoopLength
				  (LoopOverlap 0) "\n"
				  (LoopDensityCutoffNo 100) "\n"
				  (MinimalDistance 0.3) "\n"
				  (ForceMinCAsKept 1) "\n"
				  (LoopRMS 0.7) "\n"
				  (LoopStructureThreshold -1.5) "\n"
				  (LoopStructureCutoffNo 50) "\n"
				  (LoopStructureMinNo 15) "\n"
				  (LoopMainChainDensNo 5) "\n"
				  (LoopsToBuild "A39&lt;5&gt;A43") "\n"
				  (LoopSequence NRESV) "\n"
				  (UseLoopSequence 1) "\n"
				  (InputPDB ,InputPDB) "\n"
				  ;;(IncludeChains)
				  (InputMTZ ,InputMTZ) "\n"
				  (OutputMapFromMTZ ,OutputMapFromMTZ) "\n"
				  (FWTLabel ,FWTLabel) "\n"
				  (PHIWTLabel ,PHIWTLabel) "\n"
				  (SaveResult 0) "\n"
				  (SaveBestLoopSets 1) "\n"
				  (SaveBestNumber 3) "\n"
				  (SaveLoopsDir ,SaveLoopsDir) "\n"
				  (SaveLoopsBasename ,SaveLoopsBasename) "\n"
				  (SaveProposedPDBDir ,SaveProposedPDBDir) "\n"
				  (SaveProposedPDBBasename SaveProposedPDBBasename) "\n"
				  (SaveUseLoopDef 0) "\n"
				  (WeightMainChain 1.) "\n"
				  (WeightSideChain 0.) "\n"))))
			  (sxml->xml expr port)))))))))))
