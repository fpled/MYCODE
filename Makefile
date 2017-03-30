# Username
User_Name_mac = Op
User_NameLinux = pled
# Machine name
Machine_Name_mac = MacBookProOp
Machine_Name_pcmsme = pcmsmemeca52
# Directory
Dir_mac = /Users/$(User_Name_mac)/Documents/Recherche/GeM
Dir_pcmsme = /home/p/$(User_NameLinux)/Documents
Dir_cluster = /home/p/pled/$(User_NameLinux)/Documents
Dir_qsub = /usr/local/pbs/default/bin
# Cluster name
Cluster_Name = cluster
# Problem name
#Pb_Name = multiscaleDetLinDiffForm
Pb_Name = multiscaleDetLinDiff
#Pb_Name = multiscaleDetLinDiffCircHoles
#Pb_Name = multiscaleDetLinElas
#Pb_Name = multiscaleDetLinElasCircHoles
#Pb_Name = multiscaleDetLinElasEdgeCracks
#Pb_Name = multiscaleDetLinElasInteriorCracks
#Pb_Name = multiscaleDetNonlinDiffReac

#Pb_Name = multiscaleStoLinDiff
#Pb_Name = multiscaleStoLinElas
#Pb_Name = multiscaleStoNonlinDiffReac

#Pb_Name = multiscaleStoLinDiffAlignInclusionsIso
#Pb_Name = multiscaleStoLinDiffAlignInclusionsAniso
#Pb_Name = multiscaleStoLinDiffCircInclusionsIso
#Pb_Name = multiscaleStoLinDiffCircInclusionsAniso
#Pb_Name = multiscaleStoLinDiffSquareInclusionsIso
#Pb_Name = multiscaleStoLinDiffSquareInclusionsAniso
#Pb_Name = multiscaleStoNonlinDiffReacAlignInclusionsIso
#Pb_Name = multiscaleStoNonlinDiffReacAlignInclusionsAniso
#Pb_Name = multiscaleStoNonlinDiffReacSquareInclusionsIso
#Pb_Name = multiscaleStoNonlinDiffReacSquareInclusionsAniso

#Pb_Name = monoscaleStoLinDiff
#Pb_Name = monoscaleStoLinDiffCircInclusionsIso
#Pb_Name = monoscaleStoLinDiffCircInclusionsAniso
#Pb_Name = monoscaleStoNonlinDiffReac

#Pb_Name = sparseAnisotropicFunction
#Pb_Name = sparseIshigamiFunction
#Pb_Name = sparsePolynomialFunction
#Pb_Name = sparseSobolFunction
#Pb_Name = sparseGeometricBrownian

#Pb_Name = plateCircDetLinElas
#Pb_Name = plateCircDetLinElasCv
#Pb_Name = plateRectDetLinElas
#Pb_Name = plateRectDetLinElasCv
#Pb_Name = tableCircDetLinElas
#Pb_Name = tableRectDetLinElas
#Pb_Name = FCBATableCircDetLinElas
#Pb_Name = FCBATableCircStoLinElas
#Pb_Name = FCBADeskDetLinElas
#Pb_Name = FCBADeskStoLinElas

# Default ---------------------------
default:
#	matlab -nodesktop -nosplash -noFigureWindows -r "$(Pb_Name); exit" -logfile results/$(Pb_Name).log
	matlab -nodesktop -nosplash -noFigureWindows -r "plateCircDetLinElasCv; exit" -logfile results/plate/plateCircDetLinElasCv.log
	matlab -nodesktop -nosplash -noFigureWindows -r "plateRectDetLinElasCv; exit" -logfile results/plate/plateRectDetLinElasCv.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleStoLinDiffAlignInclusionsIso; exit" -logfile results/multiscaleSto/linDiffAlignInclusionsIso.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleStoLinDiffAlignInclusionsAniso; exit" -logfile results/multiscaleSto/linDiffAlignInclusionsAniso.log

# All ---------------------------
all: sparse monosto multidet multisto plate

# Multiscale ---------------------------
multi: multidet multisto

# Multiscale Deterministic ---------------------------
multidet:
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleDetLinDiffForm; exit" -logfile results/multiscaleDet/linDiffForm.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleDetLinDiff; exit" -logfile results/multiscaleDet/linDiff.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleDetLinDiffCircHoles; exit" -logfile results/multiscaleDet/linDiffCircHoles.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleDetLinElas; exit" -logfile results/multiscaleDet/linElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleDetLinElasCircHoles; exit" -logfile results/multiscaleDet/linElasCircHoles.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleDetLinElasEdgeCracks; exit" -logfile results/multiscaleDet/linElasEdgeCracks.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleDetLinElasInteriorCracks; exit" -logfile results/multiscaleDet/linElasInteriorCracks.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleDetNonlinDiffReac; exit" -logfile results/multiscaleDet/nonlinDiffReac.log

# Multiscale Stochastic ---------------------------
multisto:
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleStoLinDiff; exit" -logfile results/multiscaleSto/linDiff.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleStoLinElas; exit" -logfile results/multiscaleSto/linElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleStoLinDiffAlignInclusionsIso; exit" -logfile results/multiscaleSto/linDiffAlignInclusionsIso.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleStoLinDiffAlignInclusionsAniso; exit" -logfile results/multiscaleSto/linDiffAlignInclusionsAniso.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleStoLinDiffCircInclusionsIso; exit" -logfile results/multiscaleSto/linDiffCircInclusionsIso.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleStoLinDiffCircInclusionsAniso; exit" -logfile results/multiscaleSto/linDiffCircInclusionsAniso.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleStoLinDiffSquareInclusionsIso; exit" -logfile results/multiscaleSto/linDiffSquareInclusionsIso.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleStoLinDiffSquareInclusionsAniso; exit" -logfile results/multiscaleSto/linDiffSquareInclusionsAniso.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleStoNonlinDiffReac; exit" -logfile results/multiscaleSto/nonlinDiffReac.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleStoNonlinDiffReacAlignInclusionsIso; exit" -logfile results/multiscaleSto/nonlinDiffReacAlignInclusionsIso.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleStoNonlinDiffReacAlignInclusionsAniso; exit" -logfile results/multiscaleSto/nonlinDiffReacAlignInclusionsAniso.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleStoNonlinDiffReacSquareInclusionsIso; exit" -logfile results/multiscaleSto/nonlinDiffReacSquareInclusionsIso.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleStoNonlinDiffReacSquareInclusionsAniso; exit" -logfile results/multiscaleSto/nonlinDiffReacSquareInclusionsAniso.log

# Monoscale Stochastic ---------------------------
monosto:
	matlab -nodesktop -nosplash -noFigureWindows -r "monoscaleStoLinDiff; exit" -logfile results/monoscaleSto/linDiff.log
	matlab -nodesktop -nosplash -noFigureWindows -r "monoscaleStoLinDiffCircInclusionsIso; exit" -logfile results/monoscaleSto/linDiffCircInclusionsIso.log
	matlab -nodesktop -nosplash -noFigureWindows -r "monoscaleStoLinDiffCircInclusionsAniso; exit" -logfile results/monoscaleSto/linDiffCircInclusionsAniso.log
	matlab -nodesktop -nosplash -noFigureWindows -r "monoscaleStoNonlinDiffReac; exit" -logfile results/monoscaleSto/nonlinDiffReac.log

# Sparse approximation ---------------------------
sparse:
	matlab -nodesktop -nosplash -noFigureWindows -r "sparseAnisotropicFunction; exit" -logfile results/sparse/anisotropicFunction.log
	matlab -nodesktop -nosplash -noFigureWindows -r "sparseIshigamiFunction; exit" -logfile results/sparse/ishigamiFunction.log
	matlab -nodesktop -nosplash -noFigureWindows -r "sparsePolynomialFunction; exit" -logfile results/sparse/polynomialFunction.log
	matlab -nodesktop -nosplash -noFigureWindows -r "sparseSobolFunction; exit" -logfile results/sparse/sobolFunction.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "sparseGeometricBrownian; exit" -logfile results/sparse/geometricBrownian.log

# Plate Shell ---------------------------
plate:
	matlab -nodesktop -nosplash -noFigureWindows -r "plateCircDetLinElas; exit" -logfile results/plate/plateCircDetLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "plateCircDetLinElasCv; exit" -logfile results/plate/plateCircDetLinElasCv.log
	matlab -nodesktop -nosplash -noFigureWindows -r "plateRectDetLinElas; exit" -logfile results/plate/plateRectDetLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "plateRectDetLinElasCv; exit" -logfile results/plate/plateRectDetLinElasCv.log
	matlab -nodesktop -nosplash -noFigureWindows -r "tableCircDetLinElas; exit" -logfile results/plate/tableCircDetLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "tableRectDetLinElas; exit" -logfile results/plate/tableRectDetLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "FCBATableCircDetLinElas; exit" -logfile results/plate/FCBATableCircDetLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "FCBATableCircStoLinElas; exit" -logfile results/plate/FCBATableCircStoLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskDetLinElas; exit" -logfile results/plate/FCBADeskDetLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskStoLinElas; exit" -logfile results/plate/FCBADeskStoLinElas.log

# Clean ---------------------------
clean:
	-rm -r results/*

remove_old:
	-find ../ -name '*.m~' -exec rm -fr {} \;

# Sync ---------------------------
sync: sync_from_pcmsme_to_mac sync_from_mac_to_pcmsme

sync_cluster: sync_from_cluster_to_mac sync_from_mac_to_cluster

sync_all: sync_from_pcmsme_to_mac sync_from_mac_to_pcmsme sync_from_cluster_to_mac sync_from_mac_to_cluster

sync_results: sync_results_from_pcmsme_to_mac sync_results_from_mac_to_pcmsme

# excludes for rsync
exclude = --exclude 'links.m' --exclude 'initfemobjectoptions.m' --exclude '.svn' --exclude '.DS_Store' --exclude '.directory' --exclude '*.mat' --exclude '*.geo' --exclude '*.msh' --exclude '*.err' --exclude '*.out' --exclude '*.log' --exclude '*.err' --exclude '*.jpeg' --exclude '*.jpg' --exclude '*.png' --exclude '*.tiff' --exclude '*.eps' --exclude '*.dvi' --exclude '*.pdf' --exclude '*.fig' --exclude '*.tex' --exclude '*.tsv' --exclude '*.avi' --exclude 'spams-matlab'

sync_from_mac_to_pcmsme:
	rsync -auv $(Dir_mac)/FEMObject/MYCODE $(User_NameLinux)@$(Machine_Name_pcmsme):$(Dir_pcmsme)/FEMObject $(exclude)

sync_from_pcmsme_to_mac:
	rsync -auv $(User_NameLinux)@$(Machine_Name_pcmsme):$(Dir_pcmsme)/FEMObject/MYCODE $(Dir_mac)/FEMObject $(exclude)

sync_from_mac_to_cluster:
	rsync -auv $(Dir_mac)/FEMObject/MYCODE $(User_NameLinux)@$(Cluster_Name):$(Dir_cluster)/FEMObject $(exclude)

sync_from_cluster_to_mac:
	rsync -auv $(User_NameLinux)@$(Cluster_Name):$(Dir_cluster)/FEMObject/MYCODE $(Dir_mac)/FEMObject $(exclude)

sync_results_from_mac_to_pcmsme:
	rsync -auv $(Dir_mac)/FEMObject/MYCODE/results $(User_NameLinux)@$(Machine_Name_pcmsme):$(Dir_pcmsme)/FEMObject/MYCODE

sync_results_from_pcmsme_to_mac:
	rsync -auv $(User_NameLinux)@$(Machine_Name_pcmsme):$(Dir_pcmsme)/FEMObject/MYCODE/results $(Dir_mac)/FEMObject/MYCODE

# Cluster ---------------------------
cluster:
	ssh -YC $(User_NameLinux)@$(Cluster_Name) "hostname; cd $(Dir_cluster)/FEMObject/MYCODE/examples/multiscale; qsub $(Pb_Name).pbs"

interactive:
	ssh -YC $(User_NameLinux)@$(Cluster_Name) "hostname; cd $(Dir_cluster)/FEMObject/MYCODE/examples/multiscale; qsub -I -X $(Pb_Name).pbs"
