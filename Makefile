# Username
User_Name_mac = Op
User_Name_linux = pled
# Machine name
Machine_Name_mac = MacBookProOp
Machine_Name_pcmsme = pcmsmemeca52
# Directory
Dir_mac = /Users/$(User_Name_mac)/Documents/Recherche/GeM
Dir_pcmsme = /home/p/$(User_Name_linux)/Documents
Dir_cluster = /home/p/pled/$(User_Name_linux)/Documents
Dir_qsub = /usr/local/pbs/default/bin
# Cluster name
Cluster_Name = cluster
# Problem name
#Pb_Name = multiscale_det_lin_diff_form
#Pb_Name = multiscale_det_lin_diff
#Pb_Name = multiscale_det_lin_diff_circ_holes
#Pb_Name = multiscale_det_lin_elas
#Pb_Name = multiscale_det_lin_elas_circ_holes
#Pb_Name = multiscale_det_lin_elas_edge_cracks
#Pb_Name = multiscale_det_lin_elas_interior_cracks
#Pb_Name = multiscale_det_nonlin_diff_reac

#Pb_Name = multiscale_sto_lin_diff
#Pb_Name = multiscale_sto_lin_elas
#Pb_Name = multiscale_sto_nonlin_diff_reac

#Pb_Name = multiscale_sto_lin_diff_align_inclusions_iso
#Pb_Name = multiscale_sto_lin_diff_align_inclusions_aniso
#Pb_Name = multiscale_sto_lin_diff_circ_inclusions_iso
#Pb_Name = multiscale_sto_lin_diff_circ_inclusions_aniso
#Pb_Name = multiscale_sto_lin_diff_square_inclusions_iso
#Pb_Name = multiscale_sto_lin_diff_square_inclusions_aniso
#Pb_Name = multiscale_sto_nonlin_diff_reac_align_inclusions_iso
#Pb_Name = multiscale_sto_nonlin_diff_reac_align_inclusions_aniso
#Pb_Name = multiscale_sto_nonlin_diff_reac_square_inclusions_iso
#Pb_Name = multiscale_sto_nonlin_diff_reac_square_inclusions_aniso

#Pb_Name = monoscale_sto_lin_diff
#Pb_Name = monoscale_sto_lin_diff_circ_inclusions_iso
#Pb_Name = monoscale_sto_lin_diff_circ_inclusions_aniso
#Pb_Name = monoscale_sto_nonlin_diff_reac

#Pb_Name = sparse_anisotropic_function
#Pb_Name = sparse_ishigami_function
#Pb_Name = sparse_polynomial_function
#Pb_Name = sparse_sobol_function
#Pb_Name = sparse_geometric_brownian

#Pb_Name = plate_circ_det_lin_elas
#Pb_Name = plate_circ_det_lin_elas_cv
#Pb_Name = plate_rect_det_lin_elas
#Pb_Name = plate_rect_det_lin_elas_cv
#Pb_Name = table_circ_det_lin_elas
#Pb_Name = table_rect_det_lin_elas
#Pb_Name = FCBA_table_circ_det_lin_elas
Pb_Name = FCBA_table_circ_sto_lin_elas

# Default ---------------------------
default:
	matlab -nodesktop -nosplash -noFigureWindows -r "$(Pb_Name); exit" -logfile results/$(Pb_Name).log
#	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_sto_lin_diff_align_inclusions_iso; exit" -logfile results/multiscale_sto_lin_diff_align_inclusions_iso.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_sto_lin_diff_align_inclusions_aniso; exit" -logfile results/multiscale_sto_lin_diff_align_inclusions_aniso.log

# All ---------------------------
all: sparse monosto multidet multisto plate

# Multiscale ---------------------------
multi: multidet multisto

# Multiscale Deterministic ---------------------------
multidet:
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_det_lin_diff_form; exit" -logfile results/multiscale_det_lin_diff_form.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_det_lin_diff; exit" -logfile results/multiscale_det_lin_diff.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_det_lin_diff_circ_holes; exit" -logfile results/multiscale_det_lin_diff_circ_holes.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_det_lin_elas; exit" -logfile results/multiscale_det_lin_elas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_det_lin_elas_circ_holes; exit" -logfile results/multiscale_det_lin_elas_circ_holes.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_det_lin_elas_edge_cracks; exit" -logfile results/multiscale_det_lin_elas_edge_cracks.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_det_lin_elas_interior_cracks; exit" -logfile results/multiscale_det_lin_elas_interior_cracks.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_det_nonlin_diff_reac; exit" -logfile results/multiscale_det_nonlin_diff_reac.log

# Multiscale Stochastic ---------------------------
multisto:
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_sto_lin_diff; exit" -logfile results/multiscale_sto_lin_diff.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_sto_lin_elas; exit" -logfile results/multiscale_sto_lin_elas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_sto_nonlin_diff_reac; exit" -logfile results/multiscale_sto_nonlin_diff_reac.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_sto_lin_diff_align_inclusions_iso; exit" -logfile results/multiscale_sto_lin_diff_align_inclusions_iso.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_sto_lin_diff_align_inclusions_aniso; exit" -logfile results/multiscale_sto_lin_diff_align_inclusions_aniso.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_sto_nonlin_diff_reac_align_inclusions_iso; exit" -logfile results/multiscale_sto_nonlin_diff_reac_align_inclusions_iso.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_sto_nonlin_diff_reac_align_inclusions_aniso; exit" -logfile results/multiscale_sto_nonlin_diff_reac_align_inclusions_aniso.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_sto_lin_diff_circ_inclusions_iso; exit" -logfile results/multiscale_sto_lin_diff_circ_inclusions_iso.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_sto_lin_diff_circ_inclusions_aniso; exit" -logfile results/multiscale_sto_lin_diff_circ_inclusions_aniso.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_sto_lin_diff_square_inclusions_iso; exit" -logfile results/multiscale_sto_lin_diff_square_inclusions_iso.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_sto_lin_diff_square_inclusions_aniso; exit" -logfile results/multiscale_sto_lin_diff_square_inclusions_aniso.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_sto_nonlin_diff_reac_square_inclusions_iso; exit" -logfile results/multiscale_sto_nonlin_diff_reac_square_inclusions_iso.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscale_sto_nonlin_diff_reac_square_inclusions_aniso; exit" -logfile results/multiscale_sto_nonlin_diff_reac_square_inclusions_aniso.log

# Monoscale Stochastic ---------------------------
monosto:
	matlab -nodesktop -nosplash -noFigureWindows -r "monoscale_sto_lin_diff; exit" -logfile results/monoscale_sto_lin_diff.log
	matlab -nodesktop -nosplash -noFigureWindows -r "monoscale_sto_lin_diff_circ_inclusions_iso; exit" -logfile results/monoscale_sto_lin_diff_circ_inclusions_iso.log
	matlab -nodesktop -nosplash -noFigureWindows -r "monoscale_sto_lin_diff_circ_inclusions_aniso; exit" -logfile results/monoscale_sto_lin_diff_circ_inclusions_aniso.log
	matlab -nodesktop -nosplash -noFigureWindows -r "monoscale_sto_nonlin_diff_reac; exit" -logfile results/monoscale_sto_nonlin_diff_reac.log

# Sparse approximation ---------------------------
sparse:
	matlab -nodesktop -nosplash -noFigureWindows -r "sparse_anisotropic_function; exit" -logfile results/sparse_anisotropic_function.log
	matlab -nodesktop -nosplash -noFigureWindows -r "sparse_ishigami_function; exit" -logfile results/sparse_ishigami_function.log
	matlab -nodesktop -nosplash -noFigureWindows -r "sparse_polynomial_function; exit" -logfile results/sparse_polynomial_function.log
	matlab -nodesktop -nosplash -noFigureWindows -r "sparse_sobol_function; exit" -logfile results/sparse_sobol_function.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "sparse_geometric_brownian; exit" -logfile results/sparse_geometric_brownian.log

# Plate Shell ---------------------------
plate:
	matlab -nodesktop -nosplash -noFigureWindows -r "plate_circ_det_lin_elas; exit" -logfile results/plate_circ_det_lin_elas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "plate_circ_det_lin_elas_cv; exit" -logfile results/plate_circ_det_lin_elas_cv.log
	matlab -nodesktop -nosplash -noFigureWindows -r "plate_rect_det_lin_elas; exit" -logfile results/plate_rect_det_lin_elas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "plate_rect_det_lin_elas_cv; exit" -logfile results/plate_rect_det_lin_elas_cv.log
	matlab -nodesktop -nosplash -noFigureWindows -r "table_circ_det_lin_elas; exit" -logfile results/table_circ_det_lin_elas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "table_rect_det_lin_elas; exit" -logfile results/table_rect_det_lin_elas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "FCBA_table_circ_det_lin_elas; exit" -logfile results/FCBA_table_circ_det_lin_elas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "FCBA_table_circ_sto_lin_elas; exit" -logfile results/FCBA_table_circ_sto_lin_elas.log

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
	rsync -auv $(Dir_mac)/FEMObject/MYCODE $(User_Name_linux)@$(Machine_Name_pcmsme):$(Dir_pcmsme)/FEMObject $(exclude)

sync_from_pcmsme_to_mac:
	rsync -auv $(User_Name_linux)@$(Machine_Name_pcmsme):$(Dir_pcmsme)/FEMObject/MYCODE $(Dir_mac)/FEMObject $(exclude)

sync_from_mac_to_cluster:
	rsync -auv $(Dir_mac)/FEMObject/MYCODE $(User_Name_linux)@$(Cluster_Name):$(Dir_cluster)/FEMObject $(exclude)

sync_from_cluster_to_mac:
	rsync -auv $(User_Name_linux)@$(Cluster_Name):$(Dir_cluster)/FEMObject/MYCODE $(Dir_mac)/FEMObject $(exclude)

sync_results_from_mac_to_pcmsme:
	rsync -auv $(Dir_mac)/FEMObject/MYCODE/results $(User_Name_linux)@$(Machine_Name_pcmsme):$(Dir_pcmsme)/FEMObject/MYCODE

sync_results_from_pcmsme_to_mac:
	rsync -auv $(User_Name_linux)@$(Machine_Name_pcmsme):$(Dir_pcmsme)/FEMObject/MYCODE/results $(Dir_mac)/FEMObject/MYCODE

# Cluster ---------------------------
cluster:
	ssh -YC $(User_Name_linux)@$(Cluster_Name) "hostname; cd $(Dir_cluster)/FEMObject/MYCODE/examples/multiscale; qsub $(Pb_Name).pbs"

interactive:
	ssh -YC $(User_Name_linux)@$(Cluster_Name) "hostname; cd $(Dir_cluster)/FEMObject/MYCODE/examples/multiscale; qsub -I -X $(Pb_Name).pbs"
