# Username
User_Name_mac = Op
User_Name_pcmsme = pled
User_Name_cluster = pled
# Machine name
Machine_Name_mac = MacBookProOp
Machine_Name_pcmsme = pcmsmemeca22
# Directory
Dir_mac = /Users/$(User_Name_mac)/Documents/Recherche/GeM
Dir_pcmsme = /home/p/$(User_Name_pcmsme)/Documents
Dir_cluster = /home/p/pled/$(User_Name_cluster)/Documents
Dir_qsub = /usr/local/pbs/default/bin
# Cluster name
Cluster_Name = cluster
# Problem name
#Pb_Name = multiscale_det_lin_diff_form
Pb_Name = multiscale_det_lin_diff
#Pb_Name = multiscale_det_lin_diff_circ_holes
#Pb_Name = multiscale_det_lin_elas
#Pb_Name = multiscale_det_lin_elas_circ_holes
#Pb_Name = multiscale_det_lin_elas_edge_cracks
#Pb_Name = multiscale_det_lin_elas_interior_cracks
#Pb_Name = multiscale_det_nonlin_diff_reac

#Pb_Name = multiscale_sto_lin_diff
#Pb_Name = multiscale_sto_lin_elas
#Pb_Name = multiscale_sto_nonlin_diff_reac

#Pb_Name = multiscale_sto_lin_diff_align_inclusions
#Pb_Name = multiscale_sto_lin_diff_square_inclusions_iso
#Pb_Name = multiscale_sto_lin_diff_square_inclusions_aniso

#Pb_Name = multiscale_sto_lin_diff_circ_inclusions_iso
#Pb_Name = multiscale_sto_lin_diff_circ_inclusions_aniso

#Pb_Name = monoscale_sto_lin_diff
#Pb_Name = monoscale_sto_nonlin_diff_reac
#Pb_Name = monoscale_sto_lin_diff_circ_inclusions_iso
#Pb_Name = monoscale_sto_lin_diff_circ_inclusions_aniso

#Pb_Name = adaptive_sparse_approx_polynomial_function
#Pb_Name = adaptive_sparse_approx_ishigami_function
#Pb_Name = adaptive_sparse_approx_anisotropic_function
#Pb_Name = adaptive_sparse_approx_geometric_brownian

# Default ---------------------------
default:
	matlab -nodesktop -nosplash -r "$(Pb_Name); exit" -logfile RESULTS/$(Pb_Name).log

# All ---------------------------
all: multidet multisto monosto

# Multiscale ---------------------------
multi: multidet multisto

# Multiscale Deterministic ---------------------------
multidet:
	matlab -nodesktop -nosplash -r "multiscale_det_lin_diff_form; exit" -logfile RESULTS/multiscale_det_lin_diff_form.log
	matlab -nodesktop -nosplash -r "multiscale_det_lin_diff; exit" -logfile RESULTS/multiscale_det_lin_diff.log
	matlab -nodesktop -nosplash -r "multiscale_det_lin_diff_circ_holes; exit" -logfile RESULTS/multiscale_det_lin_diff_circ_holes.log
	matlab -nodesktop -nosplash -r "multiscale_det_lin_elas; exit" -logfile RESULTS/multiscale_det_lin_elas.log
	matlab -nodesktop -nosplash -r "multiscale_det_lin_elas_circ_holes; exit" -logfile RESULTS/multiscale_det_lin_elas_circ_holes.log
	matlab -nodesktop -nosplash -r "multiscale_det_lin_elas_edge_cracks; exit" -logfile RESULTS/multiscale_det_lin_elas_edge_cracks.log
	matlab -nodesktop -nosplash -r "multiscale_det_lin_elas_interior_cracks; exit" -logfile RESULTS/multiscale_det_lin_elas_interior_cracks.log
	matlab -nodesktop -nosplash -r "multiscale_det_nonlin_diff_reac; exit" -logfile RESULTS/multiscale_det_nonlin_diff_reac.log

# Multiscale Stochastic ---------------------------
multisto:
	matlab -nodesktop -nosplash -r "multiscale_sto_lin_diff; exit" -logfile RESULTS/multiscale_sto_lin_diff.log
	matlab -nodesktop -nosplash -r "multiscale_sto_lin_elas; exit" -logfile RESULTS/multiscale_sto_lin_elas.log
	matlab -nodesktop -nosplash -r "multiscale_sto_nonlin_diff_reac; exit" -logfile RESULTS/multiscale_sto_nonlin_diff_reac.log
	matlab -nodesktop -nosplash -r "multiscale_sto_lin_diff_align_inclusions; exit" -logfile RESULTS/multiscale_sto_lin_diff_align_inclusions.log
	matlab -nodesktop -nosplash -r "multiscale_sto_lin_diff_square_inclusions_iso; exit" -logfile RESULTS/multiscale_sto_lin_diff_square_inclusions_iso.log
	matlab -nodesktop -nosplash -r "multiscale_sto_lin_diff_square_inclusions_aniso; exit" -logfile RESULTS/multiscale_sto_lin_diff_square_inclusions_aniso.log
#	matlab -nodesktop -nosplash -r "multiscale_sto_lin_diff_circ_inclusions_iso; exit" -logfile RESULTS/multiscale_sto_lin_diff_circ_inclusions_iso.log
#	matlab -nodesktop -nosplash -r "multiscale_sto_lin_diff_circ_inclusions_aniso; exit" -logfile RESULTS/multiscale_sto_lin_diff_circ_inclusions_aniso.log

# Monoscale Stochastic ---------------------------
mono:
	matlab -nodesktop -nosplash -r "monoscale_sto_lin_diff; exit" -logfile RESULTS/monoscale_sto_lin_diff.log
	matlab -nodesktop -nosplash -r "monoscale_sto_nonlin_diff_reac; exit" -logfile RESULTS/monoscale_sto_nonlin_diff_reac.log
	matlab -nodesktop -nosplash -r "monoscale_sto_lin_diff_circ_inclusions_iso; exit" -logfile RESULTS/monoscale_sto_lin_diff_circ_inclusions_iso.log
	matlab -nodesktop -nosplash -r "monoscale_sto_lin_diff_circ_inclusions_aniso; exit" -logfile RESULTS/monoscale_sto_lin_diff_circ_inclusions_aniso.log

# Sparse approximation ---------------------------
sparseapprox:
	matlab -nodesktop -nosplash -r "adaptive_sparse_approx_polynomial_function; exit" -logfile RESULTS/adaptive_sparse_approx_polynomial_function.log
	matlab -nodesktop -nosplash -r "adaptive_sparse_approx_ishigami_function; exit" -logfile RESULTS/adaptive_sparse_approx_ishigami_function.log
	matlab -nodesktop -nosplash -r "adaptive_sparse_approx_anisotropic_function; exit" -logfile RESULTS/adaptive_sparse_approx_anisotropic_function.log
	matlab -nodesktop -nosplash -r "adaptive_sparse_approx_geometric_brownian; exit" -logfile RESULTS/adaptive_sparse_approx_geometric_brownian.log

# Clean ---------------------------
clean:
	-rm -r RESULTS/*

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
	rsync -auv $(Dir_mac)/FEMOBJECT/MYCODE $(User_Name_pcmsme)@$(Machine_Name_pcmsme):$(Dir_pcmsme)/FEMOBJECT $(exclude)

sync_from_pcmsme_to_mac:
	rsync -auv $(User_Name_pcmsme)@$(Machine_Name_pcmsme):$(Dir_pcmsme)/FEMOBJECT/MYCODE $(Dir_mac)/FEMOBJECT $(exclude)

sync_from_mac_to_cluster:
	rsync -auv $(Dir_mac)/FEMOBJECT/MYCODE $(User_Name_cluster)@$(Cluster_Name):$(Dir_cluster)/FEMOBJECT $(exclude)

sync_from_cluster_to_mac:
	rsync -auv $(User_Name_cluster)@$(Cluster_Name):$(Dir_cluster)/FEMOBJECT/MYCODE $(Dir_mac)/FEMOBJECT $(exclude)

sync_results_from_mac_to_pcmsme:
	rsync -auv $(Dir_mac)/FEMOBJECT/MYCODE/RESULTS $(User_Name_pcmsme)@$(Machine_Name_pcmsme):$(Dir_pcmsme)/FEMOBJECT/MYCODE

sync_results_from_pcmsme_to_mac:
	rsync -auv $(User_Name_pcmsme)@$(Machine_Name_pcmsme):$(Dir_pcmsme)/FEMOBJECT/MYCODE/RESULTS $(Dir_mac)/FEMOBJECT/MYCODE

# Cluster ---------------------------
cluster:
	ssh -YC $(User_Name_cluster)@$(Cluster_Name) "hostname; cd $(Dir_cluster)/FEMOBJECT/MYCODE/EXAMPLES/MULTISCALE; qsub $(Pb_Name).pbs"

interactive:
	ssh -YC $(User_Name_cluster)@$(Cluster_Name) "hostname; cd $(Dir_cluster)/FEMOBJECT/MYCODE/EXAMPLES/MULTISCALE; qsub -I -X $(Pb_Name).pbs"
