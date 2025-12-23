# Username
User_Name_mac = Op
User_Name_linux = pled
# Machine name
#Machine_Name_mac = MacBookProOp
Machine_Name_mac = MPLAV-MSME-24-001
Machine_Name_pcmsme = pcmsmemeca52
# Directory
Dir_mac = /Users/$(User_Name_mac)/Documents/Recherche/GeM
Dir_pcmsme = /home/p/$(User_Name_linux)/Documents
Dir_cluster = /home/p/pled/$(User_Name_linux)/Documents
Dir_qsub = /usr/local/pbs/default/bin
# Cluster name
Cluster_Name = cluster
# Problem name
#Pb_Name = multiscaleDetLinDiffForm
#Pb_Name = multiscaleDetLinDiff
#Pb_Name = multiscaleDetLinDiffCircHoles
#Pb_Name = multiscaleDetLinElas
#Pb_Name = multiscaleDetLinElasCircHoles
#Pb_Name = multiscaleDetLinElasEdgeCracks
#Pb_Name = multiscaleDetLinElasInteriorCracks
#Pb_Name = multiscaleDetNonlinDiffReac
#Pb_Name = multiscaleDetTransientLinAdvDiffReac
#Pb_Name = multiscaleDetDynLinElas

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
#Pb_Name = multiscaleStoTransientLinAdvDiffReac
#Pb_Name = multiscaleStoDynLinElas

#Pb_Name = monoscaleStoLinDiff
#Pb_Name = monoscaleStoLinDiffCircInclusionsIso
#Pb_Name = monoscaleStoLinDiffCircInclusionsAniso
#Pb_Name = monoscaleStoNonlinDiffReac
#Pb_Name = monoscaleDetTransientLinAdvDiffReac
#Pb_Name = monoscaleDetDynLinElas
#Pb_Name = monoscaleStoTransientLinAdvDiffReac
#Pb_Name = monoscaleStoDynLinElas

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
#Pb_Name = FCBADeskBeamJunctionDetLinElas
#Pb_Name = FCBADeskPlateDetLinElas
#Pb_Name = FCBADeskPlateStoLinElas
#Pb_Name = FCBADeskPlateJunctionDetLinElas
#Pb_Name = FCBADeskPlateJunctionStoLinElas
#Pb_Name = FCBADesk3DDetLinElas
#Pb_Name = FCBADesk3DStoLinElas
#Pb_Name = FCBAJunctionScrewDowelBeamDetLinElas
#Pb_Name = FCBAJunctionScrewDowelBeamStoLinElas

#Pb_Name = identification_ET_GL
#Pb_Name = identification_EL_NUL
#Pb_Name = identification_KS
#Pb_Name = identification_KD

#Pb_Name = phasefieldDetLinElasAsymmetricNotchedPlate
#Pb_Name = phasefieldDetLinElasAsymmetricNotchedPlateAdaptive
#Pb_Name = phasefieldDetLinElasDoubleEdgeCrack
#Pb_Name = phasefieldDetLinElasDoubleEdgeCrackAdaptive
#Pb_Name = phasefieldDetLinElasDoubleEdgeCrackSelfHealing
#Pb_Name = phasefieldDetLinElasLshapedPanel
#Pb_Name = phasefieldDetLinElasLshapedPanelAdaptive
#Pb_Name = phasefieldDetLinElasLshapedPanelSelfHealing
#Pb_Name = phasefieldDetLinElasPlateWithHole
#Pb_Name = phasefieldDetLinElasPlateWithHoleAdaptive
#Pb_Name = phasefieldDetLinElasPlateWithHoleSelfHealing
#Pb_Name = phasefieldDetLinElasPlateWithHoleWoodNoel
#Pb_Name = phasefieldDetLinElasSingleEdgeCrack
Pb_Name = phasefieldDetLinElasSingleEdgeCrackAdaptive
#Pb_Name = phasefieldDetLinElasSingleEdgeCrackSelfHealing

# Default ---------------------------
default:
	matlab -nodesktop -nosplash -noFigureWindows -r "$(Pb_Name); exit" -logfile results/$(Pb_Name).log
#	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleStoLinDiffAlignInclusionsIso; exit" -logfile results/multiscaleSto/linDiffAlignInclusionsIso.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleStoLinDiffAlignInclusionsAniso; exit" -logfile results/multiscaleSto/linDiffAlignInclusionsAniso.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleStoTransientLinAdvDiffReac; exit" -logfile results/multiscaleSto/transientLinAdvDiffReac.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateStoLinElasStaticHori1_100N; exit" -logfile results/plate/FCBADeskPlateStoLinElasStaticHori1_100N.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateStoLinElasStaticHori1_200N; exit" -logfile results/plate/FCBADeskPlateStoLinElasStaticHori1_200N.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateStoLinElasStaticHori2_100N; exit" -logfile results/plate/FCBADeskPlateStoLinElasStaticHori2_100N.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateStoLinElasStaticHori2_200N; exit" -logfile results/plate/FCBADeskPlateStoLinElasStaticHori2_200N.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateStoLinElasStaticVert_300N; exit" -logfile results/plate/FCBADeskPlateStoLinElasStaticVert_300N.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateStoLinElasStaticVert_400N; exit" -logfile results/plate/FCBADeskPlateStoLinElasStaticVert_400N.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateStoLinElasStaticVert_500N; exit" -logfile results/plate/FCBADeskPlateStoLinElasStaticVert_500N.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateStoLinElasStability_400N; exit" -logfile results/plate/FCBADeskPlateStoLinElasStability_400N.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateStoLinElasFatigue1_100N; exit" -logfile results/plate/FCBADeskPlateStoLinElasFatigue1_100N.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateStoLinElasFatigue2_100N; exit" -logfile results/plate/FCBADeskPlateStoLinElasFatigue2_100N.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateJunctionStoLinElasStaticHori1_100N; exit" -logfile results/plate/FCBADeskPlateJunctionStoLinElasStaticHori1_100N.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateJunctionStoLinElasStaticHori1_200N; exit" -logfile results/plate/FCBADeskPlateJunctionStoLinElasStaticHori1_200N.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateJunctionStoLinElasStaticHori2_100N; exit" -logfile results/plate/FCBADeskPlateJunctionStoLinElasStaticHori2_100N.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateJunctionStoLinElasStaticHori2_200N; exit" -logfile results/plate/FCBADeskPlateJunctionStoLinElasStaticHori2_200N.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateJunctionStoLinElasStaticVert_300N; exit" -logfile results/plate/FCBADeskPlateJunctionStoLinElasStaticVert_300N.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateJunctionStoLinElasStaticVert_400N; exit" -logfile results/plate/FCBADeskPlateJunctionStoLinElasStaticVert_400N.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateJunctionStoLinElasStaticVert_500N; exit" -logfile results/plate/FCBADeskPlateJunctionStoLinElasStaticVert_500N.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateJunctionStoLinElasStability_400N; exit" -logfile results/plate/FCBADeskPlateJunctionStoLinElasStability_400N.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateJunctionStoLinElasFatigue1_100N; exit" -logfile results/plate/FCBADeskPlateJunctionStoLinElasFatigue1_100N.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateJunctionStoLinElasFatigue2_100N; exit" -logfile results/plate/FCBADeskPlateJunctionStoLinElasFatigue2_100N.log

# All ---------------------------
all: sparse monosto multidet multisto identification plate phasefield

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
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleDetTransientLinAdvDiffReac; exit" -logfile results/multiscaleDet/transientLinAdvDiffReac.log
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleDetDynLinElas; exit" -logfile results/multiscaleDet/dynLinElas.log

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
	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleStoTransientLinAdvDiffReac; exit" -logfile results/multiscaleSto/transientLinAdvDiffReac.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "multiscaleStoDynLinElas; exit" -logfile results/multiscaleSto/dynLinElas.log

# Monoscale Stochastic ---------------------------
monosto:
	matlab -nodesktop -nosplash -noFigureWindows -r "monoscaleStoLinDiff; exit" -logfile results/monoscaleSto/linDiff.log
	matlab -nodesktop -nosplash -noFigureWindows -r "monoscaleStoLinDiffCircInclusionsIso; exit" -logfile results/monoscaleSto/linDiffCircInclusionsIso.log
	matlab -nodesktop -nosplash -noFigureWindows -r "monoscaleStoLinDiffCircInclusionsAniso; exit" -logfile results/monoscaleSto/linDiffCircInclusionsAniso.log
	matlab -nodesktop -nosplash -noFigureWindows -r "monoscaleStoNonlinDiffReac; exit" -logfile results/monoscaleSto/nonlinDiffReac.log
	matlab -nodesktop -nosplash -noFigureWindows -r "monoscaleStoTransientLinAdvDiffReac; exit" -logfile results/monoscaleSto/transientLinAdvDiffReac.log
	matlab -nodesktop -nosplash -noFigureWindows -r "monoscaleStoDynLinElas; exit" -logfile results/monoscaleSto/dynLinElas.log

# Sparse approximation ---------------------------
sparse:
	matlab -nodesktop -nosplash -noFigureWindows -r "sparseAnisotropicFunction; exit" -logfile results/sparse/anisotropicFunction.log
#	matlab -nodesktop -nosplash -noFigureWindows -r "sparseGeometricBrownian; exit" -logfile results/sparse/geometricBrownian.log
	matlab -nodesktop -nosplash -noFigureWindows -r "sparseIshigamiFunction; exit" -logfile results/sparse/ishigamiFunction.log
	matlab -nodesktop -nosplash -noFigureWindows -r "sparsePolynomialFunction; exit" -logfile results/sparse/polynomialFunction.log
	matlab -nodesktop -nosplash -noFigureWindows -r "sparseSobolFunction; exit" -logfile results/sparse/sobolFunction.log

# Plate Shell ---------------------------
plate:
	matlab -nodesktop -nosplash -noFigureWindows -r "plateCircDetLinElas; exit" -logfile results/FCBA/plateCircDetLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "plateCircDetLinElasCv; exit" -logfile results/FCBA/plateCircDetLinElasCv.log
	matlab -nodesktop -nosplash -noFigureWindows -r "plateRectDetLinElas; exit" -logfile results/FCBA/plateRectDetLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "plateRectDetLinElasCv; exit" -logfile results/FCBA/plateRectDetLinElasCv.log
	matlab -nodesktop -nosplash -noFigureWindows -r "tableCircDetLinElas; exit" -logfile results/FCBA/tableCircDetLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "tableRectDetLinElas; exit" -logfile results/FCBA/tableRectDetLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "FCBATableCircDetLinElas; exit" -logfile results/FCBA/FCBATableCircDetLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "FCBATableCircStoLinElas; exit" -logfile results/FCBA/FCBATableCircStoLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskBeamJunctionDetLinElas; exit" -logfile results/FCBA/FCBADeskBeamJunctionDetLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskDetLinElas; exit" -logfile results/FCBA/FCBADeskDetLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateDetLinElas; exit" -logfile results/FCBA/FCBADeskPlateDetLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateStoLinElas; exit" -logfile results/FCBA/FCBADeskPlateStoLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateJunctionDetLinElas; exit" -logfile results/FCBA/FCBADeskPlateJunctionDetLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADeskPlateJunctionStoLinElas; exit" -logfile results/FCBA/FCBADeskPlateJunctionStoLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADesk3DDetLinElas; exit" -logfile results/FCBA/FCBADesk3DDetLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "FCBADesk3DStoLinElas; exit" -logfile results/FCBA/FCBADesk3DStoLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "FCBAJunctionScrewDowelBeamDetLinElas; exit" -logfile results/FCBA/FCBAJunctionScrewDowelBeamDetLinElas.log
	matlab -nodesktop -nosplash -noFigureWindows -r "FCBAJunctionScrewDowelPlateDetLinElas; exit" -logfile results/FCBA/FCBAJunctionScrewDowelPlateDetLinElas.log

# Identification ---------------------------
identification:
	matlab -nodesktop -nosplash -noFigureWindows -r "identification_ET_GL; exit" -logfile results/identification/materialParticleBoard/identification_ET_GL.log
	matlab -nodesktop -nosplash -noFigureWindows -r "identification_EL_NUL; exit" -logfile results/identification/materialParticleBoard/identification_EL_NUL.log
	matlab -nodesktop -nosplash -noFigureWindows -r "identification_KS; exit" -logfile results/identification/materialParticleBoard/identification_KS.log
	matlab -nodesktop -nosplash -noFigureWindows -r "identification_KD; exit" -logfile results/identification/materialParticleBoard/identification_KD.log

# Multiscale ---------------------------
phasefield: phasefielddet phasefieldsto

# Phase-field Deterministic ---------------------------
phasefielddet:
	matlab -nodesktop -nosplash -noFigureWindows -r "phasefieldDetLinElasAsymmetricNotchedPlate; exit" -logfile results/phasefieldDet/phasefieldDetLinElasAsymmetricNotchedPlate.log
	matlab -nodesktop -nosplash -noFigureWindows -r "phasefieldDetLinElasAsymmetricNotchedPlateAdaptive; exit" -logfile results/phasefieldDet/phasefieldDetLinElasAsymmetricNotchedPlateAdaptive.log
	matlab -nodesktop -nosplash -noFigureWindows -r "phasefieldDetLinElasDoubleEdgeCrack; exit" -logfile results/phasefieldDet/phasefieldDetLinElasDoubleEdgeCrack.log
	matlab -nodesktop -nosplash -noFigureWindows -r "phasefieldDetLinElasDoubleEdgeCrackAdaptive; exit" -logfile results/phasefieldDet/phasefieldDetLinElasDoubleEdgeCrackAdaptive.log
	matlab -nodesktop -nosplash -noFigureWindows -r "phasefieldDetLinElasDoubleEdgeCrackSelfHealing; exit" -logfile results/phasefieldDet/phasefieldDetLinElasDoubleEdgeCrackSelfHealing.log
	matlab -nodesktop -nosplash -noFigureWindows -r "phasefieldDetLinElasLshapedPanel; exit" -logfile results/phasefieldDet/phasefieldDetLinElasLshapedPanel.log
	matlab -nodesktop -nosplash -noFigureWindows -r "phasefieldDetLinElasLshapedPanelAdaptive; exit" -logfile results/phasefieldDet/phasefieldDetLinElasLshapedPanelAdaptive.log
	matlab -nodesktop -nosplash -noFigureWindows -r "phasefieldDetLinElasLshapedPanelSelfHealing; exit" -logfile results/phasefieldDet/phasefieldDetLinElasLshapedPanelSelfHealing.log
	matlab -nodesktop -nosplash -noFigureWindows -r "phasefieldDetLinElasPlateWithHole; exit" -logfile results/phasefieldDet/phasefieldDetLinElasPlateWithHole.log
	matlab -nodesktop -nosplash -noFigureWindows -r "phasefieldDetLinElasPlateWithHoleAdaptive; exit" -logfile results/phasefieldDet/phasefieldDetLinElasPlateWithHoleAdaptive.log
	matlab -nodesktop -nosplash -noFigureWindows -r "phasefieldDetLinElasPlateWithHoleSelfHealing; exit" -logfile results/phasefieldDet/phasefieldDetLinElasPlateWithHoleSelfHealing.log
	matlab -nodesktop -nosplash -noFigureWindows -r "phasefieldDetLinElasPlateWithHoleWoodNoel; exit" -logfile results/phasefieldDet/phasefieldDetLinElasPlateWithHoleWoodNoel.log
	matlab -nodesktop -nosplash -noFigureWindows -r "phasefieldDetLinElasSingleEdgeCrack; exit" -logfile results/phasefieldDet/phasefieldDetLinElasSingleEdgeCrack.log
	matlab -nodesktop -nosplash -noFigureWindows -r "phasefieldDetLinElasSingleEdgeCrackAdaptive; exit" -logfile results/phasefieldDet/phasefieldDetLinElasSingleEdgeCrackAdaptive.log
	matlab -nodesktop -nosplash -noFigureWindows -r "phasefieldDetLinElasSingleEdgeCrackSelfHealing; exit" -logfile results/phasefieldDet/phasefieldDetLinElasSingleEdgeCrackSelfHealing.log

# Phase-field Stochastic ---------------------------
phasefieldsto:
	matlab -nodesktop -nosplash -noFigureWindows -r "phasefieldStoLinElasAsymmetricNotchedPlate; exit" -logfile results/phasefieldSto/phasefieldStoLinElasAsymmetricNotchedPlate.log
	matlab -nodesktop -nosplash -noFigureWindows -r "phasefieldStoLinElasAsymmetricNotchedPlateAdaptive; exit" -logfile results/phasefieldSto/phasefieldStoLinElasAsymmetricNotchedPlateAdaptive.log
	matlab -nodesktop -nosplash -noFigureWindows -r "phasefieldStoLinElasSingleEdgeCrack; exit" -logfile results/phasefieldSto/phasefieldStoLinElasSingleEdgeCrack.log
	matlab -nodesktop -nosplash -noFigureWindows -r "phasefieldStoLinElasSingleEdgeCrackAdaptive; exit" -logfile results/phasefieldSto/phasefieldStoLinElasSingleEdgeCrackAdaptive.log

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
