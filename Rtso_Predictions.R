#Library activation

library(miRBaseConverter)
library(miRNAtap)
library(miRNAtap.db)
library(topGO)
library(org.Hs.eg.db)
library(GOplot)

#AllGO2Genes

allGO2genes = annFUN.org(whichOnto='BP', feasibleGenes = NULL,
                         mapping="org.Hs.eg.db", ID = "entrez")

#miRNAs from RTSO

hsa_let_7a_5p	=	'let-7a-5p'
hsa_let_7b_5p	=	'let-7b-5p'
hsa_let_7c_5p	=	'let-7c-5p'
hsa_let_7d_5p	=	'let-7d-5p'
hsa_let_7e_5p	=	'let-7e-5p'
hsa_let_7f_5p	=	'let-7f-5p'
hsa_let_7g_5p	=	'let-7g-5p'
hsa_let_7i_5p	=	'let-7i-5p'
hsa_miR_1_3p	=	'miR-1-3p'
hsa_miR_100_5p	=	'miR-100-5p'
hsa_miR_101_3p	=	'miR-101-3p'
hsa_miR_101_5p	=	'miR-101-5p'
hsa_miR_103a_3p	=	'miR-103a-3p'
hsa_miR_106b_3p	=	'miR-106b-3p'
hsa_miR_107	=	'miR-107'
hsa_miR_10a_5p	=	'miR-10a-5p'
hsa_miR_10b_5p	=	'miR-10b-5p'
hsa_miR_1227_5p	=	'miR-1227-5p'
hsa_miR_1238_5p	=	'miR-1238-5p'
hsa_miR_1247_5p	=	'miR-1247-5p'
hsa_miR_125a_5p	=	'miR-125a-5p'
hsa_miR_125b_5p	=	'miR-125b-5p'
hsa_miR_126_3p	=	'miR-126-3p'
hsa_miR_126_5p	=	'miR-126-5p'
hsa_miR_128_3p	=	'miR-128-3p'
hsa_miR_129_2_3p	=	'miR-129-2-3p'
hsa_miR_1296_5p	=	'miR-1296-5p'
hsa_miR_1306_3p	=	'miR-1306-3p'
hsa_miR_130a_3p	=	'miR-130a-3p'
hsa_miR_133a_3p	=	'miR-133a-3p'
hsa_miR_133a_5p	=	'miR-133a-5p'
hsa_miR_133b	=	'miR-133b'
hsa_miR_135a_5p	=	'miR-135a-5p'
hsa_miR_138_5p	=	'miR-138-5p'
hsa_miR_139_3p	=	'miR-139-3p'
hsa_miR_139_5p	=	'miR-139-5p'
hsa_miR_140_3p	=	'miR-140-3p'
hsa_miR_140_5p	=	'miR-140-5p'
hsa_miR_142_3p	=	'miR-142-3p'
hsa_miR_142_5p	=	'miR-142-5p'
hsa_miR_143_3p	=	'miR-143-3p'
hsa_miR_143_5p	=	'miR-143-5p'
hsa_miR_144_3p	=	'miR-144-3p'
hsa_miR_144_5p	=	'miR-144-5p'
hsa_miR_145_3p	=	'miR-145-3p'
hsa_miR_145_5p	=	'miR-145-5p'
hsa_miR_146a_5p	=	'miR-146a-5p'
hsa_miR_146b_5p	=	'miR-146b-5p'
hsa_miR_147b	=	'miR-147b'
hsa_miR_150_5p	=	'miR-150-5p'
hsa_miR_151a_5p	=	'miR-151a-5p'
hsa_miR_151b	=	'miR-151b'
hsa_miR_152_3p	=	'miR-152-3p'
hsa_miR_1537_3p	=	'miR-1537-3p'
hsa_miR_15a_5p	=	'miR-15a-5p'
hsa_miR_15b_3p	=	'miR-15b-3p'
hsa_miR_15b_5p	=	'miR-15b-5p'
hsa_miR_16_2_3p	=	'miR-16-2-3p'
hsa_miR_16_5p	=	'miR-16-5p'
hsa_miR_181a_2_3p	=	'miR-181a-2-3p'
hsa_miR_181a_3p	=	'miR-181a-3p'
hsa_miR_182_3p	=	'miR-182-3p'
hsa_miR_183_3p	=	'miR-183-3p'
hsa_miR_184	=	'miR-184'
hsa_miR_185_5p	=	'miR-185-5p'
hsa_miR_186_5p	=	'miR-186-5p'
hsa_miR_187_5p	=	'miR-187-5p'
hsa_miR_18a_5p	=	'miR-18a-5p'
hsa_miR_18b_5p	=	'miR-18b-5p'
hsa_miR_190a_5p	=	'miR-190a-5p'
hsa_miR_191_5p	=	'miR-191-5p'
hsa_miR_1913	=	'miR-1913'
hsa_miR_193a_5p	=	'miR-193a-5p'
hsa_miR_195_5p	=	'miR-195-5p'
hsa_miR_199a_3p	=	'miR-199a-3p'
hsa_miR_199a_5p	=	'miR-199a-5p'
hsa_miR_199b_5p	=	'miR-199b-5p'
hsa_miR_203a_3p	=	'miR-203a-3p'
hsa_miR_204_5p	=	'miR-204-5p'
hsa_miR_205_3p	=	'miR-205-3p'
hsa_miR_20a_5p	=	'miR-20a-5p'
hsa_miR_20b_5p	=	'miR-20b-5p'
hsa_miR_21_3p	=	'miR-21-3p'
hsa_miR_214_3p	=	'miR-214-3p'
hsa_miR_214_5p	=	'miR-214-5p'
hsa_miR_218_5p	=	'miR-218-5p'
hsa_miR_22_5p	=	'miR-22-5p'
hsa_miR_221_3p	=	'miR-221-3p'
hsa_miR_221_5p	=	'miR-221-5p'
hsa_miR_222_3p	=	'miR-222-3p'
hsa_miR_223_3p	=	'miR-223-3p'
hsa_miR_223_5p	=	'miR-223-5p'
hsa_miR_224_3p	=	'miR-224-3p'
hsa_miR_2277_3p	=	'miR-2277-3p'
hsa_miR_23a_3p	=	'miR-23a-3p'
hsa_miR_23a_5p	=	'miR-23a-5p'
hsa_miR_23b_3p	=	'miR-23b-3p'
hsa_miR_24_3p	=	'miR-24-3p'
hsa_miR_26a_1_3p	=	'miR-26a-1-3p'
hsa_miR_26a_5p	=	'miR-26a-5p'
hsa_miR_26b_3p	=	'miR-26b-3p'
hsa_miR_26b_5p	=	'miR-26b-5p'
hsa_miR_27a_3p	=	'miR-27a-3p'
hsa_miR_27a_5p	=	'miR-27a-5p'
hsa_miR_27b_3p	=	'miR-27b-3p'
hsa_miR_28_5p	=	'miR-28-5p'
hsa_miR_2861	=	'miR-2861'
hsa_miR_29a_3p	=	'miR-29a-3p'
hsa_miR_29a_5p	=	'miR-29a-5p'
hsa_miR_29b_1_5p	=	'miR-29b-1-5p'
hsa_miR_29b_2_5p	=	'miR-29b-2-5p'
hsa_miR_29b_3p	=	'miR-29b-3p'
hsa_miR_29c_3p	=	'miR-29c-3p'
hsa_miR_29c_5p	=	'miR-29c-5p'
hsa_miR_3065_3p	=	'miR-3065-3p'
hsa_miR_3065_5p	=	'miR-3065-5p'
hsa_miR_30a_3p	=	'miR-30a-3p'
hsa_miR_30a_5p	=	'miR-30a-5p'
hsa_miR_30b_5p	=	'miR-30b-5p'
hsa_miR_30c_2_3p	=	'miR-30c-2-3p'
hsa_miR_30c_5p	=	'miR-30c-5p'
hsa_miR_30d_5p	=	'miR-30d-5p'
hsa_miR_30e_3p	=	'miR-30e-3p'
hsa_miR_3149	=	'miR-3149'
hsa_miR_3188	=	'miR-3188'
hsa_miR_32_5p	=	'miR-32-5p'
hsa_miR_324_3p	=	'miR-324-3p'
hsa_miR_326	=	'miR-326'
hsa_miR_328_3p	=	'miR-328-3p'
hsa_miR_331_3p	=	'miR-331-3p'
hsa_miR_335_5p	=	'miR-335-5p'
hsa_miR_338_3p	=	'miR-338-3p'
hsa_miR_339_5p	=	'miR-339-5p'
hsa_miR_340_5p	=	'miR-340-5p'
hsa_miR_342_3p	=	'miR-342-3p'
hsa_miR_342_5p	=	'miR-342-5p'
hsa_miR_34a_3p	=	'miR-34a-3p'
hsa_miR_34a_5p	=	'miR-34a-5p'
hsa_miR_34b_5p	=	'miR-34b-5p'
hsa_miR_34c_5p	=	'miR-34c-5p'
hsa_miR_3607_3p	=	'miR-3607-3p'
hsa_miR_361_3p	=	'miR-361-3p'
hsa_miR_362_3p	=	'miR-362-3p'
hsa_miR_362_5p	=	'miR-362-5p'
hsa_miR_3620_5p	=	'miR-3620-5p'
hsa_miR_363_3p	=	'miR-363-3p'
hsa_miR_3659	=	'miR-3659'
hsa_miR_365a_3p	=	'miR-365a-3p'
hsa_miR_3660	=	'miR-3660'
hsa_miR_3663_5p	=	'miR-3663-5p'
hsa_miR_371b_5p	=	'miR-371b-5p'
hsa_miR_374a_5p	=	'miR-374a-5p'
hsa_miR_374b_5p	=	'miR-374b-5p'
hsa_miR_374c_5p	=	'miR-374c-5p'
hsa_miR_424_5p	=	'miR-424-5p'
hsa_miR_4252	=	'miR-4252'
hsa_miR_4254	=	'miR-4254'
hsa_miR_4284	=	'miR-4284'
hsa_miR_4290	=	'miR-4290'
hsa_miR_4306	=	'miR-4306'
hsa_miR_4317	=	'miR-4317'
hsa_miR_4318	=	'miR-4318'
hsa_miR_4324	=	'miR-4324'
hsa_miR_4328	=	'miR-4328'
hsa_miR_4440	=	'miR-4440'
hsa_miR_4443	=	'miR-4443'
hsa_miR_4481	=	'miR-4481'
hsa_miR_449a	=	'miR-449a'
hsa_miR_450a_5p	=	'miR-450a-5p'
hsa_miR_4516	=	'miR-4516'
hsa_miR_451a	=	'miR-451a'
hsa_miR_452_5p	=	'miR-452-5p'
hsa_miR_4521	=	'miR-4521'
hsa_miR_4532	=	'miR-4532'
hsa_miR_454_3p	=	'miR-454-3p'
hsa_miR_455_3p	=	'miR-455-3p'
hsa_miR_455_5p	=	'miR-455-5p'
hsa_miR_4634	=	'miR-4634'
hsa_miR_4655_3p	=	'miR-4655-3p'
hsa_miR_4695_3p	=	'miR-4695-3p'
hsa_miR_4716_5p	=	'miR-4716-5p'
hsa_miR_4730	=	'miR-4730'
hsa_miR_4731_3p	=	'miR-4731-3p'
hsa_miR_4763_5p	=	'miR-4763-5p'
hsa_miR_4770	=	'miR-4770'
hsa_miR_4793_3p	=	'miR-4793-3p'
hsa_miR_483_3p	=	'miR-483-3p'
hsa_miR_486_3p	=	'miR-486-3p'
hsa_miR_486_5p	=	'miR-486-5p'
hsa_miR_489_3p	=	'miR-489-3p'
hsa_miR_490_3p	=	'miR-490-3p'
hsa_miR_491_3p	=	'miR-491-3p'
hsa_miR_497_3p	=	'miR-497-3p'
hsa_miR_497_5p	=	'miR-497-5p'
hsa_miR_5001_5p	=	'miR-5001-5p'
hsa_miR_500a_3p	=	'miR-500a-3p'
hsa_miR_500b_5p	=	'miR-500b-5p'
hsa_miR_501_5p	=	'miR-501-5p'
hsa_miR_502_3p	=	'miR-502-3p'
hsa_miR_504_3p	=	'miR-504-3p'
hsa_miR_505_5p	=	'miR-505-5p'
hsa_miR_511_3p	=	'miR-511-3p'
hsa_miR_516b_5p	=	'miR-516b-5p'
hsa_miR_517_5p	=	'miR-517-5p'
hsa_miR_517a_3p	=	'miR-517a-3p'
hsa_miR_517c_3p	=	'miR-517c-3p'
hsa_miR_521	=	'miR-521'
hsa_miR_522_3p	=	'miR-522-3p'
hsa_miR_532_3p	=	'miR-532-3p'
hsa_miR_532_5p	=	'miR-532-5p'
hsa_miR_548aa	=	'miR-548aa'
hsa_miR_548aw	=	'miR-548aw'
hsa_miR_548b_3p	=	'miR-548b-3p'
hsa_miR_548c_3p	=	'miR-548c-3p'
hsa_miR_548f_3p	=	'miR-548f-3p'
hsa_miR_548q	=	'miR-548q'
hsa_miR_548x_3p	=	'miR-548x-3p'
hsa_miR_551b_3p	=	'miR-551b-3p'
hsa_miR_5701	=	'miR-5701'
hsa_miR_582_5p	=	'miR-582-5p'
hsa_miR_585_3p	=	'miR-585-3p'
hsa_miR_590_5p	=	'miR-590-5p'
hsa_miR_595	=	'miR-595'
hsa_miR_598_3p	=	'miR-598-3p'
hsa_miR_6068	=	'miR-6068'
hsa_miR_6073	=	'miR-6073'
hsa_miR_6075	=	'miR-6075'
hsa_miR_610	=	'miR-610'
hsa_miR_624_5p	=	'miR-624-5p'
hsa_miR_628_5p	=	'miR-628-5p'
hsa_miR_642b_5p	=	'miR-642b-5p'
hsa_miR_6500_5p	=	'miR-6500-5p'
hsa_miR_6516_3p	=	'miR-6516-3p'
hsa_miR_652_3p	=	'miR-652-3p'
hsa_miR_653_3p	=	'miR-653-3p'
hsa_miR_660_5p	=	'miR-660-5p'
hsa_miR_664a_3p	=	'miR-664a-3p'
hsa_miR_664b_3p	=	'miR-664b-3p'
hsa_miR_6716_3p	=	'miR-6716-3p'
hsa_miR_6722_5p	=	'miR-6722-5p'
hsa_miR_6730_3p	=	'miR-6730-3p'
hsa_miR_6743_3p	=	'miR-6743-3p'
hsa_miR_6771_5p	=	'miR-6771-5p'
hsa_miR_6779_3p	=	'miR-6779-3p'
hsa_miR_6794_3p	=	'miR-6794-3p'
hsa_miR_6804_5p	=	'miR-6804-5p'
hsa_miR_6806_5p	=	'miR-6806-5p'
hsa_miR_6817_5p	=	'miR-6817-5p'
hsa_miR_6826_5p	=	'miR-6826-5p'
hsa_miR_6865_3p	=	'miR-6865-3p'
hsa_miR_6872_3p	=	'miR-6872-3p'
hsa_miR_6886_3p	=	'miR-6886-3p'
hsa_miR_6891_3p	=	'miR-6891-3p'
hsa_miR_6895_5p	=	'miR-6895-5p'
hsa_miR_7108_3p	=	'miR-7108-3p'
hsa_miR_7111_3p	=	'miR-7111-3p'
hsa_miR_7159_5p	=	'miR-7159-5p'
hsa_miR_744_5p	=	'miR-744-5p'
hsa_miR_770_5p	=	'miR-770-5p'
hsa_miR_7704	=	'miR-7704'
hsa_miR_8063	=	'miR-8063'
hsa_miR_8077	=	'miR-8077'
hsa_miR_874_5p	=	'miR-874-5p'
hsa_miR_887_3p	=	'miR-887-3p'
hsa_miR_940	=	'miR-940'
hsa_miR_95_3p	=	'miR-95-3p'
hsa_miR_98_5p	=	'miR-98-5p'
hsa_miR_99a_3p	=	'miR-99a-3p'
hsa_miR_99a_5p	=	'miR-99a-5p'
hsa_miR_99b_5p	=	'miR-99b-5p'

##Predictions

predictions_hsa_let_7a_5p	=	getPredictedTargets(  hsa_let_7a_5p , species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_let_7b_5p	=	getPredictedTargets(	hsa_let_7b_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_let_7c_5p	=	getPredictedTargets(	hsa_let_7c_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_let_7d_5p	=	getPredictedTargets(	hsa_let_7d_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_let_7e_5p	=	getPredictedTargets(	hsa_let_7e_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_let_7f_5p	=	getPredictedTargets(	hsa_let_7f_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_let_7g_5p	=	getPredictedTargets(	hsa_let_7g_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_let_7i_5p	=	getPredictedTargets(	hsa_let_7i_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_1_3p	=	getPredictedTargets(	hsa_miR_1_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_100_5p	=	getPredictedTargets(	hsa_miR_100_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_101_3p	=	getPredictedTargets(	hsa_miR_101_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_101_5p	=	getPredictedTargets(	hsa_miR_101_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_103a_3p	=	getPredictedTargets(	hsa_miR_103a_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_106b_3p	=	getPredictedTargets(	hsa_miR_106b_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_107	=	getPredictedTargets(	hsa_miR_107	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_10a_5p	=	getPredictedTargets(	hsa_miR_10a_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_10b_5p	=	getPredictedTargets(	hsa_miR_10b_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_1227_5p	=	getPredictedTargets(	hsa_miR_1227_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_1238_5p	=	getPredictedTargets(	hsa_miR_1238_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_1247_5p	=	getPredictedTargets(	hsa_miR_1247_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_125a_5p	=	getPredictedTargets(	hsa_miR_125a_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_125b_5p	=	getPredictedTargets(	hsa_miR_125b_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_126_3p	=	getPredictedTargets(	hsa_miR_126_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_126_5p	=	getPredictedTargets(	hsa_miR_126_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_128_3p	=	getPredictedTargets(	hsa_miR_128_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_129_2_3p	=	getPredictedTargets(	hsa_miR_129_2_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_1296_5p	=	getPredictedTargets(	hsa_miR_1296_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_1306_3p	=	getPredictedTargets(	hsa_miR_1306_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_130a_3p	=	getPredictedTargets(	hsa_miR_130a_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_133a_3p	=	getPredictedTargets(	hsa_miR_133a_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_133a_5p	=	getPredictedTargets(	hsa_miR_133a_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_133b	=	getPredictedTargets(	hsa_miR_133b	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_135a_5p	=	getPredictedTargets(	hsa_miR_135a_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_138_5p	=	getPredictedTargets(	hsa_miR_138_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_139_3p	=	getPredictedTargets(	hsa_miR_139_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_139_5p	=	getPredictedTargets(	hsa_miR_139_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_140_3p	=	getPredictedTargets(	hsa_miR_140_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_140_5p	=	getPredictedTargets(	hsa_miR_140_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_142_3p	=	getPredictedTargets(	hsa_miR_142_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_142_5p	=	getPredictedTargets(	hsa_miR_142_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_143_3p	=	getPredictedTargets(	hsa_miR_143_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_143_5p	=	getPredictedTargets(	hsa_miR_143_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_144_3p	=	getPredictedTargets(	hsa_miR_144_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_144_5p	=	getPredictedTargets(	hsa_miR_144_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_145_3p	=	getPredictedTargets(	hsa_miR_145_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_145_5p	=	getPredictedTargets(	hsa_miR_145_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_146a_5p	=	getPredictedTargets(	hsa_miR_146a_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_146b_5p	=	getPredictedTargets(	hsa_miR_146b_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_147b	=	getPredictedTargets(	hsa_miR_147b	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_150_5p	=	getPredictedTargets(	hsa_miR_150_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_151a_5p	=	getPredictedTargets(	hsa_miR_151a_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_151b	=	getPredictedTargets(	hsa_miR_151b	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_152_3p	=	getPredictedTargets(	hsa_miR_152_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_1537_3p	=	getPredictedTargets(	hsa_miR_1537_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_15a_5p	=	getPredictedTargets(	hsa_miR_15a_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_15b_3p	=	getPredictedTargets(	hsa_miR_15b_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_15b_5p	=	getPredictedTargets(	hsa_miR_15b_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_16_2_3p	=	getPredictedTargets(	hsa_miR_16_2_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_16_5p	=	getPredictedTargets(	hsa_miR_16_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_181a_2_3p	=	getPredictedTargets(	hsa_miR_181a_2_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_181a_3p	=	getPredictedTargets(	hsa_miR_181a_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_182_3p	=	getPredictedTargets(	hsa_miR_182_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_183_3p	=	getPredictedTargets(	hsa_miR_183_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_184	=	getPredictedTargets(	hsa_miR_184	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_185_5p	=	getPredictedTargets(	hsa_miR_185_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_186_5p	=	getPredictedTargets(	hsa_miR_186_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_187_5p	=	getPredictedTargets(	hsa_miR_187_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_18a_5p	=	getPredictedTargets(	hsa_miR_18a_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_18b_5p	=	getPredictedTargets(	hsa_miR_18b_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_190a_5p	=	getPredictedTargets(	hsa_miR_190a_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_191_5p	=	getPredictedTargets(	hsa_miR_191_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_1913	=	getPredictedTargets(	hsa_miR_1913	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_193a_5p	=	getPredictedTargets(	hsa_miR_193a_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_195_5p	=	getPredictedTargets(	hsa_miR_195_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_199a_3p	=	getPredictedTargets(	hsa_miR_199a_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_199a_5p	=	getPredictedTargets(	hsa_miR_199a_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_199b_5p	=	getPredictedTargets(	hsa_miR_199b_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_203a_3p	=	getPredictedTargets(	hsa_miR_203a_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_204_5p	=	getPredictedTargets(	hsa_miR_204_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_205_3p	=	getPredictedTargets(	hsa_miR_205_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_20a_5p	=	getPredictedTargets(	hsa_miR_20a_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_20b_5p	=	getPredictedTargets(	hsa_miR_20b_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_21_3p	=	getPredictedTargets(	hsa_miR_21_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_214_3p	=	getPredictedTargets(	hsa_miR_214_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_214_5p	=	getPredictedTargets(	hsa_miR_214_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_218_5p	=	getPredictedTargets(	hsa_miR_218_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_22_5p	=	getPredictedTargets(	hsa_miR_22_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_221_3p	=	getPredictedTargets(	hsa_miR_221_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_221_5p	=	getPredictedTargets(	hsa_miR_221_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_222_3p	=	getPredictedTargets(	hsa_miR_222_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_223_3p	=	getPredictedTargets(	hsa_miR_223_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_223_5p	=	getPredictedTargets(	hsa_miR_223_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_224_3p	=	getPredictedTargets(	hsa_miR_224_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_2277_3p	=	getPredictedTargets(	hsa_miR_2277_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_23a_3p	=	getPredictedTargets(	hsa_miR_23a_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_23a_5p	=	getPredictedTargets(	hsa_miR_23a_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_23b_3p	=	getPredictedTargets(	hsa_miR_23b_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_24_3p	=	getPredictedTargets(	hsa_miR_24_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_26a_1_3p	=	getPredictedTargets(	hsa_miR_26a_1_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_26a_5p	=	getPredictedTargets(	hsa_miR_26a_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_26b_3p	=	getPredictedTargets(	hsa_miR_26b_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_26b_5p	=	getPredictedTargets(	hsa_miR_26b_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_27a_3p	=	getPredictedTargets(	hsa_miR_27a_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_27a_5p	=	getPredictedTargets(	hsa_miR_27a_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_27b_3p	=	getPredictedTargets(	hsa_miR_27b_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_28_5p	=	getPredictedTargets(	hsa_miR_28_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_2861	=	getPredictedTargets(	hsa_miR_2861	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_29a_3p	=	getPredictedTargets(	hsa_miR_29a_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_29a_5p	=	getPredictedTargets(	hsa_miR_29a_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_29b_1_5p	=	getPredictedTargets(	hsa_miR_29b_1_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_29b_2_5p	=	getPredictedTargets(	hsa_miR_29b_2_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_29b_3p	=	getPredictedTargets(	hsa_miR_29b_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_29c_3p	=	getPredictedTargets(	hsa_miR_29c_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_29c_5p	=	getPredictedTargets(	hsa_miR_29c_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_3065_3p	=	getPredictedTargets(	hsa_miR_3065_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_3065_5p	=	getPredictedTargets(	hsa_miR_3065_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_30a_3p	=	getPredictedTargets(	hsa_miR_30a_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_30a_5p	=	getPredictedTargets(	hsa_miR_30a_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_30b_5p	=	getPredictedTargets(	hsa_miR_30b_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_30c_2_3p	=	getPredictedTargets(	hsa_miR_30c_2_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_30c_5p	=	getPredictedTargets(	hsa_miR_30c_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_30d_5p	=	getPredictedTargets(	hsa_miR_30d_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_30e_3p	=	getPredictedTargets(	hsa_miR_30e_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_3149	=	getPredictedTargets(	hsa_miR_3149	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_3188	=	getPredictedTargets(	hsa_miR_3188	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_32_5p	=	getPredictedTargets(	hsa_miR_32_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_324_3p	=	getPredictedTargets(	hsa_miR_324_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_326	=	getPredictedTargets(	hsa_miR_326	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_328_3p	=	getPredictedTargets(	hsa_miR_328_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_331_3p	=	getPredictedTargets(	hsa_miR_331_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_335_5p	=	getPredictedTargets(	hsa_miR_335_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_338_3p	=	getPredictedTargets(	hsa_miR_338_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_339_5p	=	getPredictedTargets(	hsa_miR_339_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_340_5p	=	getPredictedTargets(	hsa_miR_340_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_342_3p	=	getPredictedTargets(	hsa_miR_342_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_342_5p	=	getPredictedTargets(	hsa_miR_342_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_34a_3p	=	getPredictedTargets(	hsa_miR_34a_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_34a_5p	=	getPredictedTargets(	hsa_miR_34a_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_34b_5p	=	getPredictedTargets(	hsa_miR_34b_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_34c_5p	=	getPredictedTargets(	hsa_miR_34c_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_3607_3p	=	getPredictedTargets(	hsa_miR_3607_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_361_3p	=	getPredictedTargets(	hsa_miR_361_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_362_3p	=	getPredictedTargets(	hsa_miR_362_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_362_5p	=	getPredictedTargets(	hsa_miR_362_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_3620_5p	=	getPredictedTargets(	hsa_miR_3620_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_363_3p	=	getPredictedTargets(	hsa_miR_363_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_3659	=	getPredictedTargets(	hsa_miR_3659	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_365a_3p	=	getPredictedTargets(	hsa_miR_365a_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_3660	=	getPredictedTargets(	hsa_miR_3660	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_3663_5p	=	getPredictedTargets(	hsa_miR_3663_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_371b_5p	=	getPredictedTargets(	hsa_miR_371b_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_374a_5p	=	getPredictedTargets(	hsa_miR_374a_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_374b_5p	=	getPredictedTargets(	hsa_miR_374b_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_374c_5p	=	getPredictedTargets(	hsa_miR_374c_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_424_5p	=	getPredictedTargets(	hsa_miR_424_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4252	=	getPredictedTargets(	hsa_miR_4252	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4254	=	getPredictedTargets(	hsa_miR_4254	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4284	=	getPredictedTargets(	hsa_miR_4284	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4290	=	getPredictedTargets(	hsa_miR_4290	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4306	=	getPredictedTargets(	hsa_miR_4306	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4317	=	getPredictedTargets(	hsa_miR_4317	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4318	=	getPredictedTargets(	hsa_miR_4318	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4324	=	getPredictedTargets(	hsa_miR_4324	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4328	=	getPredictedTargets(	hsa_miR_4328	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4440	=	getPredictedTargets(	hsa_miR_4440	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4443	=	getPredictedTargets(	hsa_miR_4443	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4481	=	getPredictedTargets(	hsa_miR_4481	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_449a	=	getPredictedTargets(	hsa_miR_449a	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_450a_5p	=	getPredictedTargets(	hsa_miR_450a_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4516	=	getPredictedTargets(	hsa_miR_4516	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_451a	=	getPredictedTargets(	hsa_miR_451a	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_452_5p	=	getPredictedTargets(	hsa_miR_452_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4521	=	getPredictedTargets(	hsa_miR_4521	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4532	=	getPredictedTargets(	hsa_miR_4532	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_454_3p	=	getPredictedTargets(	hsa_miR_454_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_455_3p	=	getPredictedTargets(	hsa_miR_455_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_455_5p	=	getPredictedTargets(	hsa_miR_455_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4634	=	getPredictedTargets(	hsa_miR_4634	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4655_3p	=	getPredictedTargets(	hsa_miR_4655_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4695_3p	=	getPredictedTargets(	hsa_miR_4695_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4716_5p	=	getPredictedTargets(	hsa_miR_4716_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4730	=	getPredictedTargets(	hsa_miR_4730	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4731_3p	=	getPredictedTargets(	hsa_miR_4731_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4763_5p	=	getPredictedTargets(	hsa_miR_4763_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4770	=	getPredictedTargets(	hsa_miR_4770	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_4793_3p	=	getPredictedTargets(	hsa_miR_4793_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_483_3p	=	getPredictedTargets(	hsa_miR_483_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_486_3p	=	getPredictedTargets(	hsa_miR_486_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_486_5p	=	getPredictedTargets(	hsa_miR_486_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_489_3p	=	getPredictedTargets(	hsa_miR_489_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_490_3p	=	getPredictedTargets(	hsa_miR_490_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_491_3p	=	getPredictedTargets(	hsa_miR_491_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_497_3p	=	getPredictedTargets(	hsa_miR_497_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_497_5p	=	getPredictedTargets(	hsa_miR_497_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_5001_5p	=	getPredictedTargets(	hsa_miR_5001_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_500a_3p	=	getPredictedTargets(	hsa_miR_500a_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_500b_5p	=	getPredictedTargets(	hsa_miR_500b_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_501_5p	=	getPredictedTargets(	hsa_miR_501_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_502_3p	=	getPredictedTargets(	hsa_miR_502_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_504_3p	=	getPredictedTargets(	hsa_miR_504_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_505_5p	=	getPredictedTargets(	hsa_miR_505_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_511_3p	=	getPredictedTargets(	hsa_miR_511_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_516b_5p	=	getPredictedTargets(	hsa_miR_516b_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_517_5p	=	getPredictedTargets(	hsa_miR_517_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_517a_3p	=	getPredictedTargets(	hsa_miR_517a_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_517c_3p	=	getPredictedTargets(	hsa_miR_517c_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_521	=	getPredictedTargets(	hsa_miR_521	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_522_3p	=	getPredictedTargets(	hsa_miR_522_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_532_3p	=	getPredictedTargets(	hsa_miR_532_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_532_5p	=	getPredictedTargets(	hsa_miR_532_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_548aa	=	getPredictedTargets(	hsa_miR_548aa	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_548aw	=	getPredictedTargets(	hsa_miR_548aw	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_548b_3p	=	getPredictedTargets(	hsa_miR_548b_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_548c_3p	=	getPredictedTargets(	hsa_miR_548c_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_548f_3p	=	getPredictedTargets(	hsa_miR_548f_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_548q	=	getPredictedTargets(	hsa_miR_548q	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_548x_3p	=	getPredictedTargets(	hsa_miR_548x_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_551b_3p	=	getPredictedTargets(	hsa_miR_551b_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_5701	=	getPredictedTargets(	hsa_miR_5701	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_582_5p	=	getPredictedTargets(	hsa_miR_582_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_585_3p	=	getPredictedTargets(	hsa_miR_585_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_590_5p	=	getPredictedTargets(	hsa_miR_590_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_595	=	getPredictedTargets(	hsa_miR_595	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_598_3p	=	getPredictedTargets(	hsa_miR_598_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_6068	=	getPredictedTargets(	hsa_miR_6068	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_6073	=	getPredictedTargets(	hsa_miR_6073	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_6075	=	getPredictedTargets(	hsa_miR_6075	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_610	=	getPredictedTargets(	hsa_miR_610	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_624_5p	=	getPredictedTargets(	hsa_miR_624_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_628_5p	=	getPredictedTargets(	hsa_miR_628_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_642b_5p	=	getPredictedTargets(	hsa_miR_642b_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_6500_5p	=	getPredictedTargets(	hsa_miR_6500_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_6516_3p	=	getPredictedTargets(	hsa_miR_6516_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_652_3p	=	getPredictedTargets(	hsa_miR_652_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_653_3p	=	getPredictedTargets(	hsa_miR_653_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_660_5p	=	getPredictedTargets(	hsa_miR_660_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_664a_3p	=	getPredictedTargets(	hsa_miR_664a_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_664b_3p	=	getPredictedTargets(	hsa_miR_664b_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_6716_3p	=	getPredictedTargets(	hsa_miR_6716_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_6722_5p	=	getPredictedTargets(	hsa_miR_6722_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_6730_3p	=	getPredictedTargets(	hsa_miR_6730_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_6743_3p	=	getPredictedTargets(	hsa_miR_6743_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_6771_5p	=	getPredictedTargets(	hsa_miR_6771_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_6779_3p	=	getPredictedTargets(	hsa_miR_6779_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_6794_3p	=	getPredictedTargets(	hsa_miR_6794_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_6804_5p	=	getPredictedTargets(	hsa_miR_6804_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_6806_5p	=	getPredictedTargets(	hsa_miR_6806_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_6817_5p	=	getPredictedTargets(	hsa_miR_6817_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_6826_5p	=	getPredictedTargets(	hsa_miR_6826_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_6865_3p	=	getPredictedTargets(	hsa_miR_6865_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_6872_3p	=	getPredictedTargets(	hsa_miR_6872_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_6886_3p	=	getPredictedTargets(	hsa_miR_6886_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_6891_3p	=	getPredictedTargets(	hsa_miR_6891_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_6895_5p	=	getPredictedTargets(	hsa_miR_6895_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_7108_3p	=	getPredictedTargets(	hsa_miR_7108_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_7111_3p	=	getPredictedTargets(	hsa_miR_7111_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_7159_5p	=	getPredictedTargets(	hsa_miR_7159_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_744_5p	=	getPredictedTargets(	hsa_miR_744_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_770_5p	=	getPredictedTargets(	hsa_miR_770_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_7704	=	getPredictedTargets(	hsa_miR_7704	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_8063	=	getPredictedTargets(	hsa_miR_8063	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_8077	=	getPredictedTargets(	hsa_miR_8077	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_874_5p	=	getPredictedTargets(	hsa_miR_874_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_887_3p	=	getPredictedTargets(	hsa_miR_887_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_940	=	getPredictedTargets(	hsa_miR_940	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_95_3p	=	getPredictedTargets(	hsa_miR_95_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_98_5p	=	getPredictedTargets(	hsa_miR_98_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_99a_3p	=	getPredictedTargets(	hsa_miR_99a_3p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_99a_5p	=	getPredictedTargets(	hsa_miR_99a_5p	, species = 'hsa', method = 'geom', min_src = 2)
predictions_hsa_miR_99b_5p	=	getPredictedTargets(	hsa_miR_99b_5p	, species = 'hsa', method = 'geom', min_src = 2)

## RankedGenes

rankedGenes_hsa_let_7a_5p	=		predictions_hsa_let_7a_5p	[,'rank_product']
rankedGenes_hsa_let_7b_5p	=		predictions_hsa_let_7b_5p	[,'rank_product']
rankedGenes_hsa_let_7c_5p	=		predictions_hsa_let_7c_5p	[,'rank_product']
rankedGenes_hsa_let_7d_5p	=		predictions_hsa_let_7d_5p	[,'rank_product']
rankedGenes_hsa_let_7e_5p	=		predictions_hsa_let_7e_5p	[,'rank_product']
rankedGenes_hsa_let_7f_5p	=		predictions_hsa_let_7f_5p	[,'rank_product']
rankedGenes_hsa_let_7g_5p	=		predictions_hsa_let_7g_5p	[,'rank_product']
rankedGenes_hsa_let_7i_5p	=		predictions_hsa_let_7i_5p	[,'rank_product']
rankedGenes_hsa_miR_1_3p	=		predictions_hsa_miR_1_3p	[,'rank_product']
rankedGenes_hsa_miR_100_5p	=		predictions_hsa_miR_100_5p	[,'rank_product']
rankedGenes_hsa_miR_101_3p	=		predictions_hsa_miR_101_3p	[,'rank_product']
rankedGenes_hsa_miR_101_5p	=		predictions_hsa_miR_101_5p	[,'rank_product']
rankedGenes_hsa_miR_103a_3p	=		predictions_hsa_miR_103a_3p	[,'rank_product']
rankedGenes_hsa_miR_106b_3p	=		predictions_hsa_miR_106b_3p	[,'rank_product']
rankedGenes_hsa_miR_107	=		predictions_hsa_miR_107	[,'rank_product']
rankedGenes_hsa_miR_10a_5p	=		predictions_hsa_miR_10a_5p	[,'rank_product']
rankedGenes_hsa_miR_10b_5p	=		predictions_hsa_miR_10b_5p	[,'rank_product']
rankedGenes_hsa_miR_1227_5p	=		predictions_hsa_miR_1227_5p	[,'rank_product']
rankedGenes_hsa_miR_1238_5p	=		predictions_hsa_miR_1238_5p	[,'rank_product']
rankedGenes_hsa_miR_1247_5p	=		predictions_hsa_miR_1247_5p	[,'rank_product']
rankedGenes_hsa_miR_125a_5p	=		predictions_hsa_miR_125a_5p	[,'rank_product']
rankedGenes_hsa_miR_125b_5p	=		predictions_hsa_miR_125b_5p	[,'rank_product']
rankedGenes_hsa_miR_126_3p	=		predictions_hsa_miR_126_3p	[,'rank_product']
rankedGenes_hsa_miR_126_5p	=		predictions_hsa_miR_126_5p	[,'rank_product']
rankedGenes_hsa_miR_128_3p	=		predictions_hsa_miR_128_3p	[,'rank_product']
rankedGenes_hsa_miR_129_2_3p	=		predictions_hsa_miR_129_2_3p	[,'rank_product']
rankedGenes_hsa_miR_1296_5p	=		predictions_hsa_miR_1296_5p	[,'rank_product']
rankedGenes_hsa_miR_1306_3p	=		predictions_hsa_miR_1306_3p	[,'rank_product']
rankedGenes_hsa_miR_130a_3p	=		predictions_hsa_miR_130a_3p	[,'rank_product']
rankedGenes_hsa_miR_133a_3p	=		predictions_hsa_miR_133a_3p	[,'rank_product']
rankedGenes_hsa_miR_133a_5p	=		predictions_hsa_miR_133a_5p	[,'rank_product']
rankedGenes_hsa_miR_133b	=		predictions_hsa_miR_133b	[,'rank_product']
rankedGenes_hsa_miR_135a_5p	=		predictions_hsa_miR_135a_5p	[,'rank_product']
rankedGenes_hsa_miR_138_5p	=		predictions_hsa_miR_138_5p	[,'rank_product']
rankedGenes_hsa_miR_139_3p	=		predictions_hsa_miR_139_3p	[,'rank_product']
rankedGenes_hsa_miR_139_5p	=		predictions_hsa_miR_139_5p	[,'rank_product']
rankedGenes_hsa_miR_140_3p	=		predictions_hsa_miR_140_3p	[,'rank_product']
rankedGenes_hsa_miR_140_5p	=		predictions_hsa_miR_140_5p	[,'rank_product']
rankedGenes_hsa_miR_142_3p	=		predictions_hsa_miR_142_3p	[,'rank_product']
rankedGenes_hsa_miR_142_5p	=		predictions_hsa_miR_142_5p	[,'rank_product']
rankedGenes_hsa_miR_143_3p	=		predictions_hsa_miR_143_3p	[,'rank_product']
rankedGenes_hsa_miR_143_5p	=		predictions_hsa_miR_143_5p	[,'rank_product']
rankedGenes_hsa_miR_144_3p	=		predictions_hsa_miR_144_3p	[,'rank_product']
rankedGenes_hsa_miR_144_5p	=		predictions_hsa_miR_144_5p	[,'rank_product']
rankedGenes_hsa_miR_145_3p	=		predictions_hsa_miR_145_3p	[,'rank_product']
rankedGenes_hsa_miR_145_5p	=		predictions_hsa_miR_145_5p	[,'rank_product']
rankedGenes_hsa_miR_146a_5p	=		predictions_hsa_miR_146a_5p	[,'rank_product']
rankedGenes_hsa_miR_146b_5p	=		predictions_hsa_miR_146b_5p	[,'rank_product']
rankedGenes_hsa_miR_147b	=		predictions_hsa_miR_147b	[,'rank_product']
rankedGenes_hsa_miR_150_5p	=		predictions_hsa_miR_150_5p	[,'rank_product']
rankedGenes_hsa_miR_151a_5p	=		predictions_hsa_miR_151a_5p	[,'rank_product']
rankedGenes_hsa_miR_151b	=		predictions_hsa_miR_151b	[,'rank_product']
rankedGenes_hsa_miR_152_3p	=		predictions_hsa_miR_152_3p	[,'rank_product']
rankedGenes_hsa_miR_1537_3p	=		predictions_hsa_miR_1537_3p	[,'rank_product']
rankedGenes_hsa_miR_15a_5p	=		predictions_hsa_miR_15a_5p	[,'rank_product']
rankedGenes_hsa_miR_15b_3p	=		predictions_hsa_miR_15b_3p	[,'rank_product']
rankedGenes_hsa_miR_15b_5p	=		predictions_hsa_miR_15b_5p	[,'rank_product']
rankedGenes_hsa_miR_16_2_3p	=		predictions_hsa_miR_16_2_3p	[,'rank_product']
rankedGenes_hsa_miR_16_5p	=		predictions_hsa_miR_16_5p	[,'rank_product']
rankedGenes_hsa_miR_181a_2_3p	=		predictions_hsa_miR_181a_2_3p	[,'rank_product']
rankedGenes_hsa_miR_181a_3p	=		predictions_hsa_miR_181a_3p	[,'rank_product']
rankedGenes_hsa_miR_182_3p	=		predictions_hsa_miR_182_3p	[,'rank_product']
rankedGenes_hsa_miR_183_3p	=		predictions_hsa_miR_183_3p	[,'rank_product']
rankedGenes_hsa_miR_184	=		predictions_hsa_miR_184	[,'rank_product']
rankedGenes_hsa_miR_185_5p	=		predictions_hsa_miR_185_5p	[,'rank_product']
rankedGenes_hsa_miR_186_5p	=		predictions_hsa_miR_186_5p	[,'rank_product']
rankedGenes_hsa_miR_187_5p	=		predictions_hsa_miR_187_5p	[,'rank_product']
rankedGenes_hsa_miR_18a_5p	=		predictions_hsa_miR_18a_5p	[,'rank_product']
rankedGenes_hsa_miR_18b_5p	=		predictions_hsa_miR_18b_5p	[,'rank_product']
rankedGenes_hsa_miR_190a_5p	=		predictions_hsa_miR_190a_5p	[,'rank_product']
rankedGenes_hsa_miR_191_5p	=		predictions_hsa_miR_191_5p	[,'rank_product']
rankedGenes_hsa_miR_1913	=		predictions_hsa_miR_1913	[,'rank_product']
rankedGenes_hsa_miR_193a_5p	=		predictions_hsa_miR_193a_5p	[,'rank_product']
rankedGenes_hsa_miR_195_5p	=		predictions_hsa_miR_195_5p	[,'rank_product']
rankedGenes_hsa_miR_199a_3p	=		predictions_hsa_miR_199a_3p	[,'rank_product']
rankedGenes_hsa_miR_199a_5p	=		predictions_hsa_miR_199a_5p	[,'rank_product']
rankedGenes_hsa_miR_199b_5p	=		predictions_hsa_miR_199b_5p	[,'rank_product']
rankedGenes_hsa_miR_203a_3p	=		predictions_hsa_miR_203a_3p	[,'rank_product']
rankedGenes_hsa_miR_204_5p	=		predictions_hsa_miR_204_5p	[,'rank_product']
rankedGenes_hsa_miR_205_3p	=		predictions_hsa_miR_205_3p	[,'rank_product']
rankedGenes_hsa_miR_20a_5p	=		predictions_hsa_miR_20a_5p	[,'rank_product']
rankedGenes_hsa_miR_20b_5p	=		predictions_hsa_miR_20b_5p	[,'rank_product']
rankedGenes_hsa_miR_21_3p	=		predictions_hsa_miR_21_3p	[,'rank_product']
rankedGenes_hsa_miR_214_3p	=		predictions_hsa_miR_214_3p	[,'rank_product']
rankedGenes_hsa_miR_214_5p	=		predictions_hsa_miR_214_5p	[,'rank_product']
rankedGenes_hsa_miR_218_5p	=		predictions_hsa_miR_218_5p	[,'rank_product']
rankedGenes_hsa_miR_22_5p	=		predictions_hsa_miR_22_5p	[,'rank_product']
rankedGenes_hsa_miR_221_3p	=		predictions_hsa_miR_221_3p	[,'rank_product']
rankedGenes_hsa_miR_221_5p	=		predictions_hsa_miR_221_5p	[,'rank_product']
rankedGenes_hsa_miR_222_3p	=		predictions_hsa_miR_222_3p	[,'rank_product']
rankedGenes_hsa_miR_223_3p	=		predictions_hsa_miR_223_3p	[,'rank_product']
rankedGenes_hsa_miR_223_5p	=		predictions_hsa_miR_223_5p	[,'rank_product']
rankedGenes_hsa_miR_224_3p	=		predictions_hsa_miR_224_3p	[,'rank_product']
rankedGenes_hsa_miR_2277_3p	=		predictions_hsa_miR_2277_3p	[,'rank_product']
rankedGenes_hsa_miR_23a_3p	=		predictions_hsa_miR_23a_3p	[,'rank_product']
rankedGenes_hsa_miR_23a_5p	=		predictions_hsa_miR_23a_5p	[,'rank_product']
rankedGenes_hsa_miR_23b_3p	=		predictions_hsa_miR_23b_3p	[,'rank_product']
rankedGenes_hsa_miR_24_3p	=		predictions_hsa_miR_24_3p	[,'rank_product']
rankedGenes_hsa_miR_26a_1_3p	=		predictions_hsa_miR_26a_1_3p	[,'rank_product']
rankedGenes_hsa_miR_26a_5p	=		predictions_hsa_miR_26a_5p	[,'rank_product']
rankedGenes_hsa_miR_26b_3p	=		predictions_hsa_miR_26b_3p	[,'rank_product']
rankedGenes_hsa_miR_26b_5p	=		predictions_hsa_miR_26b_5p	[,'rank_product']
rankedGenes_hsa_miR_27a_3p	=		predictions_hsa_miR_27a_3p	[,'rank_product']
rankedGenes_hsa_miR_27a_5p	=		predictions_hsa_miR_27a_5p	[,'rank_product']
rankedGenes_hsa_miR_27b_3p	=		predictions_hsa_miR_27b_3p	[,'rank_product']
rankedGenes_hsa_miR_28_5p	=		predictions_hsa_miR_28_5p	[,'rank_product']
rankedGenes_hsa_miR_2861	=		predictions_hsa_miR_2861	[,'rank_product']
rankedGenes_hsa_miR_29a_3p	=		predictions_hsa_miR_29a_3p	[,'rank_product']
rankedGenes_hsa_miR_29a_5p	=		predictions_hsa_miR_29a_5p	[,'rank_product']
rankedGenes_hsa_miR_29b_1_5p	=		predictions_hsa_miR_29b_1_5p	[,'rank_product']
rankedGenes_hsa_miR_29b_2_5p	=		predictions_hsa_miR_29b_2_5p	[,'rank_product']
rankedGenes_hsa_miR_29b_3p	=		predictions_hsa_miR_29b_3p	[,'rank_product']
rankedGenes_hsa_miR_29c_3p	=		predictions_hsa_miR_29c_3p	[,'rank_product']
rankedGenes_hsa_miR_29c_5p	=		predictions_hsa_miR_29c_5p	[,'rank_product']
rankedGenes_hsa_miR_3065_3p	=		predictions_hsa_miR_3065_3p	[,'rank_product']
rankedGenes_hsa_miR_3065_5p	=		predictions_hsa_miR_3065_5p	[,'rank_product']
rankedGenes_hsa_miR_30a_3p	=		predictions_hsa_miR_30a_3p	[,'rank_product']
rankedGenes_hsa_miR_30a_5p	=		predictions_hsa_miR_30a_5p	[,'rank_product']
rankedGenes_hsa_miR_30b_5p	=		predictions_hsa_miR_30b_5p	[,'rank_product']
rankedGenes_hsa_miR_30c_2_3p	=		predictions_hsa_miR_30c_2_3p	[,'rank_product']
rankedGenes_hsa_miR_30c_5p	=		predictions_hsa_miR_30c_5p	[,'rank_product']
rankedGenes_hsa_miR_30d_5p	=		predictions_hsa_miR_30d_5p	[,'rank_product']
rankedGenes_hsa_miR_30e_3p	=		predictions_hsa_miR_30e_3p	[,'rank_product']
rankedGenes_hsa_miR_3149	=		predictions_hsa_miR_3149	[,'rank_product']
rankedGenes_hsa_miR_3188	=		predictions_hsa_miR_3188	[,'rank_product']
rankedGenes_hsa_miR_32_5p	=		predictions_hsa_miR_32_5p	[,'rank_product']
rankedGenes_hsa_miR_324_3p	=		predictions_hsa_miR_324_3p	[,'rank_product']
rankedGenes_hsa_miR_326	=		predictions_hsa_miR_326	[,'rank_product']
rankedGenes_hsa_miR_328_3p	=		predictions_hsa_miR_328_3p	[,'rank_product']
rankedGenes_hsa_miR_331_3p	=		predictions_hsa_miR_331_3p	[,'rank_product']
rankedGenes_hsa_miR_335_5p	=		predictions_hsa_miR_335_5p	[,'rank_product']
rankedGenes_hsa_miR_338_3p	=		predictions_hsa_miR_338_3p	[,'rank_product']
rankedGenes_hsa_miR_339_5p	=		predictions_hsa_miR_339_5p	[,'rank_product']
rankedGenes_hsa_miR_340_5p	=		predictions_hsa_miR_340_5p	[,'rank_product']
rankedGenes_hsa_miR_342_3p	=		predictions_hsa_miR_342_3p	[,'rank_product']
rankedGenes_hsa_miR_342_5p	=		predictions_hsa_miR_342_5p	[,'rank_product']
rankedGenes_hsa_miR_34a_3p	=		predictions_hsa_miR_34a_3p	[,'rank_product']
rankedGenes_hsa_miR_34a_5p	=		predictions_hsa_miR_34a_5p	[,'rank_product']
rankedGenes_hsa_miR_34b_5p	=		predictions_hsa_miR_34b_5p	[,'rank_product']
rankedGenes_hsa_miR_34c_5p	=		predictions_hsa_miR_34c_5p	[,'rank_product']
rankedGenes_hsa_miR_3607_3p	=		predictions_hsa_miR_3607_3p	[,'rank_product']
rankedGenes_hsa_miR_361_3p	=		predictions_hsa_miR_361_3p	[,'rank_product']
rankedGenes_hsa_miR_362_3p	=		predictions_hsa_miR_362_3p	[,'rank_product']
rankedGenes_hsa_miR_362_5p	=		predictions_hsa_miR_362_5p	[,'rank_product']
rankedGenes_hsa_miR_3620_5p	=		predictions_hsa_miR_3620_5p	[,'rank_product']
rankedGenes_hsa_miR_363_3p	=		predictions_hsa_miR_363_3p	[,'rank_product']
rankedGenes_hsa_miR_3659	=		predictions_hsa_miR_3659	[,'rank_product']
rankedGenes_hsa_miR_365a_3p	=		predictions_hsa_miR_365a_3p	[,'rank_product']
rankedGenes_hsa_miR_3660	=		predictions_hsa_miR_3660	[,'rank_product']
rankedGenes_hsa_miR_3663_5p	=		predictions_hsa_miR_3663_5p	[,'rank_product']
rankedGenes_hsa_miR_371b_5p	=		predictions_hsa_miR_371b_5p	[,'rank_product']
rankedGenes_hsa_miR_374a_5p	=		predictions_hsa_miR_374a_5p	[,'rank_product']
rankedGenes_hsa_miR_374b_5p	=		predictions_hsa_miR_374b_5p	[,'rank_product']
rankedGenes_hsa_miR_374c_5p	=		predictions_hsa_miR_374c_5p	[,'rank_product']
rankedGenes_hsa_miR_424_5p	=		predictions_hsa_miR_424_5p	[,'rank_product']
rankedGenes_hsa_miR_4252	=		predictions_hsa_miR_4252	[,'rank_product']
rankedGenes_hsa_miR_4254	=		predictions_hsa_miR_4254	[,'rank_product']
rankedGenes_hsa_miR_4284	=		predictions_hsa_miR_4284	[,'rank_product']
rankedGenes_hsa_miR_4290	=		predictions_hsa_miR_4290	[,'rank_product']
rankedGenes_hsa_miR_4306	=		predictions_hsa_miR_4306	[,'rank_product']
rankedGenes_hsa_miR_4317	=		predictions_hsa_miR_4317	[,'rank_product']
rankedGenes_hsa_miR_4318	=		predictions_hsa_miR_4318	[,'rank_product']
rankedGenes_hsa_miR_4324	=		predictions_hsa_miR_4324	[,'rank_product']
rankedGenes_hsa_miR_4328	=		predictions_hsa_miR_4328	[,'rank_product']
rankedGenes_hsa_miR_4440	=		predictions_hsa_miR_4440	[,'rank_product']
rankedGenes_hsa_miR_4443	=		predictions_hsa_miR_4443	[,'rank_product']
rankedGenes_hsa_miR_4481	=		predictions_hsa_miR_4481	[,'rank_product']
rankedGenes_hsa_miR_449a	=		predictions_hsa_miR_449a	[,'rank_product']
rankedGenes_hsa_miR_450a_5p	=		predictions_hsa_miR_450a_5p	[,'rank_product']
rankedGenes_hsa_miR_4516	=		predictions_hsa_miR_4516	[,'rank_product']
rankedGenes_hsa_miR_451a	=		predictions_hsa_miR_451a	[,'rank_product']
rankedGenes_hsa_miR_452_5p	=		predictions_hsa_miR_452_5p	[,'rank_product']
rankedGenes_hsa_miR_4521	=		predictions_hsa_miR_4521	[,'rank_product']
rankedGenes_hsa_miR_4532	=		predictions_hsa_miR_4532	[,'rank_product']
rankedGenes_hsa_miR_454_3p	=		predictions_hsa_miR_454_3p	[,'rank_product']
rankedGenes_hsa_miR_455_3p	=		predictions_hsa_miR_455_3p	[,'rank_product']
rankedGenes_hsa_miR_455_5p	=		predictions_hsa_miR_455_5p	[,'rank_product']
rankedGenes_hsa_miR_4634	=		predictions_hsa_miR_4634	[,'rank_product']
rankedGenes_hsa_miR_4655_3p	=		predictions_hsa_miR_4655_3p	[,'rank_product']
rankedGenes_hsa_miR_4695_3p	=		predictions_hsa_miR_4695_3p	[,'rank_product']
rankedGenes_hsa_miR_4716_5p	=		predictions_hsa_miR_4716_5p	[,'rank_product']
rankedGenes_hsa_miR_4730	=		predictions_hsa_miR_4730	[,'rank_product']
rankedGenes_hsa_miR_4731_3p	=		predictions_hsa_miR_4731_3p	[,'rank_product']
rankedGenes_hsa_miR_4763_5p	=		predictions_hsa_miR_4763_5p	[,'rank_product']
rankedGenes_hsa_miR_4770	=		predictions_hsa_miR_4770	[,'rank_product']
rankedGenes_hsa_miR_4793_3p	=		predictions_hsa_miR_4793_3p	[,'rank_product']
rankedGenes_hsa_miR_483_3p	=		predictions_hsa_miR_483_3p	[,'rank_product']
rankedGenes_hsa_miR_486_3p	=		predictions_hsa_miR_486_3p	[,'rank_product']
rankedGenes_hsa_miR_486_5p	=		predictions_hsa_miR_486_5p	[,'rank_product']
rankedGenes_hsa_miR_489_3p	=		predictions_hsa_miR_489_3p	[,'rank_product']
rankedGenes_hsa_miR_490_3p	=		predictions_hsa_miR_490_3p	[,'rank_product']
rankedGenes_hsa_miR_491_3p	=		predictions_hsa_miR_491_3p	[,'rank_product']
rankedGenes_hsa_miR_497_3p	=		predictions_hsa_miR_497_3p	[,'rank_product']
rankedGenes_hsa_miR_497_5p	=		predictions_hsa_miR_497_5p	[,'rank_product']
rankedGenes_hsa_miR_5001_5p	=		predictions_hsa_miR_5001_5p	[,'rank_product']
rankedGenes_hsa_miR_500a_3p	=		predictions_hsa_miR_500a_3p	[,'rank_product']
rankedGenes_hsa_miR_500b_5p	=		predictions_hsa_miR_500b_5p	[,'rank_product']
rankedGenes_hsa_miR_501_5p	=		predictions_hsa_miR_501_5p	[,'rank_product']
rankedGenes_hsa_miR_502_3p	=		predictions_hsa_miR_502_3p	[,'rank_product']
rankedGenes_hsa_miR_504_3p	=		predictions_hsa_miR_504_3p	[,'rank_product']
rankedGenes_hsa_miR_505_5p	=		predictions_hsa_miR_505_5p	[,'rank_product']
rankedGenes_hsa_miR_511_3p	=		predictions_hsa_miR_511_3p	[,'rank_product']
rankedGenes_hsa_miR_516b_5p	=		predictions_hsa_miR_516b_5p	[,'rank_product']
rankedGenes_hsa_miR_517_5p	=		predictions_hsa_miR_517_5p	[,'rank_product']
rankedGenes_hsa_miR_517a_3p	=		predictions_hsa_miR_517a_3p	[,'rank_product']
rankedGenes_hsa_miR_517c_3p	=		predictions_hsa_miR_517c_3p	[,'rank_product']
rankedGenes_hsa_miR_521	=		predictions_hsa_miR_521	[,'rank_product']
rankedGenes_hsa_miR_522_3p	=		predictions_hsa_miR_522_3p	[,'rank_product']
rankedGenes_hsa_miR_532_3p	=		predictions_hsa_miR_532_3p	[,'rank_product']
rankedGenes_hsa_miR_532_5p	=		predictions_hsa_miR_532_5p	[,'rank_product']
rankedGenes_hsa_miR_548aa	=		predictions_hsa_miR_548aa	[,'rank_product']
rankedGenes_hsa_miR_548aw	=		predictions_hsa_miR_548aw	[,'rank_product']
rankedGenes_hsa_miR_548b_3p	=		predictions_hsa_miR_548b_3p	[,'rank_product']
rankedGenes_hsa_miR_548c_3p	=		predictions_hsa_miR_548c_3p	[,'rank_product']
rankedGenes_hsa_miR_548f_3p	=		predictions_hsa_miR_548f_3p	[,'rank_product']
rankedGenes_hsa_miR_548q	=		predictions_hsa_miR_548q	[,'rank_product']
rankedGenes_hsa_miR_548x_3p	=		predictions_hsa_miR_548x_3p	[,'rank_product']
rankedGenes_hsa_miR_551b_3p	=		predictions_hsa_miR_551b_3p	[,'rank_product']
rankedGenes_hsa_miR_5701	=		predictions_hsa_miR_5701	[,'rank_product']
rankedGenes_hsa_miR_582_5p	=		predictions_hsa_miR_582_5p	[,'rank_product']
rankedGenes_hsa_miR_585_3p	=		predictions_hsa_miR_585_3p	[,'rank_product']
rankedGenes_hsa_miR_590_5p	=		predictions_hsa_miR_590_5p	[,'rank_product']
rankedGenes_hsa_miR_595	=		predictions_hsa_miR_595	[,'rank_product']
rankedGenes_hsa_miR_598_3p	=		predictions_hsa_miR_598_3p	[,'rank_product']
rankedGenes_hsa_miR_6068	=		predictions_hsa_miR_6068	[,'rank_product']
rankedGenes_hsa_miR_6073	=		predictions_hsa_miR_6073	[,'rank_product']
rankedGenes_hsa_miR_6075	=		predictions_hsa_miR_6075	[,'rank_product']
rankedGenes_hsa_miR_610	=		predictions_hsa_miR_610	[,'rank_product']
rankedGenes_hsa_miR_624_5p	=		predictions_hsa_miR_624_5p	[,'rank_product']
rankedGenes_hsa_miR_628_5p	=		predictions_hsa_miR_628_5p	[,'rank_product']
rankedGenes_hsa_miR_642b_5p	=		predictions_hsa_miR_642b_5p	[,'rank_product']
rankedGenes_hsa_miR_6500_5p	=		predictions_hsa_miR_6500_5p	[,'rank_product']
rankedGenes_hsa_miR_6516_3p	=		predictions_hsa_miR_6516_3p	[,'rank_product']
rankedGenes_hsa_miR_652_3p	=		predictions_hsa_miR_652_3p	[,'rank_product']
rankedGenes_hsa_miR_653_3p	=		predictions_hsa_miR_653_3p	[,'rank_product']
rankedGenes_hsa_miR_660_5p	=		predictions_hsa_miR_660_5p	[,'rank_product']
rankedGenes_hsa_miR_664a_3p	=		predictions_hsa_miR_664a_3p	[,'rank_product']
rankedGenes_hsa_miR_664b_3p	=		predictions_hsa_miR_664b_3p	[,'rank_product']
rankedGenes_hsa_miR_6716_3p	=		predictions_hsa_miR_6716_3p	[,'rank_product']
rankedGenes_hsa_miR_6722_5p	=		predictions_hsa_miR_6722_5p	[,'rank_product']
rankedGenes_hsa_miR_6730_3p	=		predictions_hsa_miR_6730_3p	[,'rank_product']
rankedGenes_hsa_miR_6743_3p	=		predictions_hsa_miR_6743_3p	[,'rank_product']
rankedGenes_hsa_miR_6771_5p	=		predictions_hsa_miR_6771_5p	[,'rank_product']
rankedGenes_hsa_miR_6779_3p	=		predictions_hsa_miR_6779_3p	[,'rank_product']
rankedGenes_hsa_miR_6794_3p	=		predictions_hsa_miR_6794_3p	[,'rank_product']
rankedGenes_hsa_miR_6804_5p	=		predictions_hsa_miR_6804_5p	[,'rank_product']
rankedGenes_hsa_miR_6806_5p	=		predictions_hsa_miR_6806_5p	[,'rank_product']
rankedGenes_hsa_miR_6817_5p	=		predictions_hsa_miR_6817_5p	[,'rank_product']
rankedGenes_hsa_miR_6826_5p	=		predictions_hsa_miR_6826_5p	[,'rank_product']
rankedGenes_hsa_miR_6865_3p	=		predictions_hsa_miR_6865_3p	[,'rank_product']
rankedGenes_hsa_miR_6872_3p	=		predictions_hsa_miR_6872_3p	[,'rank_product']
rankedGenes_hsa_miR_6886_3p	=		predictions_hsa_miR_6886_3p	[,'rank_product']
rankedGenes_hsa_miR_6891_3p	=		predictions_hsa_miR_6891_3p	[,'rank_product']
rankedGenes_hsa_miR_6895_5p	=		predictions_hsa_miR_6895_5p	[,'rank_product']
rankedGenes_hsa_miR_7108_3p	=		predictions_hsa_miR_7108_3p	[,'rank_product']
rankedGenes_hsa_miR_7111_3p	=		predictions_hsa_miR_7111_3p	[,'rank_product']
rankedGenes_hsa_miR_7159_5p	=		predictions_hsa_miR_7159_5p	[,'rank_product']
rankedGenes_hsa_miR_744_5p	=		predictions_hsa_miR_744_5p	[,'rank_product']
rankedGenes_hsa_miR_770_5p	=		predictions_hsa_miR_770_5p	[,'rank_product']
rankedGenes_hsa_miR_7704	=		predictions_hsa_miR_7704	[,'rank_product']
rankedGenes_hsa_miR_8063	=		predictions_hsa_miR_8063	[,'rank_product']
rankedGenes_hsa_miR_8077	=		predictions_hsa_miR_8077	[,'rank_product']
rankedGenes_hsa_miR_874_5p	=		predictions_hsa_miR_874_5p	[,'rank_product']
rankedGenes_hsa_miR_887_3p	=		predictions_hsa_miR_887_3p	[,'rank_product']
rankedGenes_hsa_miR_940	=		predictions_hsa_miR_940	[,'rank_product']
rankedGenes_hsa_miR_95_3p	=		predictions_hsa_miR_95_3p	[,'rank_product']
rankedGenes_hsa_miR_98_5p	=		predictions_hsa_miR_98_5p	[,'rank_product']
rankedGenes_hsa_miR_99a_3p	=		predictions_hsa_miR_99a_3p	[,'rank_product']
rankedGenes_hsa_miR_99a_5p	=		predictions_hsa_miR_99a_5p	[,'rank_product']
rankedGenes_hsa_miR_99b_5p	=		predictions_hsa_miR_99b_5p	[,'rank_product']

#FunctionX
selection_hsa_let_7a_5p	=	function(x) TRUE
selection_hsa_let_7b_5p	=	function(x) TRUE
selection_hsa_let_7c_5p	=	function(x) TRUE
selection_hsa_let_7d_5p	=	function(x) TRUE
selection_hsa_let_7e_5p	=	function(x) TRUE
selection_hsa_let_7f_5p	=	function(x) TRUE
selection_hsa_let_7g_5p	=	function(x) TRUE
selection_hsa_let_7i_5p	=	function(x) TRUE
selection_hsa_miR_1_3p	=	function(x) TRUE
selection_hsa_miR_100_5p	=	function(x) TRUE
selection_hsa_miR_101_3p	=	function(x) TRUE
selection_hsa_miR_101_5p	=	function(x) TRUE
selection_hsa_miR_103a_3p	=	function(x) TRUE
selection_hsa_miR_106b_3p	=	function(x) TRUE
selection_hsa_miR_107	=	function(x) TRUE
selection_hsa_miR_10a_5p	=	function(x) TRUE
selection_hsa_miR_10b_5p	=	function(x) TRUE
selection_hsa_miR_1227_5p	=	function(x) TRUE
selection_hsa_miR_1238_5p	=	function(x) TRUE
selection_hsa_miR_1247_5p	=	function(x) TRUE
selection_hsa_miR_125a_5p	=	function(x) TRUE
selection_hsa_miR_125b_5p	=	function(x) TRUE
selection_hsa_miR_126_3p	=	function(x) TRUE
selection_hsa_miR_126_5p	=	function(x) TRUE
selection_hsa_miR_128_3p	=	function(x) TRUE
selection_hsa_miR_129_2_3p	=	function(x) TRUE
selection_hsa_miR_1296_5p	=	function(x) TRUE
selection_hsa_miR_1306_3p	=	function(x) TRUE
selection_hsa_miR_130a_3p	=	function(x) TRUE
selection_hsa_miR_133a_3p	=	function(x) TRUE
selection_hsa_miR_133a_5p	=	function(x) TRUE
selection_hsa_miR_133b	=	function(x) TRUE
selection_hsa_miR_135a_5p	=	function(x) TRUE
selection_hsa_miR_138_5p	=	function(x) TRUE
selection_hsa_miR_139_3p	=	function(x) TRUE
selection_hsa_miR_139_5p	=	function(x) TRUE
selection_hsa_miR_140_3p	=	function(x) TRUE
selection_hsa_miR_140_5p	=	function(x) TRUE
selection_hsa_miR_142_3p	=	function(x) TRUE
selection_hsa_miR_142_5p	=	function(x) TRUE
selection_hsa_miR_143_3p	=	function(x) TRUE
selection_hsa_miR_143_5p	=	function(x) TRUE
selection_hsa_miR_144_3p	=	function(x) TRUE
selection_hsa_miR_144_5p	=	function(x) TRUE
selection_hsa_miR_145_3p	=	function(x) TRUE
selection_hsa_miR_145_5p	=	function(x) TRUE
selection_hsa_miR_146a_5p	=	function(x) TRUE
selection_hsa_miR_146b_5p	=	function(x) TRUE
selection_hsa_miR_147b	=	function(x) TRUE
selection_hsa_miR_150_5p	=	function(x) TRUE
selection_hsa_miR_151a_5p	=	function(x) TRUE
selection_hsa_miR_151b	=	function(x) TRUE
selection_hsa_miR_152_3p	=	function(x) TRUE
selection_hsa_miR_1537_3p	=	function(x) TRUE
selection_hsa_miR_15a_5p	=	function(x) TRUE
selection_hsa_miR_15b_3p	=	function(x) TRUE
selection_hsa_miR_15b_5p	=	function(x) TRUE
selection_hsa_miR_16_2_3p	=	function(x) TRUE
selection_hsa_miR_16_5p	=	function(x) TRUE
selection_hsa_miR_181a_2_3p	=	function(x) TRUE
selection_hsa_miR_181a_3p	=	function(x) TRUE
selection_hsa_miR_182_3p	=	function(x) TRUE
selection_hsa_miR_183_3p	=	function(x) TRUE
selection_hsa_miR_184	=	function(x) TRUE
selection_hsa_miR_185_5p	=	function(x) TRUE
selection_hsa_miR_186_5p	=	function(x) TRUE
selection_hsa_miR_187_5p	=	function(x) TRUE
selection_hsa_miR_18a_5p	=	function(x) TRUE
selection_hsa_miR_18b_5p	=	function(x) TRUE
selection_hsa_miR_190a_5p	=	function(x) TRUE
selection_hsa_miR_191_5p	=	function(x) TRUE
selection_hsa_miR_1913	=	function(x) TRUE
selection_hsa_miR_193a_5p	=	function(x) TRUE
selection_hsa_miR_195_5p	=	function(x) TRUE
selection_hsa_miR_199a_3p	=	function(x) TRUE
selection_hsa_miR_199a_5p	=	function(x) TRUE
selection_hsa_miR_199b_5p	=	function(x) TRUE
selection_hsa_miR_203a_3p	=	function(x) TRUE
selection_hsa_miR_204_5p	=	function(x) TRUE
selection_hsa_miR_205_3p	=	function(x) TRUE
selection_hsa_miR_20a_5p	=	function(x) TRUE
selection_hsa_miR_20b_5p	=	function(x) TRUE
selection_hsa_miR_21_3p	=	function(x) TRUE
selection_hsa_miR_214_3p	=	function(x) TRUE
selection_hsa_miR_214_5p	=	function(x) TRUE
selection_hsa_miR_218_5p	=	function(x) TRUE
selection_hsa_miR_22_5p	=	function(x) TRUE
selection_hsa_miR_221_3p	=	function(x) TRUE
selection_hsa_miR_221_5p	=	function(x) TRUE
selection_hsa_miR_222_3p	=	function(x) TRUE
selection_hsa_miR_223_3p	=	function(x) TRUE
selection_hsa_miR_223_5p	=	function(x) TRUE
selection_hsa_miR_224_3p	=	function(x) TRUE
selection_hsa_miR_2277_3p	=	function(x) TRUE
selection_hsa_miR_23a_3p	=	function(x) TRUE
selection_hsa_miR_23a_5p	=	function(x) TRUE
selection_hsa_miR_23b_3p	=	function(x) TRUE
selection_hsa_miR_24_3p	=	function(x) TRUE
selection_hsa_miR_26a_1_3p	=	function(x) TRUE
selection_hsa_miR_26a_5p	=	function(x) TRUE
selection_hsa_miR_26b_3p	=	function(x) TRUE
selection_hsa_miR_26b_5p	=	function(x) TRUE
selection_hsa_miR_27a_3p	=	function(x) TRUE
selection_hsa_miR_27a_5p	=	function(x) TRUE
selection_hsa_miR_27b_3p	=	function(x) TRUE
selection_hsa_miR_28_5p	=	function(x) TRUE
selection_hsa_miR_2861	=	function(x) TRUE
selection_hsa_miR_29a_3p	=	function(x) TRUE
selection_hsa_miR_29a_5p	=	function(x) TRUE
selection_hsa_miR_29b_1_5p	=	function(x) TRUE
selection_hsa_miR_29b_2_5p	=	function(x) TRUE
selection_hsa_miR_29b_3p	=	function(x) TRUE
selection_hsa_miR_29c_3p	=	function(x) TRUE
selection_hsa_miR_29c_5p	=	function(x) TRUE
selection_hsa_miR_3065_3p	=	function(x) TRUE
selection_hsa_miR_3065_5p	=	function(x) TRUE
selection_hsa_miR_30a_3p	=	function(x) TRUE
selection_hsa_miR_30a_5p	=	function(x) TRUE
selection_hsa_miR_30b_5p	=	function(x) TRUE
selection_hsa_miR_30c_2_3p	=	function(x) TRUE
selection_hsa_miR_30c_5p	=	function(x) TRUE
selection_hsa_miR_30d_5p	=	function(x) TRUE
selection_hsa_miR_30e_3p	=	function(x) TRUE
selection_hsa_miR_3149	=	function(x) TRUE
selection_hsa_miR_3188	=	function(x) TRUE
selection_hsa_miR_32_5p	=	function(x) TRUE
selection_hsa_miR_324_3p	=	function(x) TRUE
selection_hsa_miR_326	=	function(x) TRUE
selection_hsa_miR_328_3p	=	function(x) TRUE
selection_hsa_miR_331_3p	=	function(x) TRUE
selection_hsa_miR_335_5p	=	function(x) TRUE
selection_hsa_miR_338_3p	=	function(x) TRUE
selection_hsa_miR_339_5p	=	function(x) TRUE
selection_hsa_miR_340_5p	=	function(x) TRUE
selection_hsa_miR_342_3p	=	function(x) TRUE
selection_hsa_miR_342_5p	=	function(x) TRUE
selection_hsa_miR_34a_3p	=	function(x) TRUE
selection_hsa_miR_34a_5p	=	function(x) TRUE
selection_hsa_miR_34b_5p	=	function(x) TRUE
selection_hsa_miR_34c_5p	=	function(x) TRUE
selection_hsa_miR_3607_3p	=	function(x) TRUE
selection_hsa_miR_361_3p	=	function(x) TRUE
selection_hsa_miR_362_3p	=	function(x) TRUE
selection_hsa_miR_362_5p	=	function(x) TRUE
selection_hsa_miR_3620_5p	=	function(x) TRUE
selection_hsa_miR_363_3p	=	function(x) TRUE
selection_hsa_miR_3659	=	function(x) TRUE
selection_hsa_miR_365a_3p	=	function(x) TRUE
selection_hsa_miR_3660	=	function(x) TRUE
selection_hsa_miR_3663_5p	=	function(x) TRUE
selection_hsa_miR_371b_5p	=	function(x) TRUE
selection_hsa_miR_374a_5p	=	function(x) TRUE
selection_hsa_miR_374b_5p	=	function(x) TRUE
selection_hsa_miR_374c_5p	=	function(x) TRUE
selection_hsa_miR_424_5p	=	function(x) TRUE
selection_hsa_miR_4252	=	function(x) TRUE
selection_hsa_miR_4254	=	function(x) TRUE
selection_hsa_miR_4284	=	function(x) TRUE
selection_hsa_miR_4290	=	function(x) TRUE
selection_hsa_miR_4306	=	function(x) TRUE
selection_hsa_miR_4317	=	function(x) TRUE
selection_hsa_miR_4318	=	function(x) TRUE
selection_hsa_miR_4324	=	function(x) TRUE
selection_hsa_miR_4328	=	function(x) TRUE
selection_hsa_miR_4440	=	function(x) TRUE
selection_hsa_miR_4443	=	function(x) TRUE
selection_hsa_miR_4481	=	function(x) TRUE
selection_hsa_miR_449a	=	function(x) TRUE
selection_hsa_miR_450a_5p	=	function(x) TRUE
selection_hsa_miR_4516	=	function(x) TRUE
selection_hsa_miR_451a	=	function(x) TRUE
selection_hsa_miR_452_5p	=	function(x) TRUE
selection_hsa_miR_4521	=	function(x) TRUE
selection_hsa_miR_4532	=	function(x) TRUE
selection_hsa_miR_454_3p	=	function(x) TRUE
selection_hsa_miR_455_3p	=	function(x) TRUE
selection_hsa_miR_455_5p	=	function(x) TRUE
selection_hsa_miR_4634	=	function(x) TRUE
selection_hsa_miR_4655_3p	=	function(x) TRUE
selection_hsa_miR_4695_3p	=	function(x) TRUE
selection_hsa_miR_4716_5p	=	function(x) TRUE
selection_hsa_miR_4730	=	function(x) TRUE
selection_hsa_miR_4731_3p	=	function(x) TRUE
selection_hsa_miR_4763_5p	=	function(x) TRUE
selection_hsa_miR_4770	=	function(x) TRUE
selection_hsa_miR_4793_3p	=	function(x) TRUE
selection_hsa_miR_483_3p	=	function(x) TRUE
selection_hsa_miR_486_3p	=	function(x) TRUE
selection_hsa_miR_486_5p	=	function(x) TRUE
selection_hsa_miR_489_3p	=	function(x) TRUE
selection_hsa_miR_490_3p	=	function(x) TRUE
selection_hsa_miR_491_3p	=	function(x) TRUE
selection_hsa_miR_497_3p	=	function(x) TRUE
selection_hsa_miR_497_5p	=	function(x) TRUE
selection_hsa_miR_5001_5p	=	function(x) TRUE
selection_hsa_miR_500a_3p	=	function(x) TRUE
selection_hsa_miR_500b_5p	=	function(x) TRUE
selection_hsa_miR_501_5p	=	function(x) TRUE
selection_hsa_miR_502_3p	=	function(x) TRUE
selection_hsa_miR_504_3p	=	function(x) TRUE
selection_hsa_miR_505_5p	=	function(x) TRUE
selection_hsa_miR_511_3p	=	function(x) TRUE
selection_hsa_miR_516b_5p	=	function(x) TRUE
selection_hsa_miR_517_5p	=	function(x) TRUE
selection_hsa_miR_517a_3p	=	function(x) TRUE
selection_hsa_miR_517c_3p	=	function(x) TRUE
selection_hsa_miR_521	=	function(x) TRUE
selection_hsa_miR_522_3p	=	function(x) TRUE
selection_hsa_miR_532_3p	=	function(x) TRUE
selection_hsa_miR_532_5p	=	function(x) TRUE
selection_hsa_miR_548aa	=	function(x) TRUE
selection_hsa_miR_548aw	=	function(x) TRUE
selection_hsa_miR_548b_3p	=	function(x) TRUE
selection_hsa_miR_548c_3p	=	function(x) TRUE
selection_hsa_miR_548f_3p	=	function(x) TRUE
selection_hsa_miR_548q	=	function(x) TRUE
selection_hsa_miR_548x_3p	=	function(x) TRUE
selection_hsa_miR_551b_3p	=	function(x) TRUE
selection_hsa_miR_5701	=	function(x) TRUE
selection_hsa_miR_582_5p	=	function(x) TRUE
selection_hsa_miR_585_3p	=	function(x) TRUE
selection_hsa_miR_590_5p	=	function(x) TRUE
selection_hsa_miR_595	=	function(x) TRUE
selection_hsa_miR_598_3p	=	function(x) TRUE
selection_hsa_miR_6068	=	function(x) TRUE
selection_hsa_miR_6073	=	function(x) TRUE
selection_hsa_miR_6075	=	function(x) TRUE
selection_hsa_miR_610	=	function(x) TRUE
selection_hsa_miR_624_5p	=	function(x) TRUE
selection_hsa_miR_628_5p	=	function(x) TRUE
selection_hsa_miR_642b_5p	=	function(x) TRUE
selection_hsa_miR_6500_5p	=	function(x) TRUE
selection_hsa_miR_6516_3p	=	function(x) TRUE
selection_hsa_miR_652_3p	=	function(x) TRUE
selection_hsa_miR_653_3p	=	function(x) TRUE
selection_hsa_miR_660_5p	=	function(x) TRUE
selection_hsa_miR_664a_3p	=	function(x) TRUE
selection_hsa_miR_664b_3p	=	function(x) TRUE
selection_hsa_miR_6716_3p	=	function(x) TRUE
selection_hsa_miR_6722_5p	=	function(x) TRUE
selection_hsa_miR_6730_3p	=	function(x) TRUE
selection_hsa_miR_6743_3p	=	function(x) TRUE
selection_hsa_miR_6771_5p	=	function(x) TRUE
selection_hsa_miR_6779_3p	=	function(x) TRUE
selection_hsa_miR_6794_3p	=	function(x) TRUE
selection_hsa_miR_6804_5p	=	function(x) TRUE
selection_hsa_miR_6806_5p	=	function(x) TRUE
selection_hsa_miR_6817_5p	=	function(x) TRUE
selection_hsa_miR_6826_5p	=	function(x) TRUE
selection_hsa_miR_6865_3p	=	function(x) TRUE
selection_hsa_miR_6872_3p	=	function(x) TRUE
selection_hsa_miR_6886_3p	=	function(x) TRUE
selection_hsa_miR_6891_3p	=	function(x) TRUE
selection_hsa_miR_6895_5p	=	function(x) TRUE
selection_hsa_miR_7108_3p	=	function(x) TRUE
selection_hsa_miR_7111_3p	=	function(x) TRUE
selection_hsa_miR_7159_5p	=	function(x) TRUE
selection_hsa_miR_744_5p	=	function(x) TRUE
selection_hsa_miR_770_5p	=	function(x) TRUE
selection_hsa_miR_7704	=	function(x) TRUE
selection_hsa_miR_8063	=	function(x) TRUE
selection_hsa_miR_8077	=	function(x) TRUE
selection_hsa_miR_874_5p	=	function(x) TRUE
selection_hsa_miR_887_3p	=	function(x) TRUE
selection_hsa_miR_940	=	function(x) TRUE
selection_hsa_miR_95_3p	=	function(x) TRUE
selection_hsa_miR_98_5p	=	function(x) TRUE
selection_hsa_miR_99a_3p	=	function(x) TRUE
selection_hsa_miR_99a_5p	=	function(x) TRUE
selection_hsa_miR_99b_5p	=	function(x) TRUE

#GOdata

GOdata_hsa_let_7a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_let_7a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_let_7a_5p	, nodeSize=10)
GOdata_hsa_let_7b_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_let_7b_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_let_7b_5p	, nodeSize=10)
GOdata_hsa_let_7c_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_let_7c_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_let_7c_5p	, nodeSize=10)
GOdata_hsa_let_7d_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_let_7d_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_let_7d_5p	, nodeSize=10)
GOdata_hsa_let_7e_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_let_7e_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_let_7e_5p	, nodeSize=10)
GOdata_hsa_let_7f_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_let_7f_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_let_7f_5p	, nodeSize=10)
GOdata_hsa_let_7g_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_let_7g_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_let_7g_5p	, nodeSize=10)
GOdata_hsa_let_7i_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_let_7i_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_let_7i_5p	, nodeSize=10)
GOdata_hsa_miR_1_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_1_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_1_3p	, nodeSize=10)
GOdata_hsa_miR_100_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_100_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_100_5p	, nodeSize=10)
GOdata_hsa_miR_101_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_101_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_101_3p	, nodeSize=10)
GOdata_hsa_miR_101_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_101_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_101_5p	, nodeSize=10)
GOdata_hsa_miR_103a_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_103a_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_103a_3p	, nodeSize=10)
GOdata_hsa_miR_106b_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_106b_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_106b_3p	, nodeSize=10)
GOdata_hsa_miR_107	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_107	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_107	, nodeSize=10)
GOdata_hsa_miR_10a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_10a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_10a_5p	, nodeSize=10)
GOdata_hsa_miR_10b_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_10b_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_10b_5p	, nodeSize=10)
GOdata_hsa_miR_1227_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_1227_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_1227_5p	, nodeSize=10)
GOdata_hsa_miR_1238_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_1238_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_1238_5p	, nodeSize=10)
GOdata_hsa_miR_1247_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_1247_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_1247_5p	, nodeSize=10)
GOdata_hsa_miR_125a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_125a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_125a_5p	, nodeSize=10)
GOdata_hsa_miR_125b_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_125b_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_125b_5p	, nodeSize=10)
GOdata_hsa_miR_126_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_126_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_126_3p	, nodeSize=10)
GOdata_hsa_miR_126_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_126_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_126_5p	, nodeSize=10)
GOdata_hsa_miR_128_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_128_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_128_3p	, nodeSize=10)
GOdata_hsa_miR_129_2_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_129_2_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_129_2_3p	, nodeSize=10)
GOdata_hsa_miR_1296_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_1296_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_1296_5p	, nodeSize=10)
GOdata_hsa_miR_1306_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_1306_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_1306_3p	, nodeSize=10)
GOdata_hsa_miR_130a_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_130a_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_130a_3p	, nodeSize=10)
GOdata_hsa_miR_133a_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_133a_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_133a_3p	, nodeSize=10)
GOdata_hsa_miR_133a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_133a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_133a_5p	, nodeSize=10)
GOdata_hsa_miR_133b	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_133b	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_133b	, nodeSize=10)
GOdata_hsa_miR_135a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_135a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_135a_5p	, nodeSize=10)
GOdata_hsa_miR_138_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_138_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_138_5p	, nodeSize=10)
GOdata_hsa_miR_139_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_139_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_139_3p	, nodeSize=10)
GOdata_hsa_miR_139_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_139_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_139_5p	, nodeSize=10)
GOdata_hsa_miR_140_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_140_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_140_3p	, nodeSize=10)
GOdata_hsa_miR_140_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_140_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_140_5p	, nodeSize=10)
GOdata_hsa_miR_142_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_142_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_142_3p	, nodeSize=10)
GOdata_hsa_miR_142_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_142_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_142_5p	, nodeSize=10)
GOdata_hsa_miR_143_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_143_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_143_3p	, nodeSize=10)
GOdata_hsa_miR_143_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_143_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_143_5p	, nodeSize=10)
GOdata_hsa_miR_144_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_144_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_144_3p	, nodeSize=10)
GOdata_hsa_miR_144_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_144_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_144_5p	, nodeSize=10)
GOdata_hsa_miR_145_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_145_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_145_3p	, nodeSize=10)
GOdata_hsa_miR_145_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_145_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_145_5p	, nodeSize=10)
GOdata_hsa_miR_146a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_146a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_146a_5p	, nodeSize=10)
GOdata_hsa_miR_146b_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_146b_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_146b_5p	, nodeSize=10)
GOdata_hsa_miR_147b	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_147b	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_147b	, nodeSize=10)
GOdata_hsa_miR_150_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_150_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_150_5p	, nodeSize=10)
GOdata_hsa_miR_151a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_151a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_151a_5p	, nodeSize=10)
GOdata_hsa_miR_151b	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_151b	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_151b	, nodeSize=10)
GOdata_hsa_miR_152_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_152_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_152_3p	, nodeSize=10)
GOdata_hsa_miR_1537_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_1537_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_1537_3p	, nodeSize=10)
GOdata_hsa_miR_15a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_15a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_15a_5p	, nodeSize=10)
GOdata_hsa_miR_15b_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_15b_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_15b_3p	, nodeSize=10)
GOdata_hsa_miR_15b_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_15b_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_15b_5p	, nodeSize=10)
GOdata_hsa_miR_16_2_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_16_2_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_16_2_3p	, nodeSize=10)
GOdata_hsa_miR_16_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_16_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_16_5p	, nodeSize=10)
GOdata_hsa_miR_181a_2_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_181a_2_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_181a_2_3p	, nodeSize=10)
GOdata_hsa_miR_181a_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_181a_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_181a_3p	, nodeSize=10)
GOdata_hsa_miR_182_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_182_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_182_3p	, nodeSize=10)
GOdata_hsa_miR_183_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_183_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_183_3p	, nodeSize=10)
GOdata_hsa_miR_184	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_184	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_184	, nodeSize=10)
GOdata_hsa_miR_185_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_185_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_185_5p	, nodeSize=10)
GOdata_hsa_miR_186_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_186_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_186_5p	, nodeSize=10)
GOdata_hsa_miR_187_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_187_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_187_5p	, nodeSize=10)
GOdata_hsa_miR_18a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_18a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_18a_5p	, nodeSize=10)
GOdata_hsa_miR_18b_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_18b_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_18b_5p	, nodeSize=10)
GOdata_hsa_miR_190a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_190a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_190a_5p	, nodeSize=10)
GOdata_hsa_miR_191_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_191_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_191_5p	, nodeSize=10)
GOdata_hsa_miR_1913	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_1913	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_1913	, nodeSize=10)
GOdata_hsa_miR_193a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_193a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_193a_5p	, nodeSize=10)
GOdata_hsa_miR_195_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_195_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_195_5p	, nodeSize=10)
GOdata_hsa_miR_199a_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_199a_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_199a_3p	, nodeSize=10)
GOdata_hsa_miR_199a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_199a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_199a_5p	, nodeSize=10)
GOdata_hsa_miR_199b_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_199b_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_199b_5p	, nodeSize=10)
GOdata_hsa_miR_203a_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_203a_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_203a_3p	, nodeSize=10)
GOdata_hsa_miR_204_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_204_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_204_5p	, nodeSize=10)
GOdata_hsa_miR_205_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_205_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_205_3p	, nodeSize=10)
GOdata_hsa_miR_20a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_20a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_20a_5p	, nodeSize=10)
GOdata_hsa_miR_20b_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_20b_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_20b_5p	, nodeSize=10)
GOdata_hsa_miR_21_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_21_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_21_3p	, nodeSize=10)
GOdata_hsa_miR_214_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_214_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_214_3p	, nodeSize=10)
GOdata_hsa_miR_214_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_214_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_214_5p	, nodeSize=10)
GOdata_hsa_miR_218_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_218_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_218_5p	, nodeSize=10)
GOdata_hsa_miR_22_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_22_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_22_5p	, nodeSize=10)
GOdata_hsa_miR_221_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_221_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_221_3p	, nodeSize=10)
GOdata_hsa_miR_221_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_221_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_221_5p	, nodeSize=10)
GOdata_hsa_miR_222_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_222_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_222_3p	, nodeSize=10)
GOdata_hsa_miR_223_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_223_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_223_3p	, nodeSize=10)
GOdata_hsa_miR_223_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_223_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_223_5p	, nodeSize=10)
GOdata_hsa_miR_224_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_224_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_224_3p	, nodeSize=10)
GOdata_hsa_miR_2277_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_2277_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_2277_3p	, nodeSize=10)
GOdata_hsa_miR_23a_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_23a_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_23a_3p	, nodeSize=10)
GOdata_hsa_miR_23a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_23a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_23a_5p	, nodeSize=10)
GOdata_hsa_miR_23b_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_23b_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_23b_3p	, nodeSize=10)
GOdata_hsa_miR_24_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_24_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_24_3p	, nodeSize=10)
GOdata_hsa_miR_26a_1_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_26a_1_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_26a_1_3p	, nodeSize=10)
GOdata_hsa_miR_26a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_26a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_26a_5p	, nodeSize=10)
GOdata_hsa_miR_26b_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_26b_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_26b_3p	, nodeSize=10)
GOdata_hsa_miR_26b_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_26b_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_26b_5p	, nodeSize=10)
GOdata_hsa_miR_27a_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_27a_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_27a_3p	, nodeSize=10)
GOdata_hsa_miR_27a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_27a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_27a_5p	, nodeSize=10)
GOdata_hsa_miR_27b_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_27b_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_27b_3p	, nodeSize=10)
GOdata_hsa_miR_28_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_28_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_28_5p	, nodeSize=10)
GOdata_hsa_miR_2861	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_2861	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_2861	, nodeSize=10)
GOdata_hsa_miR_29a_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_29a_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_29a_3p	, nodeSize=10)
GOdata_hsa_miR_29a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_29a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_29a_5p	, nodeSize=10)
GOdata_hsa_miR_29b_1_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_29b_1_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_29b_1_5p	, nodeSize=10)
GOdata_hsa_miR_29b_2_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_29b_2_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_29b_2_5p	, nodeSize=10)
GOdata_hsa_miR_29b_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_29b_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_29b_3p	, nodeSize=10)
GOdata_hsa_miR_29c_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_29c_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_29c_3p	, nodeSize=10)
GOdata_hsa_miR_29c_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_29c_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_29c_5p	, nodeSize=10)
GOdata_hsa_miR_3065_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_3065_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_3065_3p	, nodeSize=10)
GOdata_hsa_miR_3065_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_3065_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_3065_5p	, nodeSize=10)
GOdata_hsa_miR_30a_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_30a_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_30a_3p	, nodeSize=10)
GOdata_hsa_miR_30a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_30a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_30a_5p	, nodeSize=10)
GOdata_hsa_miR_30b_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_30b_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_30b_5p	, nodeSize=10)
GOdata_hsa_miR_30c_2_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_30c_2_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_30c_2_3p	, nodeSize=10)
GOdata_hsa_miR_30c_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_30c_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_30c_5p	, nodeSize=10)
GOdata_hsa_miR_30d_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_30d_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_30d_5p	, nodeSize=10)
GOdata_hsa_miR_30e_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_30e_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_30e_3p	, nodeSize=10)
GOdata_hsa_miR_3149	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_3149	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_3149	, nodeSize=10)
GOdata_hsa_miR_3188	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_3188	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_3188	, nodeSize=10)
GOdata_hsa_miR_32_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_32_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_32_5p	, nodeSize=10)
GOdata_hsa_miR_324_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_324_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_324_3p	, nodeSize=10)
GOdata_hsa_miR_326	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_326	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_326	, nodeSize=10)
GOdata_hsa_miR_328_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_328_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_328_3p	, nodeSize=10)
GOdata_hsa_miR_331_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_331_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_331_3p	, nodeSize=10)
GOdata_hsa_miR_335_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_335_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_335_5p	, nodeSize=10)
GOdata_hsa_miR_338_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_338_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_338_3p	, nodeSize=10)
GOdata_hsa_miR_339_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_339_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_339_5p	, nodeSize=10)
GOdata_hsa_miR_340_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_340_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_340_5p	, nodeSize=10)
GOdata_hsa_miR_342_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_342_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_342_3p	, nodeSize=10)
GOdata_hsa_miR_342_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_342_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_342_5p	, nodeSize=10)
GOdata_hsa_miR_34a_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_34a_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_34a_3p	, nodeSize=10)
GOdata_hsa_miR_34a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_34a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_34a_5p	, nodeSize=10)
GOdata_hsa_miR_34b_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_34b_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_34b_5p	, nodeSize=10)
GOdata_hsa_miR_34c_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_34c_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_34c_5p	, nodeSize=10)
GOdata_hsa_miR_3607_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_3607_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_3607_3p	, nodeSize=10)
GOdata_hsa_miR_361_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_361_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_361_3p	, nodeSize=10)
GOdata_hsa_miR_362_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_362_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_362_3p	, nodeSize=10)
GOdata_hsa_miR_362_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_362_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_362_5p	, nodeSize=10)
GOdata_hsa_miR_3620_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_3620_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_3620_5p	, nodeSize=10)
GOdata_hsa_miR_363_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_363_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_363_3p	, nodeSize=10)
GOdata_hsa_miR_3659	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_3659	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_3659	, nodeSize=10)
GOdata_hsa_miR_365a_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_365a_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_365a_3p	, nodeSize=10)
GOdata_hsa_miR_3660	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_3660	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_3660	, nodeSize=10)
GOdata_hsa_miR_3663_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_3663_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_3663_5p	, nodeSize=10)
GOdata_hsa_miR_371b_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_371b_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_371b_5p	, nodeSize=10)
GOdata_hsa_miR_374a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_374a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_374a_5p	, nodeSize=10)
GOdata_hsa_miR_374b_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_374b_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_374b_5p	, nodeSize=10)
GOdata_hsa_miR_374c_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_374c_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_374c_5p	, nodeSize=10)
GOdata_hsa_miR_424_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_424_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_424_5p	, nodeSize=10)
GOdata_hsa_miR_4252	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4252	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4252	, nodeSize=10)
GOdata_hsa_miR_4254	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4254	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4254	, nodeSize=10)
GOdata_hsa_miR_4284	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4284	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4284	, nodeSize=10)
GOdata_hsa_miR_4290	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4290	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4290	, nodeSize=10)
GOdata_hsa_miR_4306	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4306	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4306	, nodeSize=10)
GOdata_hsa_miR_4317	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4317	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4317	, nodeSize=10)
GOdata_hsa_miR_4318	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4318	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4318	, nodeSize=10)
GOdata_hsa_miR_4324	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4324	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4324	, nodeSize=10)
GOdata_hsa_miR_4328	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4328	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4328	, nodeSize=10)
GOdata_hsa_miR_4440	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4440	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4440	, nodeSize=10)
GOdata_hsa_miR_4443	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4443	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4443	, nodeSize=10)
GOdata_hsa_miR_4481	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4481	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4481	, nodeSize=10)
GOdata_hsa_miR_449a	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_449a	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_449a	, nodeSize=10)
GOdata_hsa_miR_450a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_450a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_450a_5p	, nodeSize=10)
GOdata_hsa_miR_4516	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4516	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4516	, nodeSize=10)
GOdata_hsa_miR_451a	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_451a	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_451a	, nodeSize=10)
GOdata_hsa_miR_452_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_452_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_452_5p	, nodeSize=10)
GOdata_hsa_miR_4521	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4521	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4521	, nodeSize=10)
GOdata_hsa_miR_4532	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4532	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4532	, nodeSize=10)
GOdata_hsa_miR_454_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_454_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_454_3p	, nodeSize=10)
GOdata_hsa_miR_455_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_455_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_455_3p	, nodeSize=10)
GOdata_hsa_miR_455_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_455_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_455_5p	, nodeSize=10)
GOdata_hsa_miR_4634	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4634	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4634	, nodeSize=10)
GOdata_hsa_miR_4655_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4655_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4655_3p	, nodeSize=10)
GOdata_hsa_miR_4695_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4695_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4695_3p	, nodeSize=10)
GOdata_hsa_miR_4716_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4716_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4716_5p	, nodeSize=10)
GOdata_hsa_miR_4730	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4730	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4730	, nodeSize=10)
GOdata_hsa_miR_4731_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4731_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4731_3p	, nodeSize=10)
GOdata_hsa_miR_4763_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4763_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4763_5p	, nodeSize=10)
GOdata_hsa_miR_4770	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4770	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4770	, nodeSize=10)
GOdata_hsa_miR_4793_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_4793_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_4793_3p	, nodeSize=10)
GOdata_hsa_miR_483_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_483_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_483_3p	, nodeSize=10)
GOdata_hsa_miR_486_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_486_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_486_3p	, nodeSize=10)
GOdata_hsa_miR_486_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_486_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_486_5p	, nodeSize=10)
GOdata_hsa_miR_489_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_489_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_489_3p	, nodeSize=10)
GOdata_hsa_miR_490_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_490_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_490_3p	, nodeSize=10)
GOdata_hsa_miR_491_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_491_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_491_3p	, nodeSize=10)
GOdata_hsa_miR_497_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_497_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_497_3p	, nodeSize=10)
GOdata_hsa_miR_497_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_497_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_497_5p	, nodeSize=10)
GOdata_hsa_miR_5001_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_5001_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_5001_5p	, nodeSize=10)
GOdata_hsa_miR_500a_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_500a_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_500a_3p	, nodeSize=10)
GOdata_hsa_miR_500b_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_500b_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_500b_5p	, nodeSize=10)
GOdata_hsa_miR_501_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_501_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_501_5p	, nodeSize=10)
GOdata_hsa_miR_502_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_502_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_502_3p	, nodeSize=10)
GOdata_hsa_miR_504_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_504_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_504_3p	, nodeSize=10)
GOdata_hsa_miR_505_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_505_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_505_5p	, nodeSize=10)
GOdata_hsa_miR_511_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_511_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_511_3p	, nodeSize=10)
GOdata_hsa_miR_516b_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_516b_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_516b_5p	, nodeSize=10)
GOdata_hsa_miR_517_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_517_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_517_5p	, nodeSize=10)
GOdata_hsa_miR_517a_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_517a_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_517a_3p	, nodeSize=10)
GOdata_hsa_miR_517c_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_517c_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_517c_3p	, nodeSize=10)
GOdata_hsa_miR_521	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_521	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_521	, nodeSize=10)
GOdata_hsa_miR_522_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_522_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_522_3p	, nodeSize=10)
GOdata_hsa_miR_532_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_532_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_532_3p	, nodeSize=10)
GOdata_hsa_miR_532_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_532_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_532_5p	, nodeSize=10)
GOdata_hsa_miR_548aa	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_548aa	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_548aa	, nodeSize=10)
GOdata_hsa_miR_548aw	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_548aw	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_548aw	, nodeSize=10)
GOdata_hsa_miR_548b_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_548b_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_548b_3p	, nodeSize=10)
GOdata_hsa_miR_548c_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_548c_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_548c_3p	, nodeSize=10)
GOdata_hsa_miR_548f_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_548f_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_548f_3p	, nodeSize=10)
GOdata_hsa_miR_548q	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_548q	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_548q	, nodeSize=10)
GOdata_hsa_miR_548x_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_548x_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_548x_3p	, nodeSize=10)
GOdata_hsa_miR_551b_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_551b_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_551b_3p	, nodeSize=10)
GOdata_hsa_miR_5701	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_5701	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_5701	, nodeSize=10)
GOdata_hsa_miR_582_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_582_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_582_5p	, nodeSize=10)
GOdata_hsa_miR_585_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_585_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_585_3p	, nodeSize=10)
GOdata_hsa_miR_590_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_590_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_590_5p	, nodeSize=10)
GOdata_hsa_miR_595	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_595	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_595	, nodeSize=10)
GOdata_hsa_miR_598_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_598_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_598_3p	, nodeSize=10)
GOdata_hsa_miR_6068	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_6068	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_6068	, nodeSize=10)
GOdata_hsa_miR_6073	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_6073	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_6073	, nodeSize=10)
GOdata_hsa_miR_6075	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_6075	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_6075	, nodeSize=10)
GOdata_hsa_miR_610	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_610	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_610	, nodeSize=10)
GOdata_hsa_miR_624_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_624_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_624_5p	, nodeSize=10)
GOdata_hsa_miR_628_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_628_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_628_5p	, nodeSize=10)
GOdata_hsa_miR_642b_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_642b_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_642b_5p	, nodeSize=10)
GOdata_hsa_miR_6500_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_6500_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_6500_5p	, nodeSize=10)
GOdata_hsa_miR_6516_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_6516_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_6516_3p	, nodeSize=10)
GOdata_hsa_miR_652_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_652_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_652_3p	, nodeSize=10)
GOdata_hsa_miR_653_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_653_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_653_3p	, nodeSize=10)
GOdata_hsa_miR_660_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_660_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_660_5p	, nodeSize=10)
GOdata_hsa_miR_664a_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_664a_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_664a_3p	, nodeSize=10)
GOdata_hsa_miR_664b_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_664b_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_664b_3p	, nodeSize=10)
GOdata_hsa_miR_6716_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_6716_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_6716_3p	, nodeSize=10)
GOdata_hsa_miR_6722_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_6722_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_6722_5p	, nodeSize=10)
GOdata_hsa_miR_6730_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_6730_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_6730_3p	, nodeSize=10)
GOdata_hsa_miR_6743_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_6743_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_6743_3p	, nodeSize=10)
GOdata_hsa_miR_6771_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_6771_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_6771_5p	, nodeSize=10)
GOdata_hsa_miR_6779_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_6779_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_6779_3p	, nodeSize=10)
GOdata_hsa_miR_6794_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_6794_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_6794_3p	, nodeSize=10)
GOdata_hsa_miR_6804_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_6804_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_6804_5p	, nodeSize=10)
GOdata_hsa_miR_6806_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_6806_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_6806_5p	, nodeSize=10)
GOdata_hsa_miR_6817_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_6817_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_6817_5p	, nodeSize=10)
GOdata_hsa_miR_6826_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_6826_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_6826_5p	, nodeSize=10)
GOdata_hsa_miR_6865_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_6865_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_6865_3p	, nodeSize=10)
GOdata_hsa_miR_6872_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_6872_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_6872_3p	, nodeSize=10)
GOdata_hsa_miR_6886_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_6886_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_6886_3p	, nodeSize=10)
GOdata_hsa_miR_6891_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_6891_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_6891_3p	, nodeSize=10)
GOdata_hsa_miR_6895_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_6895_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_6895_5p	, nodeSize=10)
GOdata_hsa_miR_7108_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_7108_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_7108_3p	, nodeSize=10)
GOdata_hsa_miR_7111_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_7111_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_7111_3p	, nodeSize=10)
GOdata_hsa_miR_7159_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_7159_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_7159_5p	, nodeSize=10)
GOdata_hsa_miR_744_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_744_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_744_5p	, nodeSize=10)
GOdata_hsa_miR_770_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_770_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_770_5p	, nodeSize=10)
GOdata_hsa_miR_7704	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_7704	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_7704	, nodeSize=10)
GOdata_hsa_miR_8063	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_8063	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_8063	, nodeSize=10)
GOdata_hsa_miR_8077	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_8077	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_8077	, nodeSize=10)
GOdata_hsa_miR_874_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_874_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_874_5p	, nodeSize=10)
GOdata_hsa_miR_887_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_887_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_887_3p	, nodeSize=10)
GOdata_hsa_miR_940	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_940	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_940	, nodeSize=10)
GOdata_hsa_miR_95_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_95_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_95_3p	, nodeSize=10)
GOdata_hsa_miR_98_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_98_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_98_5p	, nodeSize=10)
GOdata_hsa_miR_99a_3p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_99a_3p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_99a_3p	, nodeSize=10)
GOdata_hsa_miR_99a_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_99a_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_99a_5p	, nodeSize=10)
GOdata_hsa_miR_99b_5p	=	new('topGOdata', ontology = 'BP', allGenes = rankedGenes_hsa_miR_99b_5p	, annot = annFUN.GO2genes, GO2genes = allGO2genes,  geneSel = selection_hsa_miR_99b_5p	, nodeSize=10)

#RuntestKs

results.ks_hsa_let_7a_5p	=	runTest(	GOdata_hsa_let_7a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_let_7b_5p	=	runTest(	GOdata_hsa_let_7b_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_let_7c_5p	=	runTest(	GOdata_hsa_let_7c_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_let_7d_5p	=	runTest(	GOdata_hsa_let_7d_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_let_7e_5p	=	runTest(	GOdata_hsa_let_7e_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_let_7f_5p	=	runTest(	GOdata_hsa_let_7f_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_let_7g_5p	=	runTest(	GOdata_hsa_let_7g_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_let_7i_5p	=	runTest(	GOdata_hsa_let_7i_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_1_3p	=	runTest(	GOdata_hsa_miR_1_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_100_5p	=	runTest(	GOdata_hsa_miR_100_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_101_3p	=	runTest(	GOdata_hsa_miR_101_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_101_5p	=	runTest(	GOdata_hsa_miR_101_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_103a_3p	=	runTest(	GOdata_hsa_miR_103a_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_106b_3p	=	runTest(	GOdata_hsa_miR_106b_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_107	=	runTest(	GOdata_hsa_miR_107	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_10a_5p	=	runTest(	GOdata_hsa_miR_10a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_10b_5p	=	runTest(	GOdata_hsa_miR_10b_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_1227_5p	=	runTest(	GOdata_hsa_miR_1227_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_1238_5p	=	runTest(	GOdata_hsa_miR_1238_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_1247_5p	=	runTest(	GOdata_hsa_miR_1247_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_125a_5p	=	runTest(	GOdata_hsa_miR_125a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_125b_5p	=	runTest(	GOdata_hsa_miR_125b_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_126_3p	=	runTest(	GOdata_hsa_miR_126_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_126_5p	=	runTest(	GOdata_hsa_miR_126_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_128_3p	=	runTest(	GOdata_hsa_miR_128_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_129_2_3p	=	runTest(	GOdata_hsa_miR_129_2_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_1296_5p	=	runTest(	GOdata_hsa_miR_1296_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_1306_3p	=	runTest(	GOdata_hsa_miR_1306_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_130a_3p	=	runTest(	GOdata_hsa_miR_130a_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_133a_3p	=	runTest(	GOdata_hsa_miR_133a_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_133a_5p	=	runTest(	GOdata_hsa_miR_133a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_133b	=	runTest(	GOdata_hsa_miR_133b	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_135a_5p	=	runTest(	GOdata_hsa_miR_135a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_138_5p	=	runTest(	GOdata_hsa_miR_138_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_139_3p	=	runTest(	GOdata_hsa_miR_139_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_139_5p	=	runTest(	GOdata_hsa_miR_139_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_140_3p	=	runTest(	GOdata_hsa_miR_140_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_140_5p	=	runTest(	GOdata_hsa_miR_140_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_142_3p	=	runTest(	GOdata_hsa_miR_142_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_142_5p	=	runTest(	GOdata_hsa_miR_142_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_143_3p	=	runTest(	GOdata_hsa_miR_143_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_143_5p	=	runTest(	GOdata_hsa_miR_143_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_144_3p	=	runTest(	GOdata_hsa_miR_144_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_144_5p	=	runTest(	GOdata_hsa_miR_144_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_145_3p	=	runTest(	GOdata_hsa_miR_145_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_145_5p	=	runTest(	GOdata_hsa_miR_145_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_146a_5p	=	runTest(	GOdata_hsa_miR_146a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_146b_5p	=	runTest(	GOdata_hsa_miR_146b_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_147b	=	runTest(	GOdata_hsa_miR_147b	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_150_5p	=	runTest(	GOdata_hsa_miR_150_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_151a_5p	=	runTest(	GOdata_hsa_miR_151a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_151b	=	runTest(	GOdata_hsa_miR_151b	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_152_3p	=	runTest(	GOdata_hsa_miR_152_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_1537_3p	=	runTest(	GOdata_hsa_miR_1537_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_15a_5p	=	runTest(	GOdata_hsa_miR_15a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_15b_3p	=	runTest(	GOdata_hsa_miR_15b_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_15b_5p	=	runTest(	GOdata_hsa_miR_15b_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_16_2_3p	=	runTest(	GOdata_hsa_miR_16_2_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_16_5p	=	runTest(	GOdata_hsa_miR_16_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_181a_2_3p	=	runTest(	GOdata_hsa_miR_181a_2_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_181a_3p	=	runTest(	GOdata_hsa_miR_181a_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_182_3p	=	runTest(	GOdata_hsa_miR_182_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_183_3p	=	runTest(	GOdata_hsa_miR_183_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_184	=	runTest(	GOdata_hsa_miR_184	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_185_5p	=	runTest(	GOdata_hsa_miR_185_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_186_5p	=	runTest(	GOdata_hsa_miR_186_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_187_5p	=	runTest(	GOdata_hsa_miR_187_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_18a_5p	=	runTest(	GOdata_hsa_miR_18a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_18b_5p	=	runTest(	GOdata_hsa_miR_18b_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_190a_5p	=	runTest(	GOdata_hsa_miR_190a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_191_5p	=	runTest(	GOdata_hsa_miR_191_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_1913	=	runTest(	GOdata_hsa_miR_1913	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_193a_5p	=	runTest(	GOdata_hsa_miR_193a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_195_5p	=	runTest(	GOdata_hsa_miR_195_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_199a_3p	=	runTest(	GOdata_hsa_miR_199a_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_199a_5p	=	runTest(	GOdata_hsa_miR_199a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_199b_5p	=	runTest(	GOdata_hsa_miR_199b_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_203a_3p	=	runTest(	GOdata_hsa_miR_203a_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_204_5p	=	runTest(	GOdata_hsa_miR_204_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_205_3p	=	runTest(	GOdata_hsa_miR_205_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_20a_5p	=	runTest(	GOdata_hsa_miR_20a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_20b_5p	=	runTest(	GOdata_hsa_miR_20b_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_21_3p	=	runTest(	GOdata_hsa_miR_21_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_214_3p	=	runTest(	GOdata_hsa_miR_214_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_214_5p	=	runTest(	GOdata_hsa_miR_214_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_218_5p	=	runTest(	GOdata_hsa_miR_218_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_22_5p	=	runTest(	GOdata_hsa_miR_22_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_221_3p	=	runTest(	GOdata_hsa_miR_221_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_221_5p	=	runTest(	GOdata_hsa_miR_221_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_222_3p	=	runTest(	GOdata_hsa_miR_222_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_223_3p	=	runTest(	GOdata_hsa_miR_223_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_223_5p	=	runTest(	GOdata_hsa_miR_223_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_224_3p	=	runTest(	GOdata_hsa_miR_224_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_2277_3p	=	runTest(	GOdata_hsa_miR_2277_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_23a_3p	=	runTest(	GOdata_hsa_miR_23a_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_23a_5p	=	runTest(	GOdata_hsa_miR_23a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_23b_3p	=	runTest(	GOdata_hsa_miR_23b_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_24_3p	=	runTest(	GOdata_hsa_miR_24_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_26a_1_3p	=	runTest(	GOdata_hsa_miR_26a_1_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_26a_5p	=	runTest(	GOdata_hsa_miR_26a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_26b_3p	=	runTest(	GOdata_hsa_miR_26b_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_26b_5p	=	runTest(	GOdata_hsa_miR_26b_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_27a_3p	=	runTest(	GOdata_hsa_miR_27a_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_27a_5p	=	runTest(	GOdata_hsa_miR_27a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_27b_3p	=	runTest(	GOdata_hsa_miR_27b_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_28_5p	=	runTest(	GOdata_hsa_miR_28_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_2861	=	runTest(	GOdata_hsa_miR_2861	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_29a_3p	=	runTest(	GOdata_hsa_miR_29a_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_29a_5p	=	runTest(	GOdata_hsa_miR_29a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_29b_1_5p	=	runTest(	GOdata_hsa_miR_29b_1_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_29b_2_5p	=	runTest(	GOdata_hsa_miR_29b_2_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_29b_3p	=	runTest(	GOdata_hsa_miR_29b_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_29c_3p	=	runTest(	GOdata_hsa_miR_29c_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_29c_5p	=	runTest(	GOdata_hsa_miR_29c_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_3065_3p	=	runTest(	GOdata_hsa_miR_3065_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_3065_5p	=	runTest(	GOdata_hsa_miR_3065_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_30a_3p	=	runTest(	GOdata_hsa_miR_30a_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_30a_5p	=	runTest(	GOdata_hsa_miR_30a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_30b_5p	=	runTest(	GOdata_hsa_miR_30b_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_30c_2_3p	=	runTest(	GOdata_hsa_miR_30c_2_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_30c_5p	=	runTest(	GOdata_hsa_miR_30c_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_30d_5p	=	runTest(	GOdata_hsa_miR_30d_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_30e_3p	=	runTest(	GOdata_hsa_miR_30e_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_3149	=	runTest(	GOdata_hsa_miR_3149	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_3188	=	runTest(	GOdata_hsa_miR_3188	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_32_5p	=	runTest(	GOdata_hsa_miR_32_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_324_3p	=	runTest(	GOdata_hsa_miR_324_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_326	=	runTest(	GOdata_hsa_miR_326	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_328_3p	=	runTest(	GOdata_hsa_miR_328_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_331_3p	=	runTest(	GOdata_hsa_miR_331_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_335_5p	=	runTest(	GOdata_hsa_miR_335_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_338_3p	=	runTest(	GOdata_hsa_miR_338_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_339_5p	=	runTest(	GOdata_hsa_miR_339_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_340_5p	=	runTest(	GOdata_hsa_miR_340_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_342_3p	=	runTest(	GOdata_hsa_miR_342_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_342_5p	=	runTest(	GOdata_hsa_miR_342_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_34a_3p	=	runTest(	GOdata_hsa_miR_34a_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_34a_5p	=	runTest(	GOdata_hsa_miR_34a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_34b_5p	=	runTest(	GOdata_hsa_miR_34b_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_34c_5p	=	runTest(	GOdata_hsa_miR_34c_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_3607_3p	=	runTest(	GOdata_hsa_miR_3607_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_361_3p	=	runTest(	GOdata_hsa_miR_361_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_362_3p	=	runTest(	GOdata_hsa_miR_362_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_362_5p	=	runTest(	GOdata_hsa_miR_362_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_3620_5p	=	runTest(	GOdata_hsa_miR_3620_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_363_3p	=	runTest(	GOdata_hsa_miR_363_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_3659	=	runTest(	GOdata_hsa_miR_3659	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_365a_3p	=	runTest(	GOdata_hsa_miR_365a_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_3660	=	runTest(	GOdata_hsa_miR_3660	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_3663_5p	=	runTest(	GOdata_hsa_miR_3663_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_371b_5p	=	runTest(	GOdata_hsa_miR_371b_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_374a_5p	=	runTest(	GOdata_hsa_miR_374a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_374b_5p	=	runTest(	GOdata_hsa_miR_374b_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_374c_5p	=	runTest(	GOdata_hsa_miR_374c_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_424_5p	=	runTest(	GOdata_hsa_miR_424_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4252	=	runTest(	GOdata_hsa_miR_4252	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4254	=	runTest(	GOdata_hsa_miR_4254	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4284	=	runTest(	GOdata_hsa_miR_4284	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4290	=	runTest(	GOdata_hsa_miR_4290	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4306	=	runTest(	GOdata_hsa_miR_4306	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4317	=	runTest(	GOdata_hsa_miR_4317	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4318	=	runTest(	GOdata_hsa_miR_4318	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4324	=	runTest(	GOdata_hsa_miR_4324	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4328	=	runTest(	GOdata_hsa_miR_4328	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4440	=	runTest(	GOdata_hsa_miR_4440	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4443	=	runTest(	GOdata_hsa_miR_4443	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4481	=	runTest(	GOdata_hsa_miR_4481	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_449a	=	runTest(	GOdata_hsa_miR_449a	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_450a_5p	=	runTest(	GOdata_hsa_miR_450a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4516	=	runTest(	GOdata_hsa_miR_4516	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_451a	=	runTest(	GOdata_hsa_miR_451a	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_452_5p	=	runTest(	GOdata_hsa_miR_452_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4521	=	runTest(	GOdata_hsa_miR_4521	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4532	=	runTest(	GOdata_hsa_miR_4532	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_454_3p	=	runTest(	GOdata_hsa_miR_454_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_455_3p	=	runTest(	GOdata_hsa_miR_455_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_455_5p	=	runTest(	GOdata_hsa_miR_455_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4634	=	runTest(	GOdata_hsa_miR_4634	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4655_3p	=	runTest(	GOdata_hsa_miR_4655_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4695_3p	=	runTest(	GOdata_hsa_miR_4695_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4716_5p	=	runTest(	GOdata_hsa_miR_4716_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4730	=	runTest(	GOdata_hsa_miR_4730	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4731_3p	=	runTest(	GOdata_hsa_miR_4731_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4763_5p	=	runTest(	GOdata_hsa_miR_4763_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4770	=	runTest(	GOdata_hsa_miR_4770	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_4793_3p	=	runTest(	GOdata_hsa_miR_4793_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_483_3p	=	runTest(	GOdata_hsa_miR_483_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_486_3p	=	runTest(	GOdata_hsa_miR_486_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_486_5p	=	runTest(	GOdata_hsa_miR_486_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_489_3p	=	runTest(	GOdata_hsa_miR_489_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_490_3p	=	runTest(	GOdata_hsa_miR_490_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_491_3p	=	runTest(	GOdata_hsa_miR_491_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_497_3p	=	runTest(	GOdata_hsa_miR_497_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_497_5p	=	runTest(	GOdata_hsa_miR_497_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_5001_5p	=	runTest(	GOdata_hsa_miR_5001_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_500a_3p	=	runTest(	GOdata_hsa_miR_500a_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_500b_5p	=	runTest(	GOdata_hsa_miR_500b_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_501_5p	=	runTest(	GOdata_hsa_miR_501_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_502_3p	=	runTest(	GOdata_hsa_miR_502_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_504_3p	=	runTest(	GOdata_hsa_miR_504_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_505_5p	=	runTest(	GOdata_hsa_miR_505_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_511_3p	=	runTest(	GOdata_hsa_miR_511_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_516b_5p	=	runTest(	GOdata_hsa_miR_516b_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_517_5p	=	runTest(	GOdata_hsa_miR_517_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_517a_3p	=	runTest(	GOdata_hsa_miR_517a_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_517c_3p	=	runTest(	GOdata_hsa_miR_517c_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_521	=	runTest(	GOdata_hsa_miR_521	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_522_3p	=	runTest(	GOdata_hsa_miR_522_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_532_3p	=	runTest(	GOdata_hsa_miR_532_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_532_5p	=	runTest(	GOdata_hsa_miR_532_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_548aa	=	runTest(	GOdata_hsa_miR_548aa	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_548aw	=	runTest(	GOdata_hsa_miR_548aw	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_548b_3p	=	runTest(	GOdata_hsa_miR_548b_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_548c_3p	=	runTest(	GOdata_hsa_miR_548c_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_548f_3p	=	runTest(	GOdata_hsa_miR_548f_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_548q	=	runTest(	GOdata_hsa_miR_548q	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_548x_3p	=	runTest(	GOdata_hsa_miR_548x_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_551b_3p	=	runTest(	GOdata_hsa_miR_551b_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_5701	=	runTest(	GOdata_hsa_miR_5701	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_582_5p	=	runTest(	GOdata_hsa_miR_582_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_585_3p	=	runTest(	GOdata_hsa_miR_585_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_590_5p	=	runTest(	GOdata_hsa_miR_590_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_595	=	runTest(	GOdata_hsa_miR_595	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_598_3p	=	runTest(	GOdata_hsa_miR_598_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_6068	=	runTest(	GOdata_hsa_miR_6068	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_6073	=	runTest(	GOdata_hsa_miR_6073	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_6075	=	runTest(	GOdata_hsa_miR_6075	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_610	=	runTest(	GOdata_hsa_miR_610	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_624_5p	=	runTest(	GOdata_hsa_miR_624_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_628_5p	=	runTest(	GOdata_hsa_miR_628_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_642b_5p	=	runTest(	GOdata_hsa_miR_642b_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_6500_5p	=	runTest(	GOdata_hsa_miR_6500_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_6516_3p	=	runTest(	GOdata_hsa_miR_6516_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_652_3p	=	runTest(	GOdata_hsa_miR_652_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_653_3p	=	runTest(	GOdata_hsa_miR_653_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_660_5p	=	runTest(	GOdata_hsa_miR_660_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_664a_3p	=	runTest(	GOdata_hsa_miR_664a_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_664b_3p	=	runTest(	GOdata_hsa_miR_664b_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_6716_3p	=	runTest(	GOdata_hsa_miR_6716_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_6722_5p	=	runTest(	GOdata_hsa_miR_6722_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_6730_3p	=	runTest(	GOdata_hsa_miR_6730_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_6743_3p	=	runTest(	GOdata_hsa_miR_6743_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_6771_5p	=	runTest(	GOdata_hsa_miR_6771_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_6779_3p	=	runTest(	GOdata_hsa_miR_6779_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_6794_3p	=	runTest(	GOdata_hsa_miR_6794_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_6804_5p	=	runTest(	GOdata_hsa_miR_6804_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_6806_5p	=	runTest(	GOdata_hsa_miR_6806_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_6817_5p	=	runTest(	GOdata_hsa_miR_6817_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_6826_5p	=	runTest(	GOdata_hsa_miR_6826_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_6865_3p	=	runTest(	GOdata_hsa_miR_6865_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_6872_3p	=	runTest(	GOdata_hsa_miR_6872_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_6886_3p	=	runTest(	GOdata_hsa_miR_6886_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_6891_3p	=	runTest(	GOdata_hsa_miR_6891_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_6895_5p	=	runTest(	GOdata_hsa_miR_6895_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_7108_3p	=	runTest(	GOdata_hsa_miR_7108_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_7111_3p	=	runTest(	GOdata_hsa_miR_7111_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_7159_5p	=	runTest(	GOdata_hsa_miR_7159_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_744_5p	=	runTest(	GOdata_hsa_miR_744_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_770_5p	=	runTest(	GOdata_hsa_miR_770_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_7704	=	runTest(	GOdata_hsa_miR_7704	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_8063	=	runTest(	GOdata_hsa_miR_8063	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_8077	=	runTest(	GOdata_hsa_miR_8077	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_874_5p	=	runTest(	GOdata_hsa_miR_874_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_887_3p	=	runTest(	GOdata_hsa_miR_887_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_940	=	runTest(	GOdata_hsa_miR_940	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_95_3p	=	runTest(	GOdata_hsa_miR_95_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_98_5p	=	runTest(	GOdata_hsa_miR_98_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_99a_3p	=	runTest(	GOdata_hsa_miR_99a_3p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_99a_5p	=	runTest(	GOdata_hsa_miR_99a_5p	, algorithm = "classic", statistic = "ks")
results.ks_hsa_miR_99b_5p	=	runTest(	GOdata_hsa_miR_99b_5p	, algorithm = "classic", statistic = "ks")

##AllRes

allRes_hsa_let_7a_5p	=	GenTable(	GOdata_hsa_let_7a_5p	, KS = results.ks_hsa_let_7a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_let_7b_5p	=	GenTable(	GOdata_hsa_let_7b_5p	, KS = results.ks_hsa_let_7b_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_let_7c_5p	=	GenTable(	GOdata_hsa_let_7c_5p	, KS = results.ks_hsa_let_7c_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_let_7d_5p	=	GenTable(	GOdata_hsa_let_7d_5p	, KS = results.ks_hsa_let_7d_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_let_7e_5p	=	GenTable(	GOdata_hsa_let_7e_5p	, KS = results.ks_hsa_let_7e_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_let_7f_5p	=	GenTable(	GOdata_hsa_let_7f_5p	, KS = results.ks_hsa_let_7f_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_let_7g_5p	=	GenTable(	GOdata_hsa_let_7g_5p	, KS = results.ks_hsa_let_7g_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_let_7i_5p	=	GenTable(	GOdata_hsa_let_7i_5p	, KS = results.ks_hsa_let_7i_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_1_3p	=	GenTable(	GOdata_hsa_miR_1_3p	, KS = results.ks_hsa_miR_1_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_100_5p	=	GenTable(	GOdata_hsa_miR_100_5p	, KS = results.ks_hsa_miR_100_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_101_3p	=	GenTable(	GOdata_hsa_miR_101_3p	, KS = results.ks_hsa_miR_101_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_101_5p	=	GenTable(	GOdata_hsa_miR_101_5p	, KS = results.ks_hsa_miR_101_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_103a_3p	=	GenTable(	GOdata_hsa_miR_103a_3p	, KS = results.ks_hsa_miR_103a_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_106b_3p	=	GenTable(	GOdata_hsa_miR_106b_3p	, KS = results.ks_hsa_miR_106b_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_107	=	GenTable(	GOdata_hsa_miR_107	, KS = results.ks_hsa_miR_107	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_10a_5p	=	GenTable(	GOdata_hsa_miR_10a_5p	, KS = results.ks_hsa_miR_10a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_10b_5p	=	GenTable(	GOdata_hsa_miR_10b_5p	, KS = results.ks_hsa_miR_10b_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_1227_5p	=	GenTable(	GOdata_hsa_miR_1227_5p	, KS = results.ks_hsa_miR_1227_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_1238_5p	=	GenTable(	GOdata_hsa_miR_1238_5p	, KS = results.ks_hsa_miR_1238_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_1247_5p	=	GenTable(	GOdata_hsa_miR_1247_5p	, KS = results.ks_hsa_miR_1247_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_125a_5p	=	GenTable(	GOdata_hsa_miR_125a_5p	, KS = results.ks_hsa_miR_125a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_125b_5p	=	GenTable(	GOdata_hsa_miR_125b_5p	, KS = results.ks_hsa_miR_125b_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_126_3p	=	GenTable(	GOdata_hsa_miR_126_3p	, KS = results.ks_hsa_miR_126_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_126_5p	=	GenTable(	GOdata_hsa_miR_126_5p	, KS = results.ks_hsa_miR_126_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_128_3p	=	GenTable(	GOdata_hsa_miR_128_3p	, KS = results.ks_hsa_miR_128_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_129_2_3p	=	GenTable(	GOdata_hsa_miR_129_2_3p	, KS = results.ks_hsa_miR_129_2_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_1296_5p	=	GenTable(	GOdata_hsa_miR_1296_5p	, KS = results.ks_hsa_miR_1296_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_1306_3p	=	GenTable(	GOdata_hsa_miR_1306_3p	, KS = results.ks_hsa_miR_1306_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_130a_3p	=	GenTable(	GOdata_hsa_miR_130a_3p	, KS = results.ks_hsa_miR_130a_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_133a_3p	=	GenTable(	GOdata_hsa_miR_133a_3p	, KS = results.ks_hsa_miR_133a_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_133a_5p	=	GenTable(	GOdata_hsa_miR_133a_5p	, KS = results.ks_hsa_miR_133a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_133b	=	GenTable(	GOdata_hsa_miR_133b	, KS = results.ks_hsa_miR_133b	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_135a_5p	=	GenTable(	GOdata_hsa_miR_135a_5p	, KS = results.ks_hsa_miR_135a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_138_5p	=	GenTable(	GOdata_hsa_miR_138_5p	, KS = results.ks_hsa_miR_138_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_139_3p	=	GenTable(	GOdata_hsa_miR_139_3p	, KS = results.ks_hsa_miR_139_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_139_5p	=	GenTable(	GOdata_hsa_miR_139_5p	, KS = results.ks_hsa_miR_139_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_140_3p	=	GenTable(	GOdata_hsa_miR_140_3p	, KS = results.ks_hsa_miR_140_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_140_5p	=	GenTable(	GOdata_hsa_miR_140_5p	, KS = results.ks_hsa_miR_140_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_142_3p	=	GenTable(	GOdata_hsa_miR_142_3p	, KS = results.ks_hsa_miR_142_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_142_5p	=	GenTable(	GOdata_hsa_miR_142_5p	, KS = results.ks_hsa_miR_142_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_143_3p	=	GenTable(	GOdata_hsa_miR_143_3p	, KS = results.ks_hsa_miR_143_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_143_5p	=	GenTable(	GOdata_hsa_miR_143_5p	, KS = results.ks_hsa_miR_143_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_144_3p	=	GenTable(	GOdata_hsa_miR_144_3p	, KS = results.ks_hsa_miR_144_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_144_5p	=	GenTable(	GOdata_hsa_miR_144_5p	, KS = results.ks_hsa_miR_144_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_145_3p	=	GenTable(	GOdata_hsa_miR_145_3p	, KS = results.ks_hsa_miR_145_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_145_5p	=	GenTable(	GOdata_hsa_miR_145_5p	, KS = results.ks_hsa_miR_145_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_146a_5p	=	GenTable(	GOdata_hsa_miR_146a_5p	, KS = results.ks_hsa_miR_146a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_146b_5p	=	GenTable(	GOdata_hsa_miR_146b_5p	, KS = results.ks_hsa_miR_146b_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_147b	=	GenTable(	GOdata_hsa_miR_147b	, KS = results.ks_hsa_miR_147b	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_150_5p	=	GenTable(	GOdata_hsa_miR_150_5p	, KS = results.ks_hsa_miR_150_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_151a_5p	=	GenTable(	GOdata_hsa_miR_151a_5p	, KS = results.ks_hsa_miR_151a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_151b	=	GenTable(	GOdata_hsa_miR_151b	, KS = results.ks_hsa_miR_151b	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_152_3p	=	GenTable(	GOdata_hsa_miR_152_3p	, KS = results.ks_hsa_miR_152_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_1537_3p	=	GenTable(	GOdata_hsa_miR_1537_3p	, KS = results.ks_hsa_miR_1537_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_15a_5p	=	GenTable(	GOdata_hsa_miR_15a_5p	, KS = results.ks_hsa_miR_15a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_15b_3p	=	GenTable(	GOdata_hsa_miR_15b_3p	, KS = results.ks_hsa_miR_15b_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_15b_5p	=	GenTable(	GOdata_hsa_miR_15b_5p	, KS = results.ks_hsa_miR_15b_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_16_2_3p	=	GenTable(	GOdata_hsa_miR_16_2_3p	, KS = results.ks_hsa_miR_16_2_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_16_5p	=	GenTable(	GOdata_hsa_miR_16_5p	, KS = results.ks_hsa_miR_16_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_181a_2_3p	=	GenTable(	GOdata_hsa_miR_181a_2_3p	, KS = results.ks_hsa_miR_181a_2_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_181a_3p	=	GenTable(	GOdata_hsa_miR_181a_3p	, KS = results.ks_hsa_miR_181a_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_182_3p	=	GenTable(	GOdata_hsa_miR_182_3p	, KS = results.ks_hsa_miR_182_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_183_3p	=	GenTable(	GOdata_hsa_miR_183_3p	, KS = results.ks_hsa_miR_183_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_184	=	GenTable(	GOdata_hsa_miR_184	, KS = results.ks_hsa_miR_184	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_185_5p	=	GenTable(	GOdata_hsa_miR_185_5p	, KS = results.ks_hsa_miR_185_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_186_5p	=	GenTable(	GOdata_hsa_miR_186_5p	, KS = results.ks_hsa_miR_186_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_187_5p	=	GenTable(	GOdata_hsa_miR_187_5p	, KS = results.ks_hsa_miR_187_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_18a_5p	=	GenTable(	GOdata_hsa_miR_18a_5p	, KS = results.ks_hsa_miR_18a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_18b_5p	=	GenTable(	GOdata_hsa_miR_18b_5p	, KS = results.ks_hsa_miR_18b_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_190a_5p	=	GenTable(	GOdata_hsa_miR_190a_5p	, KS = results.ks_hsa_miR_190a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_191_5p	=	GenTable(	GOdata_hsa_miR_191_5p	, KS = results.ks_hsa_miR_191_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_1913	=	GenTable(	GOdata_hsa_miR_1913	, KS = results.ks_hsa_miR_1913	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_193a_5p	=	GenTable(	GOdata_hsa_miR_193a_5p	, KS = results.ks_hsa_miR_193a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_195_5p	=	GenTable(	GOdata_hsa_miR_195_5p	, KS = results.ks_hsa_miR_195_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_199a_3p	=	GenTable(	GOdata_hsa_miR_199a_3p	, KS = results.ks_hsa_miR_199a_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_199a_5p	=	GenTable(	GOdata_hsa_miR_199a_5p	, KS = results.ks_hsa_miR_199a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_199b_5p	=	GenTable(	GOdata_hsa_miR_199b_5p	, KS = results.ks_hsa_miR_199b_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_203a_3p	=	GenTable(	GOdata_hsa_miR_203a_3p	, KS = results.ks_hsa_miR_203a_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_204_5p	=	GenTable(	GOdata_hsa_miR_204_5p	, KS = results.ks_hsa_miR_204_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_205_3p	=	GenTable(	GOdata_hsa_miR_205_3p	, KS = results.ks_hsa_miR_205_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_20a_5p	=	GenTable(	GOdata_hsa_miR_20a_5p	, KS = results.ks_hsa_miR_20a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_20b_5p	=	GenTable(	GOdata_hsa_miR_20b_5p	, KS = results.ks_hsa_miR_20b_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_21_3p	=	GenTable(	GOdata_hsa_miR_21_3p	, KS = results.ks_hsa_miR_21_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_214_3p	=	GenTable(	GOdata_hsa_miR_214_3p	, KS = results.ks_hsa_miR_214_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_214_5p	=	GenTable(	GOdata_hsa_miR_214_5p	, KS = results.ks_hsa_miR_214_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_218_5p	=	GenTable(	GOdata_hsa_miR_218_5p	, KS = results.ks_hsa_miR_218_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_22_5p	=	GenTable(	GOdata_hsa_miR_22_5p	, KS = results.ks_hsa_miR_22_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_221_3p	=	GenTable(	GOdata_hsa_miR_221_3p	, KS = results.ks_hsa_miR_221_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_221_5p	=	GenTable(	GOdata_hsa_miR_221_5p	, KS = results.ks_hsa_miR_221_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_222_3p	=	GenTable(	GOdata_hsa_miR_222_3p	, KS = results.ks_hsa_miR_222_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_223_3p	=	GenTable(	GOdata_hsa_miR_223_3p	, KS = results.ks_hsa_miR_223_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_223_5p	=	GenTable(	GOdata_hsa_miR_223_5p	, KS = results.ks_hsa_miR_223_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_224_3p	=	GenTable(	GOdata_hsa_miR_224_3p	, KS = results.ks_hsa_miR_224_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_2277_3p	=	GenTable(	GOdata_hsa_miR_2277_3p	, KS = results.ks_hsa_miR_2277_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_23a_3p	=	GenTable(	GOdata_hsa_miR_23a_3p	, KS = results.ks_hsa_miR_23a_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_23a_5p	=	GenTable(	GOdata_hsa_miR_23a_5p	, KS = results.ks_hsa_miR_23a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_23b_3p	=	GenTable(	GOdata_hsa_miR_23b_3p	, KS = results.ks_hsa_miR_23b_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_24_3p	=	GenTable(	GOdata_hsa_miR_24_3p	, KS = results.ks_hsa_miR_24_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_26a_1_3p	=	GenTable(	GOdata_hsa_miR_26a_1_3p	, KS = results.ks_hsa_miR_26a_1_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_26a_5p	=	GenTable(	GOdata_hsa_miR_26a_5p	, KS = results.ks_hsa_miR_26a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_26b_3p	=	GenTable(	GOdata_hsa_miR_26b_3p	, KS = results.ks_hsa_miR_26b_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_26b_5p	=	GenTable(	GOdata_hsa_miR_26b_5p	, KS = results.ks_hsa_miR_26b_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_27a_3p	=	GenTable(	GOdata_hsa_miR_27a_3p	, KS = results.ks_hsa_miR_27a_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_27a_5p	=	GenTable(	GOdata_hsa_miR_27a_5p	, KS = results.ks_hsa_miR_27a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_27b_3p	=	GenTable(	GOdata_hsa_miR_27b_3p	, KS = results.ks_hsa_miR_27b_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_28_5p	=	GenTable(	GOdata_hsa_miR_28_5p	, KS = results.ks_hsa_miR_28_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_2861	=	GenTable(	GOdata_hsa_miR_2861	, KS = results.ks_hsa_miR_2861	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_29a_3p	=	GenTable(	GOdata_hsa_miR_29a_3p	, KS = results.ks_hsa_miR_29a_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_29a_5p	=	GenTable(	GOdata_hsa_miR_29a_5p	, KS = results.ks_hsa_miR_29a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_29b_1_5p	=	GenTable(	GOdata_hsa_miR_29b_1_5p	, KS = results.ks_hsa_miR_29b_1_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_29b_2_5p	=	GenTable(	GOdata_hsa_miR_29b_2_5p	, KS = results.ks_hsa_miR_29b_2_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_29b_3p	=	GenTable(	GOdata_hsa_miR_29b_3p	, KS = results.ks_hsa_miR_29b_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_29c_3p	=	GenTable(	GOdata_hsa_miR_29c_3p	, KS = results.ks_hsa_miR_29c_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_29c_5p	=	GenTable(	GOdata_hsa_miR_29c_5p	, KS = results.ks_hsa_miR_29c_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_3065_3p	=	GenTable(	GOdata_hsa_miR_3065_3p	, KS = results.ks_hsa_miR_3065_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_3065_5p	=	GenTable(	GOdata_hsa_miR_3065_5p	, KS = results.ks_hsa_miR_3065_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_30a_3p	=	GenTable(	GOdata_hsa_miR_30a_3p	, KS = results.ks_hsa_miR_30a_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_30a_5p	=	GenTable(	GOdata_hsa_miR_30a_5p	, KS = results.ks_hsa_miR_30a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_30b_5p	=	GenTable(	GOdata_hsa_miR_30b_5p	, KS = results.ks_hsa_miR_30b_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_30c_2_3p	=	GenTable(	GOdata_hsa_miR_30c_2_3p	, KS = results.ks_hsa_miR_30c_2_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_30c_5p	=	GenTable(	GOdata_hsa_miR_30c_5p	, KS = results.ks_hsa_miR_30c_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_30d_5p	=	GenTable(	GOdata_hsa_miR_30d_5p	, KS = results.ks_hsa_miR_30d_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_30e_3p	=	GenTable(	GOdata_hsa_miR_30e_3p	, KS = results.ks_hsa_miR_30e_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_3149	=	GenTable(	GOdata_hsa_miR_3149	, KS = results.ks_hsa_miR_3149	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_3188	=	GenTable(	GOdata_hsa_miR_3188	, KS = results.ks_hsa_miR_3188	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_32_5p	=	GenTable(	GOdata_hsa_miR_32_5p	, KS = results.ks_hsa_miR_32_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_324_3p	=	GenTable(	GOdata_hsa_miR_324_3p	, KS = results.ks_hsa_miR_324_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_326	=	GenTable(	GOdata_hsa_miR_326	, KS = results.ks_hsa_miR_326	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_328_3p	=	GenTable(	GOdata_hsa_miR_328_3p	, KS = results.ks_hsa_miR_328_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_331_3p	=	GenTable(	GOdata_hsa_miR_331_3p	, KS = results.ks_hsa_miR_331_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_335_5p	=	GenTable(	GOdata_hsa_miR_335_5p	, KS = results.ks_hsa_miR_335_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_338_3p	=	GenTable(	GOdata_hsa_miR_338_3p	, KS = results.ks_hsa_miR_338_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_339_5p	=	GenTable(	GOdata_hsa_miR_339_5p	, KS = results.ks_hsa_miR_339_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_340_5p	=	GenTable(	GOdata_hsa_miR_340_5p	, KS = results.ks_hsa_miR_340_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_342_3p	=	GenTable(	GOdata_hsa_miR_342_3p	, KS = results.ks_hsa_miR_342_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_342_5p	=	GenTable(	GOdata_hsa_miR_342_5p	, KS = results.ks_hsa_miR_342_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_34a_3p	=	GenTable(	GOdata_hsa_miR_34a_3p	, KS = results.ks_hsa_miR_34a_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_34a_5p	=	GenTable(	GOdata_hsa_miR_34a_5p	, KS = results.ks_hsa_miR_34a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_34b_5p	=	GenTable(	GOdata_hsa_miR_34b_5p	, KS = results.ks_hsa_miR_34b_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_34c_5p	=	GenTable(	GOdata_hsa_miR_34c_5p	, KS = results.ks_hsa_miR_34c_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_3607_3p	=	GenTable(	GOdata_hsa_miR_3607_3p	, KS = results.ks_hsa_miR_3607_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_361_3p	=	GenTable(	GOdata_hsa_miR_361_3p	, KS = results.ks_hsa_miR_361_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_362_3p	=	GenTable(	GOdata_hsa_miR_362_3p	, KS = results.ks_hsa_miR_362_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_362_5p	=	GenTable(	GOdata_hsa_miR_362_5p	, KS = results.ks_hsa_miR_362_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_3620_5p	=	GenTable(	GOdata_hsa_miR_3620_5p	, KS = results.ks_hsa_miR_3620_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_363_3p	=	GenTable(	GOdata_hsa_miR_363_3p	, KS = results.ks_hsa_miR_363_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_3659	=	GenTable(	GOdata_hsa_miR_3659	, KS = results.ks_hsa_miR_3659	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_365a_3p	=	GenTable(	GOdata_hsa_miR_365a_3p	, KS = results.ks_hsa_miR_365a_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_3660	=	GenTable(	GOdata_hsa_miR_3660	, KS = results.ks_hsa_miR_3660	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_3663_5p	=	GenTable(	GOdata_hsa_miR_3663_5p	, KS = results.ks_hsa_miR_3663_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_371b_5p	=	GenTable(	GOdata_hsa_miR_371b_5p	, KS = results.ks_hsa_miR_371b_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_374a_5p	=	GenTable(	GOdata_hsa_miR_374a_5p	, KS = results.ks_hsa_miR_374a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_374b_5p	=	GenTable(	GOdata_hsa_miR_374b_5p	, KS = results.ks_hsa_miR_374b_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_374c_5p	=	GenTable(	GOdata_hsa_miR_374c_5p	, KS = results.ks_hsa_miR_374c_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_424_5p	=	GenTable(	GOdata_hsa_miR_424_5p	, KS = results.ks_hsa_miR_424_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4252	=	GenTable(	GOdata_hsa_miR_4252	, KS = results.ks_hsa_miR_4252	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4254	=	GenTable(	GOdata_hsa_miR_4254	, KS = results.ks_hsa_miR_4254	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4284	=	GenTable(	GOdata_hsa_miR_4284	, KS = results.ks_hsa_miR_4284	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4290	=	GenTable(	GOdata_hsa_miR_4290	, KS = results.ks_hsa_miR_4290	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4306	=	GenTable(	GOdata_hsa_miR_4306	, KS = results.ks_hsa_miR_4306	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4317	=	GenTable(	GOdata_hsa_miR_4317	, KS = results.ks_hsa_miR_4317	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4318	=	GenTable(	GOdata_hsa_miR_4318	, KS = results.ks_hsa_miR_4318	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4324	=	GenTable(	GOdata_hsa_miR_4324	, KS = results.ks_hsa_miR_4324	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4328	=	GenTable(	GOdata_hsa_miR_4328	, KS = results.ks_hsa_miR_4328	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4440	=	GenTable(	GOdata_hsa_miR_4440	, KS = results.ks_hsa_miR_4440	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4443	=	GenTable(	GOdata_hsa_miR_4443	, KS = results.ks_hsa_miR_4443	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4481	=	GenTable(	GOdata_hsa_miR_4481	, KS = results.ks_hsa_miR_4481	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_449a	=	GenTable(	GOdata_hsa_miR_449a	, KS = results.ks_hsa_miR_449a	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_450a_5p	=	GenTable(	GOdata_hsa_miR_450a_5p	, KS = results.ks_hsa_miR_450a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4516	=	GenTable(	GOdata_hsa_miR_4516	, KS = results.ks_hsa_miR_4516	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_451a	=	GenTable(	GOdata_hsa_miR_451a	, KS = results.ks_hsa_miR_451a	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_452_5p	=	GenTable(	GOdata_hsa_miR_452_5p	, KS = results.ks_hsa_miR_452_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4521	=	GenTable(	GOdata_hsa_miR_4521	, KS = results.ks_hsa_miR_4521	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4532	=	GenTable(	GOdata_hsa_miR_4532	, KS = results.ks_hsa_miR_4532	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_454_3p	=	GenTable(	GOdata_hsa_miR_454_3p	, KS = results.ks_hsa_miR_454_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_455_3p	=	GenTable(	GOdata_hsa_miR_455_3p	, KS = results.ks_hsa_miR_455_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_455_5p	=	GenTable(	GOdata_hsa_miR_455_5p	, KS = results.ks_hsa_miR_455_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4634	=	GenTable(	GOdata_hsa_miR_4634	, KS = results.ks_hsa_miR_4634	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4655_3p	=	GenTable(	GOdata_hsa_miR_4655_3p	, KS = results.ks_hsa_miR_4655_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4695_3p	=	GenTable(	GOdata_hsa_miR_4695_3p	, KS = results.ks_hsa_miR_4695_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4716_5p	=	GenTable(	GOdata_hsa_miR_4716_5p	, KS = results.ks_hsa_miR_4716_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4730	=	GenTable(	GOdata_hsa_miR_4730	, KS = results.ks_hsa_miR_4730	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4731_3p	=	GenTable(	GOdata_hsa_miR_4731_3p	, KS = results.ks_hsa_miR_4731_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4763_5p	=	GenTable(	GOdata_hsa_miR_4763_5p	, KS = results.ks_hsa_miR_4763_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4770	=	GenTable(	GOdata_hsa_miR_4770	, KS = results.ks_hsa_miR_4770	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_4793_3p	=	GenTable(	GOdata_hsa_miR_4793_3p	, KS = results.ks_hsa_miR_4793_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_483_3p	=	GenTable(	GOdata_hsa_miR_483_3p	, KS = results.ks_hsa_miR_483_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_486_3p	=	GenTable(	GOdata_hsa_miR_486_3p	, KS = results.ks_hsa_miR_486_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_486_5p	=	GenTable(	GOdata_hsa_miR_486_5p	, KS = results.ks_hsa_miR_486_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_489_3p	=	GenTable(	GOdata_hsa_miR_489_3p	, KS = results.ks_hsa_miR_489_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_490_3p	=	GenTable(	GOdata_hsa_miR_490_3p	, KS = results.ks_hsa_miR_490_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_491_3p	=	GenTable(	GOdata_hsa_miR_491_3p	, KS = results.ks_hsa_miR_491_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_497_3p	=	GenTable(	GOdata_hsa_miR_497_3p	, KS = results.ks_hsa_miR_497_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_497_5p	=	GenTable(	GOdata_hsa_miR_497_5p	, KS = results.ks_hsa_miR_497_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_5001_5p	=	GenTable(	GOdata_hsa_miR_5001_5p	, KS = results.ks_hsa_miR_5001_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_500a_3p	=	GenTable(	GOdata_hsa_miR_500a_3p	, KS = results.ks_hsa_miR_500a_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_500b_5p	=	GenTable(	GOdata_hsa_miR_500b_5p	, KS = results.ks_hsa_miR_500b_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_501_5p	=	GenTable(	GOdata_hsa_miR_501_5p	, KS = results.ks_hsa_miR_501_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_502_3p	=	GenTable(	GOdata_hsa_miR_502_3p	, KS = results.ks_hsa_miR_502_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_504_3p	=	GenTable(	GOdata_hsa_miR_504_3p	, KS = results.ks_hsa_miR_504_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_505_5p	=	GenTable(	GOdata_hsa_miR_505_5p	, KS = results.ks_hsa_miR_505_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_511_3p	=	GenTable(	GOdata_hsa_miR_511_3p	, KS = results.ks_hsa_miR_511_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_516b_5p	=	GenTable(	GOdata_hsa_miR_516b_5p	, KS = results.ks_hsa_miR_516b_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_517_5p	=	GenTable(	GOdata_hsa_miR_517_5p	, KS = results.ks_hsa_miR_517_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_517a_3p	=	GenTable(	GOdata_hsa_miR_517a_3p	, KS = results.ks_hsa_miR_517a_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_517c_3p	=	GenTable(	GOdata_hsa_miR_517c_3p	, KS = results.ks_hsa_miR_517c_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_521	=	GenTable(	GOdata_hsa_miR_521	, KS = results.ks_hsa_miR_521	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_522_3p	=	GenTable(	GOdata_hsa_miR_522_3p	, KS = results.ks_hsa_miR_522_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_532_3p	=	GenTable(	GOdata_hsa_miR_532_3p	, KS = results.ks_hsa_miR_532_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_532_5p	=	GenTable(	GOdata_hsa_miR_532_5p	, KS = results.ks_hsa_miR_532_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_548aa	=	GenTable(	GOdata_hsa_miR_548aa	, KS = results.ks_hsa_miR_548aa	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_548aw	=	GenTable(	GOdata_hsa_miR_548aw	, KS = results.ks_hsa_miR_548aw	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_548b_3p	=	GenTable(	GOdata_hsa_miR_548b_3p	, KS = results.ks_hsa_miR_548b_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_548c_3p	=	GenTable(	GOdata_hsa_miR_548c_3p	, KS = results.ks_hsa_miR_548c_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_548f_3p	=	GenTable(	GOdata_hsa_miR_548f_3p	, KS = results.ks_hsa_miR_548f_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_548q	=	GenTable(	GOdata_hsa_miR_548q	, KS = results.ks_hsa_miR_548q	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_548x_3p	=	GenTable(	GOdata_hsa_miR_548x_3p	, KS = results.ks_hsa_miR_548x_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_551b_3p	=	GenTable(	GOdata_hsa_miR_551b_3p	, KS = results.ks_hsa_miR_551b_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_5701	=	GenTable(	GOdata_hsa_miR_5701	, KS = results.ks_hsa_miR_5701	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_582_5p	=	GenTable(	GOdata_hsa_miR_582_5p	, KS = results.ks_hsa_miR_582_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_585_3p	=	GenTable(	GOdata_hsa_miR_585_3p	, KS = results.ks_hsa_miR_585_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_590_5p	=	GenTable(	GOdata_hsa_miR_590_5p	, KS = results.ks_hsa_miR_590_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_595	=	GenTable(	GOdata_hsa_miR_595	, KS = results.ks_hsa_miR_595	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_598_3p	=	GenTable(	GOdata_hsa_miR_598_3p	, KS = results.ks_hsa_miR_598_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_6068	=	GenTable(	GOdata_hsa_miR_6068	, KS = results.ks_hsa_miR_6068	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_6073	=	GenTable(	GOdata_hsa_miR_6073	, KS = results.ks_hsa_miR_6073	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_6075	=	GenTable(	GOdata_hsa_miR_6075	, KS = results.ks_hsa_miR_6075	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_610	=	GenTable(	GOdata_hsa_miR_610	, KS = results.ks_hsa_miR_610	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_624_5p	=	GenTable(	GOdata_hsa_miR_624_5p	, KS = results.ks_hsa_miR_624_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_628_5p	=	GenTable(	GOdata_hsa_miR_628_5p	, KS = results.ks_hsa_miR_628_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_642b_5p	=	GenTable(	GOdata_hsa_miR_642b_5p	, KS = results.ks_hsa_miR_642b_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_6500_5p	=	GenTable(	GOdata_hsa_miR_6500_5p	, KS = results.ks_hsa_miR_6500_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_6516_3p	=	GenTable(	GOdata_hsa_miR_6516_3p	, KS = results.ks_hsa_miR_6516_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_652_3p	=	GenTable(	GOdata_hsa_miR_652_3p	, KS = results.ks_hsa_miR_652_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_653_3p	=	GenTable(	GOdata_hsa_miR_653_3p	, KS = results.ks_hsa_miR_653_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_660_5p	=	GenTable(	GOdata_hsa_miR_660_5p	, KS = results.ks_hsa_miR_660_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_664a_3p	=	GenTable(	GOdata_hsa_miR_664a_3p	, KS = results.ks_hsa_miR_664a_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_664b_3p	=	GenTable(	GOdata_hsa_miR_664b_3p	, KS = results.ks_hsa_miR_664b_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_6716_3p	=	GenTable(	GOdata_hsa_miR_6716_3p	, KS = results.ks_hsa_miR_6716_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_6722_5p	=	GenTable(	GOdata_hsa_miR_6722_5p	, KS = results.ks_hsa_miR_6722_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_6730_3p	=	GenTable(	GOdata_hsa_miR_6730_3p	, KS = results.ks_hsa_miR_6730_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_6743_3p	=	GenTable(	GOdata_hsa_miR_6743_3p	, KS = results.ks_hsa_miR_6743_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_6771_5p	=	GenTable(	GOdata_hsa_miR_6771_5p	, KS = results.ks_hsa_miR_6771_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_6779_3p	=	GenTable(	GOdata_hsa_miR_6779_3p	, KS = results.ks_hsa_miR_6779_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_6794_3p	=	GenTable(	GOdata_hsa_miR_6794_3p	, KS = results.ks_hsa_miR_6794_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_6804_5p	=	GenTable(	GOdata_hsa_miR_6804_5p	, KS = results.ks_hsa_miR_6804_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_6806_5p	=	GenTable(	GOdata_hsa_miR_6806_5p	, KS = results.ks_hsa_miR_6806_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_6817_5p	=	GenTable(	GOdata_hsa_miR_6817_5p	, KS = results.ks_hsa_miR_6817_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_6826_5p	=	GenTable(	GOdata_hsa_miR_6826_5p	, KS = results.ks_hsa_miR_6826_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_6865_3p	=	GenTable(	GOdata_hsa_miR_6865_3p	, KS = results.ks_hsa_miR_6865_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_6872_3p	=	GenTable(	GOdata_hsa_miR_6872_3p	, KS = results.ks_hsa_miR_6872_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_6886_3p	=	GenTable(	GOdata_hsa_miR_6886_3p	, KS = results.ks_hsa_miR_6886_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_6891_3p	=	GenTable(	GOdata_hsa_miR_6891_3p	, KS = results.ks_hsa_miR_6891_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_6895_5p	=	GenTable(	GOdata_hsa_miR_6895_5p	, KS = results.ks_hsa_miR_6895_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_7108_3p	=	GenTable(	GOdata_hsa_miR_7108_3p	, KS = results.ks_hsa_miR_7108_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_7111_3p	=	GenTable(	GOdata_hsa_miR_7111_3p	, KS = results.ks_hsa_miR_7111_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_7159_5p	=	GenTable(	GOdata_hsa_miR_7159_5p	, KS = results.ks_hsa_miR_7159_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_744_5p	=	GenTable(	GOdata_hsa_miR_744_5p	, KS = results.ks_hsa_miR_744_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_770_5p	=	GenTable(	GOdata_hsa_miR_770_5p	, KS = results.ks_hsa_miR_770_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_7704	=	GenTable(	GOdata_hsa_miR_7704	, KS = results.ks_hsa_miR_7704	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_8063	=	GenTable(	GOdata_hsa_miR_8063	, KS = results.ks_hsa_miR_8063	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_8077	=	GenTable(	GOdata_hsa_miR_8077	, KS = results.ks_hsa_miR_8077	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_874_5p	=	GenTable(	GOdata_hsa_miR_874_5p	, KS = results.ks_hsa_miR_874_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_887_3p	=	GenTable(	GOdata_hsa_miR_887_3p	, KS = results.ks_hsa_miR_887_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_940	=	GenTable(	GOdata_hsa_miR_940	, KS = results.ks_hsa_miR_940	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_95_3p	=	GenTable(	GOdata_hsa_miR_95_3p	, KS = results.ks_hsa_miR_95_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_98_5p	=	GenTable(	GOdata_hsa_miR_98_5p	, KS = results.ks_hsa_miR_98_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_99a_3p	=	GenTable(	GOdata_hsa_miR_99a_3p	, KS = results.ks_hsa_miR_99a_3p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_99a_5p	=	GenTable(	GOdata_hsa_miR_99a_5p	, KS = results.ks_hsa_miR_99a_5p	, orderBy = "KS", topNodes = 10)
allRes_hsa_miR_99b_5p	=	GenTable(	GOdata_hsa_miR_99b_5p	, KS = results.ks_hsa_miR_99b_5p	, orderBy = "KS", topNodes = 10)

#AllResData

allRes_data_hsa_let_7a_5p	=	allRes_hsa_let_7a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_let_7b_5p	=	allRes_hsa_let_7b_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_let_7c_5p	=	allRes_hsa_let_7c_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_let_7d_5p	=	allRes_hsa_let_7d_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_let_7e_5p	=	allRes_hsa_let_7e_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_let_7f_5p	=	allRes_hsa_let_7f_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_let_7g_5p	=	allRes_hsa_let_7g_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_let_7i_5p	=	allRes_hsa_let_7i_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_1_3p	=	allRes_hsa_miR_1_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_100_5p	=	allRes_hsa_miR_100_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_101_3p	=	allRes_hsa_miR_101_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_101_5p	=	allRes_hsa_miR_101_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_103a_3p	=	allRes_hsa_miR_103a_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_106b_3p	=	allRes_hsa_miR_106b_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_107	=	allRes_hsa_miR_107	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_10a_5p	=	allRes_hsa_miR_10a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_10b_5p	=	allRes_hsa_miR_10b_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_1227_5p	=	allRes_hsa_miR_1227_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_1238_5p	=	allRes_hsa_miR_1238_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_1247_5p	=	allRes_hsa_miR_1247_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_125a_5p	=	allRes_hsa_miR_125a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_125b_5p	=	allRes_hsa_miR_125b_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_126_3p	=	allRes_hsa_miR_126_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_126_5p	=	allRes_hsa_miR_126_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_128_3p	=	allRes_hsa_miR_128_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_129_2_3p	=	allRes_hsa_miR_129_2_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_1296_5p	=	allRes_hsa_miR_1296_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_1306_3p	=	allRes_hsa_miR_1306_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_130a_3p	=	allRes_hsa_miR_130a_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_133a_3p	=	allRes_hsa_miR_133a_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_133a_5p	=	allRes_hsa_miR_133a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_133b	=	allRes_hsa_miR_133b	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_135a_5p	=	allRes_hsa_miR_135a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_138_5p	=	allRes_hsa_miR_138_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_139_3p	=	allRes_hsa_miR_139_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_139_5p	=	allRes_hsa_miR_139_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_140_3p	=	allRes_hsa_miR_140_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_140_5p	=	allRes_hsa_miR_140_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_142_3p	=	allRes_hsa_miR_142_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_142_5p	=	allRes_hsa_miR_142_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_143_3p	=	allRes_hsa_miR_143_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_143_5p	=	allRes_hsa_miR_143_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_144_3p	=	allRes_hsa_miR_144_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_144_5p	=	allRes_hsa_miR_144_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_145_3p	=	allRes_hsa_miR_145_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_145_5p	=	allRes_hsa_miR_145_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_146a_5p	=	allRes_hsa_miR_146a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_146b_5p	=	allRes_hsa_miR_146b_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_147b	=	allRes_hsa_miR_147b	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_150_5p	=	allRes_hsa_miR_150_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_151a_5p	=	allRes_hsa_miR_151a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_151b	=	allRes_hsa_miR_151b	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_152_3p	=	allRes_hsa_miR_152_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_1537_3p	=	allRes_hsa_miR_1537_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_15a_5p	=	allRes_hsa_miR_15a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_15b_3p	=	allRes_hsa_miR_15b_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_15b_5p	=	allRes_hsa_miR_15b_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_16_2_3p	=	allRes_hsa_miR_16_2_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_16_5p	=	allRes_hsa_miR_16_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_181a_2_3p	=	allRes_hsa_miR_181a_2_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_181a_3p	=	allRes_hsa_miR_181a_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_182_3p	=	allRes_hsa_miR_182_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_183_3p	=	allRes_hsa_miR_183_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_184	=	allRes_hsa_miR_184	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_185_5p	=	allRes_hsa_miR_185_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_186_5p	=	allRes_hsa_miR_186_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_187_5p	=	allRes_hsa_miR_187_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_18a_5p	=	allRes_hsa_miR_18a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_18b_5p	=	allRes_hsa_miR_18b_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_190a_5p	=	allRes_hsa_miR_190a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_191_5p	=	allRes_hsa_miR_191_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_1913	=	allRes_hsa_miR_1913	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_193a_5p	=	allRes_hsa_miR_193a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_195_5p	=	allRes_hsa_miR_195_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_199a_3p	=	allRes_hsa_miR_199a_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_199a_5p	=	allRes_hsa_miR_199a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_199b_5p	=	allRes_hsa_miR_199b_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_203a_3p	=	allRes_hsa_miR_203a_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_204_5p	=	allRes_hsa_miR_204_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_205_3p	=	allRes_hsa_miR_205_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_20a_5p	=	allRes_hsa_miR_20a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_20b_5p	=	allRes_hsa_miR_20b_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_21_3p	=	allRes_hsa_miR_21_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_214_3p	=	allRes_hsa_miR_214_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_214_5p	=	allRes_hsa_miR_214_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_218_5p	=	allRes_hsa_miR_218_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_22_5p	=	allRes_hsa_miR_22_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_221_3p	=	allRes_hsa_miR_221_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_221_5p	=	allRes_hsa_miR_221_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_222_3p	=	allRes_hsa_miR_222_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_223_3p	=	allRes_hsa_miR_223_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_223_5p	=	allRes_hsa_miR_223_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_224_3p	=	allRes_hsa_miR_224_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_2277_3p	=	allRes_hsa_miR_2277_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_23a_3p	=	allRes_hsa_miR_23a_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_23a_5p	=	allRes_hsa_miR_23a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_23b_3p	=	allRes_hsa_miR_23b_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_24_3p	=	allRes_hsa_miR_24_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_26a_1_3p	=	allRes_hsa_miR_26a_1_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_26a_5p	=	allRes_hsa_miR_26a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_26b_3p	=	allRes_hsa_miR_26b_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_26b_5p	=	allRes_hsa_miR_26b_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_27a_3p	=	allRes_hsa_miR_27a_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_27a_5p	=	allRes_hsa_miR_27a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_27b_3p	=	allRes_hsa_miR_27b_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_28_5p	=	allRes_hsa_miR_28_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_2861	=	allRes_hsa_miR_2861	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_29a_3p	=	allRes_hsa_miR_29a_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_29a_5p	=	allRes_hsa_miR_29a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_29b_1_5p	=	allRes_hsa_miR_29b_1_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_29b_2_5p	=	allRes_hsa_miR_29b_2_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_29b_3p	=	allRes_hsa_miR_29b_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_29c_3p	=	allRes_hsa_miR_29c_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_29c_5p	=	allRes_hsa_miR_29c_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_3065_3p	=	allRes_hsa_miR_3065_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_3065_5p	=	allRes_hsa_miR_3065_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_30a_3p	=	allRes_hsa_miR_30a_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_30a_5p	=	allRes_hsa_miR_30a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_30b_5p	=	allRes_hsa_miR_30b_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_30c_2_3p	=	allRes_hsa_miR_30c_2_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_30c_5p	=	allRes_hsa_miR_30c_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_30d_5p	=	allRes_hsa_miR_30d_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_30e_3p	=	allRes_hsa_miR_30e_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_3149	=	allRes_hsa_miR_3149	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_3188	=	allRes_hsa_miR_3188	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_32_5p	=	allRes_hsa_miR_32_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_324_3p	=	allRes_hsa_miR_324_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_326	=	allRes_hsa_miR_326	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_328_3p	=	allRes_hsa_miR_328_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_331_3p	=	allRes_hsa_miR_331_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_335_5p	=	allRes_hsa_miR_335_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_338_3p	=	allRes_hsa_miR_338_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_339_5p	=	allRes_hsa_miR_339_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_340_5p	=	allRes_hsa_miR_340_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_342_3p	=	allRes_hsa_miR_342_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_342_5p	=	allRes_hsa_miR_342_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_34a_3p	=	allRes_hsa_miR_34a_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_34a_5p	=	allRes_hsa_miR_34a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_34b_5p	=	allRes_hsa_miR_34b_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_34c_5p	=	allRes_hsa_miR_34c_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_3607_3p	=	allRes_hsa_miR_3607_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_361_3p	=	allRes_hsa_miR_361_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_362_3p	=	allRes_hsa_miR_362_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_362_5p	=	allRes_hsa_miR_362_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_3620_5p	=	allRes_hsa_miR_3620_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_363_3p	=	allRes_hsa_miR_363_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_3659	=	allRes_hsa_miR_3659	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_365a_3p	=	allRes_hsa_miR_365a_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_3660	=	allRes_hsa_miR_3660	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_3663_5p	=	allRes_hsa_miR_3663_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_371b_5p	=	allRes_hsa_miR_371b_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_374a_5p	=	allRes_hsa_miR_374a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_374b_5p	=	allRes_hsa_miR_374b_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_374c_5p	=	allRes_hsa_miR_374c_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_424_5p	=	allRes_hsa_miR_424_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4252	=	allRes_hsa_miR_4252	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4254	=	allRes_hsa_miR_4254	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4284	=	allRes_hsa_miR_4284	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4290	=	allRes_hsa_miR_4290	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4306	=	allRes_hsa_miR_4306	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4317	=	allRes_hsa_miR_4317	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4318	=	allRes_hsa_miR_4318	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4324	=	allRes_hsa_miR_4324	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4328	=	allRes_hsa_miR_4328	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4440	=	allRes_hsa_miR_4440	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4443	=	allRes_hsa_miR_4443	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4481	=	allRes_hsa_miR_4481	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_449a	=	allRes_hsa_miR_449a	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_450a_5p	=	allRes_hsa_miR_450a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4516	=	allRes_hsa_miR_4516	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_451a	=	allRes_hsa_miR_451a	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_452_5p	=	allRes_hsa_miR_452_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4521	=	allRes_hsa_miR_4521	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4532	=	allRes_hsa_miR_4532	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_454_3p	=	allRes_hsa_miR_454_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_455_3p	=	allRes_hsa_miR_455_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_455_5p	=	allRes_hsa_miR_455_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4634	=	allRes_hsa_miR_4634	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4655_3p	=	allRes_hsa_miR_4655_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4695_3p	=	allRes_hsa_miR_4695_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4716_5p	=	allRes_hsa_miR_4716_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4730	=	allRes_hsa_miR_4730	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4731_3p	=	allRes_hsa_miR_4731_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4763_5p	=	allRes_hsa_miR_4763_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4770	=	allRes_hsa_miR_4770	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_4793_3p	=	allRes_hsa_miR_4793_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_483_3p	=	allRes_hsa_miR_483_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_486_3p	=	allRes_hsa_miR_486_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_486_5p	=	allRes_hsa_miR_486_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_489_3p	=	allRes_hsa_miR_489_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_490_3p	=	allRes_hsa_miR_490_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_491_3p	=	allRes_hsa_miR_491_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_497_3p	=	allRes_hsa_miR_497_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_497_5p	=	allRes_hsa_miR_497_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_5001_5p	=	allRes_hsa_miR_5001_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_500a_3p	=	allRes_hsa_miR_500a_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_500b_5p	=	allRes_hsa_miR_500b_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_501_5p	=	allRes_hsa_miR_501_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_502_3p	=	allRes_hsa_miR_502_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_504_3p	=	allRes_hsa_miR_504_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_505_5p	=	allRes_hsa_miR_505_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_511_3p	=	allRes_hsa_miR_511_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_516b_5p	=	allRes_hsa_miR_516b_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_517_5p	=	allRes_hsa_miR_517_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_517a_3p	=	allRes_hsa_miR_517a_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_517c_3p	=	allRes_hsa_miR_517c_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_521	=	allRes_hsa_miR_521	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_522_3p	=	allRes_hsa_miR_522_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_532_3p	=	allRes_hsa_miR_532_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_532_5p	=	allRes_hsa_miR_532_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_548aa	=	allRes_hsa_miR_548aa	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_548aw	=	allRes_hsa_miR_548aw	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_548b_3p	=	allRes_hsa_miR_548b_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_548c_3p	=	allRes_hsa_miR_548c_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_548f_3p	=	allRes_hsa_miR_548f_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_548q	=	allRes_hsa_miR_548q	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_548x_3p	=	allRes_hsa_miR_548x_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_551b_3p	=	allRes_hsa_miR_551b_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_5701	=	allRes_hsa_miR_5701	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_582_5p	=	allRes_hsa_miR_582_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_585_3p	=	allRes_hsa_miR_585_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_590_5p	=	allRes_hsa_miR_590_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_595	=	allRes_hsa_miR_595	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_598_3p	=	allRes_hsa_miR_598_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_6068	=	allRes_hsa_miR_6068	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_6073	=	allRes_hsa_miR_6073	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_6075	=	allRes_hsa_miR_6075	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_610	=	allRes_hsa_miR_610	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_624_5p	=	allRes_hsa_miR_624_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_628_5p	=	allRes_hsa_miR_628_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_642b_5p	=	allRes_hsa_miR_642b_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_6500_5p	=	allRes_hsa_miR_6500_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_6516_3p	=	allRes_hsa_miR_6516_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_652_3p	=	allRes_hsa_miR_652_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_653_3p	=	allRes_hsa_miR_653_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_660_5p	=	allRes_hsa_miR_660_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_664a_3p	=	allRes_hsa_miR_664a_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_664b_3p	=	allRes_hsa_miR_664b_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_6716_3p	=	allRes_hsa_miR_6716_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_6722_5p	=	allRes_hsa_miR_6722_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_6730_3p	=	allRes_hsa_miR_6730_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_6743_3p	=	allRes_hsa_miR_6743_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_6771_5p	=	allRes_hsa_miR_6771_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_6779_3p	=	allRes_hsa_miR_6779_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_6794_3p	=	allRes_hsa_miR_6794_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_6804_5p	=	allRes_hsa_miR_6804_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_6806_5p	=	allRes_hsa_miR_6806_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_6817_5p	=	allRes_hsa_miR_6817_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_6826_5p	=	allRes_hsa_miR_6826_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_6865_3p	=	allRes_hsa_miR_6865_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_6872_3p	=	allRes_hsa_miR_6872_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_6886_3p	=	allRes_hsa_miR_6886_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_6891_3p	=	allRes_hsa_miR_6891_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_6895_5p	=	allRes_hsa_miR_6895_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_7108_3p	=	allRes_hsa_miR_7108_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_7111_3p	=	allRes_hsa_miR_7111_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_7159_5p	=	allRes_hsa_miR_7159_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_744_5p	=	allRes_hsa_miR_744_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_770_5p	=	allRes_hsa_miR_770_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_7704	=	allRes_hsa_miR_7704	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_8063	=	allRes_hsa_miR_8063	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_8077	=	allRes_hsa_miR_8077	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_874_5p	=	allRes_hsa_miR_874_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_887_3p	=	allRes_hsa_miR_887_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_940	=	allRes_hsa_miR_940	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_95_3p	=	allRes_hsa_miR_95_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_98_5p	=	allRes_hsa_miR_98_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_99a_3p	=	allRes_hsa_miR_99a_3p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_99a_5p	=	allRes_hsa_miR_99a_5p	[,c('GO.ID','Term','KS')]
allRes_data_hsa_miR_99b_5p	=	allRes_hsa_miR_99b_5p	[,c('GO.ID','Term','KS')]

Top_10_miRNAs_ranked_GOdata = rbind.data.frame(allRes_data_hsa_let_7a_5p	,
                                               allRes_data_hsa_let_7b_5p	,
                                               allRes_data_hsa_let_7c_5p	,
                                               allRes_data_hsa_let_7d_5p	,
                                               allRes_data_hsa_let_7e_5p	,
                                               allRes_data_hsa_let_7f_5p	,
                                               allRes_data_hsa_let_7g_5p	,
                                               allRes_data_hsa_let_7i_5p ,
                                               allRes_data_hsa_miR_1_3p	,
                                               allRes_data_hsa_miR_100_5p	,
                                               allRes_data_hsa_miR_101_3p	,
                                               allRes_data_hsa_miR_101_5p	,
                                               allRes_data_hsa_miR_103a_3p	,
                                               allRes_data_hsa_miR_106b_3p	,
                                               allRes_data_hsa_miR_107	,
                                               allRes_data_hsa_miR_10a_5p	,
                                               allRes_data_hsa_miR_10b_5p	,
                                               allRes_data_hsa_miR_1227_5p	,
                                               allRes_data_hsa_miR_1238_5p	,
                                               allRes_data_hsa_miR_1247_5p	,
                                               allRes_data_hsa_miR_125a_5p	,
                                               allRes_data_hsa_miR_125b_5p	,
                                               allRes_data_hsa_miR_126_3p	,
                                               allRes_data_hsa_miR_126_5p	,
                                               allRes_data_hsa_miR_128_3p	,
                                               allRes_data_hsa_miR_129_2_3p	,
                                               allRes_data_hsa_miR_1296_5p	,
                                               allRes_data_hsa_miR_1306_3p	,
                                               allRes_data_hsa_miR_130a_3p	,
                                               allRes_data_hsa_miR_133a_3p	,
                                               allRes_data_hsa_miR_133a_5p	,
                                               allRes_data_hsa_miR_133b	,
                                               allRes_data_hsa_miR_135a_5p	,
                                               allRes_data_hsa_miR_138_5p	,
                                               allRes_data_hsa_miR_139_3p	,
                                               allRes_data_hsa_miR_139_5p	,
                                               allRes_data_hsa_miR_140_3p	,
                                               allRes_data_hsa_miR_140_5p	,
                                               allRes_data_hsa_miR_142_3p	,
                                               allRes_data_hsa_miR_142_5p	,
                                               allRes_data_hsa_miR_143_3p	,
                                               allRes_data_hsa_miR_143_5p	,
                                               allRes_data_hsa_miR_144_3p	,
                                               allRes_data_hsa_miR_144_5p	,
                                               allRes_data_hsa_miR_145_3p	,
                                               allRes_data_hsa_miR_145_5p	,
                                               allRes_data_hsa_miR_146a_5p	,
                                               allRes_data_hsa_miR_146b_5p	,
                                               allRes_data_hsa_miR_147b	,
                                               allRes_data_hsa_miR_150_5p	,
                                               allRes_data_hsa_miR_151a_5p	,
                                               allRes_data_hsa_miR_151b	,
                                               allRes_data_hsa_miR_152_3p	,
                                               allRes_data_hsa_miR_1537_3p	,
                                               allRes_data_hsa_miR_15a_5p	,
                                               allRes_data_hsa_miR_15b_3p	,
                                               allRes_data_hsa_miR_15b_5p	,
                                               allRes_data_hsa_miR_16_2_3p	,
                                               allRes_data_hsa_miR_16_5p	,
                                               allRes_data_hsa_miR_181a_2_3p	,
                                               allRes_data_hsa_miR_181a_3p	,
                                               allRes_data_hsa_miR_182_3p	,
                                               allRes_data_hsa_miR_183_3p	,
                                               allRes_data_hsa_miR_184	,
                                               allRes_data_hsa_miR_185_5p	,
                                               allRes_data_hsa_miR_186_5p	,
                                               allRes_data_hsa_miR_187_5p	,
                                               allRes_data_hsa_miR_18a_5p	,
                                               allRes_data_hsa_miR_18b_5p	,
                                               allRes_data_hsa_miR_190a_5p	,
                                               allRes_data_hsa_miR_191_5p	,
                                               allRes_data_hsa_miR_1913	,
                                               allRes_data_hsa_miR_193a_5p	,
                                               allRes_data_hsa_miR_195_5p	,
                                               allRes_data_hsa_miR_199a_3p	,
                                               allRes_data_hsa_miR_199a_5p	,
                                               allRes_data_hsa_miR_199b_5p	,
                                               allRes_data_hsa_miR_203a_3p	,
                                               allRes_data_hsa_miR_204_5p	,
                                               allRes_data_hsa_miR_205_3p	,
                                               allRes_data_hsa_miR_20a_5p	,
                                               allRes_data_hsa_miR_20b_5p	,
                                               allRes_data_hsa_miR_21_3p	,
                                               allRes_data_hsa_miR_214_3p	,
                                               allRes_data_hsa_miR_214_5p	,
                                               allRes_data_hsa_miR_218_5p	,
                                               allRes_data_hsa_miR_22_5p	,
                                               allRes_data_hsa_miR_221_3p	,
                                               allRes_data_hsa_miR_221_5p	,
                                               allRes_data_hsa_miR_222_3p	,
                                               allRes_data_hsa_miR_223_3p	,
                                               allRes_data_hsa_miR_223_5p	,
                                               allRes_data_hsa_miR_224_3p	,
                                               allRes_data_hsa_miR_2277_3p	,
                                               allRes_data_hsa_miR_23a_3p	,
                                               allRes_data_hsa_miR_23a_5p	,
                                               allRes_data_hsa_miR_23b_3p	,
                                               allRes_data_hsa_miR_24_3p	,
                                               allRes_data_hsa_miR_26a_1_3p	,
                                               allRes_data_hsa_miR_26a_5p	,
                                               allRes_data_hsa_miR_26b_3p	,
                                               allRes_data_hsa_miR_26b_5p	,
                                               allRes_data_hsa_miR_27a_3p	,
                                               allRes_data_hsa_miR_27a_5p	,
                                               allRes_data_hsa_miR_27b_3p	,
                                               allRes_data_hsa_miR_28_5p	,
                                               allRes_data_hsa_miR_2861	,
                                               allRes_data_hsa_miR_29a_3p	,
                                               allRes_data_hsa_miR_29a_5p	,
                                               allRes_data_hsa_miR_29b_1_5p	,
                                               allRes_data_hsa_miR_29b_2_5p	,
                                               allRes_data_hsa_miR_29b_3p	,
                                               allRes_data_hsa_miR_29c_3p	,
                                               allRes_data_hsa_miR_29c_5p	,
                                               allRes_data_hsa_miR_3065_3p	,
                                               allRes_data_hsa_miR_3065_5p	,
                                               allRes_data_hsa_miR_30a_3p	,
                                               allRes_data_hsa_miR_30a_5p	,
                                               allRes_data_hsa_miR_30b_5p	,
                                               allRes_data_hsa_miR_30c_2_3p	,
                                               allRes_data_hsa_miR_30c_5p	,
                                               allRes_data_hsa_miR_30d_5p	,
                                               allRes_data_hsa_miR_30e_3p	,
                                               allRes_data_hsa_miR_3149	,
                                               allRes_data_hsa_miR_3188	,
                                               allRes_data_hsa_miR_32_5p	,
                                               allRes_data_hsa_miR_324_3p	,
                                               allRes_data_hsa_miR_326	,
                                               allRes_data_hsa_miR_328_3p	,
                                               allRes_data_hsa_miR_331_3p	,
                                               allRes_data_hsa_miR_335_5p	,
                                               allRes_data_hsa_miR_338_3p	,
                                               allRes_data_hsa_miR_339_5p	,
                                               allRes_data_hsa_miR_340_5p	,
                                               allRes_data_hsa_miR_342_3p	,
                                               allRes_data_hsa_miR_342_5p	,
                                               allRes_data_hsa_miR_34a_3p	,
                                               allRes_data_hsa_miR_34a_5p	,
                                               allRes_data_hsa_miR_34b_5p	,
                                               allRes_data_hsa_miR_34c_5p	,
                                               allRes_data_hsa_miR_3607_3p	,
                                               allRes_data_hsa_miR_361_3p	,
                                               allRes_data_hsa_miR_362_3p	,
                                               allRes_data_hsa_miR_362_5p	,
                                               allRes_data_hsa_miR_3620_5p	,
                                               allRes_data_hsa_miR_363_3p	,
                                               allRes_data_hsa_miR_3659	,
                                               allRes_data_hsa_miR_365a_3p	,
                                               allRes_data_hsa_miR_3660	,
                                               allRes_data_hsa_miR_3663_5p	,
                                               allRes_data_hsa_miR_371b_5p	,
                                               allRes_data_hsa_miR_374a_5p	,
                                               allRes_data_hsa_miR_374b_5p	,
                                               allRes_data_hsa_miR_374c_5p	,
                                               allRes_data_hsa_miR_424_5p	,
                                               allRes_data_hsa_miR_4252	,
                                               allRes_data_hsa_miR_4254	,
                                               allRes_data_hsa_miR_4284	,
                                               allRes_data_hsa_miR_4290	,
                                               allRes_data_hsa_miR_4306	,
                                               allRes_data_hsa_miR_4317	,
                                               allRes_data_hsa_miR_4318	,
                                               allRes_data_hsa_miR_4324	,
                                               allRes_data_hsa_miR_4328	,
                                               allRes_data_hsa_miR_4440	,
                                               allRes_data_hsa_miR_4443	,
                                               allRes_data_hsa_miR_4481	,
                                               allRes_data_hsa_miR_449a	,
                                               allRes_data_hsa_miR_450a_5p	,
                                               allRes_data_hsa_miR_4516	,
                                               allRes_data_hsa_miR_451a	,
                                               allRes_data_hsa_miR_452_5p	,
                                               allRes_data_hsa_miR_4521	,
                                               allRes_data_hsa_miR_4532	,
                                               allRes_data_hsa_miR_454_3p	,
                                               allRes_data_hsa_miR_455_3p	,
                                               allRes_data_hsa_miR_455_5p	,
                                               
                                               allRes_data_hsa_miR_4655_3p	,
                                               allRes_data_hsa_miR_4695_3p	,
                                               allRes_data_hsa_miR_4716_5p	,
                                               allRes_data_hsa_miR_4730	,
                                               allRes_data_hsa_miR_4731_3p	,
                                               allRes_data_hsa_miR_4763_5p	,
                                               allRes_data_hsa_miR_4770	,
                                               allRes_data_hsa_miR_4793_3p	,
                                               allRes_data_hsa_miR_483_3p	,
                                               allRes_data_hsa_miR_486_3p	,
                                               allRes_data_hsa_miR_486_5p	,
                                               allRes_data_hsa_miR_489_3p	,
                                               allRes_data_hsa_miR_490_3p	,
                                               allRes_data_hsa_miR_491_3p	,
                                               allRes_data_hsa_miR_497_3p	,
                                               allRes_data_hsa_miR_497_5p	,
                                               allRes_data_hsa_miR_5001_5p	,
                                               allRes_data_hsa_miR_500a_3p	,
                                               allRes_data_hsa_miR_500b_5p	,
                                               allRes_data_hsa_miR_501_5p	,
                                               allRes_data_hsa_miR_502_3p	,
                                               allRes_data_hsa_miR_504_3p	,
                                               allRes_data_hsa_miR_505_5p	,
                                               allRes_data_hsa_miR_511_3p	,
                                               allRes_data_hsa_miR_516b_5p	,
                                               allRes_data_hsa_miR_517_5p	,
                                               allRes_data_hsa_miR_517a_3p	,
                                               allRes_data_hsa_miR_517c_3p	,
                                               allRes_data_hsa_miR_521	,
                                               allRes_data_hsa_miR_522_3p	,
                                               allRes_data_hsa_miR_532_3p	,
                                               allRes_data_hsa_miR_532_5p	,
                                               allRes_data_hsa_miR_548aa	,
                                               allRes_data_hsa_miR_548aw	,
                                               allRes_data_hsa_miR_548b_3p	,
                                               allRes_data_hsa_miR_548c_3p	,
                                               allRes_data_hsa_miR_548f_3p	,
                                               allRes_data_hsa_miR_548q	,
                                               allRes_data_hsa_miR_548x_3p	,
                                               allRes_data_hsa_miR_551b_3p	,
                                               allRes_data_hsa_miR_5701	,
                                               allRes_data_hsa_miR_582_5p	,
                                               allRes_data_hsa_miR_585_3p	,
                                               allRes_data_hsa_miR_590_5p	,
                                               allRes_data_hsa_miR_595	,
                                               allRes_data_hsa_miR_598_3p	,
                                               
                                               allRes_data_hsa_miR_6073	,
                                               allRes_data_hsa_miR_6075	,
                                               allRes_data_hsa_miR_610	,
                                               allRes_data_hsa_miR_624_5p	,
                                               allRes_data_hsa_miR_628_5p	,
                                               allRes_data_hsa_miR_642b_5p	,
                                               allRes_data_hsa_miR_6500_5p	,
                                               allRes_data_hsa_miR_6516_3p	,
                                               allRes_data_hsa_miR_652_3p	,
                                               allRes_data_hsa_miR_653_3p	,
                                               allRes_data_hsa_miR_660_5p	,
                                               allRes_data_hsa_miR_664a_3p	,
                                               allRes_data_hsa_miR_664b_3p	,
                                              
                                               allRes_data_hsa_miR_6722_5p	,
                                               allRes_data_hsa_miR_6730_3p	,
                                               allRes_data_hsa_miR_6743_3p	,
                                               
                                               allRes_data_hsa_miR_6779_3p	,
                                               allRes_data_hsa_miR_6794_3p	,
                                               allRes_data_hsa_miR_6804_5p	,
                                               allRes_data_hsa_miR_6806_5p	,
                                               allRes_data_hsa_miR_6817_5p	,
                                               allRes_data_hsa_miR_6826_5p	,
                                               allRes_data_hsa_miR_6865_3p	,
                                               allRes_data_hsa_miR_6872_3p	,
                                               allRes_data_hsa_miR_6886_3p	,
                                               allRes_data_hsa_miR_6891_3p	,
                                               allRes_data_hsa_miR_6895_5p	,
                                               allRes_data_hsa_miR_7108_3p	,
                                               allRes_data_hsa_miR_7111_3p	,
                                               allRes_data_hsa_miR_7159_5p	,
                                               allRes_data_hsa_miR_744_5p	,
                                               allRes_data_hsa_miR_770_5p	,
                                               allRes_data_hsa_miR_7704	,
                                               allRes_data_hsa_miR_8063	,
                                               allRes_data_hsa_miR_8077	,
                                               allRes_data_hsa_miR_874_5p	,
                                               allRes_data_hsa_miR_887_3p	,
                                               allRes_data_hsa_miR_940	,
                                               allRes_data_hsa_miR_95_3p	,
                                               allRes_data_hsa_miR_98_5p	,
                                               allRes_data_hsa_miR_99a_3p	,
                                               allRes_data_hsa_miR_99a_5p	,
                                               allRes_data_hsa_miR_99b_5p)
                                               
write.csv(Top_10_miRNAs_ranked_GOdata, file = 'Top_10_miRNAs_ranked_GOdata.csv')
