################################################################################
# Functions --------------------------------------------------------------------
#   mass convert ---------------------------------------------------------------
setGeneric(name = 'calculateExactMass',
           def = function(
    formula
           ){
             molecule <- Rdisop::getMolecule(formula)
             # getFormula(molecule)
             Rdisop::getMass(molecule)
           })


# exact_mass <- 180.0634
# adduct <- '[M-H]-'
# delta_mz <- -1.0073
# calculateMz(exact_mass = 180.0634,
#             adduct = '[M-H]-',
#             delta_mz = -1.0073)
#
# calculateMz(exact_mass = 180.0634,
#             adduct = '[2M-H]-',
#             delta_mz = -1.0073)

setGeneric(name = 'calculateMz',
           def = function(
    exact_mass,
    adduct,
    delta_mz
           ){

             if (stringr::str_detect(adduct, pattern = '2M')) {
               mz <- exact_mass*2 + delta_mz
             } else if (stringr::str_detect(adduct, pattern = '3M')) {
               mz <- exact_mass*3 + delta_mz
             } else {
               mz <- exact_mass + delta_mz
             }

             mz
           })


# exact_mass <- 180.0634
# adduct <- '[M-H]-'
# delta_mz <- -1.0073
# transformMz(exact_mass = 180.0634, type = 'adduct', polarity = 'positive')
# transformMz(exact_mass = 180.0634, type = 'nl', polarity = 'positive')
# transformMz(exact_mass = 180.0634, adduct_list = c('[M-H]-', '[M+H]+'))

setGeneric(name = 'transformMz',
           function(
    exact_mass,
    formula = NULL,
    adduct_list = NULL,
    type = c('adduct', 'nl'),
    polarity = c('positive', 'negative'),
    # lib_adduct_nl = lib_adduct_nl,
    ...
           ){

             if (all(is.null(exact_mass), is.null(formula))) {
               stop('Please input exact_mass or formula.')
             }

             if (!is.null(formula)) {
               exact_mass <- calculateExactMass(formula)
             }

             if (is.null(adduct_list)) {
               lib <- switch(polarity,
                             'positive' = {
                               lib_adduct_nl$positive
                             },
                             'negative' = {
                               lib_adduct_nl$negative
                             }
               )

               lib <- switch(type,
                             'adduct' = {
                               lib %>%
                                 dplyr::filter(type == 'Adduct') %>%
                                 dplyr::filter(credential == 'Yes')
                             },
                             'nl' = {
                               lib %>%
                                 dplyr::filter(type == 'NeutralLoss') %>%
                                 dplyr::filter(credential == 'Yes')
                             })
             } else {
               lib <- lib_adduct_nl$positive %>%
                 dplyr::bind_rows(lib_adduct_nl$negative)

               if (!all(adduct_list %in% lib$adduct)) {
                 stop('Sorry, not all adduct included in the adduct list\n')
               }

               lib <- lib %>%
                 dplyr::filter(adduct %in% adduct_list) %>%
                 dplyr::arrange(match(adduct, adduct_list))
             }


             result_mz <- sapply(seq_along(lib$adduct), function(i){
               calculateMz(exact_mass = exact_mass,
                           adduct = lib$adduct[i],
                           delta_mz = lib$delta_mz[i])
             })

             result <- tibble::tibble(exact_mass = exact_mass,
                                      adduct = lib$adduct,
                                      mz = result_mz)


             return(result)

           })




# convertMz2Adduct(base_mz = 181.0707,
#                  base_adduct = '[M+H]+',
#                  adduct_list = NULL,
#                  type = 'adduct',
#                  polarity = 'positive')

# convertMz2Adduct(base_mz = 181.0707,
#                  base_adduct = '[M+H]+',
#                  adduct_list = NULL,
#                  type = 'nl',
#                  polarity = 'negative')

# convertMz2Adduct(base_mz = 181.0707,
#                  base_adduct = '[M+H]+',
#                  adduct_list = c('[M-H2O+H]+', '[M+Na]+'),
#                  type = 'nl',
#                  polarity = 'negative')

# convertMz2Adduct(base_mz = 143.0347,
#                  base_adduct = '[2M-H]-',
#                  adduct_list = c('[M-H2O+H]+', '[M+Na]+'))

# convertMz2Adduct(base_mz = 143.0347,
#                  base_adduct = '[3M-H]-',
#                  adduct_list = c('[M-H2O+H]+', '[M+Na]+'))

setGeneric(name = 'convertMz2Adduct',
           def = function(
    base_mz,
    base_adduct,
    # lib_adduct_nl = lib_adduct_nl,
    ...
    # adduct_list = NULL,
    # type = c('adduct', 'nl'),
    # polarity = c('positive', 'negative')
           ){

             lib <- lib_adduct_nl$positive %>%
               dplyr::bind_rows(lib_adduct_nl$negative)

             if (!(base_adduct %in% lib$adduct)) {
               stop('Sorry, base_adduct is not included\n')
             }

             temp_delta_mz <- lib %>%
               dplyr::filter(adduct == base_adduct) %>%
               dplyr::pull(delta_mz)


             if (stringr::str_detect(base_adduct, pattern = '2M')) {
               temp_exact_mass <- (base_mz - temp_delta_mz)/2
             } else if (stringr::str_detect(base_adduct, pattern = '3M')) {
               temp_exact_mass <- (base_mz - temp_delta_mz)/3
             } else {
               temp_exact_mass <- base_mz - temp_delta_mz
             }

             # temp_exact_mass <- base_mz - temp_delta_mz

             result <- transformMz(exact_mass = temp_exact_mass,
                                   ...)

             return(result)

           })




################################################################################
# MetFrag ----------------------------------------------------------------------
#   generateMetFragLocalDB -----------------------------------------------------
#' @title generateMetFragLocalDB
#' @author Zhiwei Zhou
#' @description Generate local database for MetFrag
#' @param metfrag_candidate the metfrag candidates after random order
# #' @export



setGeneric(name = 'generateMetFragLocalDB',
           def = function(
    metfrag_candidate,
    base_dir='./02_metfrag',
    is_output = TRUE
           ){
             # browser()
             inchikeys <- splitInchiKey(inchikey = metfrag_candidate$inchikey)

             local_db <- data.frame(Identifier=metfrag_candidate$id_kegg,
                                    InChI=metfrag_candidate$inchi,
                                    MonoisotopicMass=metfrag_candidate$monoisopic_mass,
                                    MolecularFormula=metfrag_candidate$formula,
                                    InChIKey1=inchikeys$inchikey1,
                                    InChIKey2=inchikeys$inchikey2,
                                    InChIKey3=inchikeys$inchikey3,
                                    # InChI=inchikeys$inchikey,
                                    InChIKey=inchikeys$inchikey,
                                    SMILES=metfrag_candidate$smiles,
                                    Name=metfrag_candidate$name,
                                    stringsAsFactors = F)

             if(is_output) {
               readr::write_csv(local_db,
                                path = file.path(base_dir,
                                                 'local_db_metfrag.csv'))

             }

             # return(local_db)
           }
)


#   runMetFrag ----------------------------------------------------------------
# library(ReSOLUTION)

setGeneric(name = 'runMetFrag',
           def = function(
    mass,
    adduct_type,
    results_filename,
    base_dir,
    peaklist_path,
    DB = c("PubChem", "LocalCSV"),
    localDB_path = '',
    metfrag_path = 'F:/software/metfrag/MetFrag2.4.5-CL.jar',
    ppm=15,
    # s = as.numeric(peak_info$PRECURSORMZ),
    #adduct_type = p
    mzabs=0.003,
    frag_ppm=15
           ){
             # browser()
             if (DB=='PubChem') {
               localDB_path <- ''
             }

             if (DB=='LocalCSV' & localDB_path=='') {
               stop('Please input localDB_path')
             }

             config_file_test <- ReSOLUTION::MetFragConfig(mass = mass,
                                                           adduct_type = 0,
                                                           results_filename = results_filename,
                                                           base_dir = base_dir,
                                                           peaklist_path = peaklist_path,
                                                           DB=DB,
                                                           localDB_path=localDB_path,
                                                           # output="XLS",
                                                           output = 'CSV',
                                                           token="",
                                                           neutralPrecursorMass=TRUE,
                                                           ppm=ppm,
                                                           mzabs=mzabs,
                                                           frag_ppm=frag_ppm,
                                                           IsPosMode=TRUE,
                                                           tree_depth=2,
                                                           num_threads=4,
                                                           add_refs=FALSE,
                                                           minInt=0,
                                                           filter_isotopes=TRUE,
                                                           filter_by_InChIKey=FALSE,
                                                           useMoNAMetFusion = FALSE,
                                                           useMonaIndiv = FALSE)


             # metfrag_dir <- "E:/zhiwei_create/metfrag"
             # MetFragCL_name <- "MetFrag2.4.5-CL.jar"

             # runMetFrag(config_file_test, metfrag_dir, MetFragCL_name)
             # runMetFrag(config_file_test, metfrag_dir, MetFragCL_name)

             config_path <- file.path(base_dir,
                                      'config',
                                      paste0(results_filename,
                                             '_config.txt')
             )

             cmd <- paste0('java -Xmx1024m -jar ', metfrag_path, ' ',
                           config_path)

             shell(cmd)


           }
)





################################################################################
# CFMID ------------------------------------------------------------------------
#   runCFMID ------------------------------------------------------------------

#' @title runCFMID
#' @author Zhiwei Zhou
#' @description R interfaces to run CFM-ID.exe
#' @param cfm_id the path of cfm_id.exe Default: 'F:/MetIMMS_MetFrag/cfm_id/cfm-id.exe'
#' @param spec_file the path of spec_file txt file. 1st column: mz; 2nd column: intensity. It could be generated by GenerateCfmIdSpec
#' @param  candidate_file the path of  candidate_file txt file. 1st column: identifier; 2nd column: SMILES. It could be generated by GenerateCfmIdCandidates
#' @param id An identifier for the target molecule (Used to retrieve input spectrum from msp (if used). Default: 'temp'
#' @param num_highest The number of (ranked) candidates to return or -1 for all. Default: -1
#' @param ppm_mass_tol The mass tolerance in ppm to use when matching peaks within the dot product comparison - will use higher resulting tolerance of ppm and abs (if not given defaults to 10ppm). Default: 10
#' @param abs_mass_tol The mass tolerance in abs Da to use when matching peaks within the dot product comparison - will use higher resulting tolerance of ppm and abs ( if not given defaults to 0.01Da). Default: 0.01
#' @param prob_thresh The probability below which to prune unlikely fragmentations (default 0.001). Default: 0.001
#' @param param_file This file is the output of cfm-train. Pre-trained models as used in the above publication can be found in the supplementary data for that paper stored within the source tree of this project. Here, we use the se sfm model for both positive and negative modes. Positive mode: F:/MetIMMS_MetFrag/cfm_id/metab_se_cfm/param_output0.log; Negative mode: F:/MetIMMS_MetFrag/cfm_id/negative_metab_se_cfm/param_output0.log. Default:  'F:/MetIMMS_MetFrag/cfm_id/metab_se_cfm/param_output0.log'
#' @param config_file The filename where the configuration parameters of the cfm model can be found. This needs to match the file passed to cfm-train during training. See cfm-train documentation below for further details. Here, we use the se sfm model for both positive and negative modes. Positive mode: F:/MetIMMS_MetFrag/cfm_id/metab_se_cfm/param_config.txt; Negative mode: F:/MetIMMS_MetFrag/cfm_id/negative_metab_se_cfm/param_config.txt. Default:  'F:/MetIMMS_MetFrag/cfm_id/metab_se_cfm/param_config.txt'
#' @param score_type The type of scoring function to use when comparing spectra. Options: Jaccard (default), DotProduct. Default: Jaccard
#' @param apply_postprocessing Whether or not to post-process predicted spectra to take the top 80% of energy (at least 5 peaks), or the highest 30 peaks (whichever comes first) (0 = OFF (default for EI-MS), 1 = ON (default for ESI-MS/MS)). Default: 1
#' @param output_file the path of output result. It should be in ".txt" format. Default: 'F:/MetIMMS_MetFrag/20190521_retrieve_candidates/test/cfm_id_output/TEST_01_result.txt'
#' @param output_msp the path of output MSP file (MS/MS file). It should be in ".msp" format. Default: 'F:/MetIMMS_MetFrag/20190521_retrieve_candidates/test/cfm_id_output/TEST_01_spec.msp'
# #' @export
#' @examples
#' runCFMID(spec_file = 'F:/MetIMMS_MetFrag/20190521_retrieve_candidates/test/cfm_id_test_result/example_spec.txt',
#'          candidate_file = 'F:/MetIMMS_MetFrag/20190521_retrieve_candidates/test/cfm_id_test_result/example_candidates.txt',
#'          ppm_mass_tol = 10,
#'          abs_mass_tol = 0.01,
#'          score_type = 'Jaccard',
#'          output_file = 'F:/MetIMMS_MetFrag/20190521_retrieve_candidates/test/cfm_id_output/TEST_01_result.txt',
#'          output_msp = 'F:/MetIMMS_MetFrag/20190521_retrieve_candidates/test/cfm_id_output/TEST_01_spec.msp')

setGeneric(name = 'runCFMID',
           def = function(
    cfm_id = 'F:/MetIMMS_MetFrag/cfm_id/cfm-id.exe',
    spec_file = 'F:/MetIMMS_MetFrag/20190521_retrieve_candidates/test/cfm_id_test_result/example_spec.txt',
    id = 'temp',
    candidate_file = 'F:/MetIMMS_MetFrag/20190521_retrieve_candidates/test/cfm_id_test_result/example_candidates.txt',
    num_highest = -1,
    ppm_mass_tol = 10,
    abs_mass_tol = 0.01,
    prob_thresh = 0.001,
    param_file = 'F:/MetIMMS_MetFrag/cfm_id/metab_se_cfm/param_output0.log',
    config_file = 'F:/MetIMMS_MetFrag/cfm_id/metab_se_cfm/param_config.txt',
    score_type = c('Jaccard', 'DotProduct'),
    apply_postprocessing = 1,
    output_file = 'F:/MetIMMS_MetFrag/20190521_retrieve_candidates/test/cfm_id_output/TEST_01_result.txt',
    output_msp = 'F:/MetIMMS_MetFrag/20190521_retrieve_candidates/test/cfm_id_output/TEST_01_spec.msp'
           ){
             score_type <- match.arg(score_type)

             if (!stringr::str_detect(param_file, 'param_output0.log')) {stop('Please check param_file [param_output0.log]\n')}
             if (!stringr::str_detect(config_file, 'param_config.txt')) {stop('Please check config_file [param_config.txt]\n')}
             if (!stringr::str_detect(output_msp, '.msp')) {stop('Please check the file name of output msp\n')}


             cmdl <- paste(cfm_id, spec_file, id, candidate_file, num_highest, ppm_mass_tol,
                           abs_mass_tol, prob_thresh, param_file, config_file, score_type, apply_postprocessing,
                           output_file, output_msp, sep=' ')

             shell(cmdl, translate = TRUE)

             cat('CFM-ID has completed\n')
           }
)






################################################################################
# MSFINDER ---------------------------------------------------------------------
#   generateMsFinderFormulaDB --------------------------------------------------
setGeneric(name = 'generateMsFinderFormulaDB',
           def = function(
    raw_data,
    dir_path = '.',
    file_name = 'MsfinderFormulaDB-VS10.efd'
           ){
             idx <- match(raw_data$id, lib_meta$id)

             unique_formula <- unique(lib_meta$formula[idx])
             exact_mass <- sapply(unique_formula, ImmsTools::Calcu_EM)
             # exact_mass <- metfrag_candidate$MonoisotopicMass[idx]

             result <- data.frame(exact_mass=exact_mass,
                                  formula=unique_formula,
                                  pubchem_cid='N/A',
                                  records=1,
                                  HMDB='N/A',
                                  KNApSAcK='N/A',
                                  ChEBI='N/A',
                                  DrugBank='N/A',
                                  SMPDB='N/A',
                                  YMDB='N/A',
                                  T3DB='N/A',
                                  FoodDB='N/A',
                                  NANPDB='N/A',
                                  STOFF='N/A',
                                  BMDB='N/A',
                                  LipidMAPS='N/A',
                                  Urine='N/A',
                                  Saliva ='N/A',
                                  Feces='N/A',
                                  ECMDB='N/A',
                                  CSF='N/A',
                                  Serum='N/A',
                                  PubChem=TRUE,
                                  PlantCyc='N/A',
                                  UNPD='N/A',
                                  MINE='N/A',
                                  stringsAsFactors = F)

             colnames(result) <- c("Exact mass", "Formula", "PubChem CID", "Records", "HMDB", "KNApSAcK", "ChEBI",
                                   "DrugBank", "SMPDB", "YMDB", "T3DB", "FooDB", "NANPDB", "STOFF", "BMDB",
                                   "LipidMAPS", "Urine", "Saliva", "Feces", "ECMDB", "CSF", "Serum",  "PubChem",
                                   "PlantCyc", "UNPD", "MINE")

             readr::write_tsv(result,
                              path = file.path(dir_path, file_name),
                              col_names = TRUE)





             return(result)
           }
)


#   generateFormulaMass --------------------------------------------------------
setGeneric('generateFormulaMass',
           def = function(smiles) {
             molecule <- rcdk::parse.smiles(smiles)[[1]]
             rcdk::convert.implicit.to.explicit(molecule)
             formula <- rcdk::get.mol2formula(molecule, charge=0)
             formula_result <- data.frame(exact_mass=formula@mass,
                                          formula=formula@string,
                                          stringsAsFactors = F)

             return(formula_result)
           }
)


#   generateUserDefineDatabaseMsFinder -----------------------------------------
setGeneric('generateUserDefinedDatabaseMsFinder',
           def = function(
    candidate_list,
    is.output = TRUE,
    dir_path = '.',
    file_name = 'defined_DB.txt'
           ) {

             inchikey <- candidate_list$inchikey
             short_inchikey <- splitInchiKey(inchikey)$inchikey1

             final_result <- data.frame(title=candidate_list$name,
                                        inchikey=candidate_list$inchikey,
                                        short_inchikey=short_inchikey,
                                        pubchem_cid=0,
                                        # formula_mass_result,
                                        exact_mass = candidate_list$monoisopic_mass,
                                        formula = candidate_list$formula,
                                        smiles = candidate_list$smiles,
                                        database_id = candidate_list$id_kegg,
                                        stringsAsFactors = F)

             colnames(final_result) <- c('Title', 'InChIKey', 'Short InChIKey', 'PubChem CID',
                                         'Exact mass', 'Formula', 'SMILES', 'Database ID')


             if (is.output) {
               readr::write_tsv(final_result,
                                file = file.path(dir_path, file_name),
                                col_names = TRUE)
             }

           }
)

#   runMSFinderPredict ---------------------------------------------------------
setGeneric(name = 'runMsFinderPredict',
           def = function(
    bat_file = 'test.bat',
    msfinder = 'I:/software/Tools/MSFinder/MSFINDER_ver_3.24/MsfinderConsoleApp.exe',
    input_folder = 'F:/MetIMMS_MSFinder/test8',
    output_folder = 'F:/MetIMMS_MSFinder/test8/result',
    method_file = 'F:/MetIMMS_MSFinder/test8/MsfinderConsoleApp-Param.txt'
           ){

             cmdl <- paste(msfinder, 'predict', '-i', input_folder, '-o', output_folder, '-m', method_file, sep=' ')
             readr::write_lines(x = cmdl, path = bat_file, append = TRUE)

             cat('MSFinder commend has been writen, please run bat file\n')
           }
)


#   generateMsFinderConsoleAppPara ---------------------------------------------
setGeneric(name = 'generateMsFinderConsoleAppPara',
           def = function(
    file_para = 'MsfinderConsoleApp-Param-L0621.txt',
    template = NULL,
    UserDefinedDbFilePath = 'F:/MetIMMS_MSFinder/test8/L0621_db.txt'
           ){
             if (length(template) == 0) {
               template <- readr::read_lines(system.file("extdata", "MsfinderConsoleApp-Param.txt", package="MetDNA2InSilicoTool"))
             } else {
               template <- readr::read_lines(file = template)
             }

             template[54] <- paste0('UserDefinedDbFilePath=', UserDefinedDbFilePath)
             readr::write_lines(template, file = file_para)
           }
)
#   convertAdduct4MsFinder -----------------------------------------------------
setGeneric(name = 'convertAdduct4MsFinder',
           def = function(adduct,
                          source = 'MetDNA2'){

             if (source == 'MetDNA2') {
               result <- switch (adduct,
                                 'M+' = '[M]+',
                                 '[M+H]+' = '[M+H]+',
                                 '[M+NH4]+' = '[M+NH4]+',
                                 '[M+Na]+' = '[M+Na]+',
                                 '[M-H+2Na]+' = '[M+2Na-H]+',
                                 '[M+K]+' = '[M+K]+',
                                 '[M-H+2K]+' = '[M+2K-H]+',
                                 '[2M+H]+' = '[2M+H]+',
                                 '[2M+NH4]+' = '[2M+NH4]+',
                                 '[2M+Na]+' = '[2M+Na]+',
                                 '[2M+K]+' = '[2M+K]+',
                                 '[M-H2O+H]+' = '[M+H-H2O]+',
                                 '[M-2H+3Na]+' = '[M+3Na-2H]+',
                                 '[M-2H+3K]+' = '[M+3K-2H]+',
                                 '[M+CH3CN+H]+' = '[M+ACN+H]+',
                                 '[M+CH3CN+Na]+' = '[M+2ACN+H]+',
                                 '[M+CH3COO+2H]+' = '[M+CH3COO+2H]+',
                                 '[M+HCOO+2H]+' = '[M+H+HCOOH]+',
                                 '[M+HCOO+H+K]+' = '[M+K+HCOOH]+',
                                 '[M+HCOO+H+Na]+' = '[M+Na+HCOOH]+',
                                 'M-' = '[M]-',
                                 '[M-H]-' = '[M-H]-',
                                 '[M+Na-2H]-' = '[M+Na-2H]-',
                                 '[M+K-2H]-' = '[M+K-2H]-',
                                 '[M+NH4-2H]-' = '[M+NH4-2H]-',
                                 '[2M-H]-' = '[2M-H]-',
                                 '[M+CH3COO]-' = '[M+CH3COO]-',
                                 '[2M+Na-2H]-' = '[2M+Na-2H]-',
                                 '[M-H2O-H]-' = '[M-H2O-H]-'
               )
             }

             return(result)

           })


################################################################################
# splitInchiKey ----------------------------------------------------------------

#' @title splitInchiKey
#' @author Zhiwei Zhou
#' @description split inchikey to 3 parts
#' @return a data.frame with 4 columns
# #' @export
#' @example
#' SplitInchiKey(inchikey = 'VGONTNSXDCQUGY-RRKCRQDMSA-N')

setGeneric(name = 'splitInchiKey',
           def = function(
    inchikey = 'AEMRFAOFKBGASW-UHFFFAOYSA-N'
           ){
             result <- stringr::str_split_fixed(inchikey, pattern = '-', n = 3)
             result <- data.frame(result, stringsAsFactors = F)
             colnames(result) <- c('inchikey1', 'inchikey2', 'inchikey3')

             result <- data.frame(result,
                                  inchikey=inchikey,
                                  stringsAsFactors = F)

             return(result)
           }
)


  # generateMGF ----------------------------------------------------------------
#' @title generateMGF
#' @author Zhiwei Zhou
#' @param file_name Default: 'test.mgf'
#' @param title feature name
#' @param precusormz numeric. Default: NULL
#' @param charge numeric. Default: NULL
#' @param rt numeric. Default: NULL
#' @param spec matrix
# #' @export
#' @examples
#' test <- ImmsTools::readMSP(file = 'F:/01 MetIMMS/00 data processing/190920_allccs_metabolite_id_demenstration/topscience_msms_pos_20v_190913.msp',
#'                            mode = 'all')
#' GenerateMGF(file_name = 'test.mgf',
#' title = test[[1]]$info$NAME,
#' precusormz = test[[1]]$info$PRECURSORMZ,
#' mslevel = "2",
#' charge = '1+',
#' spec = test[[1]]$spec)

# test <- ImmsTools::readMSP(file = 'F:/01 MetIMMS/00 data processing/190920_allccs_metabolite_id_demenstration/topscience_msms_pos_20v_190913.msp',
#                            mode = 'all')

# generateMGF(file_name = 'test.mgf',
#             title = test[[1]]$info$NAME,
#             precusormz = test[[1]]$info$PRECURSORMZ,
#             mslevel = "2",
#             charge = '1+',
#             spec = test[[1]]$spec)

setGeneric(name = 'generateMGF',
           def = function(
    file_name = "./test.mgf",
    title = NULL,
    precusormz = NULL,
    mslevel = c('2', '1'),
    charge = NULL,
    rt = NULL,
    spec = NULL
           ){

             mslevel <- match.arg(mslevel)

             if (!stringr::str_detect(file_name, '.mgf')) {
               stop('The suffix of file name must be .mgf\n')
             }

             if (is.null(title)) {
               stop('Please input title\n')
             }

             if (is.null(precusormz)) {
               stop('Please input precusormz\n')
             }

             # if (is.null(rt)) {
             #   stop('Please input rt\n')
             # }

             if (is.null(spec)) {
               stop('Please input spec\n')
             }

             if (ncol(spec)!=2) {
               stop('Please check the format of spectrum. It should be in a dataframe, 1st: "mz", 2nd: "intensity" \n')
             }

             if (!all(colnames(spec)==c('mz', 'intensity'))) {
               stop('Please check the format of spectrum. It should be in a dataframe, 1st: "mz", 2nd: "intensity" \n')
             }

             # write into MGF
             file_result <- file(description = file_name, open = "a")
             cat('BEGIN IONS\n', file = file_result)
             cat('TITLE=', title, '\n', sep = '', file = file_result)
             cat('PEPMASS=', precusormz, '\n', sep = '', file = file_result)
             cat('MSLEVEL=', mslevel, '\n', sep = '', file = file_result)

             if (!is.null(charge)) {
               cat('CHARGE=', charge, '\n', sep = '', file = file_result)
             }

             if (!is.null(rt)) {
               cat('RTINSECONDS=', rt, '\n', sep = '', file = file_result)
             }

             for (i in 1:nrow(spec)) {
               cat(paste(as.numeric(round(spec[i,1], digits = 4)),
                         as.numeric(round(spec[i,2], digits = 2)),
                         collapse = ' '),
                   '\n', sep = '', file = file_result)
             }

             cat('END IONS\n\n', file = file_result)

             close(file_result)
           }
)



  # generateMSP ----------------------------------------------------------------
#' @title generateMSP
#' @author Zhiwei Zhou
#' @description Convert MS/MS spectra to MSP files
#' @param file_name Required. The file name of MSP. The suffix of file name must be ".msp".
#' @param cmp_name Required
#' @param precusormz Required
#' @param spec Required. It should be in a dataframe, 1st: "mz", 2nd: "intensity"
#' @param adduct Default: NULL
#' @param instrument_type Default: NULL
#' @param instrument Default: NULL
#' @param smiles Default: NULL
#' @param inchikey Default: NULL
#' @param inchikey1 Default: NULL
#' @param formula Default: NULL
#' @param polarity 'Positive' or 'Negative'
#' @param ce Default: NULL
#' @param rt Default: NULL
#' @param ccs Default: NULL
#' @param zhulib_id Default: NULL
#' @param kegg_id Default: NULL
#' @param hmdb_id Default: NULL
#' @param pubchem_id Default: NULL
#' @param links Default: ''
#' @param comment Default: ''
# #' @export
#' @example
#' setwd('F:/01 MetIMMS/00 data processing/190515 external validation msms data extraction/')
#'     GenerateMSP(file_name = 'zhumetlib_validation_pos_20v_190520.msp',
#'                 cmp_name = external_data_pos$compound_name[i],
#'                 precusormz = external_data_pos$mz[i],
#'                 adduct = external_data_pos$adducts[i],
#'                 instrument_type = 'LC-ESI-qTof',
#'                 instrument = 'LC-ESI-qTof',
#'                 smiles = external_data_pos$smiles[i],
#'                 inchikey = external_data_pos$inchikey[i],
#'                 inchikey1 = external_data_pos$inchikey1[i],
#'                 formula = external_data_pos$formula[i],
#'                 polarity = 'Positive',
#'                 ce = '20',
#'                 ccs = external_data_pos$CCS[i],
#'                 zhulib_id = external_data_pos$matched_zhulib_id[i],
#'                 pubchem_id = external_data_pos$pubchem_cid[i],
#'                 comment = paste(external_data_pos$id[i], external_data_pos$source[i], sep = ' '),
#'                 spec = temp_spec)



setGeneric(name = 'generateMSP',
           def = function(
    file_name = "./test.msp",
    cmp_name = NULL,
    precusormz = NULL,
    adduct=NULL,
    instrument_type=NULL,
    instrument=NULL,
    smiles=NULL,
    inchikey=NULL,
    inchikey1=NULL,
    formula=NULL,
    polarity=c('Positive', 'Negative'),
    ce=NULL,
    rt=NULL,
    ccs=NULL,
    zhulib_id=NULL,
    kegg_id=NULL,
    hmdb_id=NULL,
    pubchem_id=NULL,
    links='',
    comment='',
    spec=NULL
           ){

             if (!stringr::str_detect(file_name, '.msp')) {
               stop('The suffix of file name must be .msp\n')
             }

             if (is.null(cmp_name)) {
               stop('Please input cmp_name\n')
             }

             if (is.null(precusormz)) {
               stop('Please input precusormz\n')
             }

             if (is.null(spec)) {
               stop('Please input spec\n')
             }

             polarity = match.arg(polarity)


             if (ncol(spec)!=2) {
               stop('Please check the format of spectrum. It should be in a dataframe, 1st: "mz", 2nd: "intensity" \n')
             }

             if (!all(colnames(spec)==c('mz', 'intensity'))) {
               stop('Please check the format of spectrum. It should be in a dataframe, 1st: "mz", 2nd: "intensity" \n')
             }



             # write into msp
             file_result <- file(description = file_name, open = "a")
             cat('NAME: ', cmp_name, '\n', sep = '', file = file_result)
             cat('PRECURSORMZ: ', precusormz, '\n', sep = '', file = file_result)

             if (!is.null(adduct)) {
               cat('PRECURSORTYPE: ', adduct, '\n', sep = '', file = file_result)
             }

             if (!is.null(instrument_type)) {
               cat('INSTRUMENTTYPE: ', instrument_type, '\n', sep = '', file = file_result)
             }

             if (!is.null(instrument)) {
               cat('INSTRUMENT: ', instrument, '\n', sep = '', file = file_result)
             }

             if (!is.null(smiles)) {
               cat('SMILES: ', smiles, '\n', sep = '', file = file_result)
             }

             if (!is.null(inchikey)) {
               cat('INCHIKEY: ', inchikey, '\n', sep = '', file = file_result)
             }

             if (!is.null(inchikey1)) {
               cat('INCHIKEY1: ', inchikey1, '\n', sep = '', file = file_result)
             }

             if (!is.null(formula)) {
               cat('FORMULA: ', formula, '\n', sep = '', file = file_result)
             }

             if (!is.null(polarity)) {
               cat('IONMODE: ', polarity, '\n', sep = '', file = file_result)
             }

             if (!is.null(ce)) {
               cat('COLLISIONENERGY: ', ce, '\n', sep = '', file = file_result)
             }

             if (!is.null(rt)) {
               cat('RETENTIONTIME: ', rt, '\n', sep = '', file = file_result)
             }

             if (!is.null(ccs)) {
               cat('COLLISIONCROSSSECTION: ', ccs, '\n', sep = '', file = file_result)
             }

             if (!is.null(zhulib_id)) {
               cat('ZHULAB: ', zhulib_id, '\n', sep = '', file = file_result)
             }

             if (!is.null(kegg_id)) {
               cat('KEGG: ', kegg_id, '\n', sep = '', file = file_result)
             }

             if (!is.null(hmdb_id)) {
               cat('HMDB: ', hmdb_id, '\n', sep = '', file = file_result)
             }

             if (!is.null(pubchem_id)) {
               cat('PUBCHEM: ', pubchem_id, '\n', sep = '', file = file_result)
             }

             cat('Links: ', links, '\n', sep = '', file = file_result)
             cat('Comment: ', comment, "\n", sep = '', file = file_result)

             cat('Num Peaks: ',  nrow(spec),  '\n',  sep = '', file = file_result)

             for (i in 1:nrow(spec)) {
               cat(paste(as.numeric(round(spec[i,1], digits = 4)),
                         as.numeric(round(spec[i,2], digits = 2)),
                         collapse = ' '),
                   '\n', sep = '', file = file_result)
             }

             cat('\n', file = file_result)

             close(file_result)
           }
)
