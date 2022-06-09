################################################################################
# copyFiles4InsilicoTool -------------------------------------------------------
#' @title copyFiles4InsilicoTool
#' @author Zhiwei Zhou
#' @param dir_path
#' @export

# dir_path <- 'G:/00_projects/03_MetDNA2/10_project/MetDNA2_project/Data/20220608_biological_samples/NIST_urine_pos/'

# dir_path <- 'D:/project/00_zhulab/01_metdna2/00_data/20220608_insilico_ms2_tutorial/MetDNA2_pos/'
# copyFiles4InsilicoTool(dir_path = dir_path)

copyFiles4InsilicoTool <- function(dir_path = '.'){
  cat('Copy files for in-silico MS/MS tool analysis\n')
  to_path <- file.path(dir_path, '07_insilico_msms')
  dir.create(path = to_path, showWarnings = FALSE, recursive = TRUE)

  file.copy(from = file.path(dir_path, '00_annotation_table/00_intermediate_data/table_identification'),
            to = to_path, overwrite = TRUE, recursive = TRUE, copy.date = TRUE)

  file.copy(from = file.path(dir_path, '03_annotation_credential/ms2_data.RData'),
            to = to_path, overwrite = TRUE, recursive = TRUE, copy.date = TRUE)

  file.copy(from = file.path(dir_path, '03_annotation_credential/ms2_data.msp'),
            to = to_path, overwrite = TRUE, recursive = TRUE, copy.date = TRUE)

}


# reformatTable ----------------------------------------------------------------
#' @title reformatTable1
#' @description
#' @author Zhiwei Zhou
#' @export
#' @importFrom magrittr '%>%' '%$%'

# dir_path <- 'd:/project/00_zhulab/01_metdna2/00_data/20220608_insilico_ms2_tutorial/MetDNA2_pos/07_insilico_msms/'
# dir_path <- 'G:/00_projects/03_MetDNA2/10_project/MetDNA2_project/Data/20220608_biological_samples/NIST_urine_pos/07_insilico_msms/'
# reformatTable1(dir_path = dir_path)
reformatTable1 <- function(dir_path = '.'){
  load(file.path(dir_path, 'table_identification'))
  load(file.path(dir_path, 'ms2_data.RData'))
  table_identification <- table_identification %>%
    # mutate(peak_name = paste(peak_name, 'POS', sep = '_')) %>%
    dplyr::mutate(with_ms2 = case_when(peak_name %in% names(raw_msms) ~ TRUE,
                                       !(peak_name %in% names(raw_msms)) ~ FALSE)) %>%
    dplyr::mutate(id_kegg = case_when(is.na(id_kegg) ~ id_zhulab,
                                      !is.na(id_kegg) ~ id_kegg)) %>%
    dplyr::mutate(label_seed = paste(peak_name, id_kegg, sep = '_'))

  # add inchi and monoisotopic mass
  data("cpd_emrn", envir = environment())
  data("cpd_lib", envir = environment())

  cpd_emrn <- cpd_emrn %>% tidyr::separate_rows(id_kegg_synonyms, sep = ';')

  idx1 <- match(table_identification$id_kegg, cpd_emrn$id)
  mono_mass_vector <- cpd_emrn$monoisotopic_mass[idx1]
  inchi_vector <- cpd_emrn$inchi[idx1]

  idx_na <- which(is.na(mono_mass_vector))
  idx2 <- match(table_identification$id_kegg[idx_na], cpd_emrn$id_kegg_synonyms)
  mono_mass_vector[idx_na] <- cpd_emrn$monoisotopic_mass[idx2]
  inchi_vector[idx_na] <- cpd_emrn$inchi[idx2]

  idx_na <- which(is.na(mono_mass_vector))
  idx3 <- match(table_identification$id_kegg[idx_na], cpd_lib$id)
  mono_mass_vector[idx_na] <- cpd_lib$monoisotopic_mass[idx3]
  inchi_vector[idx_na] <- cpd_lib$inchi[idx3]

  idx_na <- which(is.na(mono_mass_vector))
  idx4 <- match(table_identification$id_kegg[idx_na], cpd_lib$id_kegg)
  mono_mass_vector[idx_na] <- cpd_lib$monoisotopic_mass[idx4]
  inchi_vector[idx_na] <- cpd_lib$inchi[idx4]

  idx_na <- which(is.na(mono_mass_vector))
  idx5 <- match(table_identification$id_kegg[idx_na], cpd_emrn$id_kegg)
  mono_mass_vector[idx_na] <- cpd_emrn$monoisotopic_mass[idx5]
  inchi_vector[idx_na] <- cpd_lib$inchi[idx5]

  idx_na <- which(is.na(mono_mass_vector))
  idx6 <- match(table_identification$id_zhulab[idx_na], cpd_lib$id)
  mono_mass_vector[idx_na] <- cpd_emrn$monoisotopic_mass[idx6]
  inchi_vector[idx_na] <- cpd_lib$inchi[idx6]

  table_identification <- table_identification %>%
    dplyr::mutate(monoisopic_mass = mono_mass_vector, inchi = inchi_vector)

  rm(list = c('cpd_emrn', 'cpd_lib'));gc()

  save(table_identification, file = file.path(dir_path, 'table_identification'))
}


################################################################################
# generateFiles4InsilicoMsMs ---------------------------------------------------


#' @title generateFiles4InsilicoMsMs
#' @author Zhiwei Zhou
#' @param peak_id
#' @param dir_path
#' @export
#' @examples
#' generateFiles4InsilicoMsMs
#' generateFiles4InsilicoMsMs(peak_id = 'M196T420', dir_path = 'G:/00_projects/03_MetDNA2/10_project/MetDNA2_project/Data/20220608_biological_samples/NIST_urine_pos/07_insilico_msms/')

# generateFiles4InsilicoMsMs(peak_id = 'M196T420',
#                            dir_path = 'G:/00_projects/03_MetDNA2/10_project/MetDNA2_project/Data/20220608_biological_samples/NIST_urine_pos/07_insilico_msms/')

# load('G:/00_projects/03_MetDNA2/10_project/MetDNA2_project/Data/20220608_biological_samples/NIST_urine_pos/07_insilico_msms/table_identification')
# load('g:/00_projects/03_MetDNA2/02_database//PackageFile/library/20210629_update_emrn/cpd_emrn_210629.RData')
# load('G:/00_projects/03_MetDNA2/10_project/MetDNA2_project/Data/20220608_biological_samples/NIST_urine_neg/03_annotation_credential/ms2_data.RData')

# load('D:/project/00_zhulab/01_metdna2/00_data/20220608_insilico_ms2_tutorial/MetDNA2_pos/07_insilico_msms/table_identification')
# example M196T420

# dir_path <- 'G:/00_projects/03_MetDNA2/10_project/MetDNA2_project/Data/20220608_biological_samples/NIST_urine_pos/07_insilico_msms/'
# peak_id <- 'M196T420'

generateFiles4InsilicoMsMs <- function(peak_id,
  dir_path = '.') {

  cat('Generate files for in-silico MS/MS tools:', peak_id, '\n')
  load(file.path(dir_path, 'table_identification'))
  load(file.path(dir_path, 'ms2_data.RData'))

  peaks_with_ms2 <- table_identification %>%
    dplyr::filter(with_ms2)

  if (!(peak_id %in% peaks_with_ms2$peak_name)) {
    stop('The input peak_id has no annotation or MS2\n')
  }

  # peak_info
  candidate_list <- peaks_with_ms2 %>%
    dplyr::filter(peak_name == peak_id) %>%
    dplyr::select(peak_name, mz, rt, id_kegg, name, formula, monoisopic_mass,
      smiles, inchi, inchikey, adduct, isotope)

  # ms/ms spec
  ms2 <- match(peak_id, names(raw_msms)) %>% raw_msms[[.]]

  dir.create(file.path(dir_path, peak_id), showWarnings = FALSE, recursive = TRUE)
  save(candidate_list, file = file.path(dir_path, peak_id, 'candidate_list'))
  save(ms2, file = file.path(dir_path, peak_id, 'ms2'))

  readr::write_csv(candidate_list,
                   file = file.path(dir_path, peak_id, 'candidate_list.csv'))
  generateMGF(file_name = file.path(dir_path, peak_id, 'ms2.mgf'),
              title = ms2$info$NAME,
              precusormz = ms2$info$PRECURSORMZ,
              charge = ifelse(ms2$info$IONMODE == 'positive', '1+', '1-'),
              spec = ms2$spec)

  cat('Done!\n')
}




################################################################################
# runMetFragMatch --------------------------------------------------------------

#' @title runMetFragMatch
#' @author Zhiwei Zhou
#' @description this function calls ReSOLUTION package. Please refer ReSOLUTION and MetFrag
#' @param peak_id
#' @param dir_path default: '.'
#' @param ppm relative precursor tolerance. Default: 25 ppm
#' @param mzabs absolute precursor tolerance. Default: 0.01
#' @param frag_ppm relative fragment tolerance. Default: 25 ppm
#' @export
#' @examples
#'

# runMetFragMatch(peak_id = 'M196T420',
#                 dir_path = 'G:/00_projects/03_MetDNA2/10_project/MetDNA2_project/Data/20220608_biological_samples/NIST_urine_pos/07_insilico_msms/',
#                 metfrag_path = 'F:/software/metfrag/MetFrag2.4.5-CL.jar',
#                 ppm = 25,
#                 mzabs = 0.01,
#                 frag_ppm = 25)
# dir_path <- 'G:/00_projects/03_MetDNA2/10_project/MetDNA2_project/Data/20220608_biological_samples/NIST_urine_pos/07_insilico_msms/'
# peak_id <- 'M196T420'

runMetFragMatch <- function(peak_id,
                            dir_path = '.',
                            metfrag_path = 'F:/software/metfrag/MetFrag2.4.5-CL.jar',
                            ppm = 25,
                            mzabs = 0.01,
                            frag_ppm = 25) {

  if (!(peak_id %in% list.files(dir_path))) {
    stop('Please call generateFiles4InsilicoMsMs for peak', peak_id)
  }

  cat('Start run MetFrag Match for KGMN candidates\n')

  load(file.path(dir_path, peak_id, 'candidate_list'))
  load(file.path(dir_path, peak_id, 'ms2'))

  root_path <- file.path(dir_path, peak_id, '01_metfrag')

  # preparing data ---------------------------------------------
  peak_list <- as.data.frame(ms2$spec)
  dir.create(root_path, recursive = T, showWarnings = FALSE)
  spec_file_path  <- file.path(root_path,
                               "peak_list.txt")
  readr::write_delim(x = peak_list,
                     path = spec_file_path,
                     col_names = F)

  cat('Generate local DB...\n\n')
  generateMetFragLocalDB(metfrag_candidate = candidate_list,
                         base_dir = root_path)

  temp_adduct_list <- candidate_list$adduct %>% unique()
  temp_polarity <- ms2$info$IONMODE

  for (adduct_form in temp_adduct_list) {
    cat('Start processing', adduct_form, '\n\n')
    temp_mass <- candidate_list$mz[1]
    temp_mass <- convertMz2Adduct(base_mz = temp_mass,
                                  base_adduct = adduct_form,
                                  adduct_list = ifelse(temp_polarity == 'positive', '[M+H]+', '[M-H]-'))
    temp_mass <- temp_mass$exact_mass[1]

    temp_path <- file.path(root_path, paste0(adduct_form))
    dir.create(temp_path, showWarnings = FALSE, recursive = TRUE)

    # match with locol db ------------------------------------------------------
    runMetFrag(mass = temp_mass,
               adduct_type = 0,
               results_filename = 'metfrag_rank',
               base_dir = temp_path,
               metfrag_path = metfrag_path,
               peaklist_path = file.path(root_path, 'peak_list.txt'),
               DB = 'LocalCSV',
               localDB_path = file.path(root_path, 'local_db_metfrag.csv'),
               ppm = ppm,
               mzabs = mzabs,
               frag_ppm = frag_ppm)
  }

}



################################################################################
# runCfmIdMatch ----------------------------------------------------------------

#' @title runCfmIdMatch
#' @author Zhiwei Zhou
#' @param peak_id
#' @param cfmid_path the path of cfmid program. Default: 'F:/software/cfm_id/cfm-id.exe'
#' @param config_file the path of config file
#' @param param_file the path of parameter file
#' @param score_type score type in CFMID. Default: 'Jaccard'
#' @param ppm relative precursor tolerance. Default: 25 ppm
#' @param mzabs absolute precursor tolerance. Default: 0.01
#' @export
#' @examples
#'

# runCfmIdMatch(peak_id = 'M196T420',
#               dir_path = 'G:/00_projects/03_MetDNA2/10_project/MetDNA2_project/Data/20220608_biological_samples/NIST_urine_pos/07_insilico_msms/',
#               cfmid_path = 'F:/software/cfm_id/cfm-id.exe',
#               config_file = 'F:/software/cfm_id/metab_se_cfm/param_config.txt',
#               param_file = 'F:/software/cfm_id/metab_se_cfm/param_output0.log',
#               score_type = 'Jaccard',
#               ppm = 25,
#               mzabs = 0.01)

# dir_path <- 'G:/00_projects/03_MetDNA2/10_project/MetDNA2_project/Data/20220608_biological_samples/NIST_urine_pos/07_insilico_msms/'
# peak_id <- 'M196T420'

runCfmIdMatch <- function(peak_id,
                          dir_path = '.',
                          cfmid_path = 'F:/software/cfm_id/cfm-id.exe',
                          config_file = 'F:/software/cfm_id/metab_se_cfm/param_config.txt',
                          param_file = 'F:/software/cfm_id/metab_se_cfm/param_output0.log',
                          score_type = 'Jaccard',
                          ppm = 25,
                          mzabs = 0.01) {

  if (!(peak_id %in% list.files(dir_path))) {
    stop('Please call generateFiles4InsilicoMsMs for peak', peak_id)
  }

  cat('Start run CFM-ID Match for KGMN candidates\n')
  load(file.path(dir_path, peak_id, 'candidate_list'))
  load(file.path(dir_path, peak_id, 'ms2'))

  root_path <- file.path(dir_path, peak_id, '02_cfmid')

  # modify input data
  peak_list <- as.data.frame(ms2$spec)
  dir.create(root_path, recursive = T, showWarnings = FALSE)
  spec_file_path  <- file.path(root_path, "peak_list.txt")
  readr::write_delim(x = peak_list,
                     path = spec_file_path,
                     col_names = F)

  candidate_list <- candidate_list %>%
    dplyr::select(inchikey, smiles)
  dir.create(root_path, recursive = T, showWarnings = FALSE)
  candidate_path <- file.path(root_path, 'candidate_list.txt')
  readr::write_delim(x = candidate_list,
                     path = candidate_path,
                     col_names = F)

  output_file_path <- file.path(root_path, 'cfmid_result.txt')
  output_msp_path <- file.path(root_path, 'cfmid_pred_spec.msp')


  runCFMID(cfm_id = cfmid_path,
           spec_file = spec_file_path,
           candidate_file = candidate_path,
           id = peak_id,
           ppm_mass_tol = ppm,
           abs_mass_tol = mzabs,
           score_type = score_type,
           config_file = config_file,
           param_file = param_file,
           output_file = output_file_path,
           output_msp = output_msp_path)

}



################################################################################
# runMsFinderMatch -------------------------------------------------------------

#' @title runMsFinderMatch
#' @author Zhiwei Zhou
#' @description R interface to run MS-FINDER. If you use this function, please refer MS-FINDER parper.
#' @param peak_id peak name
#' @param dir_path path of working
#' @param msfinder_path msfinder program path
#' @export
#' @example
#' runMsFinderMatch(peak_id = 'M196T420', dir_path = 'G:/00_projects/03_MetDNA2/10_project/MetDNA2_project/Data/20220608_biological_samples/NIST_urine_pos/07_insilico_msms/', msfinder_path = 'F:/software/MSFINDER/MSFINDER_ver_3.24/MsfinderConsoleApp.exe')

# runMsFinderMatch(peak_id = 'M196T420',
#                  dir_path = 'G:/00_projects/03_MetDNA2/10_project/MetDNA2_project/Data/20220608_biological_samples/NIST_urine_pos/07_insilico_msms/',
#                  msfinder_path = 'F:/software/MSFINDER/MSFINDER_ver_3.24/MsfinderConsoleApp.exe')

runMsFinderMatch <- function(peak_id,
                             dir_path = '.',
                             msfinder_path = 'F:/software/MSFINDER/MSFINDER_ver_3.24/MsfinderConsoleApp.exe') {
  if (!(peak_id %in% list.files(dir_path))) {
    stop('Please call generateFiles4InsilicoMsMs for peak', peak_id)
  }

  cat('Start run CFM-ID Match for KGMN candidates\n')
  load(file.path(dir_path, peak_id, 'candidate_list'))
  load(file.path(dir_path, peak_id, 'ms2'))

  root_path <- file.path(dir_path, peak_id, '03_msfinder')

  # modify input format for MS-Finder
  for (adduct_form in unique(candidate_list$adduct)) {
    cat('Start processing', adduct_form, '\n\n')

    temp_path <- file.path(root_path, paste0(adduct_form))
    dir.create(temp_path, recursive = TRUE, showWarnings = FALSE)

    temp_adduct_form_converted <- convertAdduct4MsFinder(adduct = adduct_form)
    generateMSP(file_name = file.path(temp_path,
                                      paste0(adduct_form, '.msp')),
                cmp_name = peak_id,
                precusormz = as.numeric(ms2$info$PRECURSORMZ),
                adduct = temp_adduct_form_converted,
                polarity = ifelse(ms2$info$IONMODE == 'positive',
                                  'Positive', 'Negative'),
                spec = ms2$spec)

    cat('Generate user defined db\n')
    temp <- generateUserDefinedDatabaseMsFinder(candidate_list = candidate_list,
                                                dir_path = temp_path,
                                                file_name = paste(peak_id, 'db.txt', sep = '_'))
    cat('Generate consoleApp param\n')
    generateMsFinderConsoleAppPara(file_para = file.path(temp_path,
                                                         paste0('MsfinderConsoleApp-Param-', peak_id, '.txt')),
                                   # template = 'D:/01_r_package/MetDNA2InSilicoTool/inst/extdata/MsfinderConsoleApp-Param.txt',
                                   UserDefinedDbFilePath = file.path(temp_path,
                                                                     paste(peak_id, 'db.txt', sep = '_')))


    runMsFinderPredict(bat_file = file.path(root_path, paste0(peak_id, '_script.bat')),
                       msfinder = msfinder_path,
                       input_folder = temp_path,
                       output_folder = file.path(temp_path, 'result'),
                       method_file = file.path(temp_path,
                                               paste0('MsfinderConsoleApp-Param-', peak_id, '.txt')))

  }

  cat('\n');cat('Exclute bat file of MSFINDER\n')
  shell.exec(file = file.path(root_path, paste0(peak_id, '_script.bat')))
  cat('Done!\n')
}



################################################################################
# startup massage --------------------------------------------------------------
.onAttach <- function(libname, pkgname){
  packageStartupMessage("
If you have any questions, please send email to zhouzw@sioc.ac.cn or jiangzhu@sioc.ac.cn.
Authors: Zhiwei Zhou and Dr. Zhengjiang Zhu (jiangzhu@sioc.ac.cn).
Maintainer: Zhiwei Zhou

Version 0.1.0 (20220609)
-------------
o initial commit. An R package to connect KGMN result with other in-silico MS/MS tools.

")
}
