library(dplyr)

rm(list = ls())
gc()

load("cdp.rds")
physio = read.csv("2015_aug_dec_2016_dec.csv") %>%
  select(id,
         sex,
         age,
         hct_vena_mean,
         cms_score,
         spo2_mean,
         hr_mean,
         sbp_mean,
         dbp_mean,
         co_ppm,
         glucose,
         insulin,
         cholesterol,
         hdl,
         ldl,
         triglycerides,
         ferritin,
         iron,
         transferrin,
         total_testosterone,
         free_testosterone,
         height,
         hvr_corrected ,
         h_hvr_corrected_fixed,
         h_hvr_corrected,
         hcvr_corrected,
         hrr_hypox,
         hrr_co2,
         tv_ra,
         f_ra,
         etco2_ra,
         sao2_ra,
         hr_ra,
         peco2_ra,
         tv_21c,
         f_21c,
         etco2_21c,
         fio2_21c,
         sao2_21c,
         hr_21c,
         vibtps_kg_21c,
         peco2_21c,
         peco2_fix_21c,
         prdi_b_sleep,
         pahi_b_sleep,
         odi_b_sleep,
         meansatuation_b_sleep,
         minsaturation_b_sleep,
         maxsaturation_b_sleep,
         totalnumberofdesaturations_b_sleep,
         satbelow85_b_sleep,
         satbelow80_b_sleep,
         sleepefficiency_b_sleep,
         numberofwakes_b_sleep,
         vt_sleep,
         frequency_sleep,
         sao2_sleep,
         hr_bpm_sleep,
         vi_lmin_sleep,
         vibtps_sleep,
         vibtpskg_sleep,
         peco2_mmhg_sleep,
         hvr_sleep,
         hypercapnic_hvr_sleep,
         hcvr_sleep,
         hvr_corrected__sleep,
         h_hvr_corrected_sleep,
         hcvr_corrected_sleep,
         hrr_hypox_sleep,
         hrr_co2_sleep,
         sleepefficiency_sleep,
         totalahi_sleep,
         totalari_sleep,
         nadir_desat_sleep_sleep,
         spo2below85p_sleep_prct_sleep,
         spo2below80p_sleep_prct_sleep,
         spo2mean_sleep_sleep,
         odi_sleep_sleep,
         abpm_sbp_24h,
         abpm_dbp_24h,
         abpm_map_24h,
         abpm_sbp_vig,
         abpm_dbp_vig,
         abpm_map_vig,
         abpm_sbp_sleep,
         abpm_dbp_sleep,
         abpm_map_sleep,
         abpm_conventional_sbp,
         abpm_conventional_dbp
  )
wgs = read.csv("2018_wgs.csv")

suppressWarnings({
  wgs = cdp[,1,drop=FALSE] %>%
    left_join(wgs,by="wgs_id") %>%
    select(wgs_id,id) %>%
    left_join(physio,by="id") %>%
    cbind.data.frame(cdp[,-1])
})

rm(cdp)
wgsm = wgs[which(wgs$sex=="M"),]
wgsf = wgs[which(wgs$sex=="F"),]
wgsm_physio = wgsm[,1:90] %>%
  mutate_if(is.factor,as.numeric)
wgsf_physio = wgsf[,1:90] %>%
  mutate_if(is.factor,as.numeric)
wgsm_geno = wgsm[,-(1:90)]
wgsf_geno = wgsf[,-(1:90)]
rm(wgs,wgsm,wgsf)

save(wgsm_physio,wgsf_physio,wgsm_geno,wgsf_geno,file = "wgs_concise_physio.RData")
