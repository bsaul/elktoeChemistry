#-----------------------------------------------------------------------------#
#   Title: QC of linking processes
#  Author: B Saul
#    Date: 2018-10-07
# Purpose: 
#-----------------------------------------------------------------------------#

library(ggplot2)

test1 <- create_wide_analysis_data(
  valve_data, 
 .reference_transition = "on_", 
 .zero_function = identity) %>%
  mutate(idref_method = "default")

test2 <- create_wide_analysis_data(
  valve_data, 
   .reference_transition = "ipx_", 
   .zero_function = ca_changepoint, 
   .zf_args = list(.use_up_to_row = 150, .method = "AMOC")) %>%
  mutate(idref_method = "ca_changepoint")

test3 <- create_wide_analysis_data(
  valve_data, 
   .reference_transition = "ipx_", 
   .zero_function = pbca_changepoint, 
   .zf_args = list(.use_up_to_row = 150, .method = "AMOC")) %>%
  mutate(idref_method = "pbca_changepoint")

test4 <- create_wide_analysis_data(
  valve_data, 
  .reference_transition = "ipx_", 
  .zero_function = pbca_maxpoint, 
  .zf_args = list(.use_up_to_row = 150)) %>%
  mutate(idref_method = "pbca_maxpoint")

tests <- bind_rows(test1, test2, test3, test4) %>%
  mutate(
    idref_method = factor(idref_method, levels = c("default", "ca_changepoint", "pbca_changepoint", "pbca_maxpoint"))
  )

valve_lengths <- valve_measurements %>%
  group_by(id, transect) %>%
  filter(sum(grepl("ipx_", layer_transition)) > 0, sum(grepl("opx_", layer_transition)) > 0) %>%
  summarise(
    on_opx_distance = distance[grepl("opx_", layer_transition)],
    on_ipx_distance = distance[grepl("ipx_", layer_transition)], 
    valve_length = on_opx_distance - on_ipx_distance
  ) 


tester <- function(data){
  data %>%
    group_by(idref_method, id, transect) %>%
    filter(sum(layer == "ipx") > 0, sum(layer == "opx") > 0) %>%
    summarise(
      chem_length = distance[min(which(layer == "opx"))] - distance[max(which(layer == "ipx"))]
    ) %>% 
    inner_join(valve_lengths, by = c("id", "transect")) %>%
    mutate(diff = chem_length - valve_length)
}

test_check <- tests %>% tester %>% dplyr::filter(diff > -250) 


### Plot checks ####


create_qc_plot <- function(.idt){
  hold <- tests %>%
    group_by(idref_method) %>%
    mutate(
      PbCa_ratio = log10(abs(Pb208_CPS)/abs(Ca43_CPS)), 
      idt = paste(id, transect, sep = "_"),
      rn  = row_number()
    ) %>%
    filter(idt == .idt) 
  
  vlines <- hold %>%
    group_by(idref_method) %>%
    summarise(
      start = rn[max(which(layer == "ipx"))],
      end   = rn[min(which(layer == "opx"))])
  
  hold %>%
    ggplot(., aes(x = rn, y = PbCa_ratio, group = id)) + 
    geom_vline(data = vlines, aes(xintercept = start), color = "blue") +
    geom_vline(data = vlines, aes(xintercept = end), color = "red") +
    geom_line(alpha = .5) +
    facet_grid(idref_method ~ .) + 
    ggtitle(.idt) + 
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text.y = element_text(angle = 0),
      axis.title.x= element_blank(),
      axis.text.x = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank()
    )
}

ids <- distinct(test_check, id, transect) %>% mutate(idt = paste(id, transect, sep = "_")) %>% pull(idt)
invisible(lapply(ids, function(id){
  p <- create_qc_plot(id) 
  ggsave(sprintf("programs/QC/%s.png", id), p, width = 5, height =  2, dpi = 150, units = "in")
}))


create_qc_plot("C472_2") 

test3 %>% filter(id == "A3", transect == "1") %>% pull(cpt)

library(pracma)
temp1 <- tests %>%
  filter(id == "C475", transect == "2", idref_method == "pbca_changepoint") %>%
  filter(row_number() <= 150) %>%
  mutate(
    rn = row_number(),
    PbCa_ratio = Pb208_CPS/Ca43_CPS
  ) 
temp1

zz <- max(findpeaks(temp1$PbCa_ratio, threshold = 1)[, 2])

zz <- nrow(temp1) - changepoint::cpt.mean(rev(cummin(temp1$PbCa_ratio)))@cpts[1]
zz <- changepoint::cpt.meanvar(cumsum(rev(abs(diff(c(rep(temp1$PbCa_ratio[1], 100), temp1$PbCa_ratio))))))@cpts[1] - 100
zz <- changepoint::cpt.var(diff(temp1$PbCa_ratio))@cpts[1]

zz
zz <- which(diff(temp1$PbCa_ratio) == min(diff(temp1$PbCa_ratio)))
zz
plot(diff(temp1$PbCa_ratio))

temp1 %>%
  ggplot(. , aes(x = rn, y = PbCa_ratio)) + 
  geom_line() + 
  geom_vline(
    xintercept = zz
  )



temp2 <- create_wide_analysis_data(
  valve_data %>% filter(id == "A3", transect == "1"), 
  .reference_transition = "ipx_", 
  .zero_function = pbca_changepoint, 
  .zf_args = list(.use_up_to_row = 150, .method = "AMOC")) 


%>%
  pull(cpt)

# # test1 %>%
#   # lengths %>%
# valve_measurements %>%
#   dplyr::filter(id == "C497", transect == "2") %>%
#   View()



hist(test1_check$diff, breaks = 30)
hist(test2_check$diff, breaks = 30)
hist(test3_check$diff, breaks = 30)
hist(test4_check$diff, breaks = 30)

mean(test1_check$diff)




