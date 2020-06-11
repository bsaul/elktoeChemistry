
all_crosses <- 
  dt %>%
  filter(sublabel == "ri[gam]") %>%
  distinct(label, sublabel, element, species, signal) %>%
  purrr::lift_dl(tidyr::crossing)()



all_crosses %>%
  anti_join(
    dt %>% distinct(label, sublabel, element, species, signal) %>%
      filter(sublabel == "ri[gam]")
  ) %>%
  distinct(element)


# remove duplicates
dt %>%
  mutate(
    mtime = purrr::map_chr(sha, ~ file.info(sprintf("data/ri/%s.rds",.x))[["mtime"]])
  ) %>%
  group_by(label, sublabel, element, species, signal) %>%
  filter(n() > 1) %>%
  filter(mtime != max(mtime)) %>%
  pull(sha) %>%
  sprintf("data/ri/%s.rds", .) %>%
  unlink()


# dt %>%
#   filter(label == "A") %>%
#   distinct(label, sublabel, species, n_valves) %>%
#   arrange(species, sublabel)
# 
# dt %>%
#   filter(label == "C", sublabel == "ri[mom]") %>%
#   pull(sha) %>%
#   sprintf("data/ri/%s.rds", .) %>%
#   unlink()

