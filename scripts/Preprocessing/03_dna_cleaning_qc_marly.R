# ***********************************************************
# Title: DNA Cleanup QC - AMPure (Marly, Nov 2023)
# Purpose: Tidy and archive DNA QC before/after AMPure cleaning
# Author: Marly Erazo
# ***********************************************************

raw <- read_excel("data/raw/20241125_marly_extra_DNA_cleaning_roedeanst_DNA.xlsx",
                  sheet = "Plate_set-up DNA Plate16s r",
                  skip = 1)  # Skip the description row


# Rename the useful columns based on the real headers
colnames(raw)[c(1, 2, 3, 10, 11, 17, 18, 19, 20, 21, 22, 23)] <- c(
  "Sample_ID", "Well", "Qubit_ng_ul_pre", "Nanodrop_ng_ul_pre", "A260_pre",
  "Qubit_ng_ul_post", "Nanodrop_ng_ul_post", "Nanodrop_unit_post",
  "A260_post", "A280_post", "Ratio_260_280_post", "Ratio_260_230_post"
)

df_cleaning <- raw %>%
  transmute(
    Sample_ID           = Sample_ID,
    Well                = Well,
    Qubit_ng_ul_pre     = as.numeric(Qubit_ng_ul_pre),
    Nanodrop_ng_ul_pre  = Nanodrop_ng_ul_pre,
    A260_pre            = A260_pre,
    Qubit_ng_ul_post    = Nanodrop_unit_post,
    Nanodrop_ng_ul_post = A260_post,
    A260_post           = as.numeric(A260_post),
    A280_post           = as.numeric(Ratio_260_280_post),
    Ratio_260_280_post  = `260/280...24`,
    Ratio_260_230_post  = `260/230...25`
  )


write_csv(df_cleaning, "data/processed/dna_cleanup_marly_2023.csv")
message("✅ Cleaned DNA cleanup log saved.")

# Plot Qubit pre vs post
p1 <- ggplot(df_cleaning, aes(x = Qubit_ng_ul_pre, y = Qubit_ng_ul_post)) +
  geom_point() +
  labs(title = "Qubit: Before vs After Cleanup", x = "Before (ng/µl)", y = "After (ng/µl)") +
  theme_minimal()
ggsave("results/figures/dna_cleaning/qubit_pre_vs_post.png", p1, width = 7, height = 6)

p1
# Plot 260/280 vs 260/230
p2 <- ggplot(df_cleaning, aes(x = Ratio_260_280_post, y = Ratio_260_230_post)) +
  geom_point() +
  labs(title = "Post-Cleanup DNA Purity Ratios", x = "260/280", y = "260/230") +
  theme_minimal()
p2
ggsave("results/figures/dna_cleaning/post_cleanup_ratios.png", p2, width = 7, height = 6)

message("✅ DNA cleanup QC plots saved.")
