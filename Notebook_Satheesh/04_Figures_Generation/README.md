workdir: `/Users/vsatheesh/Documents/ISU/GIF/gif_projects/2025/04_PEdwards_Methylation`




Loading BSseq object (this may take a minute)...
BSseq object loaded: 50971214 CpG sites x 25 samples
Sample names: P (2017) | T (2018) | AP (2019) | M (2017) | AZ (2019) | L (2017) | AA (2018) | AU (2019) | U (2018) | AR (2019) | AM (2019) | G (2017) | AO (2019) | W (2018) | AT (2019) | J (2017) | I (2017) | AV (2019) | O (2017) | A (2017) | AC (2018) | S (2018) | X (2018) | Y (2018) | AL (2019)

--- Extracting individual CpG sites (Fig 2A) ---

--- Extracting DMR region means (Figs 2B, 2C, 3A, 3B) ---
  Cyp1a1_promoter_52bp — 5 CpGs
  Cyp2j13_promoter_305bp — 4 CpGs
  Mif_promoter_73bp — 7 CpGs
  Sgk1_exon_53bp — 5 CpGs
  Sgk1_intergenic_66bp — 7 CpGs
  Igf1r_intron_72bp — 7 CpGs
  Igf1r_intron_95bp — 7 CpGs

Saved: methylation_for_figures.csv  (225 rows)

=== Summary of extracted loci ===
 figure    gene                  locus             comparison_sig n_samples
  Fig2A  Atp5g3       Atp5g3_909bp_TSS               Peak-Decline        24
  Fig2A  Atp5g3       Atp5g3_910bp_TSS Peak-Decline + Low-Decline        22
  Fig2B  Cyp1a1   Cyp1a1_promoter_52bp                Low-Decline        17
  Fig2B Cyp2j13 Cyp2j13_promoter_305bp               Peak-Decline        25
  Fig2C     Mif      Mif_promoter_73bp                Low-Decline        25
  Fig3A    Sgk1         Sgk1_exon_53bp                Low-Decline        25
  Fig3A    Sgk1   Sgk1_intergenic_66bp                   Peak-Low        25
  Fig3B   Igf1r      Igf1r_intron_72bp               Peak-Decline        25
  Fig3B   Igf1r      Igf1r_intron_95bp     Low-Decline + Peak-Low        25
 n_na mean_meth_pct n_cpgs_in_region
    1          34.7                1
    3          44.1                1
    8          77.7                5
    0          15.8                4
    0          66.2                7
    0          52.4                5
    0          54.2                7
    0          22.0                7
    0          61.0                7

Group means per locus:
 figure    gene                  locus   group mean_meth n
  Fig2A  Atp5g3       Atp5g3_909bp_TSS    Peak       6.2 8
  Fig2A  Atp5g3       Atp5g3_909bp_TSS Decline      69.0 7
  Fig2A  Atp5g3       Atp5g3_909bp_TSS     Low      33.1 9
  Fig2A  Atp5g3       Atp5g3_910bp_TSS    Peak      10.1 7
  Fig2A  Atp5g3       Atp5g3_910bp_TSS Decline      88.1 7
  Fig2A  Atp5g3       Atp5g3_910bp_TSS     Low      35.4 8
  Fig2B  Cyp1a1   Cyp1a1_promoter_52bp    Peak      77.1 5
  Fig2B  Cyp1a1   Cyp1a1_promoter_52bp Decline      98.5 5
  Fig2B  Cyp1a1   Cyp1a1_promoter_52bp     Low      63.3 7
  Fig2B Cyp2j13 Cyp2j13_promoter_305bp    Peak      38.6 8
  Fig2B Cyp2j13 Cyp2j13_promoter_305bp Decline       0.0 8
  Fig2B Cyp2j13 Cyp2j13_promoter_305bp     Low       9.7 9
  Fig2C     Mif      Mif_promoter_73bp    Peak      66.5 8
  Fig2C     Mif      Mif_promoter_73bp Decline      83.2 8
  Fig2C     Mif      Mif_promoter_73bp     Low      50.8 9
  Fig3A    Sgk1         Sgk1_exon_53bp    Peak      35.7 8
  Fig3A    Sgk1         Sgk1_exon_53bp Decline      44.1 8
  Fig3A    Sgk1         Sgk1_exon_53bp     Low      74.8 9
  Fig3A    Sgk1   Sgk1_intergenic_66bp    Peak      45.3 8
  Fig3A    Sgk1   Sgk1_intergenic_66bp Decline      51.6 8
  Fig3A    Sgk1   Sgk1_intergenic_66bp     Low      64.4 9
  Fig3B   Igf1r      Igf1r_intron_72bp    Peak      15.8 8
  Fig3B   Igf1r      Igf1r_intron_72bp Decline      30.0 8
  Fig3B   Igf1r      Igf1r_intron_72bp     Low      20.3 9
  Fig3B   Igf1r      Igf1r_intron_95bp    Peak      75.8 8
  Fig3B   Igf1r      Igf1r_intron_95bp Decline      77.2 8
  Fig3B   Igf1r      Igf1r_intron_95bp     Low      33.3 9

Done. Next step: run plot_figures_2_and_3.py