# MotrpacRatTraining6moWATData 2.0.0

## Enhancements
* Moved analyses from data-raw/ to vignettes to improve visibility.
* Added online-only articles with code to generate manuscript figures.
* Updated analyses of plasma clinical analytes to more closely match what was done for the MoTrPAC Physiology manuscript.

## Bugfixes
* PROT_MITOCARTA_FGSEA results submitted in the first round of reviews calculated set size based on human genes rather than rat genes, so "Amino acid metabolism", "CI assembly factors", "OXPHOS assembly factors", and "mt-rRNA modifications" gene sets were filtered out. This has been fixed.
