# Cefprozil's Effect on Gut Microbiota

This project analyzes the effect of Cefprozil on gut microbiota based on data collected from a previous experiment involving adults aged 21-35 with regular intestinal transit. The analysis is inspired by the findings in the article: *“The initial state of the human gut microbiome determines its reshaping by antibiotics”*.

## Introduction
The gut microbiota, which contains up to 1000 species, is essential for human body functioning and is affected by environmental conditions, physical activity, diseases, diet, and hormones. The microbiota plays a crucial role in human physiology, including immune system development, digestion, and detoxification reactions. Understanding the microbiota can help adjust individual treatments and monitor physiological changes. The 16s subunit of the ribosome is often sequenced using NGS technology to study bacterial populations, and bioinformatic tools are used for advanced analysis (1-4).

The use of antibiotics for a short duration can alter the balance of microorganisms in the gut, leading to persistent imbalances that can potentially contribute to the onset or progression of illness (5). Cefprozil, a beta-lactam antibiotic, functions by targeting specific penicillin-binding proteins (PBPs) within the bacterial cell wall, hindering the final step of bacterial cell wall synthesis (6).

This project focuses on analyzing how Cefprozil affects gut microbiota. It includes data manipulations and visualizations to represent the relative abundance of bacterial families in both control and exposed participants. The main visualizations include:
- A bar plot representing the relative abundance of bacterial families.
- Boxplots showing the changes in relative abundance of the bacterial families most affected by antibiotic exposure.

## Data
The analysis is based on the data collected from adults aged 21-35 with regular intestinal transit. The data includes:
- `txt` files containing metadata.

## Analysis
The analysis involves the following steps:
1. Loading and preprocessing the data.
2. Creating a bar plot to represent the relative abundance of bacterial families in controls and exposed participants.
3. Creating boxplots to show the changes in relative abundance of the most affected bacterial families.


## Installation
To use this program, you need to install a few libraries. Use the provided `requirements.txt` file to install all necessary dependencies:

```bash
pip install -r requirements.txt
```

## Usage
To run the code, follow these steps:
1. Clone the repository:
    ```bash
    git clone https://github.com/Mazalik/CefprozilGutImpact.git
    cd CefprozilGutBiome
    ```
2. Provide the data files (`txt`) in the appropriate directories.
3. Run the analysis script and specify the input and output directories:
    ```bash
    python analyze_cefprozil_effect.py path/to/input_files path/to/save_results
    ```

## Results
The analysis will produce several output files in the specified output directory:
- `relative_abundance_plot.png`: A bar plot representing the relative abundance of bacterial families.
- `families_exposed_plot.png`: A plot showing the abundance of bacterial families for each family at different time points for the exposed patients. 
- `combined_family_<family_name>_plot.png`: Plots showing the combined data only for the significantly changed families from both control and exposed groups.

## Tests
This project includes tests to validate the functionality of the data preprocessing, analysis, and visualization scripts. Tests cover:
- Data loading and validation.
- Correct generation of bar plots and boxplots.
- Accuracy of statistical calculations (if applicable).

To run the tests, use the following command:
```bash
pytest
```

## References
1. D’Argenio, V. & Salvatore, F. The role of the gut microbiome in the healthy adult status. Clinica Chimica Acta 451, 97–102 (2015).
2. Shreiner, A. B. et al. Dicer (E-7): sc-393328. Cell 69, 393328 (2015).
3. Konstantinidis, T. et al. Effects of antibiotics upon the gut microbiome: A review of the literature. Biomedicines vol. 8 1–15 Preprint at <https://doi.org/10.3390/biomedicines8110502> (2020).
4. Weersma, R. K., Zhernakova, A. & Fu, J. Interaction between drugs and the gut microbiome. Gut 69, 1510–1519 (2020).
5. Lange, K., Buerger, M., Stallmach, A. & Bruns, T. Effects of Antibiotics on Gut Microbiota. Digestive Diseases 34, 260–268 (2016).
6. Chin L, Cummings R, Le D, Pon A, Knox C, Wilson M. DrugBank 5.0: a major update to the DrugBank database for 2018. Nucleic Acids Res. 2017 Nov 8. doi: 10.1093/nar/gkx1037.
7. Raymond, F. et al. The initial state of the human gut microbiome determines its reshaping by antibiotics. ISME Journal 10, 707–720 (2016).
8. Dai, D. et al. ‘GMrepo v2: a curated human gut microbiome database with special focus on disease markers and cross-dataset comparison’. Nucleic Acids Res (2022). Volume 50, Issue D1, Pages D777– D784.
9. Iredell, J., Brown, J. & Tagg, K. Antibiotic resistance in Enterobacteriaceae: Mechanisms and clinical implications. BMJ (Online) 352, (2016).
10. Knecht, H. et al. Effects of b-lactam antibiotics and fluoroquinolones on human gut microbiota in relation to clostridium difficile associated diarrhea. PLoS One 9, 1–8 (2014).


## Acknowledgments
This project was originally implemented as part of the Python programming course at the Weizmann Institute of Science taught by Gabor Szabo.