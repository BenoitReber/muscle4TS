# MUSCLE4TS: MUltivariate Sparse CLustering for Extremes in Time Series

MUSCLE4TS is a new method implemented in R dedicated to analyzing and modeling extreme events in heavy-tailed time series. It presents an innovative, heuristic adaptation of the **MUSCLE** (MUltivariate Sparse CLustering for Extremes) algorithm, originally developed by Nicolas Meyer and Olivier Wintenberger (2024) for the i.i.d. setting. The core purpose of this work is to detect components of a time series that experience large values simultaneously or in close temporal succession during extreme events, which can span multiple time frames. This project is a key part of a research project by Benoît Reber, supervised by Prof. Wintenberger, Assoc. Prof. Meyer, and Assoc. Prof. Buritica.

The main features of MUSCLE4TS are the following:

* **Extremal Pattern Detection**: Identifies sparse, co-occurrent patterns of extreme values in a multivariate time series.
* **Automatic Parameter Selection**: Automatically selects the optimal block length and the proportion of data considered extreme, addressing a major challenge in extreme value analysis.
* **Sparse Representation**: Provides a sparse, interpretable model of the spectral tail process, which describes the dependence structure of extremes.
* **Versatile Application**: The method is applicable to various heavy-tailed time series datasets, including financial returns and commodity prices.

---

## Content of the repository

The project contains the following core files:

* `muscle.r`: This file contains the original MUSCLE algorithm for the i.i.d. (independent and identically distributed) setting, as developed by Nicolas Meyer and Olivier Wintenberger. It has been slightly modified so that some of its components are reusable within the `muscle_TS.r` script.
* `muscle_TS.r`: This is the main implementation of the **MUSCLE4TS** algorithm, which extends the original MUSCLE algorithm to handle heavy-tailed time series data.
* `Muscle4TS.Rmd`: This R Markdown notebook provides a comprehensive overview of the **MUSCLE4TS** method. It includes some theoretical foundations, a description of the heuristic behind the new method, and several experiments on both synthetic and real-world financial datasets.
* `Muscle4TS.pdf`: This PDF file contains the compiled version of the `Muscle4TS.Rmd` notebook, providing a more accessible format for readers.

***

### Experiments and Findings

The `muscle4ts.Rmd` notebook provides a detailed walkthrough of the method and its application to various datasets.

1.  **Synthetic Data**: Experiments on simple heavy-tailed AR(1) and MA(2) models demonstrate that **MUSCLE4TS** can effectively capture known temporal patterns of extreme events.
2.  **US Industry Sector Returns**: Analysis of daily returns from 49 US industry portfolios successfully identifies plausible extreme event patterns, such as co-movement among sectors related to consumer spending. However, some detected patterns are difficult to interpret.
3.  **Commodity Futures Prices**: The method is applied to daily commodity futures prices. While it identifies some intuitive relationships (e.g., between WTI and Brent crude oil prices), it surprisingly detects a high degree of independence in the extreme regions.
4.  **High-Frequency Data**: An analysis of 1-minute data for gold, S&P 500, and cryptocurrencies shows the method can capture intra-asset dependence but struggles to identify strong inter-asset dependencies.

The results suggest that while **MUSCLE4TS** is a promising new tool, it requires further development, including a more rigorous theoretical foundation and benchmarking against other state-of-the-art methods. (We refer to the notebook for more details.)
