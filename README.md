# The Single-Index/Cox (SIC) mixture cure model

This repository contains R codes associated with the article : 

<b> Amico, M., Van Keilegom, I. and Legrand, C. (2018). The single-index/Cox mixture cure model.</b>

You will find the following files:

<b>Functions-SIC.R</b> which contains the functions that have been used to run simulations and to perform the real data application.

The file <b>SICfunctions.pdf</b> contains the explanations of these functions, that is, a short description of their goal, usage of the functions, their arguments, and the values they take. 

Examples are given in the file <b>examples.R. </b> In these examples, the single-index/Cox mixture cure model is applied on a melanoma dataset containing data from an Eastern Cooperative Oncology Group (ECOG) phase III clinical trial on melanoma patients (e1684) (Kirkwood et al., 1996). This dataset is available in the R-package smcure. The purpose of the e1684 clinical trial was to evaluate high dose interferon alpha- 2b (IFN) regimen against placebo as postoperative adjuvant. The main outcome of interest is relapse-free survival computed as the time from randomization until first relapse or death, whatever comes first. There is 283 observations, and 3 covariates are available : the treatment (0 = control - 139 patients, 1 = treatment - 144 patients), the gender (0 = male - 170 patients, 1 = female - 113 patients), and the age of the patient ranging from 17 to 79 with a median of 47.5 years.

Do not hesitate to contact us at the address mailis.amico@kuleuven.be for more details or questions.

<b>Reference:</b>

Kirkwood, J., Strawderman, M., Ernstoff, M. ,Smith, T. ,Borden, E., and Blum, R. (1996). Interferon alfa-2b adjuvant therapy of high-risk resected cutaneous melanoma: the eastern cooperative oncology group trial est 1684, <i>Journal of Clinical Oncology</i>, <b>14</b>, 7-17.


