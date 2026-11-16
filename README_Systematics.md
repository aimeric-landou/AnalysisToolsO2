# Systematics and Final Results

## 1. Run unfolding and systematics
Before running, check:

- `Jets/JetSpectrum_inputs.h`: nominal Data, MC, background files, and workflows for the selected centrality.
- `Jets/JetSpectrum_Systematics_input.h`: tracking-efficiency, background, Vertex, and Centrality variation inputs.
- `Jets/JetSpectrum_settings.h`: centrality, binning, prior, display ranges, and other analysis settings.
- `Jets/JetSpectrum_RunMacro.C`: SVD k, Bayes iterations, and the functions to execute.

Optional systematics files are opened only when their source is requested. In `JetSpectrum_Systematics_input.h`, keep the input declarations and set an unused `fileName_...` entry to `""`. For example, background-only checks do not require tracking-efficiency, Vertex, or Centrality variation files. Calling a source with a missing input gives an error; `DrawTotalSys()` requires all sources and does not silently skip any of them.

To produce all systematic uncertainties, enable these calls in the RunMacro main function:

To run the nominal unfolding and all systematics together, use `DrawAndSave_JetSpectrumResultsAndSystematics(iDataset, iRadius, unfoldParameterInputList[0], unfoldParameterInputList[1], optionsAnalysis)`. This calls both `Draw_Pt_spectrum_unfolded_singleDataset()` and `DrawTotalSys()`. ROOT output still follows the `saveAnalysisResults` switch and `analysisResultsFileName` setting; there is no need to call the two functions again separately.

`DrawTotalSys()` calculates unfolding-method, k±1, measured/unfolded pT-range, prior, background, tracking-efficiency, Vertex, and Centrality uncertainties, and includes the configured Ncoll uncertainty. The four unfolding sources are combined using their bin-by-bin maximum; the remaining sources are added in quadrature. Individual plotting calls are also available in RunMacro

## 2. Save the ROOT results
Enable saving in `Jets/JetSpectrum_settings.h` and choose a filename for each centrality
The file contains:
- Top level: `HJetSpectra_measured`, `HJetSpectra_Svdkdefault`, and `HJetSpectra_refolded`.
- `Systematics_uncertainty/`: individual and total relative uncertainties, with names starting with `h1_syst_`.
- `Systematics_variation/`: nominal and varied unfolded spectra used for each source, retained for subsequent Rcp systematic calculations.

## 3. Plot the final jet spectra and Rcp perliminary plots
`JetSpectrum_PbPbResult.C` runs independently of RunMacro; run it only when you need the final preliminary plots, otherwise skip it.

Run the standalone plotting macro; no unfolding is repeated:
$ root -l AnalysisToolsO2/Jets/JetSpectrum_PbPbResult.C $

In the file "JetSpectrum_PbPbResult.C" main Macro:
{
Draw_JetSpectrumPbPbSpectra();
Draw_JetSpectrumPbPbRcp();
}
