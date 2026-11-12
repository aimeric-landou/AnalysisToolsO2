# Systematics and Final Results

## 1. Run unfolding and systematics
Before running, check:

- `Jets/JetSpectrum_inputs.h`: nominal Data, MC, background files, and workflows for the selected centrality.
- `Jets/JetSpectrum_Systematics_input.h`: tracking-efficiency, background, Vertex, and Centrality variation inputs.
- `Jets/JetSpectrum_settings.h`: centrality, binning, prior, display ranges, and other analysis settings.
- `Jets/JetSpectrum_RunMacro.C`: SVD k, Bayes iterations, and the functions to execute.

To produce all systematic uncertainties, enable these calls in the RunMacro main function:

`DrawTotalSys()` calculates unfolding-method, k±1, measured/unfolded pT-range, prior, background, tracking-efficiency, Vertex, and Centrality uncertainties, and includes the configured Ncoll uncertainty. The four unfolding sources are combined using their bin-by-bin maximum; the remaining sources are added in quadrature. Individual plotting calls are also available in RunMacro

## 2. Save the ROOT results
Enable saving in `Jets/JetSpectrum_settings.h` and choose a filename for each centrality
The file contains:
- Top level: `HJetSpectra_measured`, `HJetSpectra_Svdkdefault`, and `HJetSpectra_refolded`.
- `Systematics_uncertainty/`: individual and total relative uncertainties, with names starting with `h1_syst_`.
- `Systematics_variation/`: nominal and varied unfolded spectra used for each source, retained for subsequent Rcp systematic calculations.

## 3. Plot the final jet spectra and Rcp perliminary plots
Wenhui's result files are stored in `rootfile_folder/`:

Run the standalone plotting macro; no unfolding is repeated:
$ root -l AnalysisToolsO2/Jets/JetSpectrum_PbPbResult.C $

In the file "JetSpectrum_PbPbResult.C" main Macro:
{
Draw_JetSpectrumPbPbSpectra();
Draw_JetSpectrumPbPbRcp();
}