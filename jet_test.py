"""CMS Phase2 simulation analysis module """

#############################################
##                                         ##
## methods and classes for a common module ##
##                                         ##
#############################################

import logging
from bamboo.analysisutils import loadPlotIt
import os.path
from bamboo.analysismodules import AnalysisModule, HistogramsModule


class CMSPhase2SimModule(AnalysisModule):
    """ Base module for processing Phase2 flat trees """

    def __init__(self, args):
        super(CMSPhase2SimModule, self).__init__(args)

    def prepareTree(self, tree, sample=None, sampleCfg=None):
        from bamboo.treedecorators import decorateCMSPhase2SimTree
        from bamboo.dataframebackend import DataframeBackend
        t = decorateCMSPhase2SimTree(tree, isMC=True)
        be, noSel = DataframeBackend.create(t)
        return t, noSel, be, tuple()

    def mergeCounters(self, outF, infileNames, sample=None):
        from bamboo.root import gbl
        outF.cd()
        hNEvts = gbl.TH1F("nEvents", "Number of events", 1, 0., 1.)
        for fn in infileNames:
            f = gbl.TFile.Open(fn)
            tup = f.Get("myana/mytree")  # FIXME
            if not tup:
                raise RuntimeError(
                    "File {0} does not have a tree {1}".format(fn, self.args.treeName))
            hNEvts.Fill(0, tup.GetEntries())
        outF.cd()
        hNEvts.Write("nEvents")

    def readCounters(self, resultsFile):
        hNEvts = resultsFile.Get("nEvents")
        return {"nEvents": hNEvts.GetBinContent(1)}

# BEGIN cutflow reports, adapted from bamboo.analysisutils


logger = logging.getLogger(__name__)

_yieldsTexPreface = "\n".join(f"{ln}" for ln in
                              r"""\documentclass[12pt, landscape]{article}
\usepackage[margin=0.5in]{geometry}
\begin{document}
""".split("\n"))


def _texProcName(procName):
    if ">" in procName:
        procName = procName.replace(">", "$ > $")
    if "=" in procName:
        procName = procName.replace("=", "$ = $")
    if "_" in procName:
        procName = procName.replace("_", "\_")
    return procName


def _makeYieldsTexTable(report, samples, entryPlots, stretch=1.5, orientation="v", align="c", yieldPrecision=1, ratioPrecision=2):
    if orientation not in ("v", "h"):
        raise RuntimeError(
            f"Unsupported table orientation: {orientation} (valid: 'h' and 'v')")
    import plotit.plotit
    from plotit.plotit import Stack
    import numpy as np
    from itertools import repeat, count

    def colEntriesFromCFREntryHists(report, entryHists, precision=1):
        stacks_t = [(entryHists[entries[0]] if len(entries) == 1 else
                     Stack(entries=[entryHists[eName] for eName in entries]))
                    for entries in report.titles.values()]
        return stacks_t, ["& {0:.2e}".format(st_t.contents[1]) for st_t in stacks_t]

    def colEntriesFromCFREntryHists_forEff(report, entryHists, precision=1):
        stacks_t = [(entryHists[entries[0]] if len(entries) == 1 else
                     Stack(entries=[entryHists[eName] for eName in entries]))
                    for entries in report.titles.values()]
        return stacks_t, [" {0} ".format(st_t.contents[1]) for st_t in stacks_t]
    smp_signal = [smp for smp in samples if smp.cfg.type == "SIGNAL"]
    smp_mc = [smp for smp in samples if smp.cfg.type == "MC"]
    smp_data = [smp for smp in samples if smp.cfg.type == "DATA"]
    sepStr_v = "|l|"
    hdrs = ["Selection"]
    entries_smp = [[_texProcName(tName) for tName in report.titles.keys()]]
    stTotMC, stTotData = None, None
    if smp_signal:
        sepStr_v += "|"
        for sigSmp in smp_signal:
            _, colEntries = colEntriesFromCFREntryHists(report,
                                                        {eName: sigSmp.getHist(p) for eName, p in entryPlots.items()}, precision=yieldPrecision)
            sepStr_v += f"{align}|"
            hdrs.append(
                f"{_texProcName(sigSmp.cfg.yields_group)} {sigSmp.cfg.cross_section:f}pb")
            entries_smp.append(colEntries)
    if smp_mc:
        sepStr_v += "|"
        sel_list = []
        for mcSmp in smp_mc:
            _, colEntries = colEntriesFromCFREntryHists(report,
                                                        {eName: mcSmp.getHist(p) for eName, p in entryPlots.items()}, precision=yieldPrecision)
            sepStr_v += f"{align}|"
            if isinstance(mcSmp, plotit.plotit.Group):
                hdrs.append(_texProcName(mcSmp.name))
            else:
                hdrs.append(_texProcName(mcSmp.cfg.yields_group))
            entries_smp.append(_texProcName(colEntries))
            _, colEntries_forEff = colEntriesFromCFREntryHists_forEff(report,
                                                                      {eName: mcSmp.getHist(p) for eName, p in entryPlots.items()}, precision=yieldPrecision)
            colEntries_matrix = np.array(colEntries_forEff)
            sel_eff = np.array([100])
            for i in range(1, len(report.titles)):
                sel_eff = np.append(sel_eff, [float(
                    colEntries_matrix[i]) / float(colEntries_matrix[i-1]) * 100]).tolist()
            for i in range(len(report.titles)):
                sel_eff[i] = str(f"({sel_eff[i]:.2f}\%)")
            entries_smp.append(sel_eff)
            sel_list.append(colEntries_forEff)
        from bamboo.root import gbl
        sel_list_array = np.array(sel_list)
        gbl.gStyle.SetPalette(1)
        c1 = gbl.TCanvas("c1", "c1", 600, 400)
        cutflow_histo_FS = gbl.TH1F(
            "cutflow_histo", "Selection Cutflow", 6, 0, 6)
        cutflow_histo_FS.GetXaxis().SetTitle("Selection")
        cutflow_histo_FS.GetYaxis().SetTitle("Nevent")
        cutflow_histo_Delphes = gbl.TH1F("cutflow_histo", "Delphes", 6, 0, 6)
        for i in range(len(colEntries_forEff)):
            cutflow_histo_FS.Fill(i, float(sel_list_array[0, i]))
        cutflow_histo_FS.SetLineColor(2)
        cutflow_histo_FS.SetLineWidth(3)
        for i in range(len(colEntries_forEff)):
            cutflow_histo_Delphes.Fill(i, float(sel_list_array[1, i]))
        cutflow_histo_Delphes.SetLineColor(4)
        cutflow_histo_Delphes.SetLineWidth(3)
        cutflow_histo_FS.Draw("HIST")
        cutflow_histo_Delphes.Draw("SAME HIST")
        gbl.gPad.SetLogy()
        leg = gbl.TLegend(0.78, 0.695, 0.980, 0.935)
        leg.AddEntry(cutflow_histo_Delphes, "Delphes", "l")
        leg.AddEntry(cutflow_histo_FS, "FS", "l")
        leg.Draw()
        c1.SaveAs("cutflow.png")
        logger.info("Plot for selection cutflow is available")
    if smp_data:
        sepStr_v += f"|{align}|"
        hdrs.append("Data")
        stTotData, colEntries = colEntriesFromCFREntryHists(report, {eName: Stack(entries=[smp.getHist(
            p) for smp in smp_data]) for eName, p in entryPlots.items()}, precision=yieldPrecision)
        entries_smp.append(colEntries)
    if smp_data and smp_mc:
        sepStr_v += f"|{align}|"
        hdrs.append("Data/MC")
        colEntries = []
        for stData, stMC in zip(stTotData, stTotMC):
            dtCont = stData.contents
            mcCont = stMC.contents
            ratio = np.where(mcCont != 0., dtCont/mcCont,
                             np.zeros(dtCont.shape))
            ratioErr = np.where(mcCont != 0., np.sqrt(
                mcCont**2*stData.sumw2 + dtCont**2*(stMC.sumw2+stMC.syst2))/mcCont**2, np.zeros(dtCont.shape))
            colEntries.append("${{0:.{0}f}} \pm {{1:.{0}f}}$".format(
                ratioPrecision).format(ratio[1], ratioErr[1]))
        entries_smp.append(colEntries)
    if len(colEntries) < 2:
        logger.warning("No samples, so no yields.tex")
    return "\n".join(([
        f"\\begin{{tabular}}{{ {sepStr_v} }}",
        "    \\hline",
        "    {0} \\\\".format(" & ".join(hdrs)),
        "    \\hline"]+[
            "    {0} \\\\".format(
                " ".join(smpEntries[i] for smpEntries in entries_smp))
            for i in range(len(report.titles))])+[
        "    \\hline",
        "\\end{tabular}",
        "\\end{document}"
    ])


def printCutFlowReports(config, reportList, workdir=".", resultsdir=".", readCounters=lambda f: -1., eras=("all", None), verbose=False):
    """
    Print yields to the log file, and write a LaTeX yields table for each

    Samples can be grouped (only for the LaTeX table) by specifying the
    ``yields-group`` key (overriding the regular ``groups`` used for plots).
    The sample (or group) name to use in this table should be specified
    through the ``yields-title`` sample key.

    In addition, the following options in the ``plotIt`` section of
    the YAML configuration file influence the layout of the LaTeX yields table:

    - ``yields-table-stretch``: ``\\arraystretch`` value, 1.15 by default
    - ``yields-table-align``: orientation, ``h`` (default), samples in rows, or ``v``, samples in columns
    - ``yields-table-text-align``: alignment of text in table cells (default: ``c``)
    - ``yields-table-numerical-precision-yields``: number of digits after the decimal point for yields (default: 1)
    - ``yields-table-numerical-precision-ratio``: number of digits after the decimal point for ratios (default: 2)
    """
    eraMode, eras = eras
    if not eras:  # from config if not specified
        eras = list(config["eras"].keys())
    # helper: print one bamboo.plots.CutFlowReport.Entry

    def printEntry(entry, printFun=logger.info, recursive=True, genEvents=None):
        effMsg = ""
        if entry.parent:
            sumPass = entry.nominal.GetBinContent(1)
            sumTotal = entry.parent.nominal.GetBinContent(1)
            if sumTotal != 0.:
                effMsg = f", Eff={sumPass/sumTotal:.2%}"
                if genEvents:
                    effMsg += f", TotalEff={sumPass/genEvents:.2%}"
        printFun(
            f"Selection {entry.name}: N={entry.nominal.GetEntries()}, SumW={entry.nominal.GetBinContent(1)}{effMsg}")
        if recursive:
            for c in entry.children:
                printEntry(c, printFun=printFun,
                           recursive=recursive, genEvents=genEvents)
    # retrieve results files, get generated events for each sample
    from bamboo.root import gbl
    resultsFiles = dict()
    generated_events = dict()
    for smp, smpCfg in config["samples"].items():
        if "era" not in smpCfg or smpCfg["era"] in eras:
            resF = gbl.TFile.Open(os.path.join(resultsdir, f"{smp}.root"))
            resultsFiles[smp] = resF
            genEvts = None
            if "generated-events" in smpCfg:
                if isinstance(smpCfg["generated-events"], str):
                    genEvts = readCounters(resF)[smpCfg["generated-events"]]
                else:
                    genEvts = smpCfg["generated-events"]
            generated_events[smp] = genEvts
    has_plotit = None
    try:
        import plotit.plotit
        has_plotit = True
    except ImportError:
        has_plotit = False
    from bamboo.plots import EquidistantBinning as EqB

    class YieldPlot:
        def __init__(self, name):
            self.name = name
            self.plotopts = dict()
            self.axisTitles = ("Yield",)
            self.binnings = [EqB(1, 0., 1.)]
    for report in reportList:
        smpReports = {smp: report.readFromResults(
            resF) for smp, resF in resultsFiles.items()}
        # debug print
        for smp, smpRep in smpReports.items():
            if smpRep.printInLog:
                logger.info(f"Cutflow report {report.name} for sample {smp}")
                for root in smpRep.rootEntries():
                    printEntry(root, genEvents=generated_events[smp])
        # save yields.tex (if needed)
        if any(len(cb) > 1 or tt != cb[0] for tt, cb in report.titles.items()):
            if not has_plotit:
                logger.error(
                    f"Could not load plotit python library, no TeX yields tables for {report.name}")
            else:
                yield_plots = [YieldPlot(f"{report.name}_{eName}")
                               for tEntries in report.titles.values() for eName in tEntries]
                out_eras = []
                if len(eras) > 1 and eraMode in ("all", "combined"):
                    out_eras.append((f"{report.name}.tex", eras))
                if len(eras) == 1 or eraMode in ("split", "all"):
                    for era in eras:
                        out_eras.append((f"{report.name}_{era}.tex", [era]))
                for outName, iEras in out_eras:
                    pConfig, samples, plots, _, _ = loadPlotIt(
                        config, yield_plots, eras=iEras, workdir=workdir, resultsdir=resultsdir, readCounters=readCounters)
                    tabBlock = _makeYieldsTexTable(report, samples,
                                                   {p.name[len(
                                                       report.name)+1:]: p for p in plots},
                                                   stretch=pConfig.yields_table_stretch,
                                                   orientation=pConfig.yields_table_align,
                                                   align=pConfig.yields_table_text_align,
                                                   yieldPrecision=pConfig.yields_table_numerical_precision_yields,
                                                   ratioPrecision=pConfig.yields_table_numerical_precision_ratio)
                    with open(os.path.join(workdir, outName), "w") as ytf:
                        ytf.write("\n".join((_yieldsTexPreface, tabBlock)))
                    logger.info("Yields table for era(s) {0} was written to {1}".format(
                        ",".join(eras), os.path.join(workdir, outName)))

# END cutflow reports, adapted from bamboo.analysisutils


class CMSPhase2SimHistoModule(CMSPhase2SimModule, HistogramsModule):
    """ Base module for producing plots from Phase2 flat trees """

    def __init__(self, args):
        super(CMSPhase2SimHistoModule, self).__init__(args)

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        """ Customised cutflow reports and plots """
        if not self.plotList:
            self.plotList = self.getPlotList(resultsdir=resultsdir)
        from bamboo.plots import Plot, DerivedPlot, CutFlowReport
        plotList_cutflowreport = [
            ap for ap in self.plotList if isinstance(ap, CutFlowReport)]
        plotList_plotIt = [ap for ap in self.plotList if (isinstance(
            ap, Plot) or isinstance(ap, DerivedPlot)) and len(ap.binnings) == 1]
        eraMode, eras = self.args.eras
        if eras is None:
            eras = list(config["eras"].keys())
        if plotList_cutflowreport:
            printCutFlowReports(config, plotList_cutflowreport, workdir=workdir, resultsdir=resultsdir,
                                readCounters=self.readCounters, eras=(eraMode, eras), verbose=self.args.verbose)
        if plotList_plotIt:
            from bamboo.analysisutils import writePlotIt, runPlotIt
            cfgName = os.path.join(workdir, "plots.yml")
            writePlotIt(config, plotList_plotIt, cfgName, eras=eras, workdir=workdir, resultsdir=resultsdir,
                        readCounters=self.readCounters, vetoFileAttributes=self.__class__.CustomSampleAttributes, plotDefaults=self.plotDefaults)
            runPlotIt(cfgName, workdir=workdir, plotIt=self.args.plotIt,
                      eras=(eraMode, eras), verbose=self.args.verbose)


################################
##                            ##
## The actual analysis module ##
##                            ##
################################

class CMSPhase2SimTest(CMSPhase2SimHistoModule):
    """ Plotter module for Phase2 flat trees """

    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot, CutFlowReport, EquidistantBinning, VariableBinning
        from bamboo import treefunctions as op

        plots = []

        # definitions

        electrons = op.select(t.elec, lambda el: op.AND(
            el.pt > 20., op.abs(el.eta) < 2.5
        ))

        muons = op.select(t.muon, lambda mu: op.AND(
            mu.pt > 20., op.abs(mu.eta) < 2.5
        ))

        cleanedElectrons = op.select(electrons, lambda el: op.NOT(
            op.rng_any(muons, lambda mu: op.deltaR(el.p4, mu.p4) < 0.3)
        ))

        # we are taking the second isopass to be on which is equal to the medium working point
        isolatedElectrons = op.select(
            cleanedElectrons, lambda el: el.isopass & (1 << 2))

        identifiedElectrons = op.select(
            isolatedElectrons, lambda el: el.idpass & (1 << 2))

        cleanedMuons = op.select(muons, lambda mu: op.NOT(
            op.rng_any(electrons, lambda el: op.deltaR(mu.p4, el.p4) < 0.3)
        ))

        isolatedMuons = op.select(
            cleanedMuons, lambda mu: mu.isopass & (1 << 2))

        identifiedMuons = op.select(
            isolatedMuons, lambda mu: mu.idpass & (1 << 2))

        cleanedJets = op.select(t.jetpuppi, lambda j: op.AND(
            op.NOT(op.rng_any(identifiedElectrons,
                              lambda el: op.deltaR(el.p4, j.p4) < 0.3)),
            op.NOT(op.rng_any(identifiedMuons,
                              lambda mu: op.deltaR(mu.p4, j.p4) < 0.3))
        ))

        cleanedGoodJets30_0to1_5 = op.select(
            cleanedJets, lambda j: op.AND(j.pt > 30, op.abs(j.eta) < 1.5))

        cleanedGoodJets30_1_5to3 = op.select(cleanedJets, lambda j: op.AND(
            j.pt > 30, op.NOT(op.AND(op.abs(j.eta) < 1.5, op.abs(j.eta) > 3))))
        
        cleanedGoodJets30_3toInf = op.select(cleanedJets, lambda j: op.AND(
            j.pt > 30, op.abs(j.eta) > 3))
        
        cleanedGoodJets50_0to1_5 = op.select(
            cleanedJets, lambda j: op.AND(j.pt > 50, op.abs(j.eta) < 1.5))

        cleanedGoodJets50_1_5to3 = op.select(cleanedJets, lambda j: op.AND(
            j.pt > 50, op.NOT(op.AND(op.abs(j.eta) < 1.5, op.abs(j.eta) > 3))))

        cleanedGoodJets50_3toInf = op.select(cleanedJets, lambda j: op.AND(
            j.pt > 50, op.abs(j.eta) > 3))

        cleanedGoodJets100_0to1_5 = op.select(
            cleanedJets, lambda j: op.AND(j.pt > 100, op.abs(j.eta) < 1.5))

        cleanedGoodJets100_1_5to3 = op.select(cleanedJets, lambda j: op.AND(
            j.pt > 100, op.NOT(op.AND(op.abs(j.eta) < 1.5, op.abs(j.eta) > 3))))

        cleanedGoodJets100_3toInf = op.select(cleanedJets, lambda j: op.AND(
            j.pt > 100, op.abs(j.eta) > 3))

        cleanedGoodJets30 = op.select(cleanedJets, lambda j: j.pt > 30)

        cleanedGoodJets50 = op.select(cleanedJets, lambda j: j.pt > 50)

        cleanedGoodJets100 = op.select(cleanedJets, lambda j: j.pt > 100)

        met = op.select(t.metpuppi)

        sel1 = noSel.refine("nJet30", cut=[op.rng_len(cleanedGoodJets30) > 0])
        sel2 = noSel.refine("nJet50", cut=[op.rng_len(cleanedGoodJets50) > 0])
        sel3 = noSel.refine("nJet100", cut=[op.rng_len(cleanedGoodJets100) > 0])

        sel1_1 = noSel.refine(
            "nJet30_1", cut=[op.rng_len(cleanedGoodJets30) > 0])
        sel1_1_1 = noSel.refine(
            "nJet30_1_1", cut=[op.rng_len(cleanedGoodJets30_0to1_5) > 0])
        sel1_1_2 = noSel.refine(
            "nJet30_1_2", cut=[op.rng_len(cleanedGoodJets30_1_5to3) > 0])
        sel1_1_3 = noSel.refine(
            "nJet30_1_3", cut=[op.rng_len(cleanedGoodJets30_3toInf) > 0])
        sel2_1 = noSel.refine(
            "nJet50_1", cut=[op.rng_len(cleanedGoodJets50) > 0])
        sel1_2_1 = noSel.refine(
            "nJet50_1_1", cut=[op.rng_len(cleanedGoodJets50_0to1_5) > 0])
        sel1_2_2 = noSel.refine(
            "nJet50_1_2", cut=[op.rng_len(cleanedGoodJets50_1_5to3) > 0])
        sel1_2_3 = noSel.refine(
            "nJet50_1_3", cut=[op.rng_len(cleanedGoodJets50_3toInf) > 0])
        sel3_1 = noSel.refine(
            "nJet100_1", cut=[op.rng_len(cleanedGoodJets100) > 0])
        sel1_3_1 = noSel.refine(
            "nJet100_1_1", cut=[op.rng_len(cleanedGoodJets100_0to1_5) > 0])
        sel1_3_2 = noSel.refine(
            "nJet100_1_2", cut=[op.rng_len(cleanedGoodJets100_1_5to3) > 0])
        sel1_3_3 = noSel.refine(
            "nJet100_1_3", cut=[op.rng_len(cleanedGoodJets100_3toInf) > 0])

        sel1_2 = noSel.refine(
            "nJet30_2", cut=[op.rng_len(cleanedGoodJets30) > 1])
        sel2_1_1 = noSel.refine(
            "nJet30_2_1", cut=[op.rng_len(cleanedGoodJets30_0to1_5) > 1])
        sel2_1_2 = noSel.refine(
            "nJet30_2_2", cut=[op.rng_len(cleanedGoodJets30_1_5to3) > 1])
        sel2_1_3 = noSel.refine(
            "nJet30_2_3", cut=[op.rng_len(cleanedGoodJets30_3toInf) > 1])
        sel2_2 = noSel.refine(
            "nJet50_2", cut=[op.rng_len(cleanedGoodJets50) > 1])
        sel2_2_1 = noSel.refine(
            "nJet50_2_1", cut=[op.rng_len(cleanedGoodJets50_0to1_5) > 1])
        sel2_2_2 = noSel.refine(
            "nJet50_2_2", cut=[op.rng_len(cleanedGoodJets50_1_5to3) > 1])
        sel2_2_3 = noSel.refine(
            "nJet50_2_3", cut=[op.rng_len(cleanedGoodJets50_3toInf) > 1])
        sel3_2 = noSel.refine(
            "nJet100_2", cut=[op.rng_len(cleanedGoodJets100) > 1])
        sel2_3_1 = noSel.refine(
            "nJet100_2_1", cut=[op.rng_len(cleanedGoodJets100_0to1_5) > 1])
        sel2_3_2 = noSel.refine(
            "nJet100_2_2", cut=[op.rng_len(cleanedGoodJets100_1_5to3) > 1])
        sel2_3_3 = noSel.refine(
            "nJet100_2_3", cut=[op.rng_len(cleanedGoodJets100_3toInf) > 1])

        sel1_3 = noSel.refine(
            "nJet30_3", cut=[op.rng_len(cleanedGoodJets30) > 2])
        sel3_1_1 = noSel.refine(
            "nJet30_3_1", cut=[op.rng_len(cleanedGoodJets30_0to1_5) > 2])
        sel3_1_2 = noSel.refine(
            "nJet30_3_2", cut=[op.rng_len(cleanedGoodJets30_1_5to3) > 2])
        sel3_1_3 = noSel.refine(
            "nJet30_3_3", cut=[op.rng_len(cleanedGoodJets30_3toInf) > 2])
        sel2_3 = noSel.refine(
            "nJet50_3", cut=[op.rng_len(cleanedGoodJets50) > 2])
        sel3_2_1 = noSel.refine(
            "nJet50_3_1", cut=[op.rng_len(cleanedGoodJets50_0to1_5) > 2])
        sel3_2_2 = noSel.refine(
            "nJet50_3_2", cut=[op.rng_len(cleanedGoodJets50_1_5to3) > 2])
        sel3_2_3 = noSel.refine(
            "nJet50_3_3", cut=[op.rng_len(cleanedGoodJets50_3toInf) > 2])
        sel3_3 = noSel.refine(
            "nJet100_3", cut=[op.rng_len(cleanedGoodJets100) > 2])
        sel3_3_1 = noSel.refine(
            "nJet100_3_1", cut=[op.rng_len(cleanedGoodJets100_0to1_5) > 2])
        sel3_3_2 = noSel.refine(
            "nJet100_3_2", cut=[op.rng_len(cleanedGoodJets100_1_5to3) > 2])
        sel3_3_3 = noSel.refine(
            "nJet100_3_3", cut=[op.rng_len(cleanedGoodJets100_3toInf) > 2])
        
        sel1_4 = noSel.refine(
            "nJet30_4", cut=[op.rng_len(cleanedGoodJets30) > 3])
        sel4_1_1 = noSel.refine(
            "nJet30_4_1", cut=[op.rng_len(cleanedGoodJets30_0to1_5) > 3])
        sel4_1_2 = noSel.refine(
            "nJet30_4_2", cut=[op.rng_len(cleanedGoodJets30_1_5to3) > 3])
        sel4_1_3 = noSel.refine(
            "nJet30_4_3", cut=[op.rng_len(cleanedGoodJets30_3toInf) > 3])
        sel2_4 = noSel.refine(
            "nJet50_4", cut=[op.rng_len(cleanedGoodJets50) > 3])
        sel4_2_1 = noSel.refine(
            "nJet50_4_1", cut=[op.rng_len(cleanedGoodJets50_0to1_5) > 3])
        sel4_2_2 = noSel.refine(
            "nJet50_4_2", cut=[op.rng_len(cleanedGoodJets50_1_5to3) > 3])
        sel4_2_3 = noSel.refine(
            "nJet50_4_3", cut=[op.rng_len(cleanedGoodJets50_3toInf) > 3])
        sel3_4 = noSel.refine(
            "nJet100_4", cut=[op.rng_len(cleanedGoodJets100) > 3])
        sel4_3_1 = noSel.refine(
            "nJet100_4_1", cut=[op.rng_len(cleanedGoodJets100_0to1_5) > 3])
        sel4_3_2 = noSel.refine(
            "nJet100_4_2", cut=[op.rng_len(cleanedGoodJets100_1_5to3) > 3])
        sel4_3_3 = noSel.refine(
            "nJet100_4_3", cut=[op.rng_len(cleanedGoodJets100_3toInf) > 3])
        
        # plots

    #     # ### 30 GeV

        plots.append(Plot.make1D("nJets_jetPT_30GeV", op.rng_len(
            cleanedGoodJets30), noSel, EquidistantBinning(15, 0., 15.), title="nJets (jet p_{T} > 30GeV)"))

        plots.append(Plot.make1D("Jet1_pT_jetPT_30GeV", cleanedGoodJets30[0].pt, sel1_1, EquidistantBinning(
            50, 0., 4000.), title="Jet1_pT (jet p_{T} > 30GeV)"))

        plots.append(Plot.make1D("Jet1_eta_jetPT_30GeV", cleanedGoodJets30[0].eta, sel1_1, EquidistantBinning(
            30, -3, 3), title="Jet1_eta (jet p_{T} > 30GeV)"))
        
        plots.append(Plot.make1D("Jet1_pT_0to1_5_jetPT_30GeV", cleanedGoodJets30_0to1_5[0].pt, sel1_1_1, EquidistantBinning(
            50, 0, 4000), title="Jet1_pT 0 < eta < 1.5, (jet p_{T} > 30GeV)"))

        plots.append(Plot.make1D("Jet1_pT_1_5to3_jetPT_30GeV", cleanedGoodJets30_1_5to3[0].pt, sel1_1_2, EquidistantBinning(
            50, 0, 4000), title="Jet1_pT 1.5 < eta < 3, (jet p_{T} > 30GeV)"))
        
        plots.append(Plot.make1D("Jet1_pT_3toInf_jetPT_30GeV", cleanedGoodJets30_3toInf[0].pt, sel1_1_3, EquidistantBinning(
            50, 0, 1000), title="Jet1_pT eta > 3, (jet p_{T} > 30GeV)"))

        plots.append(Plot.make1D("Jet2_pT_jetPT_30GeV", cleanedGoodJets30[1].pt, sel1_2, EquidistantBinning(
            50, 0., 4000.), title="Jet2_pT (jet p_{T} > 30GeV)"))

        plots.append(Plot.make1D("Jet2_eta_jetPT_30GeV", cleanedGoodJets30[1].eta, sel1_2, EquidistantBinning(
            30, -3, 3), title="Jet2_eta (jet p_{T} > 30GeV)"))

        plots.append(Plot.make1D("Jet2_pT_0to1_5_jetPT_30GeV", cleanedGoodJets30_0to1_5[1].pt, sel2_1_1, EquidistantBinning(
            50, 0, 4000), title="Jet2_pT 0 < eta < 1.5, (jet p_{T} > 30GeV)"))

        plots.append(Plot.make1D("Jet2_pT_1_5to3_jetPT_30GeV", cleanedGoodJets30_1_5to3[1].pt, sel2_1_2, EquidistantBinning(
            50, 0, 4000), title="Jet2_pT 1.5 < eta < 3, (jet p_{T} > 30GeV)"))

        plots.append(Plot.make1D("Jet2_pT_3toInf_jetPT_30GeV", cleanedGoodJets30_3toInf[1].pt, sel2_1_3, EquidistantBinning(
            50, 0, 500), title="Jet2_pT eta > 3, (jet p_{T} > 30GeV)"))

        plots.append(Plot.make1D("Jet3_pT_jetPT_30GeV", cleanedGoodJets30[2].pt, sel1_3, EquidistantBinning(
            50, 0., 2000.), title="Jet3_pT (jet p_{T} > 30GeV)"))

        plots.append(Plot.make1D("Jet3_eta_jetPT_30GeV", cleanedGoodJets30[2].eta, sel1_3, EquidistantBinning(
            30, -3, 3), title="Jet3_eta (jet p_{T} > 30GeV)"))

        plots.append(Plot.make1D("Jet3_pT_0to1_5_jetPT_30GeV", cleanedGoodJets30_0to1_5[2].pt, sel3_1_1, EquidistantBinning(
            50, 0, 2000), title="Jet3_pT 0 < eta < 1.5, (jet p_{T} > 30GeV)"))

        plots.append(Plot.make1D("Jet3_pT_1_5to3_jetPT_30GeV", cleanedGoodJets30_1_5to3[2].pt, sel3_1_2, EquidistantBinning(
            50, 0, 2000), title="Jet3_pT 1.5 < eta < 3, (jet p_{T} > 30GeV)"))

        plots.append(Plot.make1D("Jet3_pT_3toInf_jetPT_30GeV", cleanedGoodJets30_3toInf[2].pt, sel3_1_3, EquidistantBinning(
            50, 0, 200), title="Jet3_pT eta > 3, (jet p_{T} > 30GeV)"))

        plots.append(Plot.make1D("Jet4_pT_jetPT_30GeV", cleanedGoodJets30[3].pt, sel1_4, EquidistantBinning(
            50, 0., 2000.), title="Jet4_pT (jet p_{T} > 30GeV)"))

        plots.append(Plot.make1D("Jet4_eta_jetPT_30GeV", cleanedGoodJets30[3].eta, sel1_4, EquidistantBinning(
            30, -3, 3), title="Jet4_eta (jet p_{T} > 30GeV)"))

        plots.append(Plot.make1D("Jet4_pT_0to1_5_jetPT_30GeV", cleanedGoodJets30_0to1_5[3].pt, sel4_1_1, EquidistantBinning(
            50, 0, 2000), title="Jet4_pT 0 < eta < 1.5, (jet p_{T} > 30GeV)"))

        plots.append(Plot.make1D("Jet4_pT_1_5to3_jetPT_30GeV", cleanedGoodJets30_1_5to3[3].pt, sel4_1_2, EquidistantBinning(
            50, 0, 2000), title="Jet4_pT 1.5 < eta < 3, (jet p_{T} > 30GeV)"))

        plots.append(Plot.make1D("Jet4_pT_3toInf_jetPT_30GeV", cleanedGoodJets30_3toInf[3].pt, sel4_1_3, EquidistantBinning(
            50, 0, 200), title="Jet4_pT eta > 3, (jet p_{T} > 30GeV)"))

        plots.append(Plot.make1D("MET_jetPT_30GeV", met[0].pt, sel1, EquidistantBinning(
            50, 0, 1000), title="MET (jet p_{T} > 30GeV)"))

    #     # ### 50 GeV

        plots.append(Plot.make1D("nJets_jetPT_50GeV", op.rng_len(
            cleanedGoodJets50), noSel, EquidistantBinning(15, 0., 15.), title="nJets (jet p_{T} > 50GeV)"))

        plots.append(Plot.make1D("Jet1_pT_jetPT_50GeV", cleanedGoodJets50[0].pt, sel2_1, EquidistantBinning(
            50, 0., 4000.), title="Jet1_pT (jet p_{T} > 50GeV)"))

        plots.append(Plot.make1D("Jet1_eta_jetPT_50GeV", cleanedGoodJets50[0].eta, sel2_1, EquidistantBinning(
            30, -3, 3), title="Jet1_eta (jet p_{T} > 50GeV)"))
        
        plots.append(Plot.make1D("Jet1_pT_0to1_5_jetPT_50GeV", cleanedGoodJets50_0to1_5[0].pt, sel1_2_1, EquidistantBinning(
            50, 0, 4000), title="Jet1_pT 0 < eta < 1.5, (jet p_{T} > 50GeV)"))

        plots.append(Plot.make1D("Jet1_pT_1_5to3_jetPT_50GeV", cleanedGoodJets50_1_5to3[0].pt, sel1_2_2, EquidistantBinning(
            50, 0, 4000), title="Jet1_pT 1.5 < eta < 3, (jet p_{T} > 50GeV)"))

        plots.append(Plot.make1D("Jet1_pT_3toInf_jetPT_50GeV", cleanedGoodJets50_3toInf[0].pt, sel1_2_3, EquidistantBinning(
            50, 0, 1000), title="Jet1_pT eta > 3, (jet p_{T} > 50GeV)"))

        plots.append(Plot.make1D("Jet2_pT_jetPT_50GeV", cleanedGoodJets50[1].pt, sel2_2, EquidistantBinning(
            50, 0., 4000.), title="Jet2_pT (jet p_{T} > 50GeV)"))

        plots.append(Plot.make1D("Jet2_eta_jetPT_50GeV", cleanedGoodJets50[1].eta, sel2_2, EquidistantBinning(
            30, -3, 3), title="Jet2_eta (jet p_{T} > 50GeV)"))
        
        plots.append(Plot.make1D("Jet2_pT_0to1_5_jetPT_50GeV", cleanedGoodJets50_0to1_5[1].pt, sel2_2_1, EquidistantBinning(
            50, 0, 4000), title="Jet2_pT 0 < eta < 1.5, (jet p_{T} > 50GeV)"))

        plots.append(Plot.make1D("Jet2_pT_1_5to3_jetPT_50GeV", cleanedGoodJets50_1_5to3[1].pt, sel2_2_2, EquidistantBinning(
            50, 0, 4000), title="Jet2_pT 1.5 < eta < 3, (jet p_{T} > 50GeV)"))

        plots.append(Plot.make1D("Jet2_pT_3toInf_jetPT_50GeV", cleanedGoodJets50_3toInf[1].pt, sel2_2_3, EquidistantBinning(
            50, 0, 500), title="Jet2_pT eta > 3, (jet p_{T} > 50GeV)"))

        plots.append(Plot.make1D("Jet3_pT_jetPT_50GeV", cleanedGoodJets50[2].pt, sel2_3, EquidistantBinning(
            50, 0., 2000.), title="Jet3_pT (jet p_{T} > 50GeV)"))

        plots.append(Plot.make1D("Jet3_eta_jetPT_50GeV", cleanedGoodJets50[2].eta, sel2_3, EquidistantBinning(
            30, -3, 3), title="Jet3_eta (jet p_{T} > 50GeV)"))
        
        plots.append(Plot.make1D("Jet3_pT_0to1_5_jetPT_50GeV", cleanedGoodJets50_0to1_5[2].pt, sel3_2_1, EquidistantBinning(
            50, 0, 2000), title="Jet3_pT 0 < eta < 1.5, (jet p_{T} > 50GeV)"))

        plots.append(Plot.make1D("Jet3_pT_1_5to3_jetPT_50GeV", cleanedGoodJets50_1_5to3[2].pt, sel3_2_2, EquidistantBinning(
            50, 0, 2000), title="Jet3_pT 1.5 < eta < 3, (jet p_{T} > 50GeV)"))

        plots.append(Plot.make1D("Jet3_pT_3toInf_jetPT_50GeV", cleanedGoodJets50_3toInf[2].pt, sel3_2_3, EquidistantBinning(
            50, 0, 200), title="Jet3_pT eta > 3, (jet p_{T} > 50GeV)"))

        plots.append(Plot.make1D("Jet4_pT_jetPT_50GeV", cleanedGoodJets50[3].pt, sel2_4, EquidistantBinning(
            50, 0., 2000.), title="Jet4_pT (jet p_{T} > 50GeV)"))

        plots.append(Plot.make1D("Jet4_eta_jetPT_50GeV", cleanedGoodJets50[3].eta, sel2_4, EquidistantBinning(
            30, -3, 3), title="Jet4_eta (jet p_{T} > 50GeV)"))

        plots.append(Plot.make1D("Jet4_pT_0to1_5_jetPT_50GeV", cleanedGoodJets50_0to1_5[3].pt, sel4_2_1, EquidistantBinning(
            50, 0, 2000), title="Jet4_pT 0 < eta < 1.5, (jet p_{T} > 30GeV)"))

        plots.append(Plot.make1D("Jet4_pT_1_5to3_jetPT_50GeV", cleanedGoodJets50_1_5to3[3].pt, sel4_2_2, EquidistantBinning(
            50, 0, 2000), title="Jet4_pT 1.5 < eta < 3, (jet p_{T} > 50GeV)"))

        plots.append(Plot.make1D("Jet4_pT_3toInf_jetPT_50GeV", cleanedGoodJets50_3toInf[3].pt, sel4_2_3, EquidistantBinning(
            50, 0, 200), title="Jet4_pT eta > 3, (jet p_{T} > 50GeV)"))

        plots.append(Plot.make1D("MET_jetPT_50GeV", met[0].pt, sel2, EquidistantBinning(
            50, 0, 1000), title="MET (jet p_{T} > 50GeV)"))

    #  ### 100 GeV

        plots.append(Plot.make1D("nJets_jetPT_100GeV", op.rng_len(
            cleanedGoodJets100), noSel, EquidistantBinning(15, 0., 15.), title="nJets (jet p_{T} > 100GeV)"))

        plots.append(Plot.make1D("Jet1_pT_jetPT_100GeV", cleanedGoodJets100[0].pt, sel3_1, EquidistantBinning(
            50, 0., 4000.), title="Jet1_pT (jet p_{T} > 100GeV)"))

        plots.append(Plot.make1D("Jet1_eta_jetPT_100GeV", cleanedGoodJets100[0].eta, sel3_1, EquidistantBinning(
            30, -3, 3), title="Jet1_eta (jet p_{T} > 100GeV)"))
        
        plots.append(Plot.make1D("Jet1_pT_0to1_5_jetPT_100GeV", cleanedGoodJets100_0to1_5[0].pt, sel1_3_1, EquidistantBinning(
            50, 0, 4000), title="Jet1_pT 0 < eta < 1.5, (jet p_{T} > 100GeV)"))

        plots.append(Plot.make1D("Jet1_pT_1_5to3_jetPT_100GeV", cleanedGoodJets100_1_5to3[0].pt, sel1_3_2, EquidistantBinning(
            50, 0, 4000), title="Jet1_pT 1.5 < eta < 3, (jet p_{T} > 100GeV)"))

        plots.append(Plot.make1D("Jet1_pT_3toInf_jetPT_100GeV", cleanedGoodJets100_3toInf[0].pt, sel1_3_3, EquidistantBinning(
            50, 0, 1000), title="Jet1_pT eta > 3, (jet p_{T} > 100GeV)"))

        plots.append(Plot.make1D("Jet2_pT_jetPT_100GeV", cleanedGoodJets100[1].pt, sel3_2, EquidistantBinning(
            50, 0., 4000.), title="Jet2_pT (jet p_{T} > 100GeV)"))

        plots.append(Plot.make1D("Jet2_eta_jetPT_100GeV", cleanedGoodJets100[1].eta, sel3_2, EquidistantBinning(
            30, -3, 3), title="Jet2_eta (jet p_{T} > 100GeV)"))
        
        plots.append(Plot.make1D("Jet2_pT_0to1_5_jetPT_100GeV", cleanedGoodJets100_0to1_5[1].pt, sel2_3_1, EquidistantBinning(
            50, 0, 4000), title="Jet2_pT 0 < eta < 1.5, (jet p_{T} > 100GeV)"))

        plots.append(Plot.make1D("Jet2_pT_1_5to3_jetPT_100GeV", cleanedGoodJets100_1_5to3[1].pt, sel2_3_2, EquidistantBinning(
            50, 0, 4000), title="Jet2_pT 1.5 < eta < 3, (jet p_{T} > 100GeV)"))

        plots.append(Plot.make1D("Jet2_pT_3toInf_jetPT_100GeV", cleanedGoodJets100_3toInf[1].pt, sel2_3_3, EquidistantBinning(
            50, 0, 500), title="Jet2_pT eta > 3, (jet p_{T} > 100GeV)"))

        plots.append(Plot.make1D("Jet3_pT_jetPT_100GeV", cleanedGoodJets100[2].pt, sel3_3, EquidistantBinning(
            50, 0., 2000.), title="Jet3_pT (jet p_{T} > 100GeV)"))

        plots.append(Plot.make1D("Jet3_eta_jetPT_100GeV", cleanedGoodJets100[2].eta, sel3_3, EquidistantBinning(
            30, -3, 3), title="Jet3_eta (jet p_{T} > 100GeV)"))
        
        plots.append(Plot.make1D("Jet3_pT_0to1_5_jetPT_100GeV", cleanedGoodJets100_0to1_5[2].pt, sel3_3_1, EquidistantBinning(
            50, 0, 2000), title="Jet3_pT 0 < eta < 1.5, (jet p_{T} > 100GeV)"))

        plots.append(Plot.make1D("Jet3_pT_1_5to3_jetPT_100GeV", cleanedGoodJets100_1_5to3[2].pt, sel3_3_2, EquidistantBinning(
            50, 0, 2000), title="Jet3_pT 1.5 < eta < 3, (jet p_{T} > 100GeV)"))

        plots.append(Plot.make1D("Jet3_pT_3toInf_jetPT_100GeV", cleanedGoodJets100_3toInf[2].pt, sel3_3_3, EquidistantBinning(
            50, 0, 200), title="Jet3_pT eta > 3, (jet p_{T} > 100GeV)"))

        plots.append(Plot.make1D("Jet4_pT_jetPT_100GeV", cleanedGoodJets100[3].pt, sel3_4, EquidistantBinning(
            50, 0., 2000.), title="Jet4_pT (jet p_{T} > 100GeV)"))

        plots.append(Plot.make1D("Jet4_eta_jetPT_100GeV", cleanedGoodJets100[3].eta, sel3_4, EquidistantBinning(
            30, -3, 3), title="Jet4_eta (jet p_{T} > 100GeV)"))
        
        plots.append(Plot.make1D("Jet4_pT_0to1_5_jetPT_100GeV", cleanedGoodJets100_0to1_5[3].pt, sel4_3_1, EquidistantBinning(
            50, 0, 2000), title="Jet4_pT 0 < eta < 1.5, (jet p_{T} > 100GeV)"))

        plots.append(Plot.make1D("Jet4_pT_1_5to3_jetPT_100GeV", cleanedGoodJets100_1_5to3[3].pt, sel4_3_2, EquidistantBinning(
            50, 0, 2000), title="Jet4_pT 1.5 < eta < 3, (jet p_{T} > 100GeV)"))

        plots.append(Plot.make1D("Jet4_pT_3toInf_jetPT_100GeV", cleanedGoodJets100_3toInf[3].pt, sel4_3_3, EquidistantBinning(
            50, 0, 200), title="Jet4_pT eta > 3, (jet p_{T} > 100GeV)"))

        plots.append(Plot.make1D("MET_jetPT_100GeV", met[0].pt, sel3, EquidistantBinning(
            50, 0, 1000), title="MET (jet p_{T} > 100GeV)"))

        # Efficiency Report on terminal and the .tex output

        cfr = CutFlowReport("yields")
        cfr.add(noSel, "None")
        cfr.add(sel1, "30GeV")
        cfr.add(sel2, "50GeV")
        cfr.add(sel3, "100GeV")

        plots.append(cfr)

        return plots
