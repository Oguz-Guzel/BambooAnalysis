"""CMS Phase2 simulation analysis module """

#############################################
##                                         ##
## methods and classes for a common module ##
##                                         ##
#############################################

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
            tup = f.Get("myana/mytree") ## FIXME
            if not tup:
                raise RuntimeError("File {0} does not have a tree {1}".format(fn, self.args.treeName))
            hNEvts.Fill(0, tup.GetEntries())
        outF.cd()
        hNEvts.Write("nEvents")
    def readCounters(self, resultsFile):
        hNEvts = resultsFile.Get("nEvents")
        return {"nEvents": hNEvts.GetBinContent(1)}

## BEGIN cutflow reports, adapted from bamboo.analysisutils

import logging
logger = logging.getLogger(__name__)
import os.path
from bamboo.analysisutils import loadPlotIt

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
        raise RuntimeError(f"Unsupported table orientation: {orientation} (valid: 'h' and 'v')")
    import plotit.plotit
    from plotit.plotit import Stack
    import numpy as np
    from itertools import repeat, count
    def colEntriesFromCFREntryHists(report, entryHists, precision=1):
        stacks_t = [ (entryHists[entries[0]] if len(entries) == 1 else
            Stack(entries=[entryHists[eName] for eName in entries]))
            for entries in report.titles.values() ]
        return stacks_t, [ "& {0:.2e}".format(st_t.contents[1]) for st_t in stacks_t ]
    def colEntriesFromCFREntryHists_forEff(report, entryHists, precision=1):
        stacks_t = [ (entryHists[entries[0]] if len(entries) == 1 else
            Stack(entries=[entryHists[eName] for eName in entries]))
            for entries in report.titles.values() ]
        return stacks_t, [ " {0} ".format(st_t.contents[1]) for st_t in stacks_t ]
    smp_signal = [smp for smp in samples if smp.cfg.type == "SIGNAL"]
    smp_mc = [smp for smp in samples if smp.cfg.type == "MC"]
    smp_data = [smp for smp in samples if smp.cfg.type == "DATA"]
    sepStr_v = "|l|"
    hdrs = ["Selection"]
    entries_smp = [ [_texProcName(tName) for tName in report.titles.keys()] ]
    stTotMC, stTotData = None, None
    if smp_signal:
        sepStr_v += "|"
        for sigSmp in smp_signal:
            _, colEntries = colEntriesFromCFREntryHists(report,
                { eName : sigSmp.getHist(p) for eName, p in entryPlots.items() }, precision=yieldPrecision)
            sepStr_v += f"{align}|"
            hdrs.append(f"{_texProcName(sigSmp.cfg.yields_group)} {sigSmp.cfg.cross_section:f}pb")
            entries_smp.append(colEntries)
    if smp_mc:
        sepStr_v += "|"
        sel_list = []
        for mcSmp in smp_mc:
            _, colEntries = colEntriesFromCFREntryHists(report,
                { eName : mcSmp.getHist(p) for eName, p in entryPlots.items() }, precision=yieldPrecision)
            sepStr_v += f"{align}|"
            if isinstance(mcSmp, plotit.plotit.Group):
                hdrs.append(_texProcName(mcSmp.name))
            else:
                hdrs.append(_texProcName(mcSmp.cfg.yields_group))
            entries_smp.append(_texProcName(colEntries))
            _, colEntries_forEff = colEntriesFromCFREntryHists_forEff(report,
                { eName : mcSmp.getHist(p) for eName, p in entryPlots.items() }, precision=yieldPrecision)
            colEntries_matrix = np.array(colEntries_forEff)
            sel_eff = np.array([100])
            for i in range(1, len(report.titles)):
                sel_eff = np.append(sel_eff, [ float(colEntries_matrix[i]) / float(colEntries_matrix[i-1]) *100 ]).tolist()
            for i in range(len(report.titles)):
                sel_eff[i] = str(f"({sel_eff[i]:.2f}\%)")
            entries_smp.append(sel_eff)
            sel_list.append(colEntries_forEff)
        from bamboo.root import gbl
        sel_list_array = np.array(sel_list)
        gbl.gStyle.SetPalette(1)
        c1 = gbl.TCanvas("c1", "c1", 600, 400)
        cutflow_histo_FS = gbl.TH1F("cutflow_histo", "Selection Cutflow", 6, 0, 6)
        cutflow_histo_FS.GetXaxis().SetTitle("Selection")
        cutflow_histo_FS.GetYaxis().SetTitle("Nevent")
        cutflow_histo_Delphes = gbl.TH1F("cutflow_histo", "Delphes", 6, 0, 6)
        for i in range(len(colEntries_forEff)):
            cutflow_histo_FS.Fill(i, float(sel_list_array[0,i]))
        cutflow_histo_FS.SetLineColor(2)
        cutflow_histo_FS.SetLineWidth(3)
        for i in range(len(colEntries_forEff)):
            cutflow_histo_Delphes.Fill(i, float(sel_list_array[1,i]))
        cutflow_histo_Delphes.SetLineColor(4)
        cutflow_histo_Delphes.SetLineWidth(3)
        cutflow_histo_FS.Draw("HIST")
        cutflow_histo_Delphes.Draw("SAME HIST")
        gbl.gPad.SetLogy()
        leg = gbl.TLegend(0.78,0.695,0.980,0.935)
        leg.AddEntry(cutflow_histo_Delphes, "Delphes", "l")
        leg.AddEntry(cutflow_histo_FS, "FS", "l")
        leg.Draw()
        c1.SaveAs("cutflow.gif")
        logger.info("Plot for selection cutflow is available")
    if smp_data:
        sepStr_v += f"|{align}|"
        hdrs.append("Data")
        stTotData, colEntries = colEntriesFromCFREntryHists(report, { eName : Stack(entries=[smp.getHist(p) for smp in smp_data]) for eName, p in entryPlots.items() }, precision=yieldPrecision)
        entries_smp.append(colEntries)
    if smp_data and smp_mc:
        sepStr_v += f"|{align}|"
        hdrs.append("Data/MC")
        colEntries = []
        for stData,stMC in zip(stTotData, stTotMC):
            dtCont = stData.contents
            mcCont = stMC.contents
            ratio = np.where(mcCont != 0., dtCont/mcCont, np.zeros(dtCont.shape))
            ratioErr = np.where(mcCont != 0., np.sqrt(mcCont**2*stData.sumw2 + dtCont**2*(stMC.sumw2+stMC.syst2))/mcCont**2, np.zeros(dtCont.shape))
            colEntries.append("${{0:.{0}f}} \pm {{1:.{0}f}}$".format(ratioPrecision).format( ratio[1], ratioErr[1]))
        entries_smp.append(colEntries)
    if len(colEntries) < 2:
        logger.warning("No samples, so no yields.tex")
    return "\n".join(([
        f"\\begin{{tabular}}{{ {sepStr_v} }}",
        "    \\hline",
        "    {0} \\\\".format(" & ".join(hdrs)),
        "    \\hline"]+[
            "    {0} \\\\".format(" ".join(smpEntries[i] for smpEntries in entries_smp))
            for i in range(len(report.titles)) ] )+[
        "    \\hline",
        "\\end{tabular}",
        "\\end{document}"
        ])

def printCutFlowReports(config, reportList, workdir=".", resultsdir=".", readCounters=lambda f : -1., eras=("all", None), verbose=False):
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
    if not eras: ## from config if not specified
        eras = list(config["eras"].keys())
    ## helper: print one bamboo.plots.CutFlowReport.Entry
    def printEntry(entry, printFun=logger.info, recursive=True, genEvents=None):
        effMsg = ""
        if entry.parent:
            sumPass = entry.nominal.GetBinContent(1)
            sumTotal = entry.parent.nominal.GetBinContent(1)
            if sumTotal != 0.:
                effMsg = f", Eff={sumPass/sumTotal:.2%}"
                if genEvents:
                    effMsg += f", TotalEff={sumPass/genEvents:.2%}"
        printFun(f"Selection {entry.name}: N={entry.nominal.GetEntries()}, SumW={entry.nominal.GetBinContent(1)}{effMsg}")
        if recursive:
            for c in entry.children:
                printEntry(c, printFun=printFun, recursive=recursive, genEvents=genEvents)
    ## retrieve results files, get generated events for each sample
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
            self.binnings = [EqB(1, 0.,1.)]
    for report in reportList:
        smpReports = { smp: report.readFromResults(resF) for smp, resF in resultsFiles.items() }
        ## debug print
        for smp, smpRep in smpReports.items():
            if smpRep.printInLog:
                logger.info(f"Cutflow report {report.name} for sample {smp}")
                for root in smpRep.rootEntries():
                    printEntry(root, genEvents=generated_events[smp])
        ## save yields.tex (if needed)
        if any(len(cb) > 1 or tt != cb[0] for tt,cb in report.titles.items()):
            if not has_plotit:
                logger.error(f"Could not load plotit python library, no TeX yields tables for {report.name}")
            else:
                yield_plots = [ YieldPlot(f"{report.name}_{eName}") for tEntries in report.titles.values() for eName in tEntries ]
                out_eras = []
                if len(eras) > 1 and eraMode in ("all", "combined"):
                    out_eras.append((f"{report.name}.tex", eras))
                if len(eras) == 1 or eraMode in ("split", "all"):
                    for era in eras:
                        out_eras.append((f"{report.name}_{era}.tex", [era]))
                for outName, iEras in out_eras:
                    pConfig, samples, plots, _, _ = loadPlotIt(config, yield_plots, eras=iEras, workdir=workdir, resultsdir=resultsdir, readCounters=readCounters)
                    tabBlock = _makeYieldsTexTable(report, samples,
                            { p.name[len(report.name)+1:]: p for p in plots },
                            stretch=pConfig.yields_table_stretch,
                            orientation=pConfig.yields_table_align,
                            align=pConfig.yields_table_text_align,
                            yieldPrecision=pConfig.yields_table_numerical_precision_yields,
                            ratioPrecision=pConfig.yields_table_numerical_precision_ratio)
                    with open(os.path.join(workdir, outName), "w") as ytf:
                        ytf.write("\n".join((_yieldsTexPreface, tabBlock)))
                    logger.info("Yields table for era(s) {0} was written to {1}".format(",".join(eras), os.path.join(workdir, outName)))

## END cutflow reports, adapted from bamboo.analysisutils

class CMSPhase2SimHistoModule(CMSPhase2SimModule, HistogramsModule):
    """ Base module for producing plots from Phase2 flat trees """
    def __init__(self, args):
        super(CMSPhase2SimHistoModule, self).__init__(args)
    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        """ Customised cutflow reports and plots """
        if not self.plotList:
            self.plotList = self.getPlotList(resultsdir=resultsdir)
        from bamboo.plots import Plot, DerivedPlot, CutFlowReport
        plotList_cutflowreport = [ ap for ap in self.plotList if isinstance(ap, CutFlowReport) ]
        plotList_plotIt = [ ap for ap in self.plotList if ( isinstance(ap, Plot) or isinstance(ap, DerivedPlot) ) and len(ap.binnings) == 1 ]
        eraMode, eras = self.args.eras
        if eras is None:
            eras = list(config["eras"].keys())
        if plotList_cutflowreport:
            printCutFlowReports(config, plotList_cutflowreport, workdir=workdir, resultsdir=resultsdir, readCounters=self.readCounters, eras=(eraMode, eras), verbose=self.args.verbose)
        if plotList_plotIt:
            from bamboo.analysisutils import writePlotIt, runPlotIt
            cfgName = os.path.join(workdir, "plots.yml")
            writePlotIt(config, plotList_plotIt, cfgName, eras=eras, workdir=workdir, resultsdir=resultsdir, readCounters=self.readCounters, vetoFileAttributes=self.__class__.CustomSampleAttributes, plotDefaults=self.plotDefaults)
            runPlotIt(cfgName, workdir=workdir, plotIt=self.args.plotIt, eras=(eraMode, eras), verbose=self.args.verbose)


################################
##                            ##
## The actual analysis module ##
##                            ##
################################

class CMSPhase2SimTest(CMSPhase2SimHistoModule):
    """ Plotter module for Phase2 flat trees """
    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot, CutFlowReport
        from bamboo.plots import EquidistantBinning as EqB
        from bamboo import treefunctions as op

        plots = []

        #definitions

        electrons = op.select(t.elec, lambda el : op.AND(
        el.pt > 20., op.abs(el.eta) < 2.5
        ))
        
        muons = op.select(t.muon, lambda mu : op.AND(
        mu.pt > 20., op.abs(mu.eta) < 2.5
        ))
        
        cleanedElectrons = op.select(electrons, lambda el : op.NOT(
        op.rng_any(muons, lambda mu : op.deltaR(el.p4, mu.p4) < 0.3 )
        ))

        # we are taking the second isopass to be on which is equal to the medium working point
        isolatedElectrons = op.select(cleanedElectrons, lambda el : el.isopass & (1<<2) )
        
        identifiedElectrons = op.select(isolatedElectrons, lambda el : el.idpass & (1<<2) )
        
        cleanedMuons = op.select(muons, lambda mu : op.NOT(
        op.rng_any(electrons, lambda el : op.deltaR(mu.p4, el.p4) < 0.3 )
        ))
        
        isolatedMuons = op.select(cleanedMuons, lambda mu : mu.isopass & (1<<2) )
        
        identifiedMuons = op.select(isolatedMuons, lambda mu : mu.idpass & (1<<2) )
        
        InvMassMuMU = op.invariant_mass(identifiedMuons[0].p4, identifiedMuons[1].p4 )
        
        cleanedJets = op.select(t.jetpuppi, lambda j : op.AND(
        op.NOT(op.rng_any(identifiedElectrons, lambda el : op.deltaR(el.p4, j.p4) < 0.3) ),
        op.NOT(op.rng_any(identifiedMuons, lambda mu : op.deltaR(mu.p4, j.p4) < 0.3) )
        ))

        cleanedGoodJets = op.select(cleanedJets, lambda j : op.AND(
        j.pt > 30, op.abs(j.eta) < 2.5
        ))

        btaggedJets = op.select(cleanedGoodJets, lambda j : j.btag & (1<<2))

        met = op.select(t.metpuppi)

        #selections

        #selection1 : Oppositely charged MuMu selection
        sel1 = noSel.refine("nmumu", cut = [op.AND(
            (op.rng_len(identifiedMuons) > 1), (op.product(identifiedMuons[0].charge, identifiedMuons[1].charge) < 0 ))]) 

        #selection2 : Invariant mass selection
        sel2 = sel1.refine("InvM", cut = [op.NOT(op.in_range(76, InvMassMuMU, 106))])

        #selection3 : two jets selection
        sel3 = sel2.refine("njet", cut = [op.rng_len(cleanedGoodJets) > 1])
    
        #selection4 : at least 1 among two leading jets is b-tagged
        sel4 = sel3.refine("btag", cut = [op.OR(
            cleanedGoodJets[0].btag & (1<<2), cleanedGoodJets[1].btag & (1<<2))])

        #selection5 : MET > 40 GeV
        sel5 = sel4.refine("MET", cut = [met[0].pt > 40])

        #plots
            
            #noSel
        plots.append(Plot.make1D("nJetsNoSel", op.rng_len(cleanedGoodJets), noSel, EqB(10, 0., 10.), title="nJets"))

        plots.append(Plot.make1D("nbtaggedJetsNoSel", op.rng_len(btaggedJets), noSel, EqB(10, 0., 10.), title="nbtaggedJets"))
  
        plots.append(Plot.make1D("nMuNoSel", op.rng_len(identifiedMuons), noSel, EqB(15, 0., 15.), title="nMuons"))
        
        plots.append(Plot.make1D("METptNoSel", met[0].pt, noSel, EqB(50, 0., 250), title="MET_PT"))
        
            #sel1
        
        plots.append(Plot.make1D("nJetsSel1", op.rng_len(cleanedGoodJets), sel1, EqB(10, 0., 10.), title="nJets"))

        plots.append(Plot.make1D("nbtaggedJetsSel1", op.rng_len(btaggedJets), sel1, EqB(10, 0., 10.), title="nbtaggedJets"))

        plots.append(Plot.make1D("nMuSel1", op.rng_len(identifiedMuons), sel1, EqB(10, 0., 10.), title="nMuons"))

        plots.append(Plot.make1D("InvMassTwoMuonsSel1", InvMassMuMU, sel1, EqB(30, 0, 300), title="m(ll)"))
                  
        plots.append(Plot.make1D("LeadingMuonPTSel1", muons[0].pt, sel1, EqB(30, 0., 250.), title=" Leading Muon PT"))
            
        plots.append(Plot.make1D("SubLeadingMuonPTSel1", muons[1].pt, sel1, EqB(30, 0., 250.), title="SubLeading Muon PT"))

        plots.append(Plot.make1D("LeadingMuonEtaSel1", muons[0].eta, sel1, EqB(30, -3, 3), title=" Leading Muon eta"))
            
        plots.append(Plot.make1D("SubLeadingMuonEtaSel1", muons[1].eta, sel1, EqB(30, -3, 3), title="SubLeading Muon eta"))
        
        plots.append(Plot.make1D("METptSel1", met[0].pt, sel1, EqB(50, 0., 250), title="MET_PT"))
                            
            #sel2

        plots.append(Plot.make1D("nJetsSel2", op.rng_len(cleanedGoodJets), sel2, EqB(10, 0., 10.), title="nJets"))

        plots.append(Plot.make1D("nbtaggedJetsSel2", op.rng_len(btaggedJets), sel2, EqB(10, 0., 10.), title="nbtaggedJets"))

        plots.append(Plot.make1D("nMuSel2", op.rng_len(identifiedMuons), sel2, EqB(10, 0., 10.), title="nMuons"))

        plots.append(Plot.make1D("InvMassTwoMuonsSel2", InvMassMuMU, sel2, EqB(20, 20., 300.), title="m(ll)"))
        
        plots.append(Plot.make1D("LeadingMuonPTSel2", muons[0].pt, sel2, EqB(30, 0., 250.), title=" Leading Muon PT"))

        plots.append(Plot.make1D("SubLeadingMuonPTSel2", muons[1].pt, sel2, EqB(30, 0., 200.), title=" SubLeading Muon PT"))

        plots.append(Plot.make1D("LeadingMuonEtaSel2", muons[0].eta, sel2, EqB(30, -3, 3), title=" Leading Muon Eta"))
            
        plots.append(Plot.make1D("SubLeadingMuonEtaSel2", muons[1].eta, sel2, EqB(30, -3, 3), title=" SubLeading Muon Eta"))

        plots.append(Plot.make1D("METptSel2", met[0].pt, sel2, EqB(50, 0., 250), title="MET_PT"))
                        
            #sel3
            
        plots.append(Plot.make1D("nJetsSel3", op.rng_len(cleanedGoodJets), sel3, EqB(10, 0., 10.), title="nJets"))

        plots.append(Plot.make1D("nbtaggedJetsSel3", op.rng_len(btaggedJets), sel3, EqB(10, 0., 10.), title="nbtaggedJets"))

        plots.append(Plot.make1D("LeadingJetPTSel3", cleanedGoodJets[0].pt, sel3, EqB(50, 0., 350.), title="Leading jet PT"))
            
        plots.append(Plot.make1D("SubLeadingJetPTSel3", cleanedGoodJets[1].pt, sel3, EqB(50, 0., 350.), title="SubLeading jet PT"))
            
        plots.append(Plot.make1D("LeadingJetEtaSel3", cleanedGoodJets[0].eta, sel3, EqB(30, -3, 3), title="Leading jet Eta"))
            
        plots.append(Plot.make1D("SubLeadingJetEtaSel3", cleanedGoodJets[1].eta, sel3, EqB(30, -3, 3), title="SubLeading jet Eta"))

        plots.append(Plot.make1D("nMuSel3", op.rng_len(identifiedMuons), sel3, EqB(10, 0., 10.), title="nMuons"))
        
        plots.append(Plot.make1D("LeadingMuonPTSel3", muons[0].pt, sel3, EqB(30, 0., 250.), title=" Leading Muon PT"))
            
        plots.append(Plot.make1D("SubLeadingMuonPTSel3", muons[1].pt, sel3, EqB(30, 0., 200.), title=" SubLeading Muon PT"))

        plots.append(Plot.make1D("LeadingMuonEtaSel3", muons[0].eta, sel3, EqB(30, -3, 3), title=" Leading Muon Eta"))
            
        plots.append(Plot.make1D("SubLeadingMuonEtaSel3", muons[1].eta, sel3, EqB(30, -3, 3), title=" SubLeading Muon Eta"))
                
        plots.append(Plot.make1D("InvMassTwoMuonsSel3", InvMassMuMU, sel3, EqB(30, 0, 300), title="m(ll)"))
        
        plots.append(Plot.make1D("METptSel3", met[0].pt, sel3, EqB(50, 0., 250), title="MET_PT"))
        
            #sel4
             
        plots.append(Plot.make1D("nJetsSel4", op.rng_len(cleanedGoodJets), sel4, EqB(10, 0, 10), title="nJets"))

        plots.append(Plot.make1D("nbtaggedJetsSel4", op.rng_len(btaggedJets), sel4, EqB(10, 0., 10.), title="nbtaggedJets"))

        plots.append(Plot.make1D("LeadingJetPTSel4", cleanedGoodJets[0].pt, sel4, EqB(50, 0., 250.), title="Leading jet PT"))
            
        plots.append(Plot.make1D("SubLeadingJetPTSel4", cleanedGoodJets[1].pt, sel4, EqB(50, 0., 250.), title="SubLeading jet PT"))
            
        plots.append(Plot.make1D("LeadingJetEtaSel4", cleanedGoodJets[0].eta, sel4, EqB(30, -3, 3.), title="Leading jet Eta"))
            
        plots.append(Plot.make1D("SubLeadingJetEtaSel4", cleanedGoodJets[1].eta, sel4, EqB(30, -3, 3.), title="SubLeading jet Eta"))
        
        plots.append(Plot.make1D("nMuSel4", op.rng_len(identifiedMuons), sel4, EqB(10, 0., 10.), title="nMuons"))
                
        plots.append(Plot.make1D("LeadingMuonPTSel4", muons[0].pt, sel4, EqB(30, 0., 250.), title=" Leading Muon PT"))
        
        plots.append(Plot.make1D("SubLeadingMuonPTSel4", muons[1].pt, sel4, EqB(30, 0., 200.), title=" SubLeading Muon PT"))
        
        plots.append(Plot.make1D("LeadingMuonEtaSel4", muons[0].eta, sel4, EqB(30, -3, 3), title=" Leading Muon Eta"))
        
        plots.append(Plot.make1D("SubLeadingMuonEtaSel4", muons[1].eta, sel4, EqB(30, -3, 3), title=" SubLeading Muon Eta"))

        plots.append(Plot.make1D("InvMassTwoMuonsSel4", InvMassMuMU, sel4, EqB(30, 0, 300), title="m(ll)"))
        
        plots.append(Plot.make1D("METptSel4", met[0].pt, sel4, EqB(50, 0., 250), title="MET_PT"))

            #sel5
                
        plots.append(Plot.make1D("nJetsSel5", op.rng_len(cleanedGoodJets), sel5, EqB(10, 0, 10), title="nJets"))

        plots.append(Plot.make1D("nbtaggedJetsSel5", op.rng_len(btaggedJets), sel5, EqB(10, 0., 10.), title="nbtaggedJets"))

        plots.append(Plot.make1D("LeadingJetPTSel5", cleanedGoodJets[0].pt, sel5, EqB(50, 0., 250.), title="Leading jet PT"))
            
        plots.append(Plot.make1D("SubLeadingJetPTSel5", cleanedGoodJets[1].pt, sel5, EqB(50, 0., 250.), title="SubLeading jet PT"))
            
        plots.append(Plot.make1D("LeadingJetEtaSel5", cleanedGoodJets[0].eta, sel5, EqB(30, -3, 3.), title="Leading jet Eta"))
            
        plots.append(Plot.make1D("SubLeadingJetEtaSel5", cleanedGoodJets[1].eta, sel5, EqB(30, -3, 3.), title="SubLeading jet Eta"))
        
        plots.append(Plot.make1D("nMuSel5", op.rng_len(identifiedMuons), sel5, EqB(10, 0., 10.), title="nMuons"))
                
        plots.append(Plot.make1D("LeadingMuonPTSel5", muons[0].pt, sel5, EqB(30, 0., 250.), title=" Leading Muon PT"))
        
        plots.append(Plot.make1D("SubLeadingMuonPTSel5", muons[1].pt, sel5, EqB(30, 0., 200.), title=" SubLeading Muon PT"))
        
        plots.append(Plot.make1D("LeadingMuonEtaSel5", muons[0].eta, sel5, EqB(30, -3, 3), title=" Leading Muon Eta"))
        
        plots.append(Plot.make1D("SubLeadingMuonEtaSel5", muons[1].eta, sel5, EqB(30, -3, 3), title=" SubLeading Muon Eta"))

        plots.append(Plot.make1D("InvMassTwoMuonsSel5", InvMassMuMU, sel5, EqB(30, 0, 300), title="m(ll)"))

        plots.append(Plot.make1D("METptSel5", met[0].pt, sel5, EqB(50, 0., 250), title="MET_PT > 40"))

        # Efficiency Report on terminal and the .tex output

        cfr = CutFlowReport("yields")
        cfr.add(noSel, "Sel0: No selection")
        cfr.add(sel1, "Sel1: nMuMu >= 2")
        cfr.add(sel2, "Sel2: InvM")
        cfr.add(sel3, "Sel3: nJet >= 2")
        cfr.add(sel4, "Sel4: btag")
        cfr.add(sel5, "Sel5: MET")

        plots.append(cfr)
                            
        return plots
