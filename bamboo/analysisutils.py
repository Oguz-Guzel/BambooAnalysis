""" The :py:mod:`bamboo.analysisutils` module bundles a number of more
specific helper methods that use the tree decorators and integrate with
other components, connect to external services, or are factored out of the
classes in :py:mod:`bamboo.analysismodules` to facilitate reuse.
"""
import copy
import logging
logger = logging.getLogger(__name__)
import os.path
import subprocess
import urllib.parse
import yaml

bamboo_cachedir = os.path.join(os.getenv("XDG_CACHE_HOME", os.path.join(os.path.expanduser("~"), ".cache")), "bamboo")

def addLumiMask(sel, jsonName, runRange=None, runAndLS=None, name="goodlumis"):
    from . import treefunctions as op
    """ Refine selection with a luminosity block filter

    Typically applied directly to the root selection (for data).
    runAndLS should be a tuple of expressions with the run number and luminosity block ID.
    The run range is used to limit the part of the JSON file to consider,
    see the LumiMask helper class for details.
    """
    if runAndLS is None:
        raise RuntimeError("Cannot construct a filter for the good lumi blocks without accessors (backend.create(..., runAndLS=XXX)), tree->(run, LS)")
    lumiSel = op.define("LumiMask", 'const auto <<name>> = LumiMask::fromJSON("{0}"{1});'.format(
                jsonName, (", {0:d}, {1:d}".format(*runRange) if runRange is not None else "")), nameHint="bamboo_lumiMask")
    return sel.refine(name, cut=lumiSel.accept(*runAndLS))

def downloadCertifiedLumiFiles(taskArgs, workdir="."):
    """ download certified lumi files (if needed) and replace in args """
    taskArgs = copy.deepcopy(taskArgs)
    certifLumiFiles = set(kwargs["certifiedLumiFile"] for args,kwargs in taskArgs if "certifiedLumiFile" in kwargs)
    ## download if needed
    clf_downloaded = dict()
    for clfu in certifLumiFiles:
        purl = urllib.parse.urlparse(clfu)
        if purl.scheme in ("http", "https"):
            fname = os.path.join(workdir, purl.path.split("/")[-1])
            if os.path.exists(fname):
                logger.warning("File {0} exists, it will not be downloaded again from {1}".format(fname, clfu))
            else:
                subprocess.check_call(["wget", "--directory-prefix={0}".format(workdir), clfu], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            clf_downloaded[clfu] = os.path.abspath(fname)
    ## update args
    for args,kwargs in taskArgs:
        if "certifiedLumiFile" in kwargs:
            clf = kwargs["certifiedLumiFile"]
            if clf in clf_downloaded:
                kwargs["certifiedLumiFile"] = clf_downloaded[clf]

    return taskArgs, set(clf_downloaded.keys())

def parseAnalysisConfig(anaCfgName, redodbqueries=False, overwritesamplefilelists=False, envConfig=None):
    cfgDir = os.path.dirname(os.path.abspath(anaCfgName))
    with open(anaCfgName) as anaCfgF:
        analysisCfg = yaml.load(anaCfgF)
    ## finish loading samples (file lists)
    samples = dict()
    for smpName, smpCfg in analysisCfg["samples"].items():
        smp = copy.deepcopy(smpCfg)
        ## read cache, if it's there
        listfile, cachelist = None, []
        if "files" in smpCfg and str(smpCfg["files"]) == smpCfg["files"]:
            listfile = smpCfg["files"] if os.path.isabs(smpCfg["files"]) else os.path.join(cfgDir, smpCfg["files"])
            if os.path.isfile(listfile):
                with open(listfile) as smpF:
                    cachelist = [ fn for fn in [ ln.strip() for ln in smpF ] if len(fn) > 0 ]

        if "db" in smpCfg and ( "files" not in smpCfg or len(cachelist) == 0 or redodbqueries ):
            files = []
            for dbEntry in (smpCfg["db"] if str(smpCfg["db"]) != smpCfg["db"] else [smpCfg["db"]]): ## convert to list if string
                if ":" not in dbEntry:
                    raise RuntimeError("'db' entry should be of the format 'protocol:location', e.g. 'das:/SingleMuon/Run2016E-03Feb2017-v1/MINIAOD'")
                protocol, dbLoc = dbEntry.split(":")
                if protocol == "das":
                    dasConfig = envConfig["das"]
                    dasQuery = "file dataset={0}".format(dbLoc)
                    entryFiles = [ os.path.join(dasConfig["storageroot"], fn.lstrip("/")) for fn in [ ln.strip() for ln in subprocess.check_output(["dasgoclient", "-query", dasQuery]).decode().split() ] if len(fn) > 0 ]
                    files += entryFiles
                    if len(entryFiles) == 0:
                        raise RuntimeError("No files found with DAS query {0}".format(dasQuery))
                    ## TODO improve: check that files are locally available, possibly fall back to xrootd otherwise; check for grid proxy before querying; maybe do queries in parallel
                elif protocol == "samadhi":
                    logger.warning("SAMADhi queries are not implemented yet")
                else:
                    raise RuntimeError("Unsupported protocol in '{0}': {1}".format(dbEntry, protocol))
            smp["files"] = files
            if listfile and ( len(cachelist) == 0 or overwritesamplefilelists ):
                with open(listfile, "w") as listF:
                    listF.write("\n".join(files))
        elif "files" not in smpCfg:
            raise RuntimeError("Cannot load files for {0}: neither 'db' nor 'files' specified".format(smpName))
        elif listfile:
            if len(cachelist) == 0:
                raise RuntimeError("No file names read from {0}".format())
            smp["files"] = cachelist
        else: ## list in yml
            smp["files"] = [ (fn if os.path.isabs(fn) or urllib.parse.urlparse(fn).scheme != "" in fn else os.path.join(cfgDir, fn)) for fn in smpCfg["files"] ]
        samples[smpName] = smp
    analysisCfg["samples"] = samples
    return analysisCfg

def readEnvConfig(explName=None):
    """ Read computing environment config file (batch system, storage site etc.)

    For using a batch cluster, the [batch] section should have a 'backend' key,
    and there should be a section with the name of the backend (slurm, htcondor...),
    see bamboo.batch_<backend> for details.
    The storage site information needed to resolve the PFNs for datasets retrieved from DAS
    should be specified under the [das] section (sitename and storageroot).
    """
    import os
    from configparser import ConfigParser
    def readFromFile(name):
        cfgp = ConfigParser()
        cfgp.read(name)
        cfg = dict((sName, dict(cfgp[sName])) for sName in cfgp.sections())
        return cfg

    xdgCfg = os.getenv("XDG_CONFIG_HOME", os.path.join(os.path.expanduser("~"), ".config"))
    toTry = ["bamboo.ini", "bamboorc", os.path.join(xdgCfg, "bamboorc")]
    if explName:
        toTry.insert(0, explName)
    for iniName in toTry:
        if os.path.exists(iniName):
            try:
                res = readFromFile(iniName)
                logger.info("Read config from file {0}".format(iniName))
                return res
            except Exception as ex:
                logger.warning("Problem reading config file {0}: {1}".format(iniName, ex))
    raise RuntimeError("No valid config file found")

plotit_plotdefaults = {
        "x-axis"           : lambda p : "{0}".format(p.axisTitles[0]),
        "y-axis"           : "Evt",
        "y-axis-format"    : "%1% / %2$.0f",
        "normalized"       : False,
        "x-axis-format"    : lambda p : [p.binnings[0].minimum, p.binnings[0].maximum],
        "log-y"            : "both",
        "y-axis-show-zero" : True,
        "save-extensions"  : ["pdf"],
        "show-ratio"       : True,
        "sort-by-yields"   : False,
        }
def runPlotIt(config, plotList, workdir=".", resultsdir=".", plotIt="plotIt", plotDefaults=None, readCounters=lambda f : -1., era=None):
    ## TODO also pass the correct luminosity numbers
    eras = []
    if era is not None:
        if str(era) == era:
            eras.append(era)
        else:
            eras = eras["eras"]
    plotitCfg = (copy.deepcopy(config["plotIt"]) if "plotIt" in config else dict())
    plotitCfg["configuration"]["root"] = os.path.relpath(resultsdir, workdir)
    plotit_files = dict()
    for smpN, smpCfg in config["samples"].items():
        if smpCfg.get("era") in eras:
            resultsName = "{0}.root".format(smpN)
            smpOpts = dict()
            smpOpts["group"] = smpCfg["group"]
            isMC = ( smpCfg["group"] != "data" )
            smpOpts["type"] = ("mc" if isMC else "data")
            if isMC:
                smpOpts["cross-section"] = smpCfg["cross-section"]
                import ROOT
                resultsFile = ROOT.TFile.Open(os.path.join(resultsdir, resultsName))
                counters = readCounters(resultsFile)
                smpOpts["generated-events"] = counters[smpCfg["generated-events"]]
            plotit_files[resultsName] = smpOpts
    plotitCfg["files"] = plotit_files
    plotit_plots = dict()
    for plot in plotList:
        plotOpts = dict(plotit_plotdefaults)
        if plotDefaults:
            plotOpts.update(plotDefaults)
        plotOpts.update(plot.plotopts)
        plotOpts = dict((k, (v(plot) if hasattr(v, "__call__") else v)) for k,v in plotOpts.items())
        plotit_plots[plot.name] = plotOpts
    plotitCfg["plots"] = plotit_plots
    cfgName = os.path.join(workdir, "plots.yml")
    with open(cfgName, "w") as plotitFile:
        yaml.dump(plotitCfg, plotitFile)

    plotsdir = os.path.join(workdir, "plots")
    if os.path.exists(plotsdir):
        logger.warning("Directory '{0}' already exists, previous plots will be overwritten".format(plotsdir))
    else:
        os.makedirs(plotsdir)

    try:
        with open(os.path.join(plotsdir, "out.log"), "w") as logFile:
            subprocess.check_call([plotIt, "-i", workdir, "-o", plotsdir, cfgName], stdout=logFile)
        logger.info("plotIt output is available in {0}".format(plotsdir))
    except subprocess.CalledProcessError as ex:
        logger.error("Command '{0}' failed with exit code {1}\n{2}".format(" ".join(ex.cmd), ex.returncode, ex.output))

def configureJets(calc, jetType, jec=None, jecLevels="default", smear=None, useGenMatch=True, genMatchDR=0.2, genMatchDPt=3., jesUncertaintySources=None, cachedir=None):
    """ Reapply JEC, set up jet smearing, or prepare JER/JES uncertainties collections

    :param calc: jet variations calculator to configure (e.g. ``Jet.calc``)
    :param jetType: jet type, e.g. AK4PFchs
    :param smear: tag of resolution (and scalefactors) to use for smearing (no smearing is done if unspecified)
    :param jec: tag of the new JEC to apply, or for the JES uncertainties (pass an empty list to jecLevels to produce only the latter without reapplying the JEC)
    :param jesUncertaintySources: list of jet energy scale uncertainty sources (see `the JECUncertaintySources twiki page <https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECUncertaintySources>`_)

    :param useGenMatch: use matching to generator-level jets for resolution smearing
    :param genMatchDR: DeltaR for generator-level jet matching (half the cone size is recommended, default is 0.2)
    :param genMatchDPt: maximal relative PT difference (in units of the resolution) between reco and gen jet
    :param jecLevels: list of JEC levels to apply (if left out the recommendations are used: L1FastJet, L2Relative, L3Absolute, and also L2L3Residual for data)
    :param cachedir: alternative root directory to use for the txt files cache, instead of ``$XDG_CACHE_HOME/bamboo`` (usually ``~/.cache/bamboo``)
    """
    from .jetdatabasecache import JetDatabaseCache
    if smear is not None:
        with JetDatabaseCache("JRDatabase", repository="cms-jet/JRDatabase", cachedir=cachedir) as jrDBCache:
            mcPTRes = jrDBCache.getPayload(smear, "PtResolution", jetType)
            mcResSF = jrDBCache.getPayload(smear, "SF", jetType)
        calc.setSmearing(mcPTRes, mcResSF, useGenMatch, genMatchDR, genMatchDPt)
    if jec is not None:
        if jecLevels == "default":
            # "L3Absolute" left out because it is dummy according to https://twiki.cern.ch/twiki/bin/view/CMS/IntroToJEC#Mandatory_Jet_Energy_Corrections
            if jec.endswith("_DATA"):
                jecLevels = ["L1FastJet", "L2Relative", "L2L3Residual"]
            elif jec.endswith("_MC"):
                # "L2L3Residual" could be added, but it is dummy for MC according to https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetCorApplication
                jecLevels = ["L1FastJet", "L2Relative"]
            else:
                raise ValueError("JEC tag {0} does not end with '_DATA' or '_MC', so the levels cannot be guessed. Please specify the JEC levels explicitly")
        with JetDatabaseCache("JECDatabase", repository="cms-jet/JECDatabase", cachedir=cachedir) as jecDBCache:
            from cppyy import gbl
            if jecLevels:
                jecParams = getattr(gbl, "std::vector<JetCorrectorParameters>")()
                for jLev in jecLevels:
                    plf = jecDBCache.getPayload(jec, jLev, jetType)
                    params = gbl.JetCorrectorParameters(plf)
                    jecParams.push_back(params)
                calc.setJEC(jecParams)
            if jesUncertaintySources:
                plf = jecDBCache.getPayload(jec, "UncertaintySources", jetType)
                for src in jesUncertaintySources:
                    params = gbl.JetCorrectorParameters(plf, src)
                    calc.addJESUncertainty(src, params)

def forceDefine(arg, selection):
    """ Force the definition of an expression as a column at a selection stage

    Use only for really computation-intensive operations that need to be precalculated

    :param arg: expression to define as a column
    :param selection: selection for which the expression should be defined
    """
    from .treeoperations import adaptArg
    return selection._fbe.selDFs[selection.name].define(adaptArg(arg))
