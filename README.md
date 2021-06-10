# Bamboo

Bamboo is a python framework for HEP analyses. In this page, its explanation is mostly intended to be used in specific analyses done in Reconstruction Tools & Benchmarks Group under UPSG-CMS. For the original repo please check [https://gitlab.cern.ch/cp3-cms/bamboo](https://gitlab.cern.ch/cp3-cms/bamboo)

Documentation: [https://cp3.irmp.ucl.ac.be/~pdavid/bamboo/index.html](https://cp3.irmp.ucl.ac.be/~pdavid/bamboo/index.html)

### Setup instructions

The following set of commands would be sufficient to configure bamboo on lxplus.
```bash
mkdir BambooAnalysis
cd BambooAnalysis
```

Make a virtualenv:
```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc10-opt/setup.sh
python -m venv bamboovenv
source bamboovenv/bin/activate
```

Install bamboo and the [python plotit](https://gitlab.cern.ch/cp3-cms/pyplotit) module,
which is required for the .tex output:
```bash
pip install git+https://gitlab.cern.ch/cp3-cms/bamboo.git
pip install git+https://gitlab.cern.ch/cp3-cms/pyplotit.git
```

Clone and install plotIt, the C++ executable that will produce the stack plots:
```bash
git clone -o upstream https://github.com/cp3-llbb/plotIt.git
mkdir build-plotit
cd build-plotit
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV ../plotIt
make -j4 install
cd -
```

Clone this repository:
```bash
git clone -o upstream https://github.com/recotoolsbenchmarks/BambooAnalysis.git rtb
```

### Starting to use bamboo

You need to execute the following commands in the bamboo setup directory, every time logged in and want to work with bamboo:
```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc10-opt/setup.sh
source bamboovenv/bin/activate
```

### Analysis module

The analysis module of the bamboo is a python file. See an example here: /phaseII-analysis.py.
In the analysis module, the physics objects are retreived from ROOT trees and defined via simple cuts. The selections are applied and finally the plots are defined. One can also add a cutflow report on demand which prints selection efficiencies in a .tex file as well as in the terminal screen. This cutflow report comes along with a histogram in .png format available in the working directory.

### Analysis configuration file (YAML)

A YAML file is needed to specify the samples to be used in the analysis. In a configuration file, the luminosity, x-sec, sample groups and plot configurations (e.g. stacked plots, shape comparison etc.) are defined. See examples in this directory and check [documentation](https://cp3.irmp.ucl.ac.be/~pdavid/bamboo/index.html) for further details.

### RTB analysis

DelphesFlat
```bash
bambooRun -m rtb/phaseII-analysis.py:CMSPhase2SimTest rtb/phaseII-analysis-FS-Delphes.yml -o output-Delphes
```
Fullsim
```bash
bambooRun -m rtb/phaseII-analysis.py:CMSPhase2SimTest rtb/phaseII-analysis-FS.yml -o output-FS
```
Delphes vs. FS
```bash
bambooRun -m rtb/phaseII-analysis.py:CMSPhase2SimTest rtb/phaseII-analysis-FS-Delphes.yml -o output-FS-Delphes
```

The pdf version of the cutflow report is obtained via:
```bash
cd output-FS-Delphes
pdflatex yields_HL-LHC.tex
```

### Outputs

See examples: [http://aguzel.web.cern.ch/aguzel/rtb/bamboo/](http://aguzel.web.cern.ch/aguzel/rtb/bamboo/)
