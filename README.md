# HiggsExample11-12
This example guides the user to reproducing the discovery of the Higgs boson using the 2011 and 2012 datasets, in the four-lepton final states. It contains multiple levels of examples, from very simple to a full analysis, all with CMS Open Data.

# Introduction

authors: N.Z. Jomhari, A. Geiser, A. Anuar 

This research level example is a strongly simplified reimplementation of 
parts of the original CMS Higgs-to-four-lepton analysis published in
[Phys.Lett. B716 (2012) 30-61,  arXiv:1207.7235, 
https://inspirehep.net/record/1124338?ln=en](https://inspirehep.net/record/1124338?ln=en)

The published reference plot which is being approximated in this example is  
https://inspirehep.net/record/1124338/files/H4l_mass_v3.png
Other Higgs final states (e.g. Higgs to two photons), which were also part of 
the same CMS paper and strongly contributed to the Higgs discovery, are not 
covered by this example. 

The highest level, of this example addresses users who feel they have at least 
some minimal understanding of the content of this paper and of the meaning 
of this reference plot, which can be reached via (separate) educational 
exercises. 
The lower levels might also be interesting for educational applications.
The example requires a minimal acquaintance with the linux operating system 
and the Root analysis tool, which can also be obtained from corresponding 
(separate) tutorials.

The example uses legacy versions of the original CMS data sets in the 
CMS AOD format, which slightly differ from the ones used for the 
publication due to improved calibrations. It also uses legacy versions of the 
corresponding Monte Carlo simulations, which are again close to, but not 
identical to, the ones in the original publication. These legacy data and MC 
sets listed below were used in practice, exactly as they are, in many later CMS
publications.

Since according to the CMS Open Data policy the fraction of data which are 
public (and used here) is only 50% of the available LHC Run I samples, 
the statistical significance is reduced with respect to what can be 
achieved with the full dataset. However, the original paper 
[Phys.Lett. B716 (2012) 30-61,  arXiv:1207.7235, 
https://inspirehep.net/record/1124338?ln=en](https://inspirehep.net/record/1124338?ln=en)
was also obtained with only part of the Run I statistics, roughly
equivalent to the luminosity of the public set, but with only partial 
statistical overlap.

The provided analysis code recodes the spirit of the original 
analysis and recodes many of the original cuts on original data objects, 
but does not provide the original analysis code itself. Also, for the sake of
simplicity, it skips some of the more advanced analysis methods of the 
original paper. Nevertheless, it provides a qualitative insight about how the 
original result was obtained. In addition to the documented core results, 
the resulting Root files also contain many undocumented plots which grew 
as a side product from setting up this example and earlier examples.
The significance of the Higgs "excess" is about 2 standard deviations in 
this example, while it was 3.2 standard deviations in this channel alone 
in the original publication. The difference is attributed to the less 
sophisticated background suppression.
In more recent (not yet public) CMS data sets with higher statistics the 
signal is observed in an analysis with more than 5 standard 
deviations in this channel alone [https://inspirehep.net/record/1518828]. 
The most recent results on the search and observation of the Higgs boson 
in this decay mode are available here: http://cms-results.web.cern.ch/cms-results/public-results/publications/HIG/ZZ.html


## Strategy 

The analysis strategy is the following: Get the 4mu and 2mu2e final states 
from the DoubleMuParked datasets and the 4e final state from the 
DoubleElectron dataset. This avoids double counting due to trigger overlaps.
All MC contributions except top use data-driven normalization:
The DY (Z/gamma^*) contribution is scaled to the Z peak.
The ZZ contribution is scaled to describe the data in the independent 
mass range 180-600 GeV.
The Higgs contribution is scaled to describe the data in the signal region.
The (very small) top contribution remains scaled to the MC generator cross 
section.     

# Detailed instructions:
There are four levels of increasing complexity for this example:

## Level 1 (easy): Compare

1. *Compare* the provided final output plot 
   [mass4l_combine.pdf](https://raw.githubusercontent.com/cms-opendata-analyses/HiggsExample11-12/master/mass4l_combine.pdf)
   or 
   [mass4l_combine.png](https://raw.githubusercontent.com/cms-opendata-analyses/HiggsExample11-12/master/mass4l_combine.png)
   to the published one, 
   https://inspirehep.net/record/1124338/files/H4l_mass_v3.png
   keeping in mind the caveats given above 

## Level 2 (beginner): Reproduce

2. *Reproduce* the final output plot from the predefined histogram files 
   using a ROOT macro 
   (Time: ~few minutes - ~few hours, depending on setup and proficiency)
   There are two ways to do this tutorial, with either a local software or virtual installation of the ROOT framework,  very commonly used in particle physics. 
   
   If not already proficient in ROOT, consider doing a brief ROOT 
   introductory tutorial [https://root.cern.ch/introductory-tutorials] in order to understand what the ROOT 
   macro will do 

   ### With local ROOT installation
   -> - if a ROOT version on the local computer compatible with ROOT 5.32/00 
        is running, use the local version. This avoids installation of the CERN VM.
   
   ### Without a local ROOT installation
        If you do not have ROOT installed, it is also possible to run this part of the tutorial in the CERN Virtual Machine. 
        Follow the instructions of the 2nd item of level 3 
        in order to install the Virtualbox and CERNVM and run ROOT there
        (the validation step with the Demo example might be skipped) 

   ### Downloading the code
      - download the full git repository, in that case you get the directory structure automatically 

        ```
        git clone https://github.com/cms-opendata-analyses/HiggsExample11-12.git
        cd HiggsExample11-12/examplerootfiles/
        ```
        
      - OR create a new directory, e.g. examplerootfiles 
      
        ```      
        mkdir examplerootfiles
        ```
      
        switch to that directory 
        
        ```
        cd examplerootfiles
        ```
        
        and download all the preproduced *.root files given in 
        [https://github.com/cms-opendata-analyses/HiggsExample11-12/tree/master/examplerootfiles] for all relevant samples to this directory 
        and download [M4Lnormdatall.cc](https://raw.githubusercontent.com/cms-opendata-analyses/HiggsExample11-12/master/examplerootfiles/M4Lnormdatall.cc), the ROOT macro to analyse these files and keep it in the same working directory as the *.root files
      
         ```
         M4Lnormdatall.cc 
         ```
         
        from this record into the same directory
      - on the linux prompt, type 
      
         ```
         root -l M4Lnormdatall.cc
         ```
      
      -> you will get the output plot on the screen
      - either, on the ROOT canvas (picture) click 
      
         ```
         file->Quit ROOT
         ```
         
        or, on the root [] prompt, type the command to end your root session:
        
         ```
         .q     
         ```
         
        -> you will exit ROOT and find the output plot in 
           mass4l_combined_user.pdf
           
      - you can compare this plot with [the plots provided in 1.](https://github.com/cms-opendata-analyses/HiggsExample11-12/blob/master/README.md#level-1-easy-compare)
<!----
## Level 3 (intermediate): Produce

3. *Produce* a ROOT data input file from original data and MC files for one 
   Higgs signal candidate and for the simulated Higgs signal with reduced 
   statistics (for speed reasons) and reproduce the final output plot 
   containing your own input using a ROOT macro 
   (~few minutes to ~1 hour if Virtual machine is already installed, 
     depending on internet connection and computer performance, up to 
     ~few hours otherwise) 
   -> - if not already done follow instructions in 
        CMS 2011 Virtual Machines: How to install  
        [http://opendata.web.cern.ch/docs/cms-virtual-machine-2011]
        * install VirtualBox
        * install CERNVM Virtual Machine
      - in particular, if not yet done, run the Demo program from  
        * Test and Validate
        (needed by the subsequent program!)
      - replace BuildFile.xml by the version downloaded from this record
      - download HiggsDemoAnalyzer.cc from this record to the /src subdirectory
      - recompile
scram b
      - download demoanalyzer_cfg_level3data.py (data example) and 
        demoanalyzer_cfg_level3MC.py (Higgs simulation example)
      - create datasets directory if not yet existing
mkdir datasets
      - change to this directory
cd datasets
      - download the 2011 JSON validation file from [http://opendata.web.cern.ch/record/1001]
      - if not yet done at level 2, create the directory rootfiles and 
        download all the level 2 root files to this directory (see level 2)
      - run the two analysis jobs (one on data, one on MC, the input files 
        are already predefined) 
cmsRun demoanalyzer_cfg_level3data.py
        -> will produce output file DoubleMuParked2012C_10000_Higgs.root
        containing 1 Higgs candidate from the data
cmsRun demoanalyzer_cfg_level3MC.py
        -> will produce output file Higgs4L1file.root
        containing the Higgs signal distributions with reduced statistics
      - move the two .root files above to the .rootfiles directory, together
        with the predefined files
mv DoubleMuParked2012C_10000_Higgs.root rootfiles/.
mv Higgs4L1file.root rootfiles/.
      - change directory
cd rootfiles
      - download the macro M4Lnormdatall_lvl3.cc
      - on the linux prompt, type 
root -l M4Lnormdatall_lvl3.cc
        -> you will get the output plot on the screen;
        the magenta Higgs signal histogram will now be the one you produced, 
        and the one data event which you have selected will be shown as a blue 
        triangle 
      - either, on the ROOT canvas (picture) click 
file->Quit ROOT
        or, on the root [] prompt, type 
.q     
        -> you will exit ROOT and find the output plot in 
           mass4l_combined_user3.pdf 

## Level 4 (advanced): Full analysis

4. Reproduce the full example analysis 
   (up to ~1 month or more on single CPU with fast internet connection, 
    depending on internet connection speed and computer performance)
      - start by running level 3 and understand what you have done
      - download demoanalyzer_cfg_level4data.py and 
        demoanalyzer_cfg_level4MC.py 
      - at this level, instead of running over a single file, you will run 
        over so-called index files which contain chains of files
      - download all the data index files for the datasets listed in 
        List_indexfile.txt to the datasets directory
      - download the 2011 validation (JSON) file to the datasets directory
        (in which you should already have the 2012 one)
      - download all the MC index files for the MC sets listed in 
        List_indexfile.txt to the MCsets directory (after having created it)
      - edit the relevant demoanalyzer file and insert the index file you 
        want; for data, make sure to use the correct JSON validation file
        in each case; set an outputfile name of your choice for each smaple 
        which you will recognise 
      - run the analysis job (cmsRun demoanalyzer_cfg_level4...) sequentially 
        on all the input samples listed in List_indexfile.txt, i.e. produce 
        all root output files yourself.
        If you have access to a computer farm with local support for the 
        installation of the CMS software (the Open Data team can only provide 
        support for the single virtual machine mode), you may also run 
        the analysis in parallel on different CPUs, correspondingly speeding 
        up the result. 
      - merge all the files from different index files of a dataset by using
        ROOT tools   
      - go to 2., using your own Root output files instead of the predefined 
        ones
--->

# Setup instructions for the impatient

Run this code in CMS Open Data VM [http://opendata.cern.ch/VM/CMS/2011]

If you have not installed the CMSSW area do the following:
```
cmsrel CMSSW_5_3_32
```
If you already have, start directly with:

```
cd CMSSW_5_3_32/src
cmsenv
```


For this example, you need to create an additional directory, you can call it `WorkDirCMSHiggs` or choose another name.
Go to this directory, and download the example code from git. This repository contains the various examples.

```

mkdir WorkDirCMSHiggs
cd WorkDirCMSHiggs
git clone https://github.com/cms-opendata-analyses/HiggsExample11-12.git

```
And then follow the instructions above.

If you are going to do the more advanced (Level 3/4) tutorials, you also need to make a demo analyser and compile it

```
cd $CMSSW_BASE
mkdir Demo
cd Demo
mkedanlzr DemoAnalyzer
cd DemoAnalyzer 
scram b
cd $CMSSW_BASE
```
You need to move some files around for the Level 3 version of the tutorial, here again it is assumed you installed in a directory `WorkDirCMSHiggs`.
```
cd $CMSSW_BASE
cp WorkDirCMSHiggs/HiggsExample11-12/Level3/BuildFile.xml Demo/DemoAnalyzer/.
cp WorkDirCMSHiggs/HiggsExample11-12/Level3/HiggsDemoAnalyzer.cc Demo/DemoAnalyzer/src/.
scram b
```
