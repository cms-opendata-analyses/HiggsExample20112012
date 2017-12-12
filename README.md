# HiggsExample11-12
This example guides the user to reproducing the discovery of the Higgs boson using the 2011 and 2012 datasets. It contains multiple levels of examples, from very simple to a full analysis, all with CMS Open Data.

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
For this example, you need to create an additional directory, you can call it `WorkDir` or choose another name.
Go to this directory, and download the example code.

```

mkdir WorkDir
cd WorkDir
git clone https://github.com/cms-opendata-analyses/HiggsExample11-12.git
