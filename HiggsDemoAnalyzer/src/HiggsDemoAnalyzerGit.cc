//
// Original Author:
//         Created:  Fri December 7, 2017 by  N.Z. Jomhari   
//                   with contributions from  A. Geiser
//                                            A. Anuar
//                   and pieces similar to previous analysis examples
// $Id$
// ..
//
// ***************************************************************************
//  Higgs-to-four-lepton analysis example at research level                  *
// ***************************************************************************
//                                                                           *
// Built upon the DEMO setup provided by CMS open data access team,          *
// expanded/upgraded to contain a pedagocigal analysis example for the       *
// approximate reproducttion of the Higgs to four lepton mass spectrum       *
// published in CMS-HIG-12-028                                               *
// Phys.Lett. B716 (2012) 30-61,  arXiv:1207.7235                            *
//                                                                           *
// This research level example is a strongly simplified reimplementation of  *
// parts of the original CMS Higgs to four lepton analysis                   *
//                                                                           * 
// The published reference plot which is being approximated in this example  *
// (in addition to many auxiliary plots) is                                  *
// https://inspirehep.net/record/1124338/files/H4l_mass_v3.png               *
// Other Higgs final states (e.g. Higgs to two photons), which were also     *
// part of the same CMS paper and strongly contributed to the Higgs          *
// discovery, are not covered by this example.                               *
//                                                                           *
// The example addresses users who feel they have at least some minimal      *
// understanding of the content of the CMS paper and of the meaning of this  *
// reference plot, which can be reached via (separate) educational exercises.*
// A Root ntuple containing the kinematic information about the event        *
// candidates, which could be used for educational purposes, is part of the  *
// Root output file.                                                         *
//                                                                           *
// The analysis code provided here recodes the spirit of the original        *
// analysis and (approximately) recodes many of the original cuts on         *
// original data objects, but does not provide the original analysis code    *
// itself. Also, for the sake of simplicity, it skips some of the more       *
// advanced analysis methods of the original paper, and does not use any     *
// corrections beyond those already implicit in the objects used, and no     *
// systematic uncertainties.                                                 * 
// Another reason why the published spectrum is only reproduced very         *
// approximately is that the data sets only partially overlap, and that the  *
// legacy software version and corresponding calibrations differ from those  *
// of the original paper.                                                    *                                                      
// Nevertheless, the example provides a qualitative insight about how the    *
// original result was obtained.                                             *
//                                                                           *
// In addition to the documented core results, the resulting  Root files     *
// also contain many undocumented plots which grew as a side product from    *
// setting up this example and earlier examples.                             *
// And it contains an ntuple with candidate four-vectors as mentioned above. *
// ***************************************************************************

// ***************************************************************************
// Analysis strategy                                                         *
// ***************************************************************************
//
// The analysis strategy is the following: Get the histograms for the 4mu    *
// and 2mu2e final states from the DoubleMu datasets and those for 4e final  *
// state from the DoubleElectron dataset. This avoids double counting due to *
// trigger overlaps.                                                         *
// The code itself is agnostic with respect to the input data set used, and  *
// the appropriate histograms have to selected at the subsequent root        *
// analysis step.                                                            *
// ***************************************************************************

// system include files
#include <memory>
#include <vector>
#include <algorithm>
#include <utility>

// user include files, general
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//------ EXTRA HEADER FILES--------------------//
#include "math.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"

// for Root histogramming
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TLorentzVector.h"

// for tracking information
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

// for vertex information 
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// for muon information
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"

//for beamspot information
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//for electron informaton
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

// for particle flow information
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

// class declaration
class HiggsDemoAnalyzerGit: public edm::EDAnalyzer {
public:
  explicit HiggsDemoAnalyzerGit(const edm::ParameterSet&);
  ~HiggsDemoAnalyzerGit();

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  bool providesGoodLumisection(const edm::Event& iEvent);

  // Declare Root histograms or tree
  // For a description of their content see below

  // TTree *t1;
  // TTree *t2;
  // TTree *t3;
  
  TH1D *h_globalmu_size;
  TH1D *h_recomu_size;
  TH1D *h_e_size;

  TH1D *h_nggmu;
  TH1D *h_ngmu;
  TH1D *h_nge;

  TH1D *h_m1_gmu;
  TH1D *h_m2_gmu;
  TH1D *h_m3_gmu;

  TH1D *h_mZ_2mu;
  TH1D *h_mZ_2e;

  TH1D *h_mZ12_4mu;
  TH1D *h_mZ34_4mu;
  TH1D *h_mZ13_4mu;
  TH1D *h_mZ24_4mu;
  TH1D *h_mZ14_4mu;
  TH1D *h_mZ23_4mu;
  TH1D *h_mZa_4mu;
  TH1D *h_mZb_4mu;
  TH1D *h_m1_m4mu;
  TH1D *h_m2_m4mu;
  TH1D *h_m3_m4mu;
  TH1D *h_m4_m4mu;

  TH1D *h_mZ12_4e; 
  TH1D *h_mZ34_4e;
  TH1D *h_mZ13_4e;
  TH1D *h_mZ24_4e;
  TH1D *h_mZ14_4e;
  TH1D *h_mZ23_4e;
  TH1D *h_mZa_4e;
  TH1D *h_mZb_4e;
  TH1D *h_m1_m4e;
  TH1D *h_m2_m4e;
  TH1D *h_m3_m4e;
  TH1D *h_m4_m4e;

  TH1D *h_mZmu_2mu2e;
  TH1D *h_mZe_2mu2e;
  TH1D *h_mZa_2mu2e;
  TH1D *h_mZb_2mu2e;
  TH1D *h_m1_m2mu2e;
  TH1D *h_m2_m2mu2e;
  TH1D *h_m3_m2mu2e;
  TH1D *h_m4_m2mu2e;

  // Control Plot

  // Global Muon
  TH1D *h_p_gmu;
  TH1D *h_pt_gmu_b4;
  TH1D *h_eta_gmu_b4;
  TH1D *h_chi2_gmu;
  TH1D *h_phi_gmu;
  TH1D *h_ndof_gmu;
  TH1D *h_normchi2_gmu;
  TH1D *h_validhits_gmu;
  TH1D *h_pixelhits_gmu;

  TH1D *h_pt_gmu_after;
  TH1D *h_eta_gmu_after;

  // Tracker Muon
  TH1D *h_p_reco;
  TH1D *h_pt_reco_b4;
  TH1D *h_eta_reco_b4;
  TH1D *h_phi_reco;
  TH1D *h_chi2_reco;
  TH1D *h_ndof_reco;
  TH1D *h_normchi2_reco;

  // PF Muon
  TH1D *h_goodhit;
  TH1D *h_dxy_mu;
  TH1D *h_relPFIso_mu;

  TH1D *h_relPFIso_mu_after;
  TH1D *h_validhits_mu;
  TH1D *h_pixelhits_mu;

  TH1D *h_pt_after_Zto2mu;
  TH1D *h_eta_after_Zto2mu;

  TH1D *h_pt_after;
  TH1D *h_eta_after;
  TH1D *h_pt_after_2mu2e;
  TH1D *h_eta_after_2mu2e;

  // Electron
  TH1D *h_p_e;
  TH1D *h_et_e;
  TH1D *h_pt_e_b4;
  TH1D *h_eta_e_b4;
  TH1D *h_phi_e;
  TH1D *h_sc_eta;
  TH1D *h_sc_rawE;
  TH1D *h_relPFIso_e;
  TH1D *h_relPFIso_e_after;
  TH2D *h_relPFIso_pt_e;

  TH1D *h_dxy_e;

  TH1D *h_pt_e_after_Zto2e;
  TH1D *h_eta_e_after_Zto2e;

  TH1D *h_pt_e_after;
  TH1D *h_eta_e_after;

  TH1D *h_relPFIso_2mu_after;
  TH1D *h_relPFIso_2e_after;

  TH1D *h_pt_e_after_2mu2e;
  TH1D *h_eta_e_after_2mu2e;

  TH1D *h_SIP3d_mu_b4;
  TH1D *h_SIP3d_e_b4;
  TH1D *h_misshite;

  // Declare variables
  int  nGoodGlobalMuon, nGoodRecoMuon;

  double s1, s2, s3, s4, s;
  double dx,dy,dz, rap, pz;

  double mZ12, mZ34, mZ13, mZ24, mZ14, mZ23;
  double mass4mu, pt_4mu, eta_4mu, phi_4mu;
  double px4mu, py4mu, pz4mu, E4mu;
  double pt_mu1, pt_mu2, pt_mu3, pt_mu4;
  double eta_mu1, eta_mu2, eta_mu3, eta_mu4;
  double phi_mu1, phi_mu2, phi_mu3, phi_mu4;
  int cas_mu1, cas_mu2, cas_mu3, cas_mu4;

  double px_mu1, px_mu2, px_mu3, px_mu4;
  double py_mu1, py_mu2, py_mu3, py_mu4;
  double pz_mu1, pz_mu2, pz_mu3, pz_mu4;

  double E_mu1, E_mu2, E_mu3, E_mu4;

  double mZa, mZb;

  double sqm1, mZ, eZ12, eZ34, eZ13, eZ24, eZ14, eZ23; 
  double pxZ12, pxZ34, pxZ13, pxZ24,  pxZ14, pxZ23;
  double pyZ12, pyZ34, pyZ13, pyZ24, pyZ14, pyZ23;
  double pzZ12, pzZ34, pzZ13, pzZ24, pzZ14, pzZ23;
  double pZ12, pZ34, pZ13, pZ24, pZ14, pZ23;
  double pTZ12, pTZ34, pTZ13, pTZ24, pTZ14, pTZ23;

  double dZ12, dZ34, dZ13, dZ24, dZ14, dZ23;
  double dZc1, dZc2, dZc3;
  double eZa, pxZa, pyZa, pzZa, pTZa;
  double eZb, pxZb, pyZb, pzZb, pTZb;

  int nGoodElectron;
  double sqme;
  int misshits;

  double mass4e, pt_4e, eta_4e, phi_4e;
  double px4e, py4e, pz4e, E4e;
  double pt_e1, pt_e2, pt_e3, pt_e4;
  double eta_e1, eta_e2, eta_e3, eta_e4;
  double phi_e1, phi_e2, phi_e3, phi_e4;
  int cas_e1, cas_e2, cas_e3, cas_e4;

  double px_e1, px_e2, px_e3, px_e4;
  double py_e1, py_e2, py_e3, py_e4;
  double pz_e1, pz_e2, pz_e3, pz_e4;

  double E_e1, E_e2, E_e3, E_e4;

  double mass2mu2e, pt_2mu2e, eta_2mu2e, phi_2mu2e;
  double px2mu2e, py2mu2e, pz2mu2e, E2mu2e;
  double pt_2mu1, pt_2mu2, pt_2e1, pt_2e2;
  double eta_2mu1, eta_2mu2, eta_2e1, eta_2e2;
  double phi_2mu1, phi_2mu2, phi_2e1, phi_2e2;
  int cas_2mu1, cas_2mu2, cas_2e1, cas_2e2;

  double px_2mu1, px_2mu2, px_2e1, px_2e2;
  double py_2mu1, py_2mu2, py_2e1, py_2e2;
  double pz_2mu1, pz_2mu2, pz_2e1, pz_2e2;

  double E_2mu1, E_2mu2, E_2e1, E_2e2;

  double goodhit;
  double relPFIso_mu;
  double relPFIso_e;

  double IP3d_mu;
  double ErrIP3d_mu;
  double SIP3d_mu;
  double IP3d_e;
  double ErrIP3d_e;
  double SIP3d_e;

  int nRun,nEvt,nLumi;

  TLorentzVector p4Za, p4Zb, p4H;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

HiggsDemoAnalyzerGit::HiggsDemoAnalyzerGit(const edm::ParameterSet& iConfig) {

  // *****************************************************************
  // This is the main analysis routine
  // The goal is to approximately reproduce the Higgs-to-four-lepton    
  // mass spectrum from HIG-12-028
  // *****************************************************************

  // now do what ever initialization is needed
  edm::Service<TFileService> fs;

  // ************************************
  // book histograms and set axis labels
  // (called once for initialization)
  // ************************************

  // Global Muon (GM) size
  h_globalmu_size = fs->make<TH1D>("NGMuons", "GMuon size", 10, 0., 10.);
  h_globalmu_size->GetXaxis()->SetTitle("Number of GMuons");
  h_globalmu_size->GetYaxis()->SetTitle("Number of Events");

  // Reco Muon (RM) size
  h_recomu_size = fs->make<TH1D>("NMuons", "Reco muon size", 10, 0., 10.);
  h_recomu_size->GetXaxis()->SetTitle("Number of Muons");
  h_recomu_size->GetYaxis()->SetTitle("Number of Events");

  // Electron size
  h_e_size = fs->make<TH1D>("Nelectrons", "Electron size", 10, 0., 10.);
  h_e_size->GetXaxis()->SetTitle("Number of Electrons");
  h_e_size->GetYaxis()->SetTitle("Number of Events");

  // No. of good global muon 
  h_nggmu = fs->make<TH1D>("NGoodGMuons", "No. of good global muon", 10, 0., 10.);
  h_nggmu->GetXaxis()->SetTitle("Number of GMuons");
  h_nggmu->GetYaxis()->SetTitle("Number of Events");

  // No. of good reco muon 
  h_ngmu = fs->make<TH1D>("NGoodRecMuons", "No. of good reco muon", 10, 0., 10.);
  h_ngmu->GetXaxis()->SetTitle("Number of RMuons");
  h_ngmu->GetYaxis()->SetTitle("Number of Events");

  // No. of good electron
  h_nge = fs->make<TH1D>("NGoodElectron", "No. of good electron", 10, 0., 10.);
  h_nge->GetXaxis()->SetTitle("Number of Electrons");
  h_nge->GetYaxis()->SetTitle("Number of Events");

  // Dimuon mass spectrum up to 4 GeV (low mass range, rho/omega, phi, psi)with GM
  h_m1_gmu = fs->make<TH1D>("GMmass", "Global Muon mass", 400, 0., 4.);
  h_m1_gmu->GetXaxis()->SetTitle("Invariant Mass for Nmuon=2 (in GeV/c^2)");
  h_m1_gmu->GetYaxis()->SetTitle("Number of Events");

  // Dimuon mass spectrum up to 120 GeV (high mass range: upsilon, Z) with GM
  h_m2_gmu = fs->make<TH1D>("GMmass_extended", "GMmass", 240, 0., 120.);
  h_m2_gmu->GetXaxis()->SetTitle("Invariant Mass for Nmuon=2 (in GeV/c^2)");
  h_m2_gmu->GetYaxis()->SetTitle("Number of Events");

  // Dimuon mass spectrum up to 600 GeV with Global Muon
  h_m3_gmu = fs->make<TH1D>("GMmass_extended_600", "GMmass", 240, 0., 600.);
  h_m3_gmu->GetXaxis()->SetTitle("Invariant Mass for Nmuon=2 (in GeV/c^2)");
  h_m3_gmu->GetYaxis()->SetTitle("Number of Events");

  // ZTo2mu mass spectrum with reco muon
  h_mZ_2mu = fs->make<TH1D>("massZto2muon", "mass of Z to 2 muon",120, 40., 120.);
  h_mZ_2mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZ_2mu->GetYaxis()->SetTitle("Number of Events");

  // ZTo2e mass spectrum with Gsf electron
  h_mZ_2e = fs->make<TH1D>("massZto2e", "mass of Z to 2e", 120, 40., 120.);
  h_mZ_2e->GetXaxis()->SetTitle("Invariant Mass for Nelectron=2 (in GeV/c^2)");
  h_mZ_2e->GetYaxis()->SetTitle("Number of Events");

  // These histograms are for 4muon reconstruction with different combinations
  
  // First combination: 1234
  // Mass Z12
  h_mZ12_4mu = fs->make<TH1D>("mZ12_4mu", "mass of Z12", 75, 0., 150.);
  h_mZ12_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZ12_4mu->GetYaxis()->SetTitle("Number of Events");

  // Mass Z34
  h_mZ34_4mu = fs->make<TH1D>("mZ34_4mu", "mass of Z34", 75, 0., 150.);
  h_mZ34_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZ34_4mu->GetYaxis()->SetTitle("Number of Events");

  // Second combination: 1324
  // Mass Z13
  h_mZ13_4mu = fs->make<TH1D>("mZ13_4mu", "mass of Z13", 75, 0., 150.);
  h_mZ13_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZ13_4mu->GetYaxis()->SetTitle("Number of Events");

  // Mass Z24
  h_mZ24_4mu = fs->make<TH1D>("mZ24_4mu", "mass of Z24", 75, 0., 150.);
  h_mZ24_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon");
  h_mZ24_4mu->GetYaxis()->SetTitle("Number of Events");

  // Third combination: 1423
  // Mass Z14
  h_mZ14_4mu = fs->make<TH1D>("mZ14_4mu", "mass of Z14", 75, 0., 150.);
  h_mZ14_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZ14_4mu->GetYaxis()->SetTitle("Number of Events");

  // Mass Z23
  h_mZ23_4mu = fs->make<TH1D>("mZ23_4mu", "mass of Z23", 75, 0., 150.);
  h_mZ23_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZ23_4mu->GetYaxis()->SetTitle("Number of Events");

  // Mass Za: mass of ZTo2mu closest to Z mass
  h_mZa_4mu = fs->make<TH1D>("mZa_4mu", "mass Za", 120, 0., 120.);
  h_mZa_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZa_4mu->GetYaxis()->SetTitle("Number of Events");

  // Mass Zb: mass of ZTo2mu not closest to Z mass
  h_mZb_4mu = fs->make<TH1D>("mZb_4mu", "mass Zb", 120, 0., 120.);
  h_mZb_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZb_4mu->GetYaxis()->SetTitle("Number of Events");

  // 4muon mass spectrum with reco muon (paper 7 TeV)
  h_m1_m4mu = fs->make<TH1D>("mass4mu_7TeV", "mass of 4 muon", 51, 98., 608.);
  h_m1_m4mu->GetXaxis()->SetTitle("Invariant mass of 4muons (in GeV/c^2)");
  h_m1_m4mu->GetYaxis()->SetTitle("Number of Events");

  // 4muon mass spectrum with reco muon (paper 8TeV)
  h_m2_m4mu = fs->make<TH1D>("mass4mu_8TeV", "mass of 4 muon", 74, 70., 810.);
  h_m2_m4mu->GetXaxis()->SetTitle("Invariant mass of 4muons (in GeV/c^2)");
  h_m2_m4mu->GetYaxis()->SetTitle("Number of Events");

  // 4muon mass spectrum with reco muon (paper 8TeV lower range)
  h_m3_m4mu = fs->make<TH1D>("mass4mu_8TeV_low", "mass of 4 muon", 37, 70., 181.);
  h_m3_m4mu->GetXaxis()->SetTitle("Invariant mass of 4muons (in GeV/c^2)");
  h_m3_m4mu->GetYaxis()->SetTitle("Number of Events");

  // 4muon mass spectrum with reco muon (full mass range)
  h_m4_m4mu = fs->make<TH1D>("mass4mu_full", "mass of 4 muon", 300, 0., 900.);
  h_m4_m4mu->GetXaxis()->SetTitle("Invariant mass of 4muons (in GeV/c^2)");
  h_m4_m4mu->GetYaxis()->SetTitle("Number of Events");

  // These histograms are for 4electron reconstruction with different combinations
  
  // First combination: 1234
  // Mass Z12
  h_mZ12_4e = fs->make<TH1D>("mZ12_4e", "mass of Z12", 75, 0., 150.);
  h_mZ12_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZ12_4e->GetYaxis()->SetTitle("Number of Events");

  // Mass Z34
  h_mZ34_4e = fs->make<TH1D>("mZ34_4e", "mass of Z34", 75, 0., 150.);
  h_mZ34_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZ34_4e->GetYaxis()->SetTitle("Number of Events");

  // Second combination: 1324
  // Mass Z13
  h_mZ13_4e = fs->make<TH1D>("mZ13_4e", "mass of Z13", 75, 0., 150.);
  h_mZ13_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZ13_4e->GetYaxis()->SetTitle("Number of Events");

  // Mass Z24
  h_mZ24_4e = fs->make<TH1D>("mZ24_4e", "mass of Z24", 75, 0., 150.);
  h_mZ24_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZ24_4e->GetYaxis()->SetTitle("Number of Events");

  // Third combination: 1423
  // Mass Z14
  h_mZ14_4e = fs->make<TH1D>("mZ14_4e", "mass of Z14", 75, 0., 150.);
  h_mZ14_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZ14_4e->GetYaxis()->SetTitle("Number of Events");

  // Mass Z23
  h_mZ23_4e = fs->make<TH1D>("mZ23_4e", "mass of Z23", 75, 0., 150.);
  h_mZ23_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZ23_4e->GetYaxis()->SetTitle("Number of Events");

  // Mass Za: mass of Z closest to Z mass
  h_mZa_4e = fs->make<TH1D>("mZa_4e", "mass Za", 120, 0., 120.);
  h_mZa_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZa_4e->GetYaxis()->SetTitle("Number of Events");

  // Mass Zb: mass of Z not closest to Z mass
  h_mZb_4e = fs->make<TH1D>("mZb_4e", "mass Zb", 120, 0., 120.);
  h_mZb_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZb_4e->GetYaxis()->SetTitle("Number of Events");

  // 4electron mass spectrum with Gsf Electron (paper 7TeV)
  h_m1_m4e = fs->make<TH1D>("mass4e_7TeV", "mass of 4 electron", 51, 98., 608.);
  h_m1_m4e->GetXaxis()->SetTitle("Invariant mass of 4e (in GeV/c^2)");
  h_m1_m4e->GetYaxis()->SetTitle("Number of Events");

  // 4electron mass spectrum with Gsf Electron (paper 8TeV)
  h_m2_m4e = fs->make<TH1D>("mass4e_8TeV", "mass of 4 electron", 74, 70., 810.);
  h_m2_m4e->GetXaxis()->SetTitle("Invariant mass of 4e (in GeV/c^2)");
  h_m2_m4e->GetYaxis()->SetTitle("Number of Events");

  // 4electron mass spectrum with Gsf Electron (paper 8TeV lower range)
  h_m3_m4e = fs->make<TH1D>("mass4e_8TeV_low", "mass 4 electron", 37, 70., 181.);
  h_m3_m4e->GetXaxis()->SetTitle("Invariant mass of 4e (in GeV/c^2)");
  h_m3_m4e->GetYaxis()->SetTitle("Number of Events");

  // 4electron mass spectrum with Gsf Electron (full mass range)
  h_m4_m4e = fs->make<TH1D>("mass4e_full", "mass of 4 electron", 300, 0., 900.);
  h_m4_m4e->GetXaxis()->SetTitle("Invariant mass of 4e (in GeV/c^2)");
  h_m4_m4e->GetYaxis()->SetTitle("Number of Events");

  
  // These histograms are for 2mu2e reconstruction with different combinations
  
  // Mass of Z to 2mu from 2mu2e
  h_mZmu_2mu2e = fs->make<TH1D>("massZmu_2mu2e", "mass Z2mu:2mu2e", 75, 0., 150.);
  h_mZmu_2mu2e->GetXaxis()->SetTitle("Invariant mass of Z1 (in GeV/c^2)");
  h_mZmu_2mu2e->GetYaxis()->SetTitle("Number of Events");

  // Mass of Z to 2e from 2mu2e
  h_mZe_2mu2e = fs->make<TH1D>("massZe_2mu2e", "mass Z2e:2mu2e", 75, 0., 150.);
  h_mZe_2mu2e->GetXaxis()->SetTitle("Invariant Mass of Z2 (in GeV/c^2)");
  h_mZe_2mu2e->GetYaxis()->SetTitle("Number of Events");

  // Mass Za: mass of Z1 closest to Z mass
  h_mZa_2mu2e = fs->make<TH1D>("mZa_2mu2e", "mass Z higher", 120, 0., 120.);
  h_mZa_2mu2e->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZa_2mu2e->GetYaxis()->SetTitle("Number of Events");

  // Mass Zb: mass of Z2 not closest to Z mass
  h_mZb_2mu2e = fs->make<TH1D>("mZb_2mu2e", "mass Z lower", 120, 0., 120.);
  h_mZb_2mu2e->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZb_2mu2e->GetYaxis()->SetTitle("Number of Events");

  // 2muon 2electron mass spectrum (paper 7TeV)
  h_m1_m2mu2e = fs->make<TH1D>("mass2mu2e_7TeV", "mass of 2mu2e", 51, 98., 608.);
  h_m1_m2mu2e->GetXaxis()->SetTitle("Invariant mass of 2mu2e (in GeV/c^2)");
  h_m1_m2mu2e->GetYaxis()->SetTitle("Number of Events");

  // 2muon 2electron mass spectrum (paper 8TeV)
  h_m2_m2mu2e = fs->make<TH1D>("mass2mu2e_8TeV", "mass of 2mu2e", 74, 70., 810.);
  h_m2_m2mu2e->GetXaxis()->SetTitle("Invariant mass of 2mu2e (in GeV/c^2)");
  h_m2_m2mu2e->GetYaxis()->SetTitle("Number of Events");

  // 2muon 2electron mass spectrum (paper 8TeV lower range)
  h_m3_m2mu2e = fs->make<TH1D>("mass2mu2e_8TeV_low", "mass 2mu2e", 37, 70., 181.);
  h_m3_m2mu2e->GetXaxis()->SetTitle("Invariant mass of 2mu2e (in GeV/c^2)");
  h_m3_m2mu2e->GetYaxis()->SetTitle("Number of Events");

  // 2muons 2electrons mass spectrum full mass range
  h_m4_m2mu2e = fs->make<TH1D>("mass2mu2e_full", "mass of 2mu2e", 300, 0., 900.);
  h_m4_m2mu2e->GetXaxis()->SetTitle("Invariant mass of 2mu2e (in GeV/c^2)");
  h_m4_m2mu2e->GetYaxis()->SetTitle("Number of Events");

  
  //////////////////////////// CONTROL PLOTS /////////////////////////////////////
  // Below are the histograms for the control plots

  //---- Control plots for dimuons, similar to those of the dimuon example -----//
  
  // Momentum of Global Muon
  h_p_gmu = fs->make<TH1D>("GM_momentum", "GM momentum", 200, 0., 200.);
  h_p_gmu->GetXaxis()->SetTitle("Momentum (GeV/c)");
  h_p_gmu->GetYaxis()->SetTitle("Number of Events");

  // Transverse momentum of Global Muon
  h_pt_gmu_b4 = fs->make<TH1D>("b4_GM_pT", "GM pT", 200, 0., 200.);
  h_pt_gmu_b4->GetXaxis()->SetTitle("pT (GeV/c)");
  h_pt_gmu_b4->GetYaxis()->SetTitle("Number of Events");

  // Pseudorapidity of Global Muon
  h_eta_gmu_b4 = fs->make<TH1D>("b4_GM_eta", "GM eta", 140, -3.5, 3.5);
  h_eta_gmu_b4->GetXaxis()->SetTitle("eta");
  h_eta_gmu_b4->GetYaxis()->SetTitle("Number of Events");

  // Phi of Global Muon
  h_phi_gmu = fs->make<TH1D>("GM_phi", "GM phi", 314, -3.15, 3.15);
  h_phi_gmu->GetXaxis()->SetTitle("Phi");
  h_phi_gmu->GetYaxis()->SetTitle("Number of Events");

  // Chi square of Global Muon
  h_chi2_gmu = fs->make<TH1D>("GM_chi2", "GM chi2", 300, 0., 150.);
  h_chi2_gmu->GetXaxis()->SetTitle("Chi2");
  h_chi2_gmu->GetYaxis()->SetTitle("Number of Events");

  // Number of degrees of freedom of Global Muon
  h_ndof_gmu = fs->make<TH1D>("GM_ndof", "GM NDoF", 100, 0., 100.);
  h_ndof_gmu->GetXaxis()->SetTitle("Ndof");
  h_ndof_gmu->GetYaxis()->SetTitle("Number of Events");

  // Normalized chi square of Global Muon
  h_normchi2_gmu = fs->make<TH1D>("GM_normchi2", "GM NormChi2", 200, 0., 20.);
  h_normchi2_gmu->GetXaxis()->SetTitle("NormalizedChi2");
  h_normchi2_gmu->GetYaxis()->SetTitle("Number of Events");

  // Validhits of Global Muon
  h_validhits_gmu = fs->make<TH1D>("GM_validhits", "GM ValidHits", 100, 0., 100.);
  h_validhits_gmu->GetXaxis()->SetTitle("Number of valid hits");
  h_validhits_gmu->GetYaxis()->SetTitle("Number of Events");

  // Pixelhits of Global Muon
  h_pixelhits_gmu = fs->make<TH1D>("GM_pixelhits", "GM Pixelhits", 14, 0., 14.);
  h_pixelhits_gmu->GetXaxis()->SetTitle("Munber of pixel hits");
  h_pixelhits_gmu->GetYaxis()->SetTitle("Number of Events");

  // pT of Global Muon after cuts
  h_pt_gmu_after = fs->make<TH1D>("after_GM_pT", "GM pT", 200, 0., 200.);
  h_pt_gmu_after->GetXaxis()->SetTitle("pT of global muons (GeV/c)");
  h_pt_gmu_after->GetYaxis()->SetTitle("Number of Events");

  // eta of Global Muon after cuts
  h_eta_gmu_after = fs->make<TH1D>("after_GM_eta", "GM Eta", 140, -3.5, 3.5);
  h_eta_gmu_after->GetXaxis()->SetTitle("eta");
  h_eta_gmu_after->GetYaxis()->SetTitle("Number of Events");

  //-------------------- End of control plot dimuon example --------------------//

  //------------------------- Reco Muon control plots --------------------------//

  // Momentum of Reco Muon 
  h_p_reco = fs->make<TH1D>("RM_momentum", "RM Momentum", 200, 0., 200.);
  h_p_reco->GetXaxis()->SetTitle("Momentum (GeV/c)");
  h_p_reco->GetYaxis()->SetTitle("Number of Events");

  // pT of Reco Muon
  h_pt_reco_b4 = fs->make<TH1D>("b4_RM_pt", "RM pT", 200, 0., 200.);
  h_pt_reco_b4->GetXaxis()->SetTitle("pT (GeV/c)");
  h_pt_reco_b4->GetYaxis()->SetTitle("Number of Events");

  // Eta of Reco Muon
  h_eta_reco_b4 = fs->make<TH1D>("b4_RM_eta", "RM Eta", 140, -3.5, 3.5);
  h_eta_reco_b4->GetXaxis()->SetTitle("eta");
  h_eta_reco_b4->GetYaxis()->SetTitle("Number of Events");

  // Phi of Reco Muon
  h_phi_reco = fs->make<TH1D>("RM_phi", "RM Phi", 314, -3.15, 3.15);
  h_phi_reco->GetXaxis()->SetTitle("Phi");
  h_phi_reco->GetYaxis()->SetTitle("Number of Events");

  // Chi of Reco Muon
  h_chi2_reco = fs->make<TH1D>("RM_chi2", "RM Chi2", 500, 0., 100.);
  h_chi2_reco->GetXaxis()->SetTitle("Chi2");
  h_chi2_reco->GetYaxis()->SetTitle("Number of Events");

  // No. of degrees of freedom of Reco Muon
  h_ndof_reco = fs->make<TH1D>("RM_Ndof", "RM NDoF", 60, 0., 60.);
  h_ndof_reco->GetXaxis()->SetTitle("Ndof");
  h_ndof_reco->GetYaxis()->SetTitle("Number of Events");

  // Normalizedchi2 of Reco Muon
  h_normchi2_reco = fs->make<TH1D>("RM_NormChi2", "RM NormChi2", 200, 0., 20.);
  h_normchi2_reco->GetXaxis()->SetTitle("NormalizedChi2");
  h_normchi2_reco->GetYaxis()->SetTitle("Number of Events");

  // No. of valid muon hits of Reco Muon
  h_goodhit = fs->make<TH1D>("RM_goodMuonChamberHit", "RM Mu Hit", 40, 0., 40.);
  h_goodhit->GetXaxis()->SetTitle("Number of muon valid hits");
  h_goodhit->GetYaxis()->SetTitle("Number of Events");

  // Transverse impact parameter with respect to primary vertex of Reco Muon
  h_dxy_mu = fs->make<TH1D>("RM_dxy", "RM dxy", 100, 0., 1.);
  h_dxy_mu->GetXaxis()->SetTitle("Transverse impact parameter w.r.t. BS");
  h_dxy_mu->GetYaxis()->SetTitle("Number of Events");

  // Relative isolation of Reco Muon
  h_relPFIso_mu = fs->make<TH1D>("RM_RelPFIso", "R.PFIso", 100, 0., 5.);
  h_relPFIso_mu->GetXaxis()->SetTitle("Relative Isolation mu");
  h_relPFIso_mu->GetYaxis()->SetTitle("Number of Events");

  // Relative isolation of Reco Muon after cuts
  h_relPFIso_mu_after=fs->make<TH1D>("after_RM_RelPFIso", "R.PFIso", 100, 0., 5.);
  h_relPFIso_mu_after->GetXaxis()->SetTitle("Relative Isolation mu");
  h_relPFIso_mu_after->GetYaxis()->SetTitle("Number of Events");

  // Validhits of Reco Muon
  h_validhits_mu = fs->make<TH1D>("RM_validhits", "RM ValidHits", 100, 0., 100.);
  h_validhits_mu->GetXaxis()->SetTitle("Number of valid hits");
  h_validhits_mu->GetYaxis()->SetTitle("Number of Events");

  // Pixelhits of Reco Muon
  h_pixelhits_mu = fs->make<TH1D>("RM_pixelhits", "RM Pixelhits", 14, 0., 14.);
  h_pixelhits_mu->GetXaxis()->SetTitle("Munber of pixel hits");
  h_pixelhits_mu->GetYaxis()->SetTitle("Number of Events");

  // pT reco muon after cuts for Z to 2mu
  h_pt_after_Zto2mu = fs->make<TH1D>("after_RM_pt_Z2mu", "Mu pT", 200, 0., 200.);
  h_pt_after_Zto2mu->GetXaxis()->SetTitle("pT (GeV/c)");
  h_pt_after_Zto2mu->GetYaxis()->SetTitle("Number of Events");

  // Eta reco muon after cuts for Z to 2mu
  h_eta_after_Zto2mu=fs->make<TH1D>("after_RM_eta_Z2mu","Mu eta", 140, -3.5, 3.5);
  h_eta_after_Zto2mu->GetXaxis()->SetTitle("Eta");
  h_eta_after_Zto2mu->GetYaxis()->SetTitle("Number of Events");

  // pT reco muon after cuts
  h_pt_after = fs->make<TH1D>("after_RM_pt", "Muon pT", 200, 0., 200.); 
  h_pt_after->GetXaxis()->SetTitle("pT (GeV/c)");
  h_pt_after->GetYaxis()->SetTitle("Number of Events");

  // Eta reco muon after cuts
  h_eta_after = fs->make<TH1D>("after_RM_eta", "Muon eta", 140, -3.5, 3.5);
  h_eta_after->GetXaxis()->SetTitle("eta");
  h_eta_after->GetYaxis()->SetTitle("Number of Events");

  // 2mu2e case
  // Relative isolation of 2mu for 2mu2e after cuts
  h_relPFIso_2mu_after=fs->make<TH1D>("after_relPFIso_2mu","R.PFIso", 50, 0., 5.);
  h_relPFIso_2mu_after->GetXaxis()->SetTitle("Relative Isolation mu");
  h_relPFIso_2mu_after->GetYaxis()->SetTitle("Number of Events");

  // Relative isolation of 2e for 2mu2e after cuts
  h_relPFIso_2e_after = fs->make<TH1D>("after_relPFIso_2e","R.PFIso", 50, 0., 5.);
  h_relPFIso_2e_after->GetXaxis()->SetTitle("Relative Isolation e");
  h_relPFIso_2e_after->GetYaxis()->SetTitle("Number of Events");

  // pT muon after cuts for 2mu2e
  h_pt_after_2mu2e = fs->make<TH1D>("after_pt_2mu2e", "Muon pT", 200, 0., 200.); 
  h_pt_after_2mu2e->GetXaxis()->SetTitle("pT (GeV/c)");
  h_pt_after_2mu2e->GetYaxis()->SetTitle("Number of Events");

  // Eta muon after cuts for 2mu2e
  h_eta_after_2mu2e=fs->make<TH1D>("after_eta_2mu2e", "Muon eta", 140, -3.5, 3.5);
  h_eta_after_2mu2e->GetXaxis()->SetTitle("eta");
  h_eta_after_2mu2e->GetYaxis()->SetTitle("Number of Events");

  //---------------------- End of Reco Muon control plots ----------------------//

  //-------------------------- Electron control plots --------------------------//

  // Electron momentum
  h_p_e = fs->make<TH1D>("e_momentum", "Electron momentum", 200, 0., 200.);
  h_p_e->GetXaxis()->SetTitle("Momentum (GeV/c)");
  h_p_e->GetYaxis()->SetTitle("Number of Events");

  // Electron eT
  h_et_e = fs->make<TH1D>("e_eT", "Electron eT", 200, 0., 200.);
  h_et_e->GetXaxis()->SetTitle("eT (GeV/c)");
  h_et_e->GetYaxis()->SetTitle("Number of Events");

  // Electron pT before cuts
  h_pt_e_b4 = fs->make<TH1D>("b4_e_pT", "Electron pT", 200, 0, 200.);
  h_pt_e_b4->GetXaxis()->SetTitle("pT (GeV/c)");
  h_pt_e_b4->GetYaxis()->SetTitle("Number of Events");

  // Electron eta before cuts
  h_eta_e_b4 = fs->make<TH1D>("b4_e_eta", "Electron eta", 140, -3.5, 3.5);
  h_eta_e_b4->GetXaxis()->SetTitle("eta");
  h_eta_e_b4->GetYaxis()->SetTitle("Number of Events");

  // Electron phi
  h_phi_e = fs->make<TH1D>("e_phi", "Electron phi", 314, -3.17, 3.17);
  h_phi_e->GetXaxis()->SetTitle("Phi");
  h_phi_e->GetYaxis()->SetTitle("Number of Events");

  // Electron SuperCluster (SC) eta
  h_sc_eta = fs->make<TH1D>("e_SC_eta", "Electron SC eta", 140, -3.5, 3.5);
  h_sc_eta->GetXaxis()->SetTitle("Super Cluster Eta");
  h_sc_eta->GetYaxis()->SetTitle("Number of Events");

  // Electron super cluster rawenergy
  h_sc_rawE = fs->make<TH1D>("e_SC_rawE", "Electron SC rawE", 200, 0., 200.);
  h_sc_rawE->GetXaxis()->SetTitle("Super Cluster Energy");
  h_sc_rawE->GetYaxis()->SetTitle("Number of Events");

  // Relative isolation electron
  h_relPFIso_e = fs->make<TH1D>("e_RelPFIso", "R.PFIso", 100, 0., 5.);
  h_relPFIso_e->GetXaxis()->SetTitle("Relative Isolation e");
  h_relPFIso_e->GetYaxis()->SetTitle("Number of Events");

  // Relative isolation electron after cuts
  h_relPFIso_e_after = fs->make<TH1D>("after_e_RelPFIso", "R.PFIso", 100, 0., 5.);
  h_relPFIso_e_after->GetXaxis()->SetTitle("Relative Isolation e");
  h_relPFIso_e_after->GetYaxis()->SetTitle("Number of Events");

  // Just for checking the relation of pT and RelPFIso
  double Iso[12]={0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.5, 10.};
  h_relPFIso_pt_e=fs->make<TH2D>("e_RelPFIso_pT","R.PFIso 2D",11,Iso,100,0., 50.);
  h_relPFIso_pt_e->GetXaxis()->SetTitle("Relative Isolation e");
  h_relPFIso_pt_e->GetYaxis()->SetTitle("Number of Events");

  // Transverse impact parameter with respect to primary vertex for Electron
  h_dxy_e = fs->make<TH1D>("e_dxy", "Electron dxy", 100, 0., 1.);
  h_dxy_e->GetXaxis()->SetTitle("Transverse impact parameter w.r.t. primvtx");
  h_dxy_e->GetYaxis()->SetTitle("Number of Events");

  // Electron pT after cuts for Zto2e
  h_pt_e_after_Zto2e = fs->make<TH1D>("after_e_pT_Zto2e", "e pT", 240, 0., 120.);
  h_pt_e_after_Zto2e->GetXaxis()->SetTitle("Electron pT (GeV/c)");
  h_pt_e_after_Zto2e->GetYaxis()->SetTitle("Number of Events");

  // Electron eta after cuts for Zto2e
  h_eta_e_after_Zto2e=fs->make<TH1D>("after_e_eta_Zto2e","e eta", 140, -3.5, 3.5);
  h_eta_e_after_Zto2e->GetXaxis()->SetTitle("Eta");
  h_eta_e_after_Zto2e->GetYaxis()->SetTitle("Number of Events");

  // Electron pT after cuts
  h_pt_e_after = fs->make<TH1D>("after_e_pT", "Electron pT", 240, 0., 120.);
  h_pt_e_after->GetXaxis()->SetTitle("Electron pT (GeV/c)");
  h_pt_e_after->GetYaxis()->SetTitle("Number of Events");

  // Electron eta after cuts
  h_eta_e_after = fs->make<TH1D>("after_e_eta", "Electron eta", 140, -3.5, 3.5);
  h_eta_e_after->GetXaxis()->SetTitle("eta");
  h_eta_e_after->GetYaxis()->SetTitle("Number of Events");

  // 2mu2e
  // Electron pT after cuts for 2mu2e
  h_pt_e_after_2mu2e = fs->make<TH1D>("after_e_pT_2mu2e", "e pT", 240, 0., 120.);
  h_pt_e_after_2mu2e->GetXaxis()->SetTitle("Electron pT (GeV/c)");
  h_pt_e_after_2mu2e->GetYaxis()->SetTitle("Number of Events");

  // Electron eta after cuts for 2mu2e
  h_eta_e_after_2mu2e=fs->make<TH1D>("after_e_eta_2mu2e","e eta", 140, -3.5, 3.5);
  h_eta_e_after_2mu2e->GetXaxis()->SetTitle("eta");
  h_eta_e_after_2mu2e->GetYaxis()->SetTitle("Number of Events");

  // SIP for muon
  h_SIP3d_mu_b4 = fs->make<TH1D>("SIP3d_mu", "SIP_3D for Muon", 100, 0., 10.);
  h_SIP3d_mu_b4->GetXaxis()->SetTitle("SIP_3D");
  h_SIP3d_mu_b4->GetYaxis()->SetTitle("Number of Events");

  // SIP for electron
  h_SIP3d_e_b4 = fs->make<TH1D>("SIP3d_e", "SIP_3D for Electron", 100, 0., 10.);
  h_SIP3d_e_b4->GetXaxis()->SetTitle("SIP_3D");
  h_SIP3d_e_b4->GetYaxis()->SetTitle("Number of Events");

  h_misshite = fs->make<TH1D>("e_misshit", "e track missing hits", 5, 0., 5.);
  h_misshite->GetXaxis()->SetTitle("gsfTrack Hit type");
  h_misshite->GetYaxis()->SetTitle("Number of Events");

}

HiggsDemoAnalyzerGit::~HiggsDemoAnalyzerGit() {
  //do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------//
void HiggsDemoAnalyzerGit::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

// **********************************************
// here each relevant event will get analyzed 
// **********************************************

  nRun  = iEvent.run();
  nEvt  = (iEvent.id()).event(); // iEvent: no class named event()
  nLumi = iEvent.luminosityBlock();
  
#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif

  // Event is to be analyzed

  edm::LogInfo("Demo")
  << "Starting to analyze \n"
  << "Event number: " << (iEvent.id()).event()
  << ", Run number: " << iEvent.run()
  << ", Lumisection: " << iEvent.luminosityBlock();

  //------------------Load (relevant) Event information------------------------//
  // INFO: Getting Data From an Event
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookChapter4#GetData
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEDMGetDataFromEvent#get_ByLabel
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAodDataTable
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideRecoDataTable

  // INFO: globalMuons
  // NB: note that when using keyword "globalMuons" getByLabel-function returns 
  // reco::TrackCollection
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel("generalTracks", tracks);

  edm::Handle<reco::TrackCollection> gmuons;
  iEvent.getByLabel("globalMuons", gmuons);

  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByLabel("muons", muons);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);

  edm::Handle<reco::VertexCollection> primvtxHandle;
  iEvent.getByLabel("offlinePrimaryVertices", primvtxHandle);

  edm::Handle<reco::GsfElectronCollection> electrons;
  iEvent.getByLabel("gsfElectrons",electrons);

  reco::BeamSpot beamSpot;
  if ( beamSpotHandle.isValid() )
    {
      beamSpot = *beamSpotHandle;

    } else
    {
      edm::LogInfo("Demo")
	<< "No beam spot available from EventSetup \n";
    }

  reco::VertexCollection primvtx;
  if ( primvtxHandle.isValid() )
    {
      primvtx = *primvtxHandle;

    } else
    {
      edm::LogInfo("Demo")
	<< "No primary vertex available from EventSetup \n";
    }


  // Declare vector that contain a pair of variables u wanna save
  // In this case: size and pT
  std::vector< std::pair<int, double> > vIdPt;
  std::vector< std::pair<int, double> > vIdPtmu;
  std::vector< std::pair<int, double> > vIdPte;

  // Initialize variables
  eZ12 = -9999.; eZ34 = -9999.; eZ13 = -9999.; eZ24 = -9999.; eZ14 = -9999.; eZ23 = -9999.; // select largest, init -
  pxZ12 = -9999.; pxZ34 = -9999.; pxZ13 = -9999.; pxZ24 = -9999.; pxZ14 = -9999.; pxZ23 = -9999.;
  pyZ12 = -9999.; pyZ34 = -9999.; pyZ13 = -9999.; pyZ24 = -9999.; pyZ14 = -9999.; pyZ23 = -9999.;
  pzZ12 = -9999.; pzZ34 = -9999.; pzZ13 = -9999.; pzZ24 = -9999.; pzZ14 = -9999.; pzZ23 = -9999.; 
  pZ12 = -9999.; pZ34 = -9999.; pZ13 = -9999.; pZ24 = -9999.; pZ14 = -9999.; pZ23 = -9999.; 
  pTZ12 = -9999.; pTZ34 = -9999.; pTZ13 = -9999.; pTZ24 = -9999.; pTZ14 = -9999.; pTZ23 = -9999.; 

  mZ12 = -9999.; mZ34 = -9999.; mZ13 = -9999.; mZ24 = -9999.; mZ14 = -9999.; mZ23 = -9999.;

  dZ12 = 9999.; dZ34 = 9999.; dZ13 = 9999.; dZ24 = 9999.; dZ14 = 9999.; dZ23 = 9999.; // select smallest, init +
  dZc1 = 9999.; dZc2 = 9999.; dZc3 = 9999.;

  eZa = -9999.; pxZa = -9999.; pyZa = -9999.; pzZa = -9999.; pTZa = -9999.; mZa = -9999.;
  eZb = -9999.; pxZb = -9999.; pyZb = -9999.; pzZb = -9999.; pTZb = -9999.; mZb = -9999.;
  mass4mu = -9999.; pt_4mu = -9999.; eta_4mu = -9999.; phi_4mu = -9999.;
  px4mu = -9999.; py4mu = -9999.; pz4mu = -9999.; E4mu = -9999.;
  
  pt_mu1 = -9999.; pt_mu2 = -9999.; pt_mu3 = -9999.; pt_mu4 = -9999.;
  eta_mu1 = -9999.; eta_mu2 = -9999.; eta_mu3 = -9999.; eta_mu4 = -9999.;
  phi_mu1 = -9999.; phi_mu2 = -9999.; phi_mu3 = -9999.; phi_mu4 = -9999.;
  cas_mu1 = -999; cas_mu2 = -999; cas_mu3 = -999; cas_mu4 = -999;

  px_mu1 = -9999.; px_mu2 = -9999.; px_mu3 = -9999.; px_mu4 = -9999.;
  py_mu1 = -9999.; py_mu2 = -9999.; py_mu3 = -9999.; py_mu4 = -9999.;
  pz_mu1 = -9999.; pz_mu2 = -9999.; pz_mu3 = -9999.; pz_mu4 = -9999.;

  E_mu1 = -9999.; E_mu2 = -9999.; E_mu3 = -9999.; E_mu4 = -9999.;

  mass4e = -9999.; pt_4e = -9999.; eta_4e = -9999.; phi_4e = -9999.;
  px4e = -9999.; py4e = -9999.; pz4e = -9999.; E4e = -9999.;
  
  pt_e1 = -9999.; pt_e2 = -9999.; pt_e3 = -9999.; pt_e4 = -9999.;
  eta_e1 = -9999.; eta_e2 = -9999.; eta_e3 = -9999.; eta_e4 = -9999.;
  phi_e1 = -9999.; phi_e2 = -9999.; phi_e3 = -9999.; phi_e4 = -9999.;
  cas_e1 = -999; cas_e2 = -999; cas_e3 = -999; cas_e4 = -999;

  px_e1 = -9999.; px_e2 = -9999.; px_e3 = -9999.; px_e4 = -9999.;
  py_e1 = -9999.; py_e2 = -9999.; py_e3 = -9999.; py_e4 = -9999.;
  pz_e1 = -9999.; pz_e2 = -9999.; pz_e3 = -9999.; pz_e4 = -9999.;

  E_e1 = -9999.; E_e2 = -9999.; E_e3 = -9999.; E_e4 = -9999.;
  
  s = -9999.;
  s1 = -9999.; s2 = -9999.; s3 = -9999.; s4 = -9999.;
  dx = 9999.; dy = 9999.; dz = 9999.; rap = -9999.; pz = -9999.;

  mass2mu2e = -9999.; pt_2mu2e = -9999.; eta_2mu2e = -9999.;  phi_2mu2e = -9999.;
  px2mu2e = -9999.; py2mu2e = -9999.; pz2mu2e = -9999.; E2mu2e = -9999.;
  
  pt_2mu1 = -9999.; pt_2mu2 = -9999.; pt_2e1 = -9999.; pt_2e2 = -9999.;
  eta_2mu1 = -9999.; eta_2mu2 = -9999.; eta_2e1 = -9999.; eta_2e2 = -9999.;
  phi_2mu1 = -9999.; phi_2mu2 = -9999.; phi_2e1 = -9999.; phi_2e2 = -9999.;
  
  cas_2mu1 = -999; cas_2mu2 = -999; cas_2e1 = -999; cas_2e2 = -999;

  px_2mu1 = -9999.; px_2mu2 = -9999.; px_2e1 = -9999.; px_2e2 = -9999.;
  py_2mu1 = -9999.; py_2mu2 = -9999.; py_2e1 = -9999.; py_2e2 = -9999.;
  pz_2mu1 = -9999.; pz_2mu2 = -9999.; pz_2e1 = -9999.; pz_2e2 = -9999.;

  E_2mu1 = -9999.; E_2mu2 = -9999.; E_2e1 = -9999.; E_2e2 = -9999.;
  
  goodhit = -9999.;
  relPFIso_mu = -9999.;
  relPFIso_e = -9999.;

  IP3d_mu = -9999;
  ErrIP3d_mu = -9999;
  SIP3d_mu = -9999.;
  IP3d_e = -9999;
  ErrIP3d_e = -9999;
  SIP3d_e = -9999.;

  p4Za.SetPxPyPzE(0., 0., 0., 0.);
  p4Zb.SetPxPyPzE(0., 0., 0., 0.);
  p4H.SetPxPyPzE(0., 0., 0., 0.);

  // constant value squared muon mass, squared electron mass and Z mass
  sqm1 = (0.105658) * (0.105658);
  sqme = (0.0005109989) * (0.0005109989);
  mZ = 91.1876;

  h_globalmu_size->Fill(gmuons->size());
  h_recomu_size->Fill(muons->size());
  h_e_size->Fill(electrons->size());

  ///////////////////////////////////////////////////////////////////////////////
  /////////////////////// Global Muon Collection Start //////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  //******************************************************************************
  // This is some code similar to the one of the 'dimuon analysis example',      *
  // illustrating the use of global muons, but not directly used for the Higgs   *
  // analysis                                                                    *
  //******************************************************************************

  // Loop over global muons size and select good muons

  for (unsigned t = 0; t < gmuons->size(); t++)
    {
      const reco::Track &iMuon = (*gmuons)[t];
      const reco::HitPattern& p = iMuon.hitPattern();

      h_p_gmu->Fill(iMuon.p());
      h_pt_gmu_b4->Fill(iMuon.pt());
      h_eta_gmu_b4->Fill(iMuon.eta());
      h_phi_gmu->Fill(iMuon.phi());

      h_chi2_gmu->Fill(iMuon.chi2());
      h_ndof_gmu->Fill(iMuon.ndof());
      h_normchi2_gmu->Fill(iMuon.normalizedChi2());

      // Counter
      int GM_ValidHits = 0;
      int GM_PixelHits = 0;
   
      // loop over the hits of the track
      for (int i = 0; i < p.numberOfHits(); i++) 
	{
	  uint32_t hit = p.getHitPattern(i);

	  // if the hit is valid and in pixel
	  if (p.validHitFilter(hit) && p.pixelHitFilter(hit)) {GM_PixelHits++;}
	  if (p.validHitFilter(hit)) {GM_ValidHits++;}
	}

      h_validhits_gmu->Fill(GM_ValidHits);
      h_pixelhits_gmu->Fill(GM_PixelHits);

      if (GM_ValidHits >= 12 && GM_PixelHits >= 2 && iMuon.normalizedChi2() < 4.0)
	{
	  // Save a vector that contains 2 variables: gmuon size (.first) and the pT (.second)
	  vIdPt.push_back( std::make_pair(t, iMuon.pt()) );
	}
    }

  // Sort the pT (.second) to decending order (from highest pT to lowest pT)
  std::sort(vIdPt.begin(), vIdPt.end(), [](const std::pair<int, double> &idPt1, const std::pair<int, double> &idPt2) {return (idPt1.second > idPt2.second);});

  ///////////////////////////////////////////////////////////////////////////////
  //////////////////////// Global Muon Collection End ///////////////////////////
  ///////////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////////////
  ///////////////////////// Reco Muon Collection Start ///////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  
  //******************************************************************************
  // This muon collection is being used here for the Higgs->4 lepton analysis    *
  //******************************************************************************

  // Loop over muons size and select good muons

  for (unsigned u = 0; u < muons->size(); u++)
    {
      const reco::Muon &itMuon = (*muons)[u];

      // math::XYZPoint point(beamSpot.position());
      math::XYZPoint point(primvtx[0].position());
      
      // select global particle flow muons
      // some muons might not have valid track references
      if (itMuon.isPFMuon() && itMuon.isPFIsolationValid() && (itMuon.globalTrack()).isNonnull())
	{
	  h_p_reco->Fill(itMuon.p());
	  h_pt_reco_b4->Fill(itMuon.pt());
	  h_eta_reco_b4->Fill(itMuon.eta());
	  h_phi_reco->Fill(itMuon.phi());

	  h_chi2_reco->Fill((itMuon.globalTrack())->chi2());
	  h_ndof_reco->Fill((itMuon.globalTrack())->ndof());
	  h_normchi2_reco->Fill((itMuon.globalTrack())->normalizedChi2());

	  //======= Use Particle Flow (PF) Muon & PF relative isolation ========//

	  // PF Relative isolation for muons
	  relPFIso_mu = ((itMuon.pfIsolationR04()).sumChargedHadronPt +
			 (itMuon.pfIsolationR04()).sumNeutralHadronEt + 
			 (itMuon.pfIsolationR04()).sumPhotonEt) / itMuon.pt(); 

	  h_relPFIso_mu->Fill(relPFIso_mu);
	  
	  // Checking hit pattern info
	  const reco::HitPattern& RM_p = (itMuon.globalTrack())->hitPattern();

	  goodhit = RM_p.numberOfValidMuonHits();
	  h_goodhit->Fill(goodhit);

	  h_dxy_mu->Fill((itMuon.globalTrack())->dxy(point));

	  IP3d_mu = sqrt((itMuon.globalTrack()->dxy(point) * itMuon.globalTrack()->dxy(point)) + (itMuon.globalTrack()->dz(point) * itMuon.globalTrack()->dz(point)));
	  
	  ErrIP3d_mu = sqrt((itMuon.globalTrack()->d0Error() * itMuon.globalTrack()->d0Error()) + (itMuon.globalTrack()->dzError() * itMuon.globalTrack()->dzError()));
	  
	  SIP3d_mu = IP3d_mu / ErrIP3d_mu;

	  h_SIP3d_mu_b4->Fill(SIP3d_mu);

	  int RM_ValidHits = 0;
	  int RM_PixelHits = 0;

	  for (int i = 0; i < RM_p.numberOfHits(); i++) {
	    uint32_t hit = RM_p.getHitPattern(i);

	    // If the hit is valid and in pixel
	    if (RM_p.validHitFilter(hit) && RM_p.pixelHitFilter(hit))
	      {RM_PixelHits++;}
	    
	      if (RM_p.validHitFilter(hit)) {RM_ValidHits++;}
	    }

	  h_validhits_mu->Fill(RM_ValidHits);
	  h_pixelhits_mu->Fill(RM_PixelHits);

	  if (std::abs(SIP3d_mu) < 4. && std::abs((itMuon.globalTrack())->dxy(point)) < 0.5 && std::abs((itMuon.globalTrack())->dz(point)) < 1. && relPFIso_mu < 0.4)
	    {
	      if (itMuon.pt() > 5. && std::abs(itMuon.eta()) < 2.4)
		{
		  vIdPtmu.push_back( std::make_pair(u, itMuon.pt()) );
		}
	    }
	} // end of if (itMuon.isPFMuon().........
    } // end muons loop

  // Sort the pT (.second) to decending order (from highest pT to lowest pT)
  std::sort(vIdPtmu.begin(), vIdPtmu.end(), [] (const std::pair<int, double> &idPtmu1, const std::pair<int, double> &idPtmu2) {return (idPtmu1.second > idPtmu2.second);});

  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////// Muon Collection end /////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  
  ///////////////////////////////////////////////////////////////////////////////
  ////////////////////// Electron Collection Start //////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  for (unsigned te = 0; te < electrons->size(); te++)
    {
      const reco::GsfElectron &iElectron = (*electrons)[te];

      // math::XYZPoint point(beamSpot.position());
       math::XYZPoint point(primvtx[0].position());

       if (iElectron.passingPflowPreselection())
	 {
	   misshits = ((iElectron.gsfTrack())->trackerExpectedHitsInner()).numberOfHits();

	   IP3d_e = sqrt ( (iElectron.gsfTrack()->dxy(point) * iElectron.gsfTrack()->dxy(point)) + (iElectron.gsfTrack()->dz(point) * iElectron.gsfTrack()->dz(point)) );
	   ErrIP3d_e = sqrt ( (iElectron.gsfTrack()->d0Error() * iElectron.gsfTrack()->d0Error()) + (iElectron.gsfTrack()->dzError() * iElectron.gsfTrack()->dzError()) );
	   SIP3d_e = IP3d_e / ErrIP3d_e;

	   h_SIP3d_e_b4->Fill(SIP3d_e);

	   h_p_e->Fill(iElectron.p());
	   h_et_e->Fill(iElectron.et());
	   h_pt_e_b4->Fill(iElectron.pt());
	   h_eta_e_b4->Fill(iElectron.eta());
	   h_phi_e->Fill(iElectron.phi());
	   h_sc_eta->Fill((iElectron.superCluster())->eta());
	   h_sc_rawE->Fill(std::abs((iElectron.superCluster())->rawEnergy()));
	   h_misshite->Fill(misshits);
      
	   // Relative isolation for electron
	   relPFIso_e = ((iElectron.pfIsolationVariables()).chargedHadronIso +
			 (iElectron.pfIsolationVariables()).neutralHadronIso +
			 (iElectron.pfIsolationVariables()).photonIso) /iElectron.pt();
	
	   h_relPFIso_e->Fill(relPFIso_e);

	   h_relPFIso_pt_e->Fill(relPFIso_e, iElectron.pt());

	   h_dxy_e->Fill((iElectron.gsfTrack())->dxy(point));

	   // Electron selection
	   if (iElectron.pt() > 7.)
	     {
	       if ((std::abs((iElectron.superCluster())->eta())) < 2.5)
		 {
		   if (misshits <= 1 && std::abs(SIP3d_e) < 4.)
		     {
		       if (std::abs(iElectron.gsfTrack()->dxy(point)) < 0.5 && std::abs(iElectron.gsfTrack()->dz(point)) < 1.)
			 {
			   if (iElectron.isEB())
			     {
			       if (relPFIso_e < 0.4)
				 {
				   vIdPte.push_back( std::make_pair(te, iElectron.pt()) );
				 }
			     }
			   else if (iElectron.isEE())
			     {
			       if (relPFIso_e < 0.4)
				 {
				   vIdPte.push_back( std::make_pair(te, iElectron.pt()) );
				 }
			     }
			 }
		     }
		 }
	     }
	 }
    } // for (unsigned te = 0; te < electrons->size(); te++)
  
  // Sort the pT (.second) to decending order (from highest pT to lowest pT)
  std::sort(vIdPte.begin(), vIdPte.end(), [](const std::pair<int, double> &idPte1, const std::pair<int, double> &idPte2) {return (idPte1.second > idPte2.second);});

  ////////////////////////////////////////////////////////////////////////////////
  ///////////////////////// Electron Collection end //////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  // these index (vIdPt...) is declare as good candidates in integer
  nGoodGlobalMuon = vIdPt.size(); 
  nGoodRecoMuon = vIdPtmu.size(); 
  nGoodElectron = vIdPte.size(); 

  h_nggmu->Fill(nGoodGlobalMuon);
  h_ngmu->Fill(nGoodRecoMuon);
  h_nge->Fill(nGoodElectron);

  ////////////////////////////////////////////////////////////////////////////////
  /////////////////////// All calculation start here /////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  
  //====================== Dimuon using Global Muon start ======================//

  // For the case of nGoodGlobalMuon > 2, the selected two muons are always have the highest pT as we have sorted it before
  if (nGoodGlobalMuon >= 2)
    {
      const reco::Track &gmuon1 = (*gmuons)[vIdPt.at(0).first];
      const reco::Track &gmuon2 = (*gmuons)[vIdPt.at(1).first];
      
      // The sum of 2 charge global muon should be 0 (neutral)
      if (gmuon1.charge() + gmuon2.charge() == 0)
	{
	  for (unsigned i = 0; i < vIdPt.size(); i++)
	    {
	      // These pT and eta are filled after all the cuts
	      h_pt_gmu_after->Fill(vIdPt.at(i).second); // access directly .second as the second pair is already pT
	      h_eta_gmu_after->Fill(((*gmuons)[vIdPt.at(i).first]).eta());
	    }

	  s1 = sqrt(((gmuon1.p()) * (gmuon1.p()) + sqm1) * ((gmuon2.p()) * (gmuon2.p()) + sqm1));
	  s2 = gmuon1.px() * gmuon2.px() + gmuon1.py() * gmuon2.py() + gmuon1.pz() * gmuon2.pz();
	  s = sqrt(2.0 * (sqm1 + (s1 - s2)));

	  h_m1_gmu->Fill(s);
	  h_m2_gmu->Fill(s);
	  h_m3_gmu->Fill(s);

	}
    }

  //======================== Dimuon using Global Muon end ======================//

  //===================== ZTo2Muon using Reco Muon start =======================//

  if (nGoodRecoMuon >= 2)
    {
      const reco::Muon &muon1 = (*muons)[vIdPtmu.at(0).first];
      const reco::Muon &muon2 = (*muons)[vIdPtmu.at(1).first];

      if (muon1.charge() + muon2.charge() == 0)
	{
	  for (unsigned i = 0; i < vIdPtmu.size(); i++)
	    {
	      // These pT and eta are filled after all the cuts
	      // access directly .second as the second pair is already pT
	      h_pt_after_Zto2mu->Fill(vIdPtmu.at(i).second); 
	      h_eta_after_Zto2mu->Fill(((*muons)[vIdPtmu.at(i).first]).eta());
	    }

	  s1 = sqrt(((muon1.p()) * (muon1.p()) + sqm1) * ((muon2.p()) * (muon2.p()) + sqm1));
	  s2 = muon1.px() * muon2.px() + muon1.py() * muon2.py() + muon1.pz() * muon2.pz();
	  s = sqrt(2.0 * (sqm1 + (s1 - s2)));

	  h_mZ_2mu->Fill(s);
	}
    }

  //===================== ZTo2Muon using Reco Muon end =========================//


  //============================ ZTo2Electron start ============================//

  if (nGoodElectron >= 2)
    {
      const reco::GsfElectron &elec1 = (*electrons)[vIdPte.at(0).first];
      const reco::GsfElectron &elec2 = (*electrons)[vIdPte.at(1).first];

      if (elec1.charge() + elec2.charge() == 0)
	{
	  for (unsigned i = 0; i < vIdPte.size(); i++)
	    {
	      // These pT and eta are filled after all the cuts
	      // access directly .second as the second pair is already pT
	      h_pt_e_after_Zto2e->Fill(vIdPte.at(i).second); 
	      h_eta_e_after_Zto2e->Fill((((*electrons)[vIdPte.at(i).first]).superCluster())->eta());
	    }

	  s1 = sqrt(((elec1.p()) * (elec1.p()) + sqme) * ((elec2.p()) * (elec2.p()) + sqme));
	  s2 = elec1.px() * elec2.px() + elec1.py() * elec2.py() + elec1.pz() * elec2.pz();
	  s = sqrt(2.0 * (sqme + (s1 - s2)));

	  h_mZ_2e->Fill(s);
	}
    }

  //========================== ZTo2Electron end ================================//


  //======================== ZZ/ZZ*To4Muon start ===============================//

  // Now, for these goodmuons, pair up and calculate mass
  if (nGoodRecoMuon >= 4)
    {
      const reco::Muon &muon1 = (*muons)[vIdPtmu.at(0).first];
      const reco::Muon &muon2 = (*muons)[vIdPtmu.at(1).first];
      const reco::Muon &muon3 = (*muons)[vIdPtmu.at(2).first];
      const reco::Muon &muon4 = (*muons)[vIdPtmu.at(3).first];

      if (muon1.charge() + muon2.charge() + muon3.charge() + muon4.charge() == 0)
	{
	  // First combination: Combine muon 1234
	  if (muon1.charge() + muon2.charge() == 0) // each lepton pair cas = 0
	    {
	      eZ12 = (sqrt(muon1.p() * muon1.p() + sqm1)) +
		     (sqrt(muon2.p() * muon2.p() + sqm1));
	
	      pxZ12 = muon1.px() + muon2.px();
	      pyZ12 = muon1.py() + muon2.py();
	      pzZ12 = muon1.pz() + muon2.pz();

	      if (muon3.charge() + muon4.charge() == 0)
		{
		  eZ34 = (sqrt(muon3.p() * muon3.p() + sqm1)) +
		         (sqrt(muon4.p() * muon4.p() + sqm1));

		  pxZ34 = muon3.px() + muon4.px();
		  pyZ34 = muon3.py() + muon4.py();
		  pzZ34 = muon3.pz() + muon4.pz();

		  // Calculate p4
		  pZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12) + (pzZ12 * pzZ12));
		  pZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34) + (pzZ34 * pzZ34));

		  pTZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12));
		  pTZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34));

		  mZ12 = sqrt((eZ12 * eZ12) - (pZ12 * pZ12));
		  mZ34 = sqrt((eZ34 * eZ34) - (pZ34 * pZ34));

		  if (mZ12 > 0.) h_mZ12_4mu->Fill(mZ12);
		  if (mZ34 > 0.) h_mZ34_4mu->Fill(mZ34);
		}
	    }

	  dZ12 = std::abs( mZ12 - mZ );
	  dZ34 = std::abs( mZ34 - mZ );

	  // take the smallest difference between mass
	  // to use for 4muon combination
	  dZc1 = (dZ12 < dZ34) ? dZ12 : dZ34; 
 
	  // Second combination: Combine muon 1324
	  if (muon1.charge() + muon3.charge() == 0)
	    {
	      eZ13 = (sqrt(muon1.p() * muon1.p() + sqm1)) +
		     (sqrt(muon3.p() * muon3.p() + sqm1));

	      pxZ13 = muon1.px() + muon3.px();
	      pyZ13 = muon1.py() + muon3.py();
	      pzZ13 = muon1.pz() + muon3.pz();

	      if (muon2.charge() + muon4.charge() == 0)
		{
		  eZ24 = (sqrt(muon2.p() * muon2.p() + sqm1)) +
		         (sqrt(muon4.p() * muon4.p() + sqm1));

		  pxZ24 = muon2.px() + muon4.px();
		  pyZ24 = muon2.py() + muon4.py();
		  pzZ24 = muon2.pz() + muon4.pz();

		  // Calculate p4
		  pZ13 = sqrt((pxZ13 * pxZ13) + (pyZ13 * pyZ13) + (pzZ13 * pzZ13));
		  pZ24 = sqrt((pxZ24 * pxZ24) + (pyZ24 * pyZ24) + (pzZ24 * pzZ24));

		  pTZ13 = sqrt((pxZ13 * pxZ13) + (pyZ13 * pyZ13));
		  pTZ24 = sqrt((pxZ24 * pxZ24) + (pyZ24 * pyZ24));

		  mZ13 = sqrt((eZ13 * eZ13) - (pZ13 * pZ13));
		  mZ24 = sqrt((eZ24 * eZ24) - (pZ24 * pZ24));

		  if (mZ13 > 0.) h_mZ13_4mu->Fill(mZ13);
		  if (mZ24 > 0.) h_mZ24_4mu->Fill(mZ24);
		}
	    }

	  dZ13 = std::abs( mZ13 - mZ );
	  dZ24 = std::abs( mZ24 - mZ );

	  dZc2 = (dZ13 < dZ24) ? dZ13 : dZ24; 

	  // Third combination: Combine muon 1423
	  if (muon1.charge() + muon4.charge() == 0)
	    {
	      eZ14 = (sqrt(muon1.p() * muon1.p() + sqm1)) +
		     (sqrt(muon4.p() * muon4.p() + sqm1));

	      pxZ14 = muon1.px() + muon4.px();
	      pyZ14 = muon1.py() + muon4.py();
	      pzZ14 = muon1.pz() + muon4.pz();

	      if (muon2.charge() + muon3.charge() == 0)
		{
		  eZ23 = sqrt((muon2.p() * muon2.p() + sqm1)) +
		         (sqrt(muon3.p() * muon3.p() + sqm1));
       
		  pxZ23 = muon2.px() + muon3.px();
		  pyZ23 = muon2.py() + muon3.py();
		  pzZ23 = muon2.pz() + muon3.pz();

		  // Calculate p4
		  pZ14 = sqrt((pxZ14 * pxZ14) + (pyZ14 * pyZ14) + (pzZ14 * pzZ14));
		  pZ23 = sqrt((pxZ23 * pxZ23) + (pyZ23 * pyZ23) + (pzZ23 * pzZ23));

		  pTZ14 = sqrt((pxZ14 * pxZ14) + (pyZ14 * pyZ14));
		  pTZ23 = sqrt((pxZ23 * pxZ23) + (pyZ23 * pyZ23));

		  mZ14 = sqrt((eZ14 * eZ14) - (pZ14 * pZ14));
		  mZ23 = sqrt((eZ23 * eZ23) - (pZ23 * pZ23));

		  if (mZ14 > 0.) h_mZ14_4mu->Fill(mZ14);
		  if (mZ23 > 0.) h_mZ23_4mu->Fill(mZ23);
		}
	    }

	  dZ14 = std::abs( mZ14 - mZ );
	  dZ23 = std::abs( mZ23 - mZ );

	  dZc3 = (dZ14 < dZ23) ? dZ14 : dZ23;

	  bool ptZadaug = false;

	  if (dZc1 < dZc2 && dZc1 < dZc3)
	    {
	      if (dZ12 < dZ34)
		{
		  eZa  = eZ12;     
		  pxZa = pxZ12;
		  pyZa = pyZ12;
		  pzZa = pzZ12;
		  pTZa = pTZ12;
		  mZa  = mZ12;

		  if (muon1.pt() > 20. and muon2.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ34;
		  pxZb = pxZ34;
		  pyZb = pyZ34;
		  pzZb = pzZ34;
		  pTZb = pTZ34;
		  mZb  = mZ34;
		}
	      else
		{
		  eZa  = eZ34;
		  pxZa = pxZ34;
		  pyZa = pyZ34;
		  pzZa = pzZ34;
		  pTZa = pTZ34;
		  mZa  = mZ34;

		  if (muon3.pt() > 20. and muon4.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ12;
		  pxZb = pxZ12;
		  pyZb = pyZ12;
		  pzZb = pzZ12;
		  pTZb = pTZ12;
		  mZb  = mZ12;
		}
	    }

	  else if (dZc2 < dZc1 && dZc2 < dZc3)
	    {
	      if (dZ13 < dZ24)
		{
		  eZa  = eZ13;
		  pxZa = pxZ13;
		  pyZa = pyZ13;
		  pzZa = pzZ13;
		  pTZa = pTZ13;
		  mZa  = mZ13;

		  if (muon1.pt() > 20. and muon3.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ24;
		  pxZb = pxZ24;
		  pyZb = pyZ24;
		  pzZb = pzZ24;
		  pTZb = pTZ24;
		  mZb  = mZ24;
		}
	      else
		{
		  eZa  = eZ24;
		  pxZa = pxZ24;
		  pyZa = pyZ24;
		  pzZa = pzZ24;
		  pTZa = pTZ24;
		  mZa  = mZ24;

		  if (muon2.pt() > 20. and muon4.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ13;
		  pxZb = pxZ13;
		  pyZb = pyZ13;
		  pzZb = pzZ13;
		  pTZb = pTZ13;
		  mZb  = mZ13;
		}
	    }

	  else if (dZc3 < dZc1 && dZc3 < dZc2)
	    {
	      if (dZ14 < dZ23)
		{
		  eZa  = eZ14;
		  pxZa = pxZ14;
		  pyZa = pyZ14;
		  pzZa = pzZ14;
		  pTZa = pTZ14;
		  mZa  = mZ14;

		  if (muon1.pt() > 20. and muon4.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ23;
		  pxZb = pxZ23;
		  pyZb = pyZ23;
		  pzZb = pzZ23;
		  pTZb = pTZ23;
		  mZb  = mZ23;
		}
	      else
		{
		  eZa  = eZ23;
		  pxZa = pxZ23;
		  pyZa = pyZ23;
		  pzZa = pzZ23;
		  pTZa = pTZ23;
		  mZa  = mZ23;

		  if (muon2.pt() > 20. and muon3.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ14;
		  pxZb = pxZ14;
		  pyZb = pyZ14;
		  pzZb = pzZ14;
		  pTZb = pTZ14;
		  mZb  = mZ14;
		}
	    }

	   if (ptZadaug) {
	    if (mZa > 40. && mZa < 120.) {
	      if (mZb > 12. && mZb < 120.) {

		h_mZa_4mu->Fill(mZa);
		h_mZb_4mu->Fill(mZb);

		// 4 vector
		p4Za.SetPxPyPzE(pxZa, pyZa, pzZa, eZa);
		p4Zb.SetPxPyPzE(pxZb, pyZb, pzZb, eZb);

		p4H = p4Za + p4Zb;

		mass4mu = p4H.M();
		pt_4mu = p4H.Pt();
		eta_4mu = p4H.Eta();
		phi_4mu = p4H.Phi();

		px4mu = p4H.Px();
		py4mu = p4H.Py();
		pz4mu = p4H.Pz();
		E4mu = p4H.E();

		pt_mu1 = muon1.pt();
		pt_mu2 = muon2.pt();
		pt_mu3 = muon3.pt();
		pt_mu4 = muon4.pt();

		eta_mu1 = muon1.eta();
		eta_mu2 = muon2.eta();
		eta_mu3 = muon3.eta();
		eta_mu4 = muon4.eta();

		phi_mu1 = muon1.phi();
		phi_mu2 = muon2.phi();
		phi_mu3 = muon3.phi();
		phi_mu4 = muon4.phi();

		cas_mu1 = muon1.charge();
		cas_mu2 = muon2.charge();
		cas_mu3 = muon3.charge();
		cas_mu4 = muon4.charge();

		px_mu1 = muon1.px();
		px_mu2 = muon2.px();
		px_mu3 = muon3.px();
		px_mu4 = muon4.px();

		py_mu1 = muon1.py();
		py_mu2 = muon2.py();
		py_mu3 = muon3.py();
		py_mu4 = muon4.py();
		
		pz_mu1 = muon1.pz();
		pz_mu2 = muon2.pz();
		pz_mu3 = muon3.pz();
		pz_mu4 = muon4.pz();

		E_mu1 = sqrt((muon1.p() * muon1.p()) + sqm1);
		E_mu2 = sqrt((muon2.p() * muon2.p()) + sqm1);
		E_mu3 = sqrt((muon3.p() * muon3.p()) + sqm1);
		E_mu4 = sqrt((muon4.p() * muon4.p()) + sqm1);

		if (mass4mu > 70.)
		  {
		    h_m1_m4mu->Fill(mass4mu);
		    h_m2_m4mu->Fill(mass4mu);
		    h_m3_m4mu->Fill(mass4mu);
		    h_m4_m4mu->Fill(mass4mu);

		    for (unsigned i = 0; i < vIdPtmu.size(); i++)
		      {
			relPFIso_mu = ((((*muons)[vIdPtmu.at(i).first]).pfIsolationR04()).sumChargedHadronPt +
				       (((*muons)[vIdPtmu.at(i).first]).pfIsolationR04()).sumNeutralHadronEt +
				       (((*muons)[vIdPtmu.at(i).first]).pfIsolationR04()).sumPhotonEt) / (((*muons)[vIdPtmu.at(i).first]).pt()); 

			h_relPFIso_mu_after->Fill(relPFIso_mu);

			h_pt_after->Fill(vIdPtmu.at(i).second);
			h_eta_after->Fill(((*muons)[vIdPtmu.at(i).first]).eta());
		      }
		    // t1->Fill();
		  }
	      } // end of mZb
	    } // end of mZa
	   } // end of ptZadaug
	} // end of total charge
    } // end of nGoodRecoMuon

  //============================ ZZ/ZZ*To4Muon end =============================//

  //============================ ZZ/ZZ*To4e start ==============================//

  // Now, for these goodelectrons, pair up and calculate mass
  if (nGoodElectron >= 4)
    {
      const reco::GsfElectron &elec1 = (*electrons)[vIdPte.at(0).first];
      const reco::GsfElectron &elec2 = (*electrons)[vIdPte.at(1).first];
      const reco::GsfElectron &elec3 = (*electrons)[vIdPte.at(2).first];
      const reco::GsfElectron &elec4 = (*electrons)[vIdPte.at(3).first];

      if (elec1.charge() + elec2.charge() + elec3.charge() + elec4.charge() == 0)
	{
	  // First combination: Combine elec 1234
	  if (elec1.charge() + elec2.charge() == 0)
	    {
	      eZ12 = (sqrt(elec1.p() * elec1.p() + sqme)) + (sqrt(elec2.p() * elec2.p() + sqme));
	
	      pxZ12 = elec1.px() + elec2.px();
	      pyZ12 = elec1.py() + elec2.py();
	      pzZ12 = elec1.pz() + elec2.pz();

	      if (elec3.charge() + elec4.charge() == 0)
		{
		  eZ34 = (sqrt(elec3.p() * elec3.p() + sqme)) + (sqrt(elec4.p() * elec4.p() + sqme));

		  pxZ34 = elec3.px() + elec4.px();
		  pyZ34 = elec3.py() + elec4.py();
		  pzZ34 = elec3.pz() + elec4.pz();

		  // Calculate the momentum and invariant mass of elec 1 and 2, elec 3 and 4
		  pZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12) + (pzZ12 * pzZ12));
		  pZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34) + (pzZ34 * pzZ34));

		  pTZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12));
		  pTZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34));

		  mZ12 = sqrt((eZ12 * eZ12) - (pZ12 * pZ12));
		  mZ34 = sqrt((eZ34 * eZ34) - (pZ34 * pZ34));

		  if (mZ12 > 0.) h_mZ12_4e->Fill(mZ12);
		  if (mZ34 > 0.) h_mZ34_4e->Fill(mZ34);
		}
	    }

	  dZ12 = std::abs( mZ12 - mZ );
	  dZ34 = std::abs( mZ34 - mZ );

	  // take the smallest diff between mass to use for 4electron combination
	  dZc1 = (dZ12 < dZ34) ? dZ12 : dZ34; 

	  // Second combination: Combine elec 1324
	  if (elec1.charge() + elec3.charge() == 0)
	    {
	      eZ13 = (sqrt(elec1.p() * elec1.p() + sqme)) + (sqrt(elec3.p() * elec3.p() + sqme));

	      pxZ13 = elec1.px() + elec3.px();
	      pyZ13 = elec1.py() + elec3.py();
	      pzZ13 = elec1.pz() + elec3.pz();

	      if (elec2.charge() + elec4.charge() == 0)
		{

		  eZ24 = (sqrt(elec2.p() * elec2.p() + sqme)) + (sqrt(elec4.p() * elec4.p() + sqme));
		  
		  pxZ24 = elec2.px() + elec4.px();
		  pyZ24 = elec2.py() + elec4.py();
		  pzZ24 = elec2.pz() + elec4.pz();

		  // Calculate the momentum and invariant mass of elec 1 and 3, elec 2 and 4
		  pZ13 = sqrt((pxZ13 * pxZ13) + (pyZ13 * pyZ13) + (pzZ13 * pzZ13));
		  pZ24 = sqrt((pxZ24 * pxZ24) + (pyZ24 * pyZ24) + (pzZ24 * pzZ24));

		  pTZ13 = sqrt((pxZ13 * pxZ13) + (pyZ13 * pyZ13) + (pzZ13 * pzZ13));
		  pTZ24 = sqrt((pxZ24 * pxZ24) + (pyZ24 * pyZ24) + (pzZ24 * pzZ24));

		  mZ13 = sqrt((eZ13 * eZ13) - (pZ13 * pZ13));
		  mZ24 = sqrt((eZ24 * eZ24) - (pZ24 * pZ24));

		  if (mZ13 > 0.) h_mZ13_4e->Fill(mZ13);
		  if (mZ24 > 0.) h_mZ24_4e->Fill(mZ24);
		}
	    }

	  dZ13 = std::abs( mZ13 - mZ );
	  dZ24 = std::abs( mZ24 - mZ );

	  dZc2 = (dZ13 < dZ24) ? dZ13 : dZ24; 

	  // Third combination: Combine elec 1423
	  if (elec1.charge() + elec4.charge() == 0)
	    {
	      eZ14 = (sqrt(elec1.p() * elec1.p() + sqme)) + (sqrt(elec4.p() * elec4.p() + sqme));

	      pxZ14 = elec1.px() + elec4.px();
	      pyZ14 = elec1.py() + elec4.py();
	      pzZ14 = elec1.pz() + elec4.pz();

	      if (elec2.charge() + elec3.charge() == 0)
		{
		  eZ23 = sqrt((elec2.p() * elec2.p() + sqme)) + (sqrt(elec3.p() * elec3.p() + sqme));
       
		  pxZ23 = elec2.px() + elec3.px();
		  pyZ23 = elec2.py() + elec3.py();
		  pzZ23 = elec2.pz() + elec3.pz();

		  // Calculate the momentum and invariant mass of elec 1 and 4, elec 2 and 3
		  pZ14 = sqrt((pxZ14 * pxZ14) + (pyZ14 * pyZ14) + (pzZ14 * pzZ14));
		  pZ23 = sqrt((pxZ23 * pxZ23) + (pyZ23 * pyZ23) + (pzZ23 * pzZ23));

		  pTZ14 = sqrt((pxZ14 * pxZ14) + (pyZ14 * pyZ14) + (pzZ14 * pzZ14));
		  pTZ23 = sqrt((pxZ23 * pxZ23) + (pyZ23 * pyZ23) + (pzZ23 * pzZ23));

		  mZ14 = sqrt((eZ14 * eZ14) - (pZ14 * pZ14));
		  mZ23 = sqrt((eZ23 * eZ23) - (pZ23 * pZ23));

		  if (mZ14 > 0.) h_mZ14_4e->Fill(mZ14);
		  if (mZ23 > 0.) h_mZ23_4e->Fill(mZ23);
		}
	    }

	  dZ14 = std::abs( mZ14 - mZ );
	  dZ23 = std::abs( mZ23 - mZ );

	  dZc3 = (dZ14 < dZ23) ? dZ14 : dZ23; 

	  bool ptZadaug = false;

	  // Now whichever have the smallest diff is considered the best comb. 
	  if (dZc1 < dZc2 && dZc1 < dZc3)
	    {
	      if (dZ12 < dZ34)
		{
		  eZa  = eZ12;     
		  pxZa = pxZ12;
		  pyZa = pyZ12;
		  pzZa = pzZ12;
		  pTZa = pTZ12;
		  mZa  = mZ12;

		  if (elec1.pt() > 20. and elec2.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ34;
		  pxZb = pxZ34;
		  pyZb = pyZ34;
		  pzZb = pzZ34;
		  pTZb = pTZ34;
		  mZb  = mZ34;
		}
	      else
		{
		  eZa  = eZ34;  
		  pxZa = pxZ34;
		  pyZa = pyZ34;
		  pzZa = pzZ34;
		  pTZa = pTZ34;
		  mZa  = mZ34;

		  if (elec3.pt() > 20. and elec4.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ12;
		  pxZb = pxZ12;
		  pyZb = pyZ12;
		  pzZb = pzZ12;
		  pTZb = pTZ12;
		  mZb  = mZ12;
		}
	    }

	  else if (dZc2 < dZc1 && dZc2 < dZc3)
	    {
	      if (dZ13 < dZ24)
		{
		  eZa  = eZ13;
		  pxZa = pxZ13;
		  pyZa = pyZ13;
		  pzZa = pzZ13;
		  pTZa = pTZ13;
		  mZa  = mZ13;

		  if (elec1.pt() > 20. and elec3.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ24;
		  pxZb = pxZ24;
		  pyZb = pyZ24;
		  pzZb = pzZ24;
		  pTZb = pTZ24;
		  mZb  = mZ24;
		}
	      else
		{
		  eZa  = eZ24;
		  pxZa = pxZ24;
		  pyZa = pyZ24;
		  pzZa = pzZ24;
		  pTZa = pTZ24;
		  mZa  = mZ24;

		  if (elec2.pt() > 20. and elec4.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ13;
		  pxZb = pxZ13;
		  pyZb = pyZ13;
		  pzZb = pzZ13;
		  pTZb = pTZ13;
		  mZb  = mZ13;
		}
	    }

	  else if (dZc3 < dZc1 && dZc3 < dZc2)
	    {
	      if (dZ14 < dZ23)
		{
		  eZa  = eZ14;
		  pxZa = pxZ14;
		  pyZa = pyZ14;
		  pzZa = pzZ14;
		  pTZa = pTZ14;
		  mZa  = mZ14;

		  if (elec1.pt() > 20. and elec4.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ23;
		  pxZb = pxZ23;
		  pyZb = pyZ23;
		  pzZb = pzZ23;
		  pTZb = pTZ23;
		  mZb  = mZ23;
		}
	      else
		{
		  eZa  = eZ23;
		  pxZa = pxZ23;
		  pyZa = pyZ23;
		  pzZa = pzZ23;
		  pTZa = pTZ23;
		  mZa  = mZ23;

		  if (elec2.pt() > 20. and elec3.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ14;
		  pxZb = pxZ14;
		  pyZb = pyZ14;
		  pzZb = pzZ14;
		  pTZb = pTZ14;
		  mZb  = mZ14;
		}
	    }

	  if (ptZadaug) {
	    if (mZa > 40. && mZa < 120.) {
	      if (mZb > 12. && mZb < 120.) {
		h_mZa_4e->Fill(mZa);
		h_mZb_4e->Fill(mZb);
		
		// Calculate 4 elec
		p4Za.SetPxPyPzE(pxZa, pyZa, pzZa, eZa);
		p4Zb.SetPxPyPzE(pxZb, pyZb, pzZb, eZb);

		p4H = p4Za + p4Zb;

		mass4e = p4H.M();
		pt_4e = p4H.Pt();
		eta_4e = p4H.Eta();
		phi_4e = p4H.Phi();
		
		px4e = p4H.Px();
		py4e = p4H.Py();
		pz4e = p4H.Pz();
		E4e = p4H.E();		

		pt_e1 = elec1.pt();
		pt_e2 = elec2.pt();
		pt_e3 = elec3.pt();
		pt_e4 = elec4.pt();

		eta_e1 = elec1.eta();
		eta_e2 = elec2.eta();
		eta_e3 = elec3.eta();
		eta_e4 = elec4.eta();

		phi_e1 = elec1.phi();
		phi_e2 = elec2.phi();
		phi_e3 = elec3.phi();
		phi_e4 = elec4.phi();

		cas_e1 = elec1.charge();
		cas_e2 = elec2.charge();
		cas_e3 = elec3.charge();
		cas_e4 = elec4.charge();

		px_e1 = elec1.px();
		px_e2 = elec2.px();
		px_e3 = elec3.px();
		px_e4 = elec4.px();

		py_e1 = elec1.py();
		py_e2 = elec2.py();
		py_e3 = elec3.py();
		py_e4 = elec4.py();
		
		pz_e1 = elec1.pz();
		pz_e2 = elec2.pz();
		pz_e3 = elec3.pz();
		pz_e4 = elec4.pz();

		E_e1 = elec1.energy();
		E_e2 = elec2.energy();
		E_e3 = elec3.energy();
		E_e4 = elec4.energy();

		if (mass4e > 70.)
		  {
		    h_m1_m4e->Fill(mass4e);
		    h_m2_m4e->Fill(mass4e);
		    h_m3_m4e->Fill(mass4e);
		    h_m4_m4e->Fill(mass4e);

		    for (unsigned i = 0; i < vIdPte.size(); i++)
		      {
			relPFIso_e = ((((*electrons)[vIdPte.at(i).first]).pfIsolationVariables()).chargedHadronIso +
				      (((*electrons)[vIdPte.at(i).first]).pfIsolationVariables()).neutralHadronIso +
				      (((*electrons)[vIdPte.at(i).first]).pfIsolationVariables()).photonIso) / (((*electrons)[vIdPte.at(i).first]).pt());

			h_relPFIso_e_after->Fill(relPFIso_e);

			h_pt_e_after->Fill(vIdPte.at(i).second);
			h_eta_e_after->Fill((((*electrons)[vIdPte.at(i).first]).superCluster())->eta());
		      }
		    // t2->Fill();
		  }
	      } // end of mZb 
	    } // end of mZa
	  } // end of ptZdaug
	} // end of total charge
    } // end of nGoodElectron

  //=============================== ZZ/ZZ*To4e end =============================//


  //=========================== ZZ/ZZ*To2mu2e start ============================//

  if (nGoodRecoMuon >= 2 && nGoodElectron >= 2)
    {
      const reco::Muon &muon1 = (*muons)[vIdPtmu.at(0).first];
      const reco::Muon &muon2 = (*muons)[vIdPtmu.at(1).first];
      const reco::GsfElectron &elec1 = (*electrons)[vIdPte.at(0).first];
      const reco::GsfElectron &elec2 = (*electrons)[vIdPte.at(1).first];

      if (muon1.charge() + muon2.charge() + elec1.charge() + elec2.charge() == 0)
	{
	  // For case 2mu2e, there is only 1 combination
	  if (muon1.charge() + muon2.charge() == 0)
	    {
	      eZ12 = (sqrt(muon1.p() * muon1.p() + sqm1)) + (sqrt(muon2.p() * muon2.p() + sqm1));
	
	      pxZ12 = muon1.px() + muon2.px();
	      pyZ12 = muon1.py() + muon2.py();
	      pzZ12 = muon1.pz() + muon2.pz();

	      if (elec1.charge() + elec2.charge() == 0)
		{
		  eZ34 = (sqrt(elec1.p() * elec1.p() + sqme)) + (sqrt(elec2.p() * elec2.p() + sqme));

		  pxZ34 = elec1.px() + elec2.px();
		  pyZ34 = elec1.py() + elec2.py();
		  pzZ34 = elec1.pz() + elec2.pz();

		  // Calculate the momentum and invariant mass of muon 12, elec 34
		  pZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12) + (pzZ12 * pzZ12));
		  pZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34) + (pzZ34 * pzZ34));

		  pTZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12));
		  pTZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34));

		  mZ12 = sqrt((eZ12 * eZ12) - (pZ12 * pZ12));
		  mZ34 = sqrt((eZ34 * eZ34) - (pZ34 * pZ34));

		  if (mZ12 > 0.) h_mZmu_2mu2e->Fill(mZ12);
		  if (mZ34 > 0.) h_mZe_2mu2e->Fill(mZ34);

		}
	    }

	  dZ12 = std::abs(mZ12 - mZ); // mu
	  dZ34 = std::abs(mZ34 - mZ); // e

	  bool ptZadaug = false;

	  if (dZ12 < dZ34)
	    {
	      eZa  = eZ12;
	      pxZa = pxZ12;
	      pyZa = pyZ12;
	      pzZa = pzZ12;
	      pTZa = pTZ12;
	      mZa  = mZ12;

	      if (muon1.pt() > 20. and muon2.pt() > 10.)
		ptZadaug = true;

	      eZb  = eZ34;
	      pxZb = pxZ34;
	      pyZb = pyZ34;
	      pzZb = pzZ34;
	      pTZb = pTZ34;
	      mZb  = mZ34;
	    }
	  else
	    {
	      eZa  = eZ34;
	      pxZa = pxZ34;
	      pyZa = pyZ34;
	      pzZa = pzZ34;
	      pTZa = pTZ34;
	      mZa  = mZ34;

	      if (elec1.pt() > 20. and elec2.pt() > 10.)
		ptZadaug = true;

	      eZb  = eZ12;
	      pxZb = pxZ12;
	      pyZb = pyZ12;
	      pzZb = pzZ12;
	      pTZb = pTZ12;
	      mZb  = mZ12;
	  }

	  if (ptZadaug) {
	    if (mZa > 40. && mZa < 120.) {
	      if (mZb > 12. && mZb < 120.) {
		h_mZa_2mu2e->Fill(mZa);
		h_mZb_2mu2e->Fill(mZb);

		// Now combine these 2 muons and 2 electrons
 
		// Calculate 4 lepton: 2muons 2electrons

		p4Za.SetPxPyPzE(pxZa, pyZa, pzZa, eZa);
		p4Zb.SetPxPyPzE(pxZb, pyZb, pzZb, eZb);

		p4H = p4Za + p4Zb;

		mass2mu2e = p4H.M();
		pt_2mu2e = p4H.Pt();
		eta_2mu2e = p4H.Eta();
		phi_2mu2e = p4H.Phi();
		
		px2mu2e = p4H.Px();
		py2mu2e = p4H.Py();
		pz2mu2e = p4H.Pz();
		E2mu2e = p4H.E();

		pt_2mu1 = muon1.pt();
		pt_2mu2 = muon2.pt();
		pt_2e1 = elec1.pt();
		pt_2e2 = elec2.pt();

		eta_2mu1 = muon1.eta();
		eta_2mu2 = muon2.eta();
		eta_2e1 = elec1.eta();
		eta_2e2 = elec2.eta();

		phi_2mu1 = muon1.phi();
		phi_2mu2 = muon2.phi();
		phi_2e1 = elec1.phi();
		phi_2e2 = elec2.phi();

		cas_2mu1 = muon1.charge();
		cas_2mu2 = muon2.charge();
		cas_2e1 = elec1.charge();
		cas_2e2 = elec2.charge();

		px_2mu1 = muon1.px();
		px_2mu2 = muon2.px();
		px_2e1 = elec1.px();
		px_2e2 = elec2.px();

		py_2mu1 = muon1.py();
		py_2mu2 = muon2.py();
		py_2e1 = elec1.py();
		py_2e2 = elec2.py();
		
		pz_2mu1 = muon1.pz();
		pz_2mu2 = muon2.pz();
		pz_2e1 = elec1.pz();
		pz_2e2 = elec2.pz();

		E_2mu1 = sqrt((muon1.p() * muon1.p()) + sqm1);
		E_2mu2 = sqrt((muon2.p() * muon2.p()) + sqm1);
		E_2e1 = elec1.energy();
		E_2e2 = elec2.energy();

		if (mass2mu2e > 70.)
		  {
		    h_m1_m2mu2e->Fill(mass2mu2e);
		    h_m2_m2mu2e->Fill(mass2mu2e);
		    h_m3_m2mu2e->Fill(mass2mu2e);
		    h_m4_m2mu2e->Fill(mass2mu2e);

		    for (unsigned i = 0; i < vIdPtmu.size(); i++)
		      {
			relPFIso_mu = ( (((*muons)[vIdPtmu.at(i).first]).pfIsolationR04()).sumChargedHadronPt +
					(((*muons)[vIdPtmu.at(i).first]).pfIsolationR04()).sumNeutralHadronEt +
					(((*muons)[vIdPtmu.at(i).first]).pfIsolationR04()).sumPhotonEt ) / (((*muons)[vIdPtmu.at(i).first]).pt());

			h_relPFIso_2mu_after->Fill(relPFIso_mu);
			h_pt_after_2mu2e->Fill(vIdPtmu.at(i).second);
			h_eta_after_2mu2e->Fill(((*muons)[vIdPtmu.at(i).first]).eta());
		      }

		    for (unsigned i = 0; i < vIdPte.size(); i++)
		      {
			relPFIso_e = ((((*electrons)[vIdPte.at(i).first]).pfIsolationVariables()).chargedHadronIso +
				      (((*electrons)[vIdPte.at(i).first]).pfIsolationVariables()).neutralHadronIso +
				      (((*electrons)[vIdPte.at(i).first]).pfIsolationVariables()).photonIso) / (((*electrons)[vIdPte.at(i).first]).pt()); 
	      
			h_relPFIso_2e_after->Fill(relPFIso_e);
			h_pt_e_after_2mu2e->Fill(vIdPte.at(i).second);
			h_eta_e_after_2mu2e->Fill(((*electrons)[vIdPte.at(i).first]).eta());
		      }
		    // t3->Fill();
		  }
	      }
	    }
	  }
	} // end of total charge
    } // end of nGoodMuonNElectron

  //============================= ZZ/ZZ*To2mu2e end ============================//
  
} // HiggsDemoAnalyzerGit::analyze ends


// ------ method called once each job just before starting event loop ---------//

void HiggsDemoAnalyzerGit::beginJob() {

  // *******************************************************
  // book the ntuple for the surviving 4 lepton candidates *
  // in the mass range 70 < m4l < 181 GeV                  *
  // *******************************************************

/*
  t1 = new TTree("tree4mu", "tree4mu");
  t2 = new TTree("tree4e", "tree4e");
  t3 = new TTree("tree2mu2e", "tree2mu2e");

  // tree 4mu
  t1->Branch("nRun", &nRun, "nRun/I");
  t1->Branch("nEvt", &nEvt, "nEvt/I");
  t1->Branch("nLumi", &nLumi, "nLumi/I");
  t1->Branch("mass4mu", &mass4mu, "mass4mu/D");
  t1->Branch("pt_4mu", &pt_4mu, "pt_4mu/D");
  t1->Branch("eta_4mu", &eta_4mu, "eta_4mu/D");
  t1->Branch("phi_4mu", &phi_4mu, "phi_4mu/D");

  t1->Branch("px4mu", &px4mu, "px4mu/D");
  t1->Branch("py4mu", &py4mu, "py4mu/D");
  t1->Branch("pz4mu", &pz4mu, "pz4mu/D");
  t1->Branch("E4mu", &E4mu, "E4mu/D");
  
  t1->Branch("pt_mu1", &pt_mu1, "pt_mu1/D");
  t1->Branch("pt_mu2", &pt_mu2, "pt_mu2/D");
  t1->Branch("pt_mu3", &pt_mu3, "pt_mu3/D");
  t1->Branch("pt_mu4", &pt_mu4, "pt_mu4/D");
  t1->Branch("eta_mu1", &eta_mu1, "eta_mu1/D");
  t1->Branch("eta_mu2", &eta_mu2, "eta_mu2/D");
  t1->Branch("eta_mu3", &eta_mu3, "eta_mu3/D");
  t1->Branch("eta_mu4", &eta_mu4, "eta_mu4/D");

  t1->Branch("phi_mu1", &phi_mu1, "phi_mu1/D");
  t1->Branch("phi_mu2", &phi_mu2, "phi_mu2/D");
  t1->Branch("phi_mu3", &phi_mu3, "phi_mu3/D");
  t1->Branch("phi_mu4", &phi_mu4, "phi_mu4/D");

  t1->Branch("cas_mu1", &cas_mu1, "cas_mu1/I");
  t1->Branch("cas_mu2", &cas_mu2, "cas_mu2/I");
  t1->Branch("cas_mu3", &cas_mu3, "cas_mu3/I");
  t1->Branch("cas_mu4", &cas_mu4, "cas_mu4/I");

  t1->Branch("px_mu1", &px_mu1, "px_mu1/D");
  t1->Branch("px_mu2", &px_mu2, "px_mu2/D");
  t1->Branch("px_mu3", &px_mu3, "px_mu3/D");
  t1->Branch("px_mu4", &px_mu4, "px_mu4/D");

  t1->Branch("py_mu1", &py_mu1, "py_mu1/D");
  t1->Branch("py_mu2", &py_mu2, "py_mu2/D");
  t1->Branch("py_mu3", &py_mu3, "py_mu3/D");
  t1->Branch("py_mu4", &py_mu4, "py_mu4/D");

  t1->Branch("pz_mu1", &pz_mu1, "pz_mu1/D");
  t1->Branch("pz_mu2", &pz_mu2, "pz_mu2/D");
  t1->Branch("pz_mu3", &pz_mu3, "pz_mu3/D");
  t1->Branch("pz_mu4", &pz_mu4, "pz_mu4/D");

  t1->Branch("E_mu1", &E_mu1, "E_mu1/D");
  t1->Branch("E_mu2", &E_mu2, "E_mu2/D");
  t1->Branch("E_mu3", &E_mu3, "E_mu3/D");
  t1->Branch("E_mu4", &E_mu4, "E_mu4/D");

  t1->Branch("mZa", &mZa, "mZa/D");
  t1->Branch("mZb", &mZb, "mZb/D");
  
  // tree 4e
  t2->Branch("nRun", &nRun, "nRun/I");
  t2->Branch("nEvt", &nEvt, "nEvt/I");
  t2->Branch("nLumi", &nLumi, "nLumi/I");
  t2->Branch("mass4e", &mass4e, "mass4e/D");
  t2->Branch("pt_4e", &pt_4e, "pt_4e/D");
  t2->Branch("eta_4e", &eta_4e, "eta_4e/D");
  t2->Branch("phi_4e", &phi_4e, "phi_4e/D");

  t2->Branch("px4e", &px4e, "px4e/D");
  t2->Branch("py4e", &py4e, "py4e/D");
  t2->Branch("pz4e", &pz4e, "pz4e/D");
  t2->Branch("E4e", &E4e, "E4e/D");
  
  t2->Branch("pt_e1", &pt_e1, "pt_e1/D");
  t2->Branch("pt_e2", &pt_e2, "pt_e2/D");
  t2->Branch("pt_e3", &pt_e3, "pt_e3/D");
  t2->Branch("pt_e4", &pt_e4, "pt_e4/D");
  t2->Branch("eta_e1", &eta_e1, "eta_e1/D");
  t2->Branch("eta_e2", &eta_e2, "eta_e2/D");
  t2->Branch("eta_e3", &eta_e3, "eta_e3/D");
  t2->Branch("eta_e4", &eta_e4, "eta_e4/D");

  t2->Branch("phi_e1", &phi_e1, "phi_e1/D");
  t2->Branch("phi_e2", &phi_e2, "phi_e2/D");
  t2->Branch("phi_e3", &phi_e3, "phi_e3/D");
  t2->Branch("phi_e4", &phi_e4, "phi_e4/D");

  t2->Branch("cas_e1", &cas_e1, "cas_e1/I");
  t2->Branch("cas_e2", &cas_e2, "cas_e2/I");
  t2->Branch("cas_e3", &cas_e3, "cas_e3/I");
  t2->Branch("cas_e4", &cas_e4, "cas_e4/I");

  t2->Branch("px_e1", &px_e1, "px_e1/D");
  t2->Branch("px_e2", &px_e2, "px_e2/D");
  t2->Branch("px_e3", &px_e3, "px_e3/D");
  t2->Branch("px_e4", &px_e4, "px_e4/D");

  t2->Branch("py_e1", &py_e1, "py_e1/D");
  t2->Branch("py_e2", &py_e2, "py_e2/D");
  t2->Branch("py_e3", &py_e3, "py_e3/D");
  t2->Branch("py_e4", &py_e4, "py_e4/D");

  t2->Branch("pz_e1", &pz_e1, "pz_e1/D");
  t2->Branch("pz_e2", &pz_e2, "pz_e2/D");
  t2->Branch("pz_e3", &pz_e3, "pz_e3/D");
  t2->Branch("pz_e4", &pz_e4, "pz_e4/D");

  t2->Branch("E_e1", &E_e1, "E_e1/D");
  t2->Branch("E_e2", &E_e2, "E_e2/D");
  t2->Branch("E_e3", &E_e3, "E_e3/D");
  t2->Branch("E_e4", &E_e4, "E_e4/D");

  t2->Branch("mZa", &mZa, "mZa/D");
  t2->Branch("mZb", &mZb, "mZb/D");
  
  // tree 2mu 2e
  t3->Branch("nRun", &nRun, "nRun/I");
  t3->Branch("nEvt", &nEvt, "nEvt/I");
  t3->Branch("nLumi", &nLumi, "nLumi/I");
  t3->Branch("mass2mu2e", &mass2mu2e, "mass2mu2e/D");
  t3->Branch("pt_2mu2e", &pt_2mu2e, "pt_2mu2e/D");
  t3->Branch("eta_2mu2e", &eta_2mu2e, "eta_2mu2e/D");
  t3->Branch("phi_2mu2e", &phi_2mu2e, "phi_2mu2e/D");

  t3->Branch("px2mu2e", &px2mu2e, "px2mu2e/D");
  t3->Branch("py2mu2e", &py2mu2e, "py2mu2e/D");
  t3->Branch("pz2mu2e", &pz2mu2e, "pz2mu2e/D");
  t3->Branch("E2mu2e", &E2mu2e, "E2mu2e/D");
  
  t3->Branch("pt_2mu1", &pt_2mu1, "pt_2mu1/D");
  t3->Branch("pt_2mu2", &pt_2mu2, "pt_2mu2/D");
  t3->Branch("pt_2e1", &pt_2e1, "pt_2e1/D");
  t3->Branch("pt_2e2", &pt_2e2, "pt_2e2/D");
  t3->Branch("eta_2mu1", &eta_2mu1, "eta_2mu1/D");
  t3->Branch("eta_2mu2", &eta_2mu2, "eta_2mu2/D");
  t3->Branch("eta_2e1", &eta_2e1, "eta_2e1/D");
  t3->Branch("eta_2e2", &eta_2e2, "eta_2e2/D");

  t3->Branch("phi_2mu1", &phi_2mu1, "phi_2mu1/D");
  t3->Branch("phi_2mu2", &phi_2mu2, "phi_2mu2/D");
  t3->Branch("phi_2e1", &phi_2e1, "phi_2e1/D");
  t3->Branch("phi_2e2", &phi_2e2, "phi_2e2/D");

  t3->Branch("cas_2mu1", &cas_2mu1, "cas_2mu1/I");
  t3->Branch("cas_2mu2", &cas_2mu2, "cas_2mu2/I");
  t3->Branch("cas_2e1", &cas_2e1, "cas_2e1/I");
  t3->Branch("cas_2e2", &cas_2e2, "cas_2e2/I");

  t3->Branch("px_2mu1", &px_2mu1, "px_2mu1/D");
  t3->Branch("px_2mu2", &px_2mu2, "px_2mu2/D");
  t3->Branch("px_2e1", &px_2e1, "px_2e1/D");
  t3->Branch("px_2e2", &px_2e2, "px_2e2/D");

  t3->Branch("py_2mu1", &py_2mu1, "py_2mu1/D");
  t3->Branch("py_2mu2", &py_2mu2, "py_2mu2/D");
  t3->Branch("py_2e1", &py_2e1, "py_2e1/D");
  t3->Branch("py_2e2", &py_2e2, "py_2e2/D");

  t3->Branch("pz_2mu1", &pz_2mu1, "pz_2mu1/D");
  t3->Branch("pz_2mu2", &pz_2mu2, "pz_2mu2/D");
  t3->Branch("pz_2e1", &pz_2e1, "pz_2e1/D");
  t3->Branch("pz_2e2", &pz_2e2, "pz_2e2/D");

  t3->Branch("E_2mu1", &E_2mu1, "E_2mu1/D");
  t3->Branch("E_2mu2", &E_2mu2, "E_2mu2/D");
  t3->Branch("E_2e1", &E_2e1, "E_2e1/D");
  t3->Branch("E_2e2", &E_2e2, "E_2e2/D");

  t3->Branch("mZa", &mZa, "mZa/D");
  t3->Branch("mZb", &mZb, "mZb/D");
  */
  
}

// ------------ method called once each job just after ending the event loop  ------------
void HiggsDemoAnalyzerGit::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiggsDemoAnalyzerGit);

