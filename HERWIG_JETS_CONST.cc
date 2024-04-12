// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"

#include "cmath"
#include "iostream"
#include "fstream"
#include "random"
#include "vector"
#include "stdio.h"
#include "stdlib.h"

using namespace std;

namespace Rivet
{
  class HERWIG_JETS_CONST : public Analysis
  {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(HERWIG_JETS_CONST);

    void init()
    {
      const FinalState fs(Cuts::pt > 0.2 && Cuts::pt < 30 && Cuts::abseta < 1.0);
      declare(fs, "fs");

      // open output file
      mytxtfile.open("herwig_11pthat15_evt_particles.txt");
    }
    /// Perform the per-event analysis
    void analyze(const Event &event)
    {

      //===========================================
      // process particles
      Particles fsParticles = applyProjection<FinalState>(event, "fs").particles();
      PseudoJets parts;

      for (const Particle &p : fsParticles)
      {
        if (p.pid() == 12 || p.pid() == 14 || p.pid() == 16)
          continue;
        PseudoJet pseudojet(p.px(), p.py(), p.pz(), p.E());
        parts.push_back(pseudojet);
      }

      //===========================================
      // jet selectors
      fastjet::Selector selector = fastjet::SelectorPtMin(15.0) * fastjet::SelectorEtaRange(-0.6, 0.6);

      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4);
      fastjet::ClusterSequence cs(parts, jet_def);
      vector<fastjet::PseudoJet> jets = sorted_by_pt(selector(cs.inclusive_jets()));
     
      if (jets.size() != 0)
      {
	for (unsigned int i = 0; i < parts.size(); i++)
	{
	  mytxtfile << parts.at(i).perp() << endl;
	}
      }	
      /*for (unsigned int i = 0; i < jets.size(); i++)
      {
        for (unsigned int j = 0; j < jets[i].constituents().size(); j++)  
        {
          mytxtfile << jets[i].constituents().at(j).perp() << endl;
        }
      }*/
    }

    /// Normalise histograms etc., after the run
    void finalize()
    {
      mytxtfile.close();
    }

  private:
    std::ofstream mytxtfile;
  };

  RIVET_DECLARE_PLUGIN(HERWIG_JETS_CONST);
}
